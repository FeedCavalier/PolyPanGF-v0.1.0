"""Optional enrichment, expression and author-marker analyses."""

from collections import defaultdict
from itertools import combinations
from pathlib import Path
import math
import numpy as np
from scipy.stats import fisher_exact

from .core import InputError, read_tsv, resolve, save_json, write_tsv
from .family import resolve_member
from .pan import CLASSES, bh_adjust


def loci(cfg):
    return read_tsv(Path(cfg["output"]) / "pan/loci.tsv")


def enrichment(cfg, dest):
    settings = cfg["enrichment"]
    sample = settings["dataset"]
    rows = [row for row in loci(cfg) if row["dataset"] == sample]
    if not rows:
        raise InputError(f"No analysis candidates for enrichment dataset {sample}")
    annotation = read_tsv(
        resolve(cfg["_base"], settings["term2gene"]), ["gene_id", "term"]
    )
    permitted = None
    if settings.get("valid_terms"):
        permitted = {
            row["term"]
            for row in read_tsv(
                resolve(cfg["_base"], settings["valid_terms"]), ["term"]
            )
        }
    terms = defaultdict(set)
    names = {}
    for row in annotation:
        if permitted is not None and row["term"] not in permitted:
            continue
        terms[row["term"]].add(row["gene_id"])
        names[row["term"]] = row.get("description", row["term"])
    universe = set().union(*terms.values()) if terms else set()
    if settings.get("background"):
        universe &= {
            row["gene_id"]
            for row in read_tsv(
                resolve(cfg["_base"], settings["background"]), ["gene_id"]
            )
        }
    foregrounds = {"All": {row["gene_id"] for row in rows}}
    foregrounds.update(
        {
            category: {
                row["gene_id"] for row in rows if row["coverage_class"] == category
            }
            for category in CLASSES
        }
    )
    tests = []
    min_size = int(settings.get("min_term_size", 2))
    for category, genes in foregrounds.items():
        genes &= universe
        if not genes:
            continue
        for term, members in sorted(terms.items()):
            members &= universe
            if not min_size <= len(members) < len(universe):
                continue
            overlap = len(genes & members)
            table = [
                [overlap, len(genes) - overlap],
                [len(members) - overlap, len(universe) - len(genes | members)],
            ]
            pvalue = fisher_exact(table, alternative="greater")[1]
            tests.append(
                {
                    "scope": "all_family" if category == "All" else "coverage_classes",
                    "category": category,
                    "term": term,
                    "description": names[term],
                    "overlap": overlap,
                    "foreground": len(genes),
                    "background_term": len(members),
                    "background": len(universe),
                    "pvalue": pvalue,
                    "qvalue": 1.0,
                    "genes": ",".join(sorted(genes & members)),
                }
            )
    for scope in ["all_family", "coverage_classes"]:
        subset = [row for row in tests if row["scope"] == scope]
        for row, qvalue in zip(subset, bh_adjust([row["pvalue"] for row in subset])):
            row["qvalue"] = float(qvalue)
    fields = [
        "scope",
        "category",
        "term",
        "description",
        "overlap",
        "foreground",
        "background_term",
        "background",
        "pvalue",
        "qvalue",
        "genes",
    ]
    write_tsv(dest / "all_tests.tsv", tests, fields)
    significant = [
        row for row in tests if row["qvalue"] <= float(settings.get("fdr", 0.05))
    ]
    write_tsv(dest / "significant.tsv", significant, fields)
    result = {
        "dataset": sample,
        "database": settings.get("database", "custom"),
        "annotated_background": len(universe),
        "annotated_candidates": len(foregrounds["All"] & universe),
        "tests": len(tests),
        "significant": len(significant),
        "multiplicity": "all-family and joint coverage-class scopes separately, including zero-hit terms",
    }
    save_json(dest / "summary.json", result)
    return result


def expression(cfg, dest):
    settings = cfg["expression"]
    dataset = settings["dataset"]
    candidate_rows = [row for row in loci(cfg) if row["dataset"] == dataset]
    data = read_tsv(resolve(cfg["_base"], settings["matrix"]), ["gene_id"])
    metadata = read_tsv(resolve(cfg["_base"], settings["samples"]), ["sample", "stage"])
    sample_names = [row["sample"] for row in metadata]
    if len(set(sample_names)) != len(sample_names) or len(sample_names) < 2:
        raise InputError(
            "Expression metadata requires unique sample IDs and at least two samples"
        )
    stages = list(dict.fromkeys(row["stage"] for row in metadata))
    if len(stages) < 2 or any(not row["stage"] for row in metadata):
        raise InputError("At least two nonempty expression stages required")
    tracks = {}
    seen = set()
    for row in data:
        if row["gene_id"] in seen:
            raise InputError("Duplicate expression gene ID: " + row["gene_id"])
        seen.add(row["gene_id"])
        possible = [
            r
            for r in candidate_rows
            if row["gene_id"] in [r["uid"], r["gene_id"], r["protein_id"]]
        ]
        if not possible:
            continue
        uid = resolve_member(candidate_rows, dataset, row["gene_id"])
        if uid in tracks:
            raise InputError(f"Multiple expression IDs resolve to {uid}")
        try:
            values = np.array([float(row[sample]) for sample in sample_names])
        except (KeyError, ValueError) as exc:
            raise InputError(
                "Expression matrix must contain finite numeric values for every sample"
            ) from exc
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise InputError("Expression values must be finite and nonnegative")
        tracks[uid] = np.array(
            [
                np.mean(values[[m["stage"] == stage for m in metadata]])
                for stage in stages
            ]
        )
    summary = []
    grouped = defaultdict(list)
    within_group = []
    for row in candidate_rows:
        values = tracks.get(row["uid"])
        tau = ""
        if values is not None and np.max(values) > 0:
            tau = float(np.sum(1 - values / np.max(values)) / (len(values) - 1))
        summary.append(
            {
                "uid": row["uid"],
                "gene_id": row["gene_id"],
                "group": row["group"],
                "class": row["coverage_class"],
                "status": "matched" if values is not None else "missing",
                "tau": tau,
                **{
                    stage: float(values[i]) if values is not None else ""
                    for i, stage in enumerate(stages)
                },
            }
        )
        if values is not None and row["group"]:
            grouped[row["group"]].append(values)
    group_tracks = []
    variable_groups = {}
    for group, vectors in sorted(grouped.items()):
        mean = np.mean(vectors, axis=0)
        group_tracks.append(
            {
                "group": group,
                "available_copies": len(vectors),
                **dict(zip(stages, map(float, mean))),
            }
        )
        variable = [
            np.log1p(vector)
            for vector in vectors
            if max(vector) > float(settings.get("min_max_expression", 0))
            and np.ptp(vector) > 0
        ]
        if len(variable) >= 2:
            variable_groups[group] = variable
            correlations = [
                float(np.corrcoef(a, b)[0, 1]) for a, b in combinations(variable, 2)
            ]
            within_group.append(
                {
                    "group": group,
                    "copies": len(variable),
                    "pairs": len(correlations),
                    "mean_pearson": float(np.mean(correlations)),
                }
            )
    write_tsv(
        dest / "locus_stage_expression.tsv",
        summary,
        ["uid", "gene_id", "group", "class", "status", "tau"] + stages,
    )
    write_tsv(
        dest / "group_stage_expression.tsv",
        group_tracks,
        ["group", "available_copies"] + stages,
    )
    write_tsv(
        dest / "within_group_correlations.tsv",
        within_group,
        ["group", "copies", "pairs", "mean_pearson"],
    )
    permutations = int(settings.get("permutations", 999))
    if permutations < 1:
        raise InputError("expression.permutations must be positive")
    result = {
        "dataset": dataset,
        "matched_candidates": len(tracks),
        "total_candidates": len(candidate_rows),
        "stages": stages,
        "multicopy_groups": len(within_group),
        "observed_mean_r": None,
        "permutation_p": None,
        "permutations": permutations,
        "permutation_scope": "eligible copies in evaluable multicopy groups; group sizes fixed",
    }
    if variable_groups:
        observed = float(np.mean([row["mean_pearson"] for row in within_group]))
        vectors = np.array([v for group in variable_groups.values() for v in group])
        correlation = np.corrcoef(vectors)
        sizes = [len(group) for group in variable_groups.values()]
        rng = np.random.default_rng(cfg["seed"])
        null = []
        for _ in range(permutations):
            order = rng.permutation(len(vectors))
            offset = 0
            means = []
            for size in sizes:
                block = order[offset : offset + size]
                offset += size
                means.append(
                    float(
                        np.mean(
                            correlation[np.ix_(block, block)][np.triu_indices(size, 1)]
                        )
                    )
                )
            null.append(float(np.mean(means)))
        result.update(
            observed_mean_r=observed,
            permutation_p=(1 + sum(value >= observed for value in null))
            / (permutations + 1),
        )
        write_tsv(
            dest / "correlation_null.tsv",
            [{"permutation": i + 1, "mean_r": value} for i, value in enumerate(null)],
            ["permutation", "mean_r"],
        )
    save_json(dest / "summary.json", result)
    return result


def markers(cfg, dest):
    settings = cfg["markers"]
    records = read_tsv(
        resolve(cfg["_base"], settings["table"]),
        ["dataset", "gene_id", "assay", "stage", "cluster", "padj"],
    )
    candidate_rows = loci(cfg)
    lookup = {row["uid"]: row for row in candidate_rows}
    matches = []
    unknown = []
    significant = defaultdict(set)
    for number, record in enumerate(records, 2):
        possible = [
            row
            for row in candidate_rows
            if row["dataset"] == record["dataset"]
            and record["gene_id"] in [row["uid"], row["gene_id"], row["protein_id"]]
        ]
        if not possible:
            unknown.append(dict(record, source_row=number))
            continue
        uid = resolve_member(candidate_rows, record["dataset"], record["gene_id"])
        pvalue = float(record["padj"])
        if not math.isfinite(pvalue) or not 0 <= pvalue <= 1:
            raise InputError("Marker padj must be in [0, 1]")
        is_significant = pvalue <= float(settings.get("fdr", 0.05))
        matches.append(
            {
                "uid": uid,
                "dataset": record["dataset"],
                "gene_id": lookup[uid]["gene_id"],
                "class": lookup[uid]["coverage_class"],
                "assay": record["assay"],
                "stage": record["stage"],
                "cluster": record["cluster"],
                "padj": pvalue,
                "significant": is_significant,
                "source_row": number,
            }
        )
        if is_significant:
            significant[(uid, record["assay"])].add(
                (record["stage"], record["cluster"])
            )
    write_tsv(
        dest / "matched_records.tsv",
        matches,
        [
            "uid",
            "dataset",
            "gene_id",
            "class",
            "assay",
            "stage",
            "cluster",
            "padj",
            "significant",
            "source_row",
        ],
    )
    write_tsv(
        dest / "breadth.tsv",
        [
            {"uid": uid, "assay": assay, "stage_cluster_count": len(values)}
            for (uid, assay), values in sorted(significant.items())
        ],
        ["uid", "assay", "stage_cluster_count"],
    )
    write_tsv(
        dest / "unmatched_records.tsv",
        unknown,
        list(records[0]) + ["source_row"] if records else ["source_row"],
    )
    return {
        "matched_records": len(matches),
        "unmatched_records": len(unknown),
        "marker_genes": len({uid for uid, _ in significant}),
    }
