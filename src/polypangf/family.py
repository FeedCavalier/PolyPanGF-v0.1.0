"""Family identification, sequence properties and external orthology adapters."""

from collections import defaultdict
from pathlib import Path
import shutil
import sys

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .core import InputError, command, read_tsv, resolve, tool, write_tsv
from .inputs import MODEL_FIELDS, fasta, manifest, write_fasta


def models(cfg):
    return read_tsv(
        Path(cfg["output"]) / "prepare/models.tsv", ["uid", "dataset", "gene_id"]
    )


def candidates(cfg):
    return read_tsv(
        Path(cfg["output"]) / "identify/candidates.tsv", ["uid", "dataset", "gene_id"]
    )


def all_sequences(cfg):
    result = {}
    for path in (Path(cfg["output"]) / "prepare/proteomes").glob("*.faa"):
        result.update(fasta(path, protein=True))
    return result


def identify(cfg, dest):
    settings = cfg.get("identification", {})
    model_rows = models(cfg)
    sequences = all_sequences(cfg)
    by_protein = {(row["dataset"], row["protein_id"]): row["uid"] for row in model_rows}
    hits = defaultdict(list)
    evidence = []
    supplied_architecture = {}
    if settings.get("candidates"):
        for row in read_tsv(
            resolve(cfg["_base"], settings["candidates"]), ["dataset", "protein_id"]
        ):
            key = (row["dataset"], row["protein_id"])
            if key not in by_protein:
                raise InputError(f"Unknown supplied candidate: {key}")
            uid = by_protein[key]
            if uid in supplied_architecture:
                raise InputError(f"Duplicate supplied candidate: {key}")
            supplied_architecture[uid] = row.get("architecture") or "unclassified"
            hits[uid].append(
                {
                    "profile": "supplied",
                    "coverage": 0.0,
                    "i_evalue": 0.0,
                    "score": 0.0,
                    "primary": True,
                }
            )
    else:
        profiles = settings.get("profiles", [])
        if not profiles or not any(
            profile.get("primary", i == 0) for i, profile in enumerate(profiles)
        ):
            raise InputError(
                "Provide identification.candidates or at least one primary HMM profile."
            )
        evalue = float(settings.get("domain_evalue", 1e-5))
        coverage = float(settings.get("min_domain_coverage", 0.0))
        if evalue <= 0 or not 0 <= coverage <= 1:
            raise InputError("Invalid HMM domain E-value or coverage threshold")
        merged = dest / "all_models.faa"
        write_fasta(merged, sequences.items())
        labels = set()
        for i, profile in enumerate(profiles):
            label = profile["label"]
            if label in labels:
                raise InputError("HMM profile labels must be unique")
            labels.add(label)
            table = dest / f"profile_{i}.domtbl"
            command(
                [
                    tool("hmmsearch"),
                    "--cpu",
                    str(cfg["threads"]),
                    "--noali",
                    "--domtblout",
                    table,
                    "-E",
                    str(evalue),
                    "--domE",
                    str(evalue),
                    resolve(cfg["_base"], profile["path"]),
                    merged,
                ],
                dest,
                dest / "commands.log",
            )
            with table.open() as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    values = line.split()
                    cov = (int(values[16]) - int(values[15]) + 1) / int(values[5])
                    if float(values[12]) > evalue or cov < coverage:
                        continue
                    hit = {
                        "profile": label,
                        "coverage": cov,
                        "i_evalue": float(values[12]),
                        "score": float(values[13]),
                        "primary": profile.get("primary", i == 0),
                    }
                    hits[values[0]].append(hit)
                    evidence.append(
                        dict(
                            uid=values[0],
                            **hit,
                            protein_start=int(values[17]),
                            protein_end=int(values[18]),
                        )
                    )
    loci = defaultdict(list)
    for row in model_rows:
        if any(hit["primary"] for hit in hits[row["uid"]]):
            loci[(row["dataset"], row["gene_id"])].append(row)
    selected = []
    isoforms = []
    for key, members in sorted(loci.items()):

        def rank(row):
            primary = [hit for hit in hits[row["uid"]] if hit["primary"]]
            return (
                -max(hit["score"] for hit in primary),
                -max(hit["coverage"] for hit in primary),
                -int(row["protein_length"]),
                row["protein_id"],
            )

        representative = min(members, key=rank)
        uid = representative["uid"]
        labels = {hit["profile"] for hit in hits[uid]}
        architecture = supplied_architecture.get(uid, "unclassified")
        for rule in settings.get("architecture_rules", []):
            if set(rule["requires"]).issubset(labels):
                architecture = rule["label"]
                break
        primary = [hit for hit in hits[uid] if hit["primary"]]
        selected.append(
            dict(
                representative,
                architecture=architecture,
                best_domain_evalue=min(hit["i_evalue"] for hit in primary),
                max_domain_coverage=max(hit["coverage"] for hit in primary),
                supporting_models=len(members),
            )
        )
        isoforms.extend(
            {
                "dataset": key[0],
                "gene_id": key[1],
                "uid": row["uid"],
                "selected": row["uid"] == uid,
            }
            for row in members
        )
    fields = MODEL_FIELDS + [
        "architecture",
        "best_domain_evalue",
        "max_domain_coverage",
        "supporting_models",
    ]
    write_tsv(dest / "candidates.tsv", selected, fields)
    write_tsv(
        dest / "supporting_isoforms.tsv",
        isoforms,
        ["dataset", "gene_id", "uid", "selected"],
    )
    write_tsv(
        dest / "domains.tsv",
        evidence,
        [
            "uid",
            "profile",
            "coverage",
            "i_evalue",
            "score",
            "primary",
            "protein_start",
            "protein_end",
        ],
    )
    family_dir = dest / "family_proteomes"
    family_dir.mkdir()
    for sample in manifest(cfg):
        write_fasta(
            family_dir / (sample["dataset"] + ".faa"),
            [
                (row["uid"], sequences[row["uid"]])
                for row in selected
                if row["dataset"] == sample["dataset"]
            ],
        )
    write_fasta(
        dest / "family.faa", [(row["uid"], sequences[row["uid"]]) for row in selected]
    )
    return {"candidate_loci": len(selected), "supporting_models": len(isoforms)}


def resolve_member(rows, dataset, identifier):
    matches = {
        row["uid"]
        for row in rows
        if row["dataset"] == dataset
        and identifier in [row["uid"], row["gene_id"], row["protein_id"]]
    }
    if len(matches) != 1:
        raise InputError(
            f"Orthology member must resolve to one candidate: {dataset}/{identifier}; matches={len(matches)}"
        )
    return matches.pop()


def parse_groups(path, rows, allow_nonfamily=False, permitted_datasets=None):
    permitted_datasets = (
        set(permitted_datasets)
        if permitted_datasets is not None
        else {row["dataset"] for row in rows}
    )
    lookup = {}
    for row in rows:
        for key in {row["uid"], row["gene_id"], row["protein_id"]}:
            lookup.setdefault((row["dataset"], key), set()).add(row["uid"])
    raw = read_tsv(path)
    groups = {}
    membership = {}
    nonfamily = 0
    for record in raw:
        if {"group", "dataset", "gene_id"}.issubset(record):
            group = record["group"]
            cells = [(record["dataset"], [record["gene_id"]])]
        else:
            first = next(iter(record))
            group = record[first]
            cells = [
                (sample, [gene.strip() for gene in value.split(",") if gene.strip()])
                for sample, value in record.items()
                if sample != first
            ]
        if not group:
            raise InputError(f"{path}: empty group identifier")
        for dataset, genes in cells:
            if dataset not in permitted_datasets:
                raise InputError(f"{path}: undeclared dataset column {dataset}")
            for gene in genes:
                matches = lookup.get((dataset, gene), set())
                if not matches and allow_nonfamily:
                    nonfamily += 1
                    continue
                if len(matches) != 1:
                    raise InputError(
                        f"{path}: unknown or ambiguous member {dataset}/{gene}"
                    )
                uid = next(iter(matches))
                if uid in membership:
                    raise InputError(
                        f"{path}: candidate repeated in orthology input: {uid}"
                    )
                membership[uid] = group
                groups.setdefault(group, []).append(uid)
    return groups, membership, nonfamily


def orthology(cfg, dest):
    settings = dict(cfg.get("orthology", {}))
    backend = settings.get("backend", "import")
    if backend == "psonic" and not settings.get("prefix"):
        prepared = Path(cfg["output"]) / "synteny"
        settings.update(
            prefix=str(prepared / "genome"),
            groups=str(prepared / "Orthogroups.tsv"),
            ploidy=str(prepared / "relative_ploidy.tsv"),
            sequence_ids=str(prepared / "SequenceIDs.txt"),
            species_ids=str(prepared / "SpeciesIDs.txt"),
        )
    rows = candidates(cfg)
    allow_nonfamily = settings.get("allow_nonfamily", False)
    if backend == "import":
        path = resolve(cfg["_base"], settings["groups"])
    elif backend == "orthofinder":
        source = Path(cfg["output"]) / "identify/family_proteomes"
        proteins = dest / "proteomes"
        proteins.mkdir()
        for path in source.glob("*.faa"):
            if path.stat().st_size:
                shutil.copy2(path, proteins / path.name)
        if len(list(proteins.glob("*.faa"))) < 2:
            raise InputError("OrthoFinder requires candidates in at least two datasets")
        command(
            [
                tool("orthofinder"),
                "-f",
                proteins,
                "-t",
                str(cfg["threads"]),
                "-a",
                str(cfg["threads"]),
                "-og",
                "-S",
                settings.get("search", "diamond"),
            ],
            dest,
            dest / "commands.log",
        )
        files = list(proteins.glob("OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv"))
        if len(files) != 1:
            raise InputError("Expected one OrthoFinder Orthogroups.tsv output")
        path = files[0]
    elif backend == "psonic":
        # These are genuine whole-genome MCScanX/OrthoFinder inputs, supplied by the user.
        prefix = dest / "genome"
        for suffix in ["collinearity", "gff", "tandem"]:
            shutil.copy2(
                resolve(cfg["_base"], settings["prefix"] + "." + suffix),
                str(prefix) + "." + suffix,
            )
        script = resolve(cfg["_base"], settings["script"])
        command(
            [
                settings.get("python", sys.executable),
                script,
                prefix,
                "-og",
                resolve(cfg["_base"], settings["groups"]),
                "-p",
                resolve(cfg["_base"], settings["ploidy"]),
                "-sID",
                resolve(cfg["_base"], settings["sequence_ids"]),
                "-specID",
                resolve(cfg["_base"], settings["species_ids"]),
                "-t",
                str(cfg["threads"]),
            ],
            dest,
            dest / "commands.log",
        )
        raw_path = dest / "pSONIC.txt"
        if not raw_path.is_file():
            raise InputError("Expected official pSONIC.txt; inspect commands.log")
        with resolve(cfg["_base"], settings["groups"]).open() as handle:
            header = handle.readline().rstrip("\n").split("\t")
        # Official pSONIC output is headerless, in the OrthoFinder species-column order.
        path = dest / "psonic_with_header.tsv"
        path.write_text("\t".join(header) + "\n" + raw_path.read_text())
        allow_nonfamily = True
    else:
        raise InputError(f"Unknown orthology backend: {backend}")
    groups, membership, nonfamily = parse_groups(
        path,
        rows,
        allow_nonfamily=allow_nonfamily,
        permitted_datasets=[sample["dataset"] for sample in manifest(cfg)],
    )
    write_tsv(
        dest / "membership.tsv",
        [
            {
                "group": membership.get(row["uid"], ""),
                "uid": row["uid"],
                "dataset": row["dataset"],
                "gene_id": row["gene_id"],
                "status": "assigned" if row["uid"] in membership else "unassigned",
            }
            for row in rows
        ],
        ["group", "uid", "dataset", "gene_id", "status"],
    )
    shutil.copy2(path, dest / "source_groups.tsv")
    return {
        "backend": backend,
        "groups": len(groups),
        "assigned": len(membership),
        "unassigned": len(rows) - len(membership),
        "nonfamily_members_filtered": nonfamily,
        "scope": "whole_genome_syntenic"
        if backend == "psonic"
        else "supplied"
        if backend == "import"
        else "family_orthogroups",
    }


def properties(cfg, dest):
    sequences = all_sequences(cfg)
    rows = []
    for row in candidates(cfg):
        sequence = sequences[row["uid"]]
        valid = set(sequence).issubset(set("ACDEFGHIKLMNPQRSTVWY"))
        result = {
            "uid": row["uid"],
            "dataset": row["dataset"],
            "gene_id": row["gene_id"],
            "length_aa": len(sequence),
            "molecular_weight": "",
            "pI": "",
            "GRAVY": "",
            "status": "ambiguous_residues",
        }
        if valid:
            protein = ProteinAnalysis(sequence)
            result.update(
                molecular_weight=protein.molecular_weight(),
                pI=protein.isoelectric_point(),
                GRAVY=protein.gravy(),
                status="calculated",
            )
        rows.append(result)
    write_tsv(
        dest / "protein_properties.tsv",
        rows,
        [
            "uid",
            "dataset",
            "gene_id",
            "length_aa",
            "molecular_weight",
            "pI",
            "GRAVY",
            "status",
        ],
    )
    return {
        "candidates": len(rows),
        "calculated": sum(row["status"] == "calculated" for row in rows),
    }


def tree(cfg, dest):
    source = Path(cfg["output"]) / "identify/family.faa"
    if len(candidates(cfg)) < 3:
        raise InputError("At least three family candidates are required for a tree")
    with (dest / "alignment.faa").open("w") as output:
        command(
            [tool("mafft"), "--auto", "--thread", str(cfg["threads"]), source],
            dest,
            dest / "commands.log",
            stdout=output,
        )
    executable = shutil.which("FastTree") or shutil.which("fasttree")
    if not executable:
        raise InputError("FastTree executable not found")
    with (dest / "tree.nwk").open("w") as output:
        command(
            [executable, "-lg", dest / "alignment.faa"],
            dest,
            dest / "commands.log",
            stdout=output,
        )
    return {"method": "MAFFT_auto_FastTree_LG", "rooting": "unrooted"}
