"""Coordinate-aware TE profiles, duplication imports and codon-based Ka/Ks."""

from collections import defaultdict
from itertools import combinations
from pathlib import Path
import math
from statistics import median

from Bio.Phylo.PAML import yn00

from .core import InputError, command, read_tsv, resolve, tool, write_tsv
from .family import candidates, all_sequences, models, resolve_member
from .inputs import fasta, write_fasta
from .analysis import loci


def merge_intervals(intervals):
    result = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise InputError(
                "TE coordinates require 0-based half-open intervals with start < end"
            )
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def covered(intervals, start, end):
    return sum(
        max(0, min(end, b) - max(start, a))
        for a, b in intervals
        if a < end and b > start
    )


def te(cfg, dest):
    settings = cfg["te"]
    rows = read_tsv(
        resolve(cfg["_base"], settings["intervals"]),
        ["dataset", "chrom", "start", "end"],
    )
    interval_map = defaultdict(list)
    for row in rows:
        interval_map[(row["dataset"], row["chrom"])].append(
            (int(row["start"]), int(row["end"]))
        )
    interval_map = {
        key: merge_intervals(values) for key, values in interval_map.items()
    }
    lengths = {}
    if settings.get("chrom_sizes"):
        for row in read_tsv(
            resolve(cfg["_base"], settings["chrom_sizes"]),
            ["dataset", "chrom", "length"],
        ):
            key = row["dataset"], row["chrom"]
            if key in lengths or int(row["length"]) < 1:
                raise InputError("Duplicate or invalid chromosome size")
            lengths[key] = int(row["length"])
        for key, intervals in interval_map.items():
            if key not in lengths or any(end > lengths[key] for _, end in intervals):
                raise InputError(
                    "TE interval exceeds supplied chromosome bounds or chromosome size is missing"
                )
    candidate_rows = candidates(cfg)
    candidate_keys = {(row["dataset"], row["gene_id"]) for row in candidate_rows}
    background = {}
    for row in models(cfg):
        key = row["dataset"], row["gene_id"]
        if key not in candidate_keys and row["start"]:
            background[key] = row
    flank = int(settings.get("flank", 2000))
    bins = int(settings.get("bins_per_region", 30))
    if flank < 1 or bins < 1:
        raise InputError("TE flank and bins_per_region must be positive")
    profiles = []
    unavailable = []
    for row in candidate_rows:
        key = row["dataset"], row["chrom"]
        if (
            not row["start"]
            or key not in interval_map
            or row["strand"] not in ["+", "-"]
        ):
            unavailable.append(
                {
                    "uid": row["uid"],
                    "reason": "missing_coordinates_strand_or_TE_chromosome",
                }
            )
            continue
        matched = [
            other
            for other in background.values()
            if other["dataset"] == row["dataset"] and other["chrom"] == row["chrom"]
        ]
        length = int(row["end"]) - int(row["start"]) + 1
        match = (
            min(
                matched,
                key=lambda other: (
                    abs((int(other["end"]) - int(other["start"]) + 1) - length),
                    other["gene_id"],
                ),
            )
            if matched
            else None
        )
        for kind, target in [("family", row), ("background", match)]:
            if target is None:
                continue
            start = int(target["start"]) - 1
            end = int(target["end"])
            if key in lengths and end > lengths[key]:
                raise InputError("Gene exceeds supplied chromosome bounds")
            limits = [(start - flank, start), (start, end), (end, end + flank)]
            values = []
            for lower, upper in limits:
                for index in range(bins):
                    left = max(0, lower + (upper - lower) * index // bins)
                    right = max(0, lower + (upper - lower) * (index + 1) // bins)
                    if key in lengths:
                        left, right = min(left, lengths[key]), min(right, lengths[key])
                    values.append(
                        covered(interval_map[key], left, right) / (right - left)
                        if right > left
                        else ""
                    )
            if target["strand"] == "-":
                values.reverse()
            profiles.extend(
                {
                    "anchor_uid": row["uid"],
                    "kind": kind,
                    "dataset": target["dataset"],
                    "gene_id": target["gene_id"],
                    "bin": index + 1,
                    "coverage": value,
                }
                for index, value in enumerate(values)
            )
    write_tsv(
        dest / "profiles.tsv",
        profiles,
        ["anchor_uid", "kind", "dataset", "gene_id", "bin", "coverage"],
    )
    write_tsv(dest / "unavailable.tsv", unavailable, ["uid", "reason"])
    return {
        "candidate_profiles": len(
            {row["anchor_uid"] for row in profiles if row["kind"] == "family"}
        ),
        "unavailable_candidates": len(unavailable),
        "merged_TE_intervals": sum(map(len, interval_map.values())),
        "coordinate_system": "TE: BED 0-based half-open; genes: GFF3 1-based closed; orientation: transcription direction",
        "chromosome_end_boundary": "clipped to provided chromosome sizes"
        if lengths
        else "not clipped without chromosome size metadata",
    }


def duplication(cfg, dest):
    settings = cfg["duplication"]
    supplied = read_tsv(
        resolve(cfg["_base"], settings["table"]),
        ["dataset", "gene_id", "duplication_class"],
    )
    rows = candidates(cfg)
    known = {}
    allowed = {
        "wgd_segmental",
        "tandem",
        "proximal",
        "transposed",
        "dispersed",
        "unsupported",
    }
    for row in supplied:
        if row["duplication_class"] not in allowed:
            raise InputError(f"Unknown duplication_class: {row['duplication_class']}")
        possible = [
            r
            for r in rows
            if r["dataset"] == row["dataset"]
            and row["gene_id"] in [r["uid"], r["gene_id"], r["protein_id"]]
        ]
        if not possible:
            continue
        uid = resolve_member(rows, row["dataset"], row["gene_id"])
        if uid in known:
            raise InputError(
                f"Duplicate or conflicting duplication assignment for {uid}"
            )
        known[uid] = row["duplication_class"]
    output = [
        {
            "uid": row["uid"],
            "dataset": row["dataset"],
            "gene_id": row["gene_id"],
            "duplication_class": known.get(row["uid"], "not_evaluated"),
        }
        for row in rows
    ]
    write_tsv(
        dest / "duplication.tsv",
        output,
        ["uid", "dataset", "gene_id", "duplication_class"],
    )
    return {
        "imported_assignments": len(known),
        "source_method": settings.get("source_method", "user_supplied"),
        "candidates": len(output),
    }


def kaks(cfg, dest):
    settings = cfg["kaks"]
    sequences = all_sequences(cfg)
    coding = {}
    for path in (Path(cfg["output"]) / "prepare").glob("*.cds.fa"):
        if path.stat().st_size:
            coding.update(fasta(path))
    groups = defaultdict(lambda: defaultdict(list))
    for row in loci(cfg):
        if row["group"]:
            groups[row["group"]][row["dataset"]].append(row["uid"])
    results = []
    excluded = []
    for group, samples in sorted(groups.items()):
        for a, b in combinations(sorted(samples), 2):
            if settings.get("anchor") and settings["anchor"] not in [a, b]:
                continue
            if len(samples[a]) != 1 or len(samples[b]) != 1:
                excluded.append(
                    {
                        "group": group,
                        "dataset_a": a,
                        "dataset_b": b,
                        "reason": "not_one_to_one",
                    }
                )
                continue
            uid_a, uid_b = samples[a][0], samples[b][0]
            result = {
                "group": group,
                "uid_a": uid_a,
                "uid_b": uid_b,
                "status": "missing_validated_CDS",
                "codons": "",
                "coverage": "",
                "Ka": "",
                "Ks": "",
                "omega": "",
            }
            results.append(result)
            if uid_a not in coding or uid_b not in coding:
                continue
            if coding[uid_a] == coding[uid_b]:
                result["status"] = "identical_CDS"
                continue
            work = dest / ("pair_" + str(len(results)).zfill(6))
            work.mkdir()
            write_fasta(
                work / "proteins.fa",
                [("seq1", sequences[uid_a]), ("seq2", sequences[uid_b])],
            )
            write_fasta(
                work / "cds.fa", [("seq1", coding[uid_a]), ("seq2", coding[uid_b])]
            )
            timeout = int(settings.get("timeout", 300))
            with (work / "alignment.fa").open("w") as output:
                command(
                    [
                        tool("mafft"),
                        "--auto",
                        "--quiet",
                        "--thread",
                        "1",
                        "proteins.fa",
                    ],
                    work,
                    work / "commands.log",
                    stdout=output,
                    timeout=timeout,
                )
            with (work / "codons.fa").open("w") as output:
                command(
                    [
                        tool("pal2nal.pl"),
                        "alignment.fa",
                        "cds.fa",
                        "-output",
                        "fasta",
                        "-nogap",
                    ],
                    work,
                    work / "commands.log",
                    stdout=output,
                    timeout=timeout,
                )
            raw = fasta_alignment(work / "codons.fa")
            if (
                set(raw) != {"seq1", "seq2"}
                or len(raw["seq1"]) != len(raw["seq2"])
                or len(raw["seq1"]) % 3
            ):
                raise InputError("PAL2NAL returned an invalid codon alignment")
            starts = [
                i
                for i in range(0, len(raw["seq1"]), 3)
                if all(set(seq[i : i + 3]) <= set("ACGT") for seq in raw.values())
            ]
            clean = {
                key: "".join(seq[i : i + 3] for i in starts) for key, seq in raw.items()
            }
            codons = len(starts)
            coverage = codons / max(len(sequences[uid_a]), len(sequences[uid_b]))
            result.update(
                codons=codons, coverage=coverage, status="insufficient_alignment"
            )
            if codons < int(settings.get("min_codons", 50)) or coverage < float(
                settings.get("min_coverage", 0.7)
            ):
                continue
            (work / "codon.phy").write_text(
                f"2 {codons * 3}\nseq1\n{clean['seq1']}\nseq2\n{clean['seq2']}\n"
            )
            (work / "yn00.ctl").write_text(
                "seqfile = codon.phy\noutfile = yn00.out\nverbose = 0\nicode = 0\nweighting = 0\ncommonf3x4 = 0\n"
            )
            command(
                [tool("yn00"), "yn00.ctl"], work, work / "commands.log", timeout=timeout
            )
            estimate = yn00.read(str(work / "yn00.out"))["seq1"]["seq2"]["YN00"]
            ka, ks, omega = (
                float(estimate["dN"]),
                float(estimate["dS"]),
                float(estimate["omega"]),
            )
            if not all(
                math.isfinite(value) and value >= 0 for value in [ka, ks, omega]
            ):
                result["status"] = "invalid_estimate"
                continue
            result.update(
                Ka=ka,
                Ks=ks,
                omega=omega,
                status="pass"
                if float(settings.get("min_ks", 0.01))
                <= ks
                <= float(settings.get("max_ks", 2))
                and omega < float(settings.get("max_omega", 99))
                else "outside_filters",
            )
    write_tsv(
        dest / "pairs.tsv",
        results,
        [
            "group",
            "uid_a",
            "uid_b",
            "status",
            "codons",
            "coverage",
            "Ka",
            "Ks",
            "omega",
        ],
    )
    write_tsv(
        dest / "excluded.tsv", excluded, ["group", "dataset_a", "dataset_b", "reason"]
    )
    passing = defaultdict(list)
    classes = {row["group"]: row["coverage_class"] for row in loci(cfg) if row["group"]}
    for row in results:
        if row["status"] == "pass":
            passing[row["group"]].append(row["omega"])
    write_tsv(
        dest / "group_medians.tsv",
        [
            {
                "group": group,
                "class": classes[group],
                "passing_pairs": len(values),
                "median_omega": median(values),
            }
            for group, values in sorted(passing.items())
        ],
        ["group", "class", "passing_pairs", "median_omega"],
    )
    return {
        "evaluated_pairs": len(results),
        "passing_pairs": sum(row["status"] == "pass" for row in results),
        "excluded_multicopy_comparisons": len(excluded),
    }


def fasta_alignment(path):
    from Bio import SeqIO

    return {record.id: str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}
