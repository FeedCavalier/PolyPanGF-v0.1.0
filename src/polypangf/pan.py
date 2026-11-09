"""Occupancy and copy-number statistics with separate representation metadata."""

from collections import Counter, defaultdict
from pathlib import Path
import math
import numpy as np

from .core import InputError, read_tsv, save_json, write_tsv
from .family import candidates
from .inputs import manifest

CLASSES = ["Core", "Softcore", "Shell", "Specific", "Unassigned"]


def classify(present, total, softcore=0.9):
    if present == total:
        return "Core"
    if present == 1:
        return "Specific"
    return "Softcore" if present / total >= softcore else "Shell"


def pan(cfg, dest):
    settings = cfg.get("pan", {})
    unit = settings.get("unit", "species")
    threshold = float(settings.get("softcore_fraction", 0.9))
    if unit not in ["species", "dataset"] or not 0 < threshold <= 1:
        raise InputError(
            "pan.unit must be species or dataset, with 0 < softcore_fraction <= 1"
        )
    samples = [row for row in manifest(cfg) if row["role"] == "analysis"]
    dataset_map = {row["dataset"]: row for row in samples}
    units = list(dict.fromkeys(row[unit] for row in samples))
    if len(units) < 2:
        raise InputError(f"At least two independent {unit} units are required")
    rows = [row for row in candidates(cfg) if row["role"] == "analysis"]
    assignment = {
        row["uid"]: row["group"]
        for row in read_tsv(Path(cfg["output"]) / "orthology/membership.tsv")
    }
    counts = defaultdict(Counter)
    for row in rows:
        if assignment[row["uid"]]:
            counts[assignment[row["uid"]]][row["dataset"]] += 1
    summaries = []
    presence = []
    copies = []
    normalized = []
    group_classes = {}
    for group, copy_counts in sorted(counts.items()):
        occupied = {dataset_map[dataset][unit] for dataset in copy_counts}
        category = classify(len(occupied), len(units), threshold)
        group_classes[group] = category
        summaries.append(
            {
                "group": group,
                "class": category,
                "present_units": len(occupied),
                "total_units": len(units),
                "occupancy": len(occupied) / len(units),
                "candidate_loci": sum(copy_counts.values()),
            }
        )
        presence.append(
            dict(group=group, **{name: int(name in occupied) for name in units})
        )
        copies.append(
            dict(
                group=group,
                **{
                    sample["dataset"]: copy_counts[sample["dataset"]]
                    for sample in samples
                },
            )
        )
        for sample in samples:
            count = copy_counts[sample["dataset"]]
            normalized.append(
                {
                    "group": group,
                    "dataset": sample["dataset"],
                    "species": sample["species"],
                    "observed_copies": count,
                    "ploidy": sample["ploidy"],
                    "represented_sets": sample["represented_sets"],
                    "copies_per_represented_set": count
                    / float(sample["represented_sets"]),
                }
            )
    locus_rows = [
        dict(
            row,
            group=assignment[row["uid"]],
            coverage_class=group_classes.get(assignment[row["uid"]], "Unassigned"),
        )
        for row in rows
    ]
    write_tsv(
        dest / "groups.tsv",
        summaries,
        [
            "group",
            "class",
            "present_units",
            "total_units",
            "occupancy",
            "candidate_loci",
        ],
    )
    write_tsv(dest / "presence.tsv", presence, ["group"] + units)
    write_tsv(
        dest / "copy_counts.tsv",
        copies,
        ["group"] + [row["dataset"] for row in samples],
    )
    write_tsv(
        dest / "copy_representation.tsv",
        normalized,
        [
            "group",
            "dataset",
            "species",
            "observed_copies",
            "ploidy",
            "represented_sets",
            "copies_per_represented_set",
        ],
    )
    write_tsv(
        dest / "loci.tsv",
        locus_rows,
        list(rows[0]) + ["group", "coverage_class"]
        if rows
        else ["uid", "group", "coverage_class"],
    )
    # Exact expected richness under uniformly drawn subsets. No Monte Carlo subset bias.
    curve = []
    for size in range(1, len(units) + 1):
        denominator = math.comb(len(units), size)
        expected_pan = sum(
            1
            - (
                math.comb(len(units) - row["present_units"], size) / denominator
                if len(units) - row["present_units"] >= size
                else 0
            )
            for row in summaries
        )
        expected_core = sum(
            math.comb(row["present_units"], size) / denominator
            if row["present_units"] >= size
            else 0
            for row in summaries
        )
        curve.append(
            {
                "units": size,
                "expected_pan": expected_pan,
                "expected_core": expected_core,
            }
        )
    write_tsv(
        dest / "accumulation.tsv", curve, ["units", "expected_pan", "expected_core"]
    )
    result = {
        "unit": unit,
        "units": len(units),
        "groups": len(summaries),
        "assigned_loci": sum(bool(row["group"]) for row in locus_rows),
        "unassigned_loci": sum(not row["group"] for row in locus_rows),
        "classes": dict(Counter(row["class"] for row in summaries)),
        "softcore_fraction": threshold,
        "representation_normalization": "observed copies divided by represented sets, descriptive only",
    }
    save_json(dest / "summary.json", result)
    return result


def bh_adjust(pvalues):
    values = np.asarray(pvalues, dtype=float)
    if not len(values):
        return values
    order = np.argsort(values, kind="stable")
    adjusted = np.minimum.accumulate(
        (values[order] * len(values) / np.arange(1, len(values) + 1))[::-1]
    )[::-1]
    result = np.empty(len(values))
    result[order] = np.minimum(adjusted, 1)
    return result
