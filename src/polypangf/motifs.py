"""Optional native MEME discovery and FIMO scanning of all family candidates."""

from collections import defaultdict
from .core import InputError, command, tool, write_tsv
from .family import all_sequences, candidates
from .inputs import write_fasta


def motifs(cfg, dest):
    settings = cfg["motifs"]
    rows = candidates(cfg)
    sequences = all_sequences(cfg)
    unique = defaultdict(list)
    for row in rows:
        unique[sequences[row["uid"]]].append(row["uid"])
    if len(unique) < 2:
        raise InputError(
            "MEME discovery requires at least two unique family protein sequences"
        )
    write_fasta(
        dest / "unique.faa", [(uids[0], sequence) for sequence, uids in unique.items()]
    )
    write_fasta(
        dest / "all_candidates.faa",
        [(row["uid"], sequences[row["uid"]]) for row in rows],
    )
    write_tsv(
        dest / "sequence_members.tsv",
        [
            {"representative": uids[0], "uid": uid}
            for uids in unique.values()
            for uid in uids
        ],
        ["representative", "uid"],
    )
    number = int(settings.get("number", 10))
    minimum = int(settings.get("min_width", 6))
    maximum = int(settings.get("max_width", 50))
    if number < 1 or minimum < 2 or maximum < minimum:
        raise InputError("Invalid MEME motif count or widths")
    maximum_residues = sum(len(sequence) for sequence in unique)
    command(
        [
            tool("meme"),
            dest / "unique.faa",
            "-protein",
            "-oc",
            dest / "meme",
            "-mod",
            "zoops",
            "-nmotifs",
            str(number),
            "-minw",
            str(minimum),
            "-maxw",
            str(maximum),
            "-maxsize",
            str(maximum_residues + 1),
            "-seed",
            str(cfg["seed"]),
        ],
        dest,
        dest / "commands.log",
    )
    command(
        [
            tool("fimo"),
            "--oc",
            dest / "fimo",
            "--thresh",
            str(settings.get("pvalue", 1e-4)),
            dest / "meme/meme.xml",
            dest / "all_candidates.faa",
        ],
        dest,
        dest / "commands.log",
    )
    return {
        "unique_sequences": len(unique),
        "scanned_candidates": len(rows),
        "requested_motifs": number,
        "fimo_threshold_type": "P-value",
        "fimo_threshold": float(settings.get("pvalue", 1e-4)),
    }
