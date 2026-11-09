"""Whole-proteome OrthoFinder and MCScanX preparation for genuine pSONIC runs."""

from collections import defaultdict
import shutil

from .core import InputError, command, open_text, tool
from .family import models, candidates, all_sequences
from .inputs import manifest, write_fasta


def synteny(cfg, dest):
    samples = manifest(cfg)
    if len(samples) > 676:
        raise InputError("MCScanX adapter supports at most 676 datasets")
    by_gene = defaultdict(list)
    for row in models(cfg):
        by_gene[(row["dataset"], row["gene_id"])].append(row)
    family_models = {
        (row["dataset"], row["gene_id"]): row["uid"] for row in candidates(cfg)
    }
    selected = []
    for key, members in sorted(by_gene.items()):
        chosen = family_models.get(key)
        selected.append(
            next(row for row in members if row["uid"] == chosen)
            if chosen
            else min(members, key=lambda row: (-int(row["protein_length"]), row["uid"]))
        )
    missing = [row for row in selected if not row["chrom"] or not row["start"]]
    if missing:
        raise InputError(
            f"Whole-genome synteny requires coordinates for every representative; {len(missing)} models lack coordinates"
        )
    for sample in samples:
        relative = sample.get("relative_ploidy", "")
        if not relative.isdigit() or int(relative) < 1:
            raise InputError(
                "Automatic pSONIC preparation requires explicit positive relative_ploidy for each manifest dataset"
            )
    sequences = all_sequences(cfg)
    proteins = dest / "proteomes"
    proteins.mkdir()
    for sample in samples:
        write_fasta(
            proteins / (sample["dataset"] + ".faa"),
            [
                (row["uid"], sequences[row["uid"]])
                for row in selected
                if row["dataset"] == sample["dataset"]
            ],
        )
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
            "diamond",
        ],
        dest,
        dest / "commands.log",
    )
    results = list(proteins.glob("OrthoFinder/Results_*"))
    if len(results) != 1:
        raise InputError("Expected one whole-proteome OrthoFinder result directory")
    result = results[0]
    for name, relative in [
        ("Orthogroups.tsv", "Orthogroups/Orthogroups.tsv"),
        ("SequenceIDs.txt", "WorkingDirectory/SequenceIDs.txt"),
        ("SpeciesIDs.txt", "WorkingDirectory/SpeciesIDs.txt"),
    ]:
        shutil.copy2(result / relative, dest / name)
    coded = {}
    for line in (dest / "SequenceIDs.txt").read_text().splitlines():
        code, uid = line.split(": ", 1)
        coded[uid.split()[0]] = code
    prefixes = {
        sample["dataset"]: chr(65 + i // 26) + chr(65 + i % 26)
        for i, sample in enumerate(samples)
    }
    ordered = sorted(
        selected,
        key=lambda row: (
            prefixes[row["dataset"]],
            row["chrom"],
            int(row["start"]),
            row["uid"],
        ),
    )
    with (dest / "genome.gff").open("w") as handle:
        for row in ordered:
            handle.write(
                "\t".join(
                    [
                        prefixes[row["dataset"]] + row["chrom"],
                        coded[row["uid"]],
                        row["start"],
                        row["end"],
                    ]
                )
                + "\n"
            )
    blast_files = sorted((result / "WorkingDirectory").glob("Blast*.txt*"))
    if not blast_files:
        raise InputError("OrthoFinder BLAST-format intermediate files not found")
    with (dest / "genome.blast").open("w") as output:
        for path in blast_files:
            with open_text(path) as handle:
                shutil.copyfileobj(handle, output)
    command([tool("MCScanX"), "genome"], dest, dest / "commands.log")
    tandem_source = "MCScanX_output"
    if not (dest / "genome.tandem").exists():
        # Some MCScanX distributions do not export tandem edges. Derive a clearly
        # documented adjacent-gene BLAST-supported edge list for the adapter.
        positions = {}
        chromosome_genes = defaultdict(list)
        for row in ordered:
            chromosome = prefixes[row["dataset"]] + row["chrom"]
            chromosome_genes[chromosome].append(coded[row["uid"]])
        for chromosome, genes in chromosome_genes.items():
            for i, gene in enumerate(genes):
                positions[gene] = (chromosome, i)
        edges = set()
        with (dest / "genome.blast").open() as handle:
            for line in handle:
                values = line.split()
                if len(values) < 12:
                    raise InputError(
                        "Expected BLAST outfmt 6 for tandem-edge preparation"
                    )
                a, b = values[:2]
                if (
                    a == b
                    or a not in positions
                    or b not in positions
                    or float(values[10]) > 1e-5
                ):
                    continue
                ca, ia = positions[a]
                cb, ib = positions[b]
                if ca == cb and abs(ia - ib) == 1:
                    edges.add(tuple(sorted((a, b))))
        (dest / "genome.tandem").write_text(
            "".join(a + "," + b + "\n" for a, b in sorted(edges))
        )
        tandem_source = "adapter_adjacent_gene_BLAST_edges_E1e-5"
    for suffix in ["collinearity", "tandem"]:
        if not (dest / ("genome." + suffix)).is_file():
            raise InputError("MCScanX did not generate genome." + suffix)
    with (dest / "Orthogroups.tsv").open() as handle:
        order = handle.readline().rstrip("\n").split("\t")[1:]
    sample_map = {sample["dataset"]: sample for sample in samples}
    (dest / "relative_ploidy.tsv").write_text(
        "\t".join(sample_map[dataset]["relative_ploidy"] for dataset in order) + "\n"
    )
    return {
        "whole_genome_representatives": len(selected),
        "datasets": len(samples),
        "family_representatives": "domain-supported selected models",
        "background_representatives": "longest protein per declared locus",
        "mcscanx": "intra- and inter-species comparisons",
        "relative_ploidy_source": "explicit manifest column",
        "tandem_source": tandem_source,
    }
