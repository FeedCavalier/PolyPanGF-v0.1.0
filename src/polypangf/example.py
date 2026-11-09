"""Small explicitly synthetic example for installation and command smoke tests."""

from pathlib import Path
import json

from .core import InputError, write_tsv
from .inputs import write_fasta


def initialize(directory):
    directory = Path(directory).resolve()
    if directory.exists() and any(directory.iterdir()):
        raise InputError("init requires a new or empty directory")
    directory.mkdir(parents=True, exist_ok=True)
    sequence = "MARGKKIGKRIENKSRQVTFSKRRNGLLKKAYELSVLCDAEVAVVFSPSGKLYEYANNNINRQVT"
    samples = []
    supplied = []
    grouping = []
    annotations = []
    coordinates = []
    expression = []
    for dataset, species, ploidy, represented, genes in [
        ("DIP", "Synthetic_species_A", 2, 1, ["a1", "a2", "a3"]),
        ("TET", "Synthetic_species_B", 4, 2, ["b1", "b2", "b3", "b4"]),
        ("OCT", "Synthetic_species_C", 8, 4, ["c1", "c2", "c3", "c4"]),
    ]:
        records = [(gene, sequence + "A" * i) for i, gene in enumerate(genes)]
        records += [
            (dataset + "_background" + str(i), "M" + "A" * (65 + i))
            for i in range(1, 7)
        ]
        write_fasta(directory / (dataset + ".faa"), records)
        model_map = []
        gff = []
        for i, (gene, _) in enumerate(records):
            start = 500 + i * 500
            model_map.append(
                {"protein_id": gene, "gene_id": gene, "transcript_id": gene + ".t"}
            )
            gff += [
                f"chr1\texample\tgene\t{start}\t{start + 249}\t.\t+\t.\tID={gene}",
                f"chr1\texample\tmRNA\t{start}\t{start + 249}\t.\t+\t.\tID={gene}.t;Parent={gene};protein_id={gene}",
            ]
            annotations.append(
                {
                    "gene_id": gene,
                    "term": "DEMO_BINDING" if gene in genes else "DEMO_BACKGROUND",
                    "description": "Synthetic annotation",
                }
            )
        write_tsv(
            directory / (dataset + ".models.tsv"),
            model_map,
            ["protein_id", "gene_id", "transcript_id"],
        )
        (directory / (dataset + ".gff3")).write_text(
            "##gff-version 3\n" + "\n".join(gff) + "\n"
        )
        samples.append(
            {
                "dataset": dataset,
                "species": species,
                "ploidy": ploidy,
                "represented_sets": represented,
                "proteins": dataset + ".faa",
                "gff": dataset + ".gff3",
                "model_map": dataset + ".models.tsv",
                "role": "analysis",
            }
        )
        for gene in genes:
            supplied.append(
                {"dataset": dataset, "protein_id": gene, "architecture": "synthetic"}
            )
        first_group = genes[:1] if dataset == "DIP" else genes[:2]
        grouping += [
            {"group": "DEMO_CORE", "dataset": dataset, "gene_id": gene}
            for gene in first_group
        ]
        if dataset != "OCT":
            grouping.append(
                {"group": "DEMO_SHELL", "dataset": dataset, "gene_id": genes[-2]}
            )
        grouping.append(
            {
                "group": "DEMO_SPECIFIC_" + dataset,
                "dataset": dataset,
                "gene_id": genes[-1],
            }
        )
        coordinates += [
            {"dataset": dataset, "chrom": "chr1", "start": 450, "end": 700},
            {"dataset": dataset, "chrom": "chr1", "start": 600, "end": 800},
        ]
    write_tsv(
        directory / "datasets.tsv",
        samples,
        [
            "dataset",
            "species",
            "ploidy",
            "represented_sets",
            "proteins",
            "gff",
            "model_map",
            "role",
        ],
    )
    write_tsv(
        directory / "candidates.tsv",
        supplied,
        ["dataset", "protein_id", "architecture"],
    )
    write_tsv(directory / "groups.tsv", grouping, ["group", "dataset", "gene_id"])
    write_tsv(directory / "terms.tsv", annotations, ["gene_id", "term", "description"])
    write_tsv(directory / "te.tsv", coordinates, ["dataset", "chrom", "start", "end"])
    for i, gene in enumerate(["c1", "c2", "c3", "c4"]):
        expression.append(dict(gene_id=gene, S1=1 + i, S2=4 + i, S3=2 + i, S4=8 + i))
    write_tsv(
        directory / "expression.tsv", expression, ["gene_id", "S1", "S2", "S3", "S4"]
    )
    write_tsv(
        directory / "samples.tsv",
        [{"sample": f"S{i}", "stage": str(i)} for i in range(1, 5)],
        ["sample", "stage"],
    )
    cfg = {
        "schema_version": 1,
        "family": "Synthetic demonstration",
        "manifest": "datasets.tsv",
        "output": "results",
        "threads": 2,
        "seed": 42,
        "identification": {"candidates": "candidates.tsv"},
        "orthology": {"backend": "import", "groups": "groups.tsv"},
        "pan": {"unit": "species", "softcore_fraction": 0.9},
        "expression": {
            "dataset": "OCT",
            "matrix": "expression.tsv",
            "samples": "samples.tsv",
            "permutations": 99,
        },
        "enrichment": {
            "dataset": "OCT",
            "term2gene": "terms.tsv",
            "database": "synthetic",
        },
        "te": {"intervals": "te.tsv", "flank": 200, "bins_per_region": 10},
    }
    (directory / "config.json").write_text(json.dumps(cfg, indent=2) + "\n")
    (directory / "README.txt").write_text(
        "Synthetic installation example. It does not provide biological validation.\nRun polypangf run -c config.json\n"
    )
    return {
        "project": str(directory),
        "config": str(directory / "config.json"),
        "synthetic": True,
    }
