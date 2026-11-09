"""Proteome and gene-model preparation with explicit locus provenance."""

from collections import defaultdict
from pathlib import Path
import re
import hashlib
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Seq import Seq

from .core import InputError, open_text, read_tsv, resolve, write_tsv

MANIFEST_FIELDS = ["dataset", "species", "ploidy", "represented_sets", "proteins"]
MODEL_FIELDS = [
    "uid",
    "dataset",
    "species",
    "role",
    "ploidy",
    "represented_sets",
    "protein_id",
    "gene_id",
    "locus_basis",
    "chrom",
    "start",
    "end",
    "strand",
    "cds_id",
    "protein_length",
    "protein_sha256",
]


def manifest(cfg):
    path = Path(cfg["manifest"])
    rows = read_tsv(path, MANIFEST_FIELDS)
    if not rows:
        raise InputError("Provide at least one dataset.")
    seen = set()
    for row in rows:
        dataset = row["dataset"]
        if not re.fullmatch("[A-Za-z][A-Za-z0-9_-]*", dataset) or dataset in seen:
            raise InputError(f"Duplicate or unsafe dataset name: {dataset}")
        seen.add(dataset)
        for key in ["ploidy", "represented_sets"]:
            try:
                value = float(row[key])
            except ValueError as exc:
                raise InputError(f"{dataset}: invalid {key}") from exc
            if value <= 0 or not value.is_integer():
                raise InputError(f"{dataset}: {key} must be a positive integer")
        row["role"] = row.get("role") or "analysis"
        if row["role"] not in ["analysis", "outgroup"] or not row["species"].strip():
            raise InputError(f"{dataset}: invalid role or species")
        for key in ["proteins", "cds", "gff", "model_map"]:
            if row.get(key):
                row[key] = str(resolve(path.parent, row[key]))
                if not Path(row[key]).is_file():
                    raise InputError(f"{dataset}: file not found: {row[key]}")
    if not any(row["role"] == "analysis" for row in rows):
        raise InputError("At least one analysis dataset is required.")
    return rows


def source_files(cfg):
    return [cfg["manifest"]] + [
        row[key]
        for row in manifest(cfg)
        for key in ["proteins", "cds", "gff", "model_map"]
        if row.get(key)
    ]


def fasta(path, protein=False, invalid_policy="error", audit=None):
    if invalid_policy not in ["error", "mask", "skip"]:
        raise InputError("invalid_protein must be error, mask or skip")
    records = {}
    seen = set()
    with open_text(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in seen:
                raise InputError(f"{path}: duplicate FASTA ID {record.id}")
            seen.add(record.id)
            sequence = str(record.seq).upper()
            if protein:
                sequence = sequence.rstrip("*")
                valid = set("ABCDEFGHIKLMNPQRSTVWXYZOUJ")
            else:
                valid = set("ACGTRYSWKMBDHVN")
            invalid = set(sequence) - valid
            if not sequence or invalid:
                if not protein or invalid_policy == "error" or not sequence:
                    raise InputError(f"{path}: invalid sequence for {record.id}")
                if audit is not None:
                    audit.append(
                        {
                            "protein_id": record.id,
                            "action": invalid_policy,
                            "invalid_residues": "".join(sorted(invalid)),
                            "source_sha256": hashlib.sha256(
                                sequence.encode()
                            ).hexdigest(),
                        }
                    )
                if invalid_policy == "skip":
                    continue
                sequence = "".join(
                    residue if residue in valid else "X" for residue in sequence
                )
            records[record.id] = sequence
    if not records:
        raise InputError(f"{path}: empty FASTA")
    return records


def write_fasta(path, records):
    with Path(path).open("w") as handle:
        for uid, sequence in records:
            handle.write(">" + uid + "\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def gff_models(path):
    transcripts = {}
    links = {}
    intervals = defaultdict(list)
    genes = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) != 9:
                raise InputError(f"{path}: expected 9 GFF3 columns")
            chrom, _, feature, start, end, _, strand, _, attributes = values
            start, end = int(start), int(end)
            if start < 1 or end < start:
                raise InputError(f"{path}: invalid GFF3 coordinates")
            attrs = {
                key: unquote(value)
                for key, value in (
                    part.split("=", 1) for part in attributes.split(";") if "=" in part
                )
            }
            uid = attrs.get("ID", "")
            parents = attrs.get("Parent", "").split(",")
            info = {"chrom": chrom, "start": start, "end": end, "strand": strand}
            if feature == "gene" and uid:
                if uid in genes and genes[uid] != info:
                    raise InputError(f"{path}: conflicting gene coordinates for {uid}")
                genes[uid] = info
            if feature in ["mRNA", "transcript"] and uid:
                transcript_info = dict(info, gene_id=parents[0] or uid)
                if uid in transcripts and transcripts[uid] != transcript_info:
                    raise InputError(f"{path}: conflicting transcript model {uid}")
                transcripts[uid] = transcript_info
                for key in ["ID", "Name", "protein_id", "transcript_id"]:
                    if attrs.get(key):
                        links.setdefault(attrs[key], set()).add(uid)
            if feature == "CDS":
                for parent in parents:
                    if parent:
                        intervals[parent].append((start, end))
                for key in ["protein_id", "transcript_id"]:
                    if attrs.get(key) and parents[0]:
                        links.setdefault(attrs[key], set()).add(parents[0])
    return transcripts, links, genes, intervals


def prepare(cfg, dest):
    models = []
    cds_checks = []
    sequence_audit = []
    proteome_dir = dest / "proteomes"
    proteome_dir.mkdir()
    for sample in manifest(cfg):
        audit = []
        proteins = fasta(
            sample["proteins"],
            protein=True,
            invalid_policy=cfg.get("preparation", {}).get("invalid_protein", "error"),
            audit=audit,
        )
        sequence_audit.extend(dict(row, dataset=sample["dataset"]) for row in audit)
        cds = fasta(sample["cds"]) if sample.get("cds") else {}
        overrides = {}
        skipped = {row["protein_id"] for row in audit if row["action"] == "skip"}
        if sample.get("model_map"):
            for row in read_tsv(sample["model_map"], ["protein_id", "gene_id"]):
                if row["protein_id"] in skipped:
                    continue
                if (
                    row["protein_id"] in overrides
                    or row["protein_id"] not in proteins
                    or not row["gene_id"]
                ):
                    raise InputError(
                        f"{sample['dataset']}: duplicate, unknown or empty model-map entry"
                    )
                if row["gene_id"] in ["UNMAPPED_PROTEIN_ONLY", "NA", ".", "unmapped"]:
                    raise InputError(
                        "Use a unique protein ID and unit_type=protein for an unmapped model, rather than a shared placeholder gene ID"
                    )
                if row.get("unit_type", "gene") not in ["gene", "protein"]:
                    raise InputError("model_map.unit_type must be gene or protein")
                overrides[row["protein_id"]] = row
        transcripts, links, genes, intervals = (
            gff_models(sample["gff"]) if sample.get("gff") else ({}, {}, {}, {})
        )
        normalized = []
        normalized_cds = []
        for number, (original, sequence) in enumerate(proteins.items(), 1):
            uid = sample["dataset"] + "__P" + str(number).zfill(8)
            override = overrides.get(original, {})
            aliases = links.get(original, {original})
            if not override.get("transcript_id") and len(aliases) > 1:
                raise InputError(
                    f"{sample['dataset']}: ambiguous GFF3 protein alias {original}; provide model_map"
                )
            transcript = override.get("transcript_id") or next(iter(aliases))
            info = transcripts.get(transcript, {})
            gene = override.get("gene_id") or info.get("gene_id") or original
            location = genes.get(gene, info)
            cid = override.get("cds_id") or transcript
            if cid not in cds and original in cds:
                cid = original
            status = "not_provided"
            if sample.get("cds"):
                status = "missing"
                if cid in cds:
                    dna = cds[cid]
                    status = "frame_error" if len(dna) % 3 else "translation_mismatch"
                    translated = str(Seq(dna).translate()) if len(dna) % 3 == 0 else ""
                    without_stop = (
                        translated[:-1] if translated.endswith("*") else translated
                    )
                    if len(dna) % 3 == 0 and without_stop == sequence:
                        status = "matched"
                        normalized_cds.append(
                            (uid, dna[:-3] if translated.endswith("*") else dna)
                        )
            cds_checks.append({"uid": uid, "cds_id": cid, "status": status})
            row = {
                key: sample[key]
                for key in ["dataset", "species", "role", "ploidy", "represented_sets"]
            }
            row.update(
                uid=uid,
                protein_id=original,
                gene_id=gene,
                locus_basis="protein_id_assumption"
                if override.get("unit_type") == "protein"
                else "model_map"
                if override
                else "GFF3"
                if info
                else "protein_id_assumption",
                cds_id=cid if status == "matched" else "",
                protein_length=len(sequence),
                protein_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
            )
            row.update(
                {
                    key: location.get(key, "")
                    for key in ["chrom", "start", "end", "strand"]
                }
            )
            models.append(row)
            normalized.append((uid, sequence))
        write_fasta(proteome_dir / (sample["dataset"] + ".faa"), normalized)
        write_fasta(dest / (sample["dataset"] + ".cds.fa"), normalized_cds)
    write_tsv(dest / "models.tsv", models, MODEL_FIELDS)
    write_tsv(dest / "cds_checks.tsv", cds_checks, ["uid", "cds_id", "status"])
    write_tsv(
        dest / "sequence_audit.tsv",
        sequence_audit,
        ["dataset", "protein_id", "action", "invalid_residues", "source_sha256"],
    )
    # Full gene order is retained for TE backgrounds and external synteny preparation.
    return {
        "datasets": len(manifest(cfg)),
        "protein_models": len(models),
        "cds_matched": sum(row["status"] == "matched" for row in cds_checks),
        "protein_id_assumptions": sum(
            row["locus_basis"] == "protein_id_assumption" for row in models
        ),
        "sequence_policy_actions": len(sequence_audit),
    }
