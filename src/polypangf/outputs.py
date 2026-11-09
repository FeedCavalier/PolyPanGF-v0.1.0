"""Standard figures, self-contained HTML reporting and portable result export."""

from collections import defaultdict
from pathlib import Path
import html
import json
import zipfile

import numpy as np

from . import __version__
from .core import InputError, read_tsv, save_json, sha256
from .pan import CLASSES

COLORS = {
    "Core": "#4B9C99",
    "Softcore": "#A591C0",
    "Shell": "#D6B063",
    "Specific": "#719BC4",
    "Unassigned": "#B5B5B5",
}


def plot(cfg, dest):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": 300,
        }
    )
    root = Path(cfg["output"])
    figures = []

    def save(fig, name, title):
        fig.tight_layout()
        for suffix in ["pdf", "svg", "png"]:
            fig.savefig(dest / (name + "." + suffix), dpi=300)
        plt.close(fig)
        figures.append({"file": name, "title": title})

    summary = json.loads((root / "pan/summary.json").read_text())
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    categories = CLASSES[:-1]
    ax.bar(
        categories,
        [summary["classes"].get(category, 0) for category in categories],
        color=[COLORS[c] for c in categories],
    )
    ax.set_ylabel("Orthogroups")
    save(fig, "coverage_classes", "Pan-family composition")
    curve = read_tsv(root / "pan/accumulation.tsv")
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    for key, label, color in [
        ("expected_pan", "Pan", COLORS["Core"]),
        ("expected_core", "Shared", COLORS["Shell"]),
    ]:
        ax.plot(
            [int(row["units"]) for row in curve],
            [float(row[key]) for row in curve],
            "o-",
            color=color,
            label=label,
        )
    ax.set_xlabel("Number of " + summary["unit"] + " units")
    ax.set_ylabel("Expected orthogroups")
    ax.legend(frameon=False)
    save(fig, "accumulation", "Exact expected accumulation")
    presence = read_tsv(root / "pan/presence.tsv")
    if presence:
        group_info = {row["group"]: row for row in read_tsv(root / "pan/groups.tsv")}
        presence.sort(
            key=lambda row: (
                CLASSES.index(group_info[row["group"]]["class"]),
                row["group"],
            )
        )
        columns = list(presence[0])[1:]
        matrix = np.array(
            [[float(row[column]) for row in presence] for column in columns]
        )
        fig, ax = plt.subplots(figsize=(7.2, max(3, len(columns) * 0.2)))
        ax.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap=ListedColormap(["#F3F0EB", "#4B9C99"]),
            vmin=0,
            vmax=1,
        )
        ax.set_yticks(range(len(columns)))
        ax.set_yticklabels(columns)
        ax.set_xticks([])
        ax.set_xlabel("Orthogroups ordered by coverage class")
        save(fig, "presence_matrix", "Family occupancy")
    copy_rows = read_tsv(root / "pan/copy_representation.tsv")
    if copy_rows:
        datasets = list(dict.fromkeys(row["dataset"] for row in copy_rows))
        groups = list(dict.fromkeys(row["group"] for row in copy_rows))
        lookup = {
            (row["dataset"], row["group"]): float(row["copies_per_represented_set"])
            for row in copy_rows
        }
        matrix = np.array(
            [[lookup[(dataset, group)] for group in groups] for dataset in datasets]
        )
        fig, ax = plt.subplots(figsize=(7.2, max(3, len(datasets) * 0.2)))
        image = ax.imshow(
            matrix, aspect="auto", interpolation="nearest", cmap="YlGnBu", vmin=0
        )
        ax.set_yticks(range(len(datasets)))
        ax.set_yticklabels(datasets)
        ax.set_xticks([])
        ax.set_xlabel("Orthogroups")
        fig.colorbar(image, ax=ax, label="Copies per represented set")
        save(
            fig,
            "copy_representation",
            "Representation-adjusted descriptive copy counts",
        )
    properties = root / "properties/protein_properties.tsv"
    if properties.exists():
        rows = read_tsv(properties)
        assignments = {
            row["uid"]: row["coverage_class"] for row in read_tsv(root / "pan/loci.tsv")
        }
        fig, axes = plt.subplots(1, 3, figsize=(7.2, 3))
        for ax, (key, label) in zip(
            axes,
            [
                ("length_aa", "Protein length (aa)"),
                ("pI", "Theoretical pI"),
                ("GRAVY", "GRAVY"),
            ],
        ):
            groups = {
                category: [
                    float(row[key])
                    for row in rows
                    if assignments.get(row["uid"]) == category and row[key] != ""
                ]
                for category in CLASSES
            }
            groups = {key: values for key, values in groups.items() if values}
            if groups:
                boxes = ax.boxplot(list(groups.values()), patch_artist=True)
                for patch, category in zip(boxes["boxes"], groups):
                    patch.set_facecolor(COLORS[category])
                ax.set_xticks(range(1, len(groups) + 1))
                ax.set_xticklabels(groups, rotation=55, ha="right")
            ax.set_ylabel(label)
        save(fig, "protein_properties", "Protein properties by occupancy class")
    expression_file = root / "expression/locus_stage_expression.tsv"
    if cfg.get("expression") is not None and expression_file.exists():
        rows = [row for row in read_tsv(expression_file) if row["status"] == "matched"]
        stages = json.loads((root / "expression/summary.json").read_text())["stages"]
        if rows:
            matrix = np.log2(
                1 + np.array([[float(row[stage]) for stage in stages] for row in rows])
            )
            fig, ax = plt.subplots(figsize=(6, max(3, min(12, len(rows) * 0.13))))
            image = ax.imshow(
                matrix, aspect="auto", interpolation="nearest", cmap="YlGnBu"
            )
            ax.set_xticks(range(len(stages)))
            ax.set_xticklabels(stages, rotation=45, ha="right")
            ax.set_yticks(range(len(rows)) if len(rows) <= 40 else [])
            ax.set_yticklabels(
                [row["gene_id"] for row in rows] if len(rows) <= 40 else []
            )
            fig.colorbar(image, ax=ax, label="log2(value + 1)")
            save(fig, "expression", "Stage-average expression")
    enrich_file = root / "enrichment/significant.tsv"
    if cfg.get("enrichment") is not None and enrich_file.exists():
        rows = read_tsv(enrich_file)
        if rows:
            fig, ax = plt.subplots(figsize=(7.2, max(3, len(rows) * 0.18)))
            image = ax.scatter(
                [int(row["overlap"]) / int(row["foreground"]) for row in rows],
                range(len(rows)),
                s=[12 + 5 * int(row["overlap"]) for row in rows],
                c=[-np.log10(max(float(row["qvalue"]), 1e-300)) for row in rows],
                cmap="viridis",
            )
            ax.set_yticks(range(len(rows)))
            ax.set_yticklabels(
                [row["category"] + " / " + row["description"] for row in rows],
                fontsize=6,
            )
            ax.set_xlabel("Annotated foreground fraction")
            fig.colorbar(image, ax=ax, label="-log10(FDR)")
            save(fig, "enrichment", "Significant functional enrichment")
    te_file = root / "te/profiles.tsv"
    if cfg.get("te") is not None and te_file.exists():
        rows = read_tsv(te_file)
        values = defaultdict(list)
        for row in rows:
            if row["coverage"] != "":
                values[(row["kind"], int(row["bin"]))].append(float(row["coverage"]))
        fig, ax = plt.subplots(figsize=(5.5, 3.2))
        for kind, color in [("family", "#4B9C99"), ("background", "#A591C0")]:
            xs = sorted(index for key, index in values if key == kind)
            ax.plot(
                xs,
                [np.mean(values[(kind, index)]) for index in xs],
                color=color,
                label=kind,
            )
        ax.set_xlabel("Upstream / gene body / downstream bins")
        ax.set_ylabel("TE-covered fraction")
        ax.legend(frameon=False)
        save(fig, "te_environment", "TE environment with matched background")
    marker_file = root / "markers/breadth.tsv"
    if cfg.get("markers") is not None and marker_file.exists():
        rows = read_tsv(marker_file)
        if rows:
            fig, ax = plt.subplots(figsize=(6, max(3, len(rows) * 0.14)))
            ax.barh(
                range(len(rows)),
                [int(row["stage_cluster_count"]) for row in rows],
                color="#4B9C99",
            )
            ax.set_yticks(range(len(rows)))
            ax.set_yticklabels(
                [row["uid"] + "/" + row["assay"] for row in rows], fontsize=5
            )
            ax.set_xlabel("Significant stage-cluster memberships")
            save(fig, "marker_breadth", "Author marker breadth")
    save_json(dest / "figures.json", figures)
    return {"figures": len(figures), "formats": ["PDF", "SVG", "300 dpi PNG"]}


def report(cfg, dest):
    root = Path(cfg["output"])
    title = str(cfg.get("family", "Gene family")) + " pan-family analysis"
    sections = [f"# {title}", "", f"PolyPanGF {__version__}", ""]
    blocks = [
        '<!doctype html><html lang="en"><meta charset="utf-8"><title>'
        + html.escape(title)
        + "</title>",
        "<style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:20px}img{max-width:100%}pre{white-space:pre-wrap}table{border-collapse:collapse}td,th{padding:6px;border:1px solid #ddd}</style>",
        "<h1>" + html.escape(title) + "</h1>",
    ]
    for stage in [
        "prepare",
        "identify",
        "synteny",
        "orthology",
        "pan",
        "properties",
        "tree",
        "motifs",
        "expression",
        "enrichment",
        "te",
        "kaks",
        "duplication",
        "markers",
    ]:
        if stage == "synteny" and (
            cfg.get("orthology", {}).get("backend") != "psonic"
            or cfg["orthology"].get("prefix")
        ):
            continue
        if (
            stage
            in [
                "tree",
                "motifs",
                "expression",
                "enrichment",
                "te",
                "kaks",
                "duplication",
                "markers",
            ]
            and cfg.get(stage) is None
        ):
            continue
        state_path = root / stage / "state.json"
        if state_path.exists():
            state = json.loads(state_path.read_text())
            summary = json.dumps(state["summary"], ensure_ascii=False, indent=2)
            sections += ["## " + stage, "", "```json", summary, "```", ""]
            blocks += ["<h2>" + stage + "</h2><pre>" + html.escape(summary) + "</pre>"]
    figures_path = root / "plot/figures.json"
    if figures_path.exists():
        for figure in json.loads(figures_path.read_text()):
            blocks.append(
                "<h2>"
                + html.escape(figure["title"])
                + '</h2><img src="../plot/'
                + figure["file"]
                + '.png" alt="'
                + html.escape(figure["title"])
                + '">'
            )
    notes = [
        "Species occupancy and dataset copy counts use separate statistical units.",
        "Copies per represented set are descriptive and depend on the declared assembly representation.",
        "Unassigned candidates are retained separately and are not automatically classified as Specific.",
        "Expression and marker modules use supplied processed values. Marker membership alone does not establish gene function.",
        "The pSONIC backend calls the official program using supplied whole-genome inputs. Imported groups retain their supplied provenance.",
    ]
    sections += ["## Interpretation", ""] + ["- " + note for note in notes]
    blocks.append(
        "<h2>Interpretation</h2><ul>"
        + "".join("<li>" + html.escape(note) + "</li>" for note in notes)
        + "</ul></html>"
    )
    (dest / "report.md").write_text("\n".join(sections) + "\n", encoding="utf-8")
    (dest / "report.html").write_text("\n".join(blocks), encoding="utf-8")
    save_json(
        dest / "configuration.json",
        {key: value for key, value in cfg.items() if not key.startswith("_")},
    )
    return {"report": "report.html", "markdown": "report.md", "network_assets": False}


def export_result(root, archive):
    root = Path(root).resolve()
    archive = Path(archive).resolve()
    if archive.exists():
        raise InputError(f"Export already exists: {archive}")
    files = [
        path
        for path in root.rglob("*")
        if path.is_file()
        and path != archive
        and not any(
            part.startswith(".") or part == "failures"
            for part in path.relative_to(root).parts
        )
    ]
    if not files:
        raise InputError("No completed results to export")
    archive.parent.mkdir(parents=True, exist_ok=True)
    temporary = archive.with_name(archive.name + ".partial")
    try:
        with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED) as output:
            manifest = []
            for path in files:
                name = path.relative_to(root).as_posix()
                output.write(path, name)
                manifest.append(sha256(path) + "  " + name)
            output.writestr("EXPORT_SHA256SUMS.txt", "\n".join(manifest) + "\n")
        with zipfile.ZipFile(temporary) as output:
            if output.testzip() is not None:
                raise InputError("Export ZIP failed CRC validation")
        temporary.rename(archive)
    finally:
        temporary.unlink(missing_ok=True)
    return {"archive": str(archive), "files": len(files), "sha256": sha256(archive)}
