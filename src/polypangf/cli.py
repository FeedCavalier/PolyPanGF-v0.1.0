"""Public command-line interface and dependency-aware stage orchestration."""

import argparse
import importlib.util
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys

from . import __version__
from .core import (
    InputError,
    Runner,
    load_config,
    project_lock,
    read_tsv,
    resolve,
    sha256,
)
from .inputs import fasta, gff_models, manifest, source_files, prepare
from .family import identify, orthology, properties, tree
from .pan import pan
from .analysis import enrichment, expression, markers
from .evolution import te, kaks, duplication
from .synteny import synteny
from .motifs import motifs
from .outputs import plot, report, export_result
from .example import initialize

FUNCTIONS = {
    "prepare": prepare,
    "identify": identify,
    "synteny": synteny,
    "orthology": orthology,
    "pan": pan,
    "properties": properties,
    "tree": tree,
    "expression": expression,
    "enrichment": enrichment,
    "te": te,
    "kaks": kaks,
    "duplication": duplication,
    "markers": markers,
    "motifs": motifs,
    "plot": plot,
    "report": report,
}
OPTIONAL = [
    "tree",
    "motifs",
    "expression",
    "enrichment",
    "te",
    "kaks",
    "duplication",
    "markers",
]
PATHS = {
    "identification": ["candidates"],
    "orthology": ["groups", "script", "ploidy", "sequence_ids", "species_ids"],
    "expression": ["matrix", "samples"],
    "enrichment": ["term2gene", "background", "valid_terms"],
    "te": ["intervals", "chrom_sizes"],
    "duplication": ["table"],
    "markers": ["table"],
}
SECTION_KEYS = {
    "preparation": {"invalid_protein"},
    "identification": {
        "candidates",
        "profiles",
        "domain_evalue",
        "min_domain_coverage",
        "architecture_rules",
    },
    "orthology": {
        "backend",
        "groups",
        "script",
        "prefix",
        "ploidy",
        "sequence_ids",
        "species_ids",
        "python",
        "allow_nonfamily",
        "search",
    },
    "pan": {"unit", "softcore_fraction"},
    "tree": set(),
    "motifs": {"number", "min_width", "max_width", "pvalue"},
    "expression": {
        "dataset",
        "matrix",
        "samples",
        "permutations",
        "min_max_expression",
    },
    "enrichment": {
        "dataset",
        "database",
        "term2gene",
        "valid_terms",
        "background",
        "min_term_size",
        "fdr",
    },
    "te": {"intervals", "chrom_sizes", "flank", "bins_per_region"},
    "kaks": {
        "anchor",
        "min_codons",
        "min_coverage",
        "min_ks",
        "max_ks",
        "max_omega",
        "timeout",
    },
    "duplication": {"table", "source_method"},
    "markers": {"table", "fdr"},
}


def external_inputs(cfg, section):
    settings = cfg.get(section, {})
    result = [
        resolve(cfg["_base"], settings[key])
        for key in PATHS.get(section, [])
        if settings.get(key)
    ]
    if section == "identification":
        result += [
            resolve(cfg["_base"], profile["path"])
            for profile in settings.get("profiles", [])
        ]
    if (
        section == "orthology"
        and settings.get("backend") == "psonic"
        and settings.get("prefix")
    ):
        result += [
            resolve(cfg["_base"], settings["prefix"] + "." + suffix)
            for suffix in ["gff", "tandem", "collinearity"]
        ]
    return result


def enabled(cfg):
    automatic = cfg.get("orthology", {}).get("backend") == "psonic" and not cfg[
        "orthology"
    ].get("prefix")
    return (
        ["prepare", "identify"]
        + (["synteny"] if automatic else [])
        + ["orthology", "pan", "properties"]
        + [name for name in OPTIONAL if cfg.get(name) is not None]
        + ["plot", "report"]
    )


def dependencies(cfg, name):
    base = {
        "prepare": [],
        "identify": ["prepare"],
        "synteny": ["prepare", "identify"],
        "orthology": ["prepare", "identify"],
        "pan": ["orthology"],
        "properties": ["identify"],
        "tree": ["identify"],
        "expression": ["pan"],
        "enrichment": ["pan"],
        "te": ["identify"],
        "kaks": ["pan"],
        "duplication": ["identify"],
        "markers": ["pan"],
        "motifs": ["identify"],
    }
    if name in ["plot", "report"]:
        return [stage for stage in enabled(cfg) if stage not in ["plot", "report"]] + (
            ["plot"] if name == "report" else []
        )
    if (
        name == "orthology"
        and cfg.get("orthology", {}).get("backend") == "psonic"
        and not cfg["orthology"].get("prefix")
    ):
        return ["synteny"]
    return base[name]


def required_tools(cfg):
    tools = []
    identification = cfg.get("identification", {})
    if not identification.get("candidates"):
        tools.append("hmmsearch")
    backend = cfg.get("orthology", {}).get("backend", "import")
    if backend == "orthofinder":
        tools += [
            "orthofinder",
            "mcl",
            cfg.get("orthology", {}).get("search", "diamond"),
        ]
    if backend == "psonic" and not cfg["orthology"].get("prefix"):
        tools += ["orthofinder", "diamond", "mcl", "MCScanX"]
    if cfg.get("tree") is not None:
        tools += ["mafft", "FastTree"]
    if cfg.get("kaks") is not None:
        tools += ["mafft", "pal2nal.pl", "yn00"]
    if cfg.get("motifs") is not None:
        tools += ["meme", "fimo"]
    return sorted(set(tools))


def validate(cfg, deep=True, targets=None):
    active = set()

    def collect(name):
        if name in active:
            return
        active.add(name)
        for dependency in dependencies(cfg, name):
            collect(dependency)

    for name in targets or enabled(cfg):
        collect(name)
    for section, keys in SECTION_KEYS.items():
        value = cfg.get(section)
        if value is None:
            continue
        if not isinstance(value, dict):
            raise InputError(section + " must be a mapping")
        if set(value) - keys:
            raise InputError(f"Unknown {section} options: {sorted(set(value) - keys)}")
    if cfg.get("preparation", {}).get("invalid_protein", "error") not in [
        "error",
        "mask",
        "skip",
    ]:
        raise InputError("preparation.invalid_protein must be error, mask or skip")
    if cfg.get("pan", {}).get("unit", "species") not in ["species", "dataset"]:
        raise InputError("pan.unit must be species or dataset")
    if not 0 < float(cfg.get("pan", {}).get("softcore_fraction", 0.9)) <= 1:
        raise InputError("pan.softcore_fraction must be in (0,1]")
    samples = manifest(cfg)
    errors = []
    if (
        "identify" in active
        and not cfg.get("identification", {}).get("candidates")
        and not cfg.get("identification", {}).get("profiles")
    ):
        errors.append("identification requires candidates or profiles")
    orth = cfg.get("orthology", {})
    backend = orth.get("backend", "import")
    if backend not in ["import", "orthofinder", "psonic"]:
        errors.append("orthology.backend must be import, orthofinder or psonic")
    required = {
        "import": ["groups"],
        "orthofinder": [],
        "psonic": [
            "groups",
            "script",
            "prefix",
            "ploidy",
            "sequence_ids",
            "species_ids",
        ]
        if orth.get("prefix")
        else ["script"],
    }
    errors += [
        "orthology." + key + " is required"
        for key in required.get(backend, [])
        if "orthology" in active and not orth.get(key)
    ]
    sections = {
        "expression": ["dataset", "matrix", "samples"],
        "enrichment": ["dataset", "term2gene"],
        "te": ["intervals"],
        "duplication": ["table"],
        "markers": ["table"],
    }
    for section, keys in sections.items():
        if section in active and cfg.get(section) is not None:
            errors += [
                section + "." + key + " is required"
                for key in keys
                if not cfg[section].get(key)
            ]
    if errors:
        raise InputError("; ".join(errors))
    for section in PATHS:
        if ("identify" if section == "identification" else section) not in active:
            continue
        for path in external_inputs(cfg, section):
            if not path.is_file():
                raise InputError("Input file not found: " + str(path))
    if "orthology" in active and backend == "psonic" and orth.get("prefix"):
        header = list(read_tsv(resolve(cfg["_base"], orth["groups"]))[0])
        with resolve(cfg["_base"], orth["ploidy"]).open() as handle:
            ploidies = handle.readline().strip().split("\t")
        if len(ploidies) != len(header) - 1 or any(
            not value.isdigit() or int(value) < 1 for value in ploidies
        ):
            raise InputError(
                "pSONIC ploidy file requires one tab-delimited line of positive relative ploidies in OrthoFinder column order"
            )
        if set(header[1:]) != {row["dataset"] for row in samples}:
            raise InputError(
                "pSONIC OrthoFinder column names must equal manifest datasets"
            )
    proteins = 0
    if deep:
        for sample in samples:
            proteins += len(
                fasta(
                    sample["proteins"],
                    protein=True,
                    invalid_policy=cfg.get("preparation", {}).get(
                        "invalid_protein", "error"
                    ),
                )
            )
            if sample.get("gff"):
                gff_models(sample["gff"])
    return {
        "status": "valid",
        "datasets": len(samples),
        "protein_models": proteins if deep else "not_scanned",
        "optional_modules": [name for name in OPTIONAL if cfg.get(name) is not None],
    }


def doctor(cfg=None):
    required = (
        required_tools(cfg)
        if cfg
        else [
            "hmmsearch",
            "orthofinder",
            "diamond",
            "mafft",
            "FastTree",
            "pal2nal.pl",
            "yn00",
        ]
    )
    result = {
        "python": sys.version.split()[0],
        "platform": sys.platform,
        "packages": {},
        "executables": {},
        "required_missing": [],
    }
    for package in ["Bio", "numpy", "scipy", "matplotlib", "yaml"]:
        result["packages"][package] = importlib.util.find_spec(package) is not None
    for name in required:
        path = shutil.which(name) or (
            shutil.which("fasttree") if name == "FastTree" else None
        )
        result["executables"][name] = {
            "path": path,
            "sha256": sha256(path) if path else None,
        }
        if not path:
            result["required_missing"].append(name)
    if cfg and cfg.get("orthology", {}).get("backend") == "psonic":
        interpreter = cfg["orthology"].get("python", sys.executable)
        probe = subprocess.run(
            [interpreter, "-c", "import igraph; print(igraph.__version__)"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        result["psonic_igraph"] = (
            probe.stdout.strip() if probe.returncode == 0 else None
        )
        if probe.returncode:
            result["required_missing"].append("python-igraph for pSONIC")
    result["ready"] = (
        sys.platform.startswith("linux")
        and all(result["packages"].values())
        and not result["required_missing"]
    )
    return result


def execute(cfg, selected, args):
    validate(cfg, deep=False, targets=selected)
    runner = Runner(cfg, resume=args.resume, force=args.force, dry_run=args.dry_run)
    seen = set()
    params = {key: value for key, value in cfg.items() if key != "_config"}

    def visit(name):
        if name in seen:
            return
        if name in OPTIONAL and cfg.get(name) is None:
            raise InputError(f"{name} requires its configuration section")
        for dependency in dependencies(cfg, name):
            visit(dependency)
        seen.add(name)
        if args.dry_run:
            print(
                json.dumps(
                    {
                        "stage": name,
                        "destination": str(runner.directory(name)),
                        "dependencies": dependencies(cfg, name),
                    }
                )
            )
            return
        files = source_files(cfg) if name == "prepare" else []
        files += external_inputs(cfg, "identification" if name == "identify" else name)
        for dependency in dependencies(cfg, name):
            files += [
                path
                for path in runner.directory(dependency).rglob("*")
                if path.is_file()
            ]
        # Track executable content as well as configuration and data hashes.
        relevant = {
            "identify": ["hmmsearch"]
            if not cfg.get("identification", {}).get("candidates")
            else [],
            "orthology": [
                "orthofinder",
                cfg.get("orthology", {}).get("search", "diamond"),
            ]
            if cfg.get("orthology", {}).get("backend") == "orthofinder"
            else [],
            "tree": ["mafft", "FastTree"],
            "kaks": ["mafft", "pal2nal.pl", "yn00"],
        }.get(name, [])
        if name == "synteny":
            relevant = ["orthofinder", "diamond", "MCScanX"]
        if name == "motifs":
            relevant = ["meme", "fimo"]
        for executable in relevant:
            path = shutil.which(executable) or (
                shutil.which("fasttree") if executable == "FastTree" else None
            )
            if path:
                files.append(Path(path))
        runner.run(name, lambda work: FUNCTIONS[name](cfg, work), files, params)

    if args.dry_run:
        for stage in selected:
            visit(stage)
    else:
        with project_lock(cfg["output"]):
            log_path = Path(cfg["output"]) / "pipeline.log"
            handler = logging.FileHandler(log_path, encoding="utf-8")
            handler.setFormatter(
                logging.Formatter("%(asctime)s %(levelname)s %(message)s")
            )
            logging.getLogger().addHandler(handler)
            try:
                for stage in selected:
                    visit(stage)
            finally:
                logging.getLogger().removeHandler(handler)
                handler.close()
    return {
        "status": "planned" if args.dry_run else "complete",
        "stages": list(stage for stage in FUNCTIONS if stage in seen),
        "output": cfg["output"],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="polypangf",
        description="Mixed-ploidy pan-gene-family analysis with explicit assembly representation.",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    subs = parser.add_subparsers(dest="command", required=True)
    init = subs.add_parser(
        "init", help="Create an explicitly synthetic runnable example"
    )
    init.add_argument("directory")
    check = subs.add_parser(
        "doctor", help="Inspect Python packages and external executables"
    )
    check.add_argument("-c", "--config")
    check = subs.add_parser("validate", help="Check input schemas, FASTA and GFF3")
    check.add_argument("-c", "--config", required=True)
    check.add_argument("--quick", action="store_true")
    check = subs.add_parser(
        "status", help="Show completed stages and failed work directories"
    )
    check.add_argument("-c", "--config", required=True)
    check.add_argument(
        "--verify", action="store_true", help="Verify stored output SHA256 values"
    )
    check = subs.add_parser("export", help="Create a checksummed ZIP of results")
    check.add_argument("-c", "--config", required=True)
    check.add_argument("-o", "--archive", required=True)
    for name in ["run"] + list(FUNCTIONS):
        sub = subs.add_parser(
            name,
            help="Run configured workflow"
            if name == "run"
            else "Run " + name + " and its prerequisites",
        )
        sub.add_argument("-c", "--config", required=True)
        sub.add_argument(
            "-o",
            "--output",
            help="Override results directory, relative to configuration",
        )
        sub.add_argument("-t", "--threads", type=int)
        sub.add_argument("--seed", type=int, help="Override deterministic random seed")
        group = sub.add_mutually_exclusive_group()
        group.add_argument("--resume", action="store_true")
        group.add_argument("--force", action="store_true")
        sub.add_argument(
            "--dry-run",
            action="store_true",
            help="Print dependency plan without creating outputs",
        )
        verbosity = sub.add_mutually_exclusive_group()
        verbosity.add_argument("-v", "--verbose", action="store_true")
        verbosity.add_argument("-q", "--quiet", action="store_true")
        if name == "run":
            sub.add_argument(
                "--stages",
                help="Comma-separated target stages with automatic prerequisites",
            )
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG
        if getattr(args, "verbose", False)
        else logging.ERROR
        if getattr(args, "quiet", False)
        else logging.INFO,
        format="%(levelname)s %(message)s",
    )
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    try:
        if args.command == "init":
            result = initialize(args.directory)
        elif args.command == "doctor":
            result = doctor(load_config(args.config) if args.config else None)
            print(json.dumps(result, indent=2))
            return 0 if result["ready"] else 1
        else:
            cfg = load_config(
                args.config,
                getattr(args, "output", None),
                getattr(args, "threads", None),
            )
            if getattr(args, "seed", None) is not None:
                if args.seed < 0:
                    raise InputError("seed must be nonnegative")
                cfg["seed"] = args.seed
            if args.command == "validate":
                result = validate(cfg, deep=not args.quick)
            elif args.command == "status":
                root = Path(cfg["output"])
                result = {}
                for path in root.glob("*/state.json"):
                    state = json.loads(path.read_text())
                    result[path.parent.name] = {
                        "status": state["status"],
                        "summary": state["summary"],
                    }
                    if args.verify:
                        result[path.parent.name]["outputs_valid"] = all(
                            (path.parent / rel).is_file()
                            and sha256(path.parent / rel) == digest
                            for rel, digest in state["outputs"].items()
                        )
                result["failures"] = [
                    json.loads(path.read_text())
                    for path in (root / "failures").glob("*.json")
                ]
            elif args.command == "export":
                with project_lock(cfg["output"]):
                    result = export_result(cfg["output"], args.archive)
            else:
                selected = (
                    args.stages.split(",")
                    if args.command == "run" and args.stages
                    else enabled(cfg)
                    if args.command == "run"
                    else [args.command]
                )
                invalid = set(selected) - set(FUNCTIONS)
                if invalid:
                    raise InputError(f"Unknown stages: {sorted(invalid)}")
                result = execute(cfg, selected, args)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    except (
        InputError,
        KeyError,
        ValueError,
        OSError,
        subprocess.SubprocessError,
    ) as exc:
        logging.error("%s", exc)
        if getattr(args, "verbose", False):
            logging.exception("Details")
        return 2
    except KeyboardInterrupt:
        logging.error(
            "Interrupted; completed stages are retained. Retry with --resume."
        )
        return 130
