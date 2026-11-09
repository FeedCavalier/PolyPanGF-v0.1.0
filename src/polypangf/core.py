"""Configuration, validated tabular I/O and transactional stage execution."""

import csv
import gzip
import hashlib
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import time
import uuid
from contextlib import contextmanager

from . import __version__


class InputError(ValueError):
    pass


def open_text(path):
    path = Path(path)
    return (
        gzip.open(path, "rt", encoding="utf-8-sig")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8-sig")
    )


def read_tsv(path, required=()):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
            raise InputError(f"{path}: required columns {list(required)}")
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise InputError(f"{path}: duplicate TSV column names")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise InputError(f"{path}: malformed TSV row")
    return rows


def write_tsv(path, rows, fields):
    path = Path(path)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fields, delimiter="\t", lineterminator="\n", extrasaction="raise"
        )
        writer.writeheader()
        writer.writerows(rows)


def save_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        temporary.write_text(
            json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve(base, value):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (Path(base) / path).resolve()


def load_config(path, output=None, threads=None):
    import yaml

    path = Path(path).resolve()
    try:
        with path.open(encoding="utf-8-sig") as handle:
            cfg = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        raise InputError(f"Invalid JSON/YAML configuration: {exc}") from exc
    if not isinstance(cfg, dict) or cfg.get("schema_version") != 1:
        raise InputError("Configuration must contain schema_version: 1")
    allowed = {
        "schema_version",
        "family",
        "manifest",
        "output",
        "threads",
        "seed",
        "preparation",
        "identification",
        "orthology",
        "pan",
        "expression",
        "enrichment",
        "te",
        "markers",
        "tree",
        "kaks",
        "duplication",
        "motifs",
    }
    unknown = set(cfg) - allowed
    if unknown:
        raise InputError(f"Unknown configuration sections: {sorted(unknown)}")
    cfg["manifest"] = str(resolve(path.parent, cfg["manifest"]))
    cfg["output"] = str(resolve(path.parent, output or cfg.get("output", "results")))
    cfg["threads"] = threads or cfg.get("threads", 4)
    cfg["seed"] = cfg.get("seed", 42)
    if not isinstance(cfg["threads"], int) or cfg["threads"] < 1:
        raise InputError("threads must be a positive integer")
    if not isinstance(cfg["seed"], int) or cfg["seed"] < 0:
        raise InputError("seed must be a nonnegative integer")
    cfg["_base"] = str(path.parent)
    cfg["_config"] = str(path)
    return cfg


def tool(name):
    result = shutil.which(name)
    if not result:
        raise InputError(
            f"Required executable not found: {name}. Run polypangf doctor."
        )
    return result


def command(args, cwd, log, stdout=None, timeout=None):
    logging.info("Running %s", " ".join(map(str, args)))
    with Path(log).open("a", encoding="utf-8") as handle:
        handle.write(json.dumps([str(arg) for arg in args]) + "\n")
        handle.flush()
        subprocess.run(
            [str(arg) for arg in args],
            cwd=cwd,
            stdout=stdout or handle,
            stderr=handle,
            stdin=subprocess.DEVNULL,
            check=True,
            timeout=timeout,
        )


@contextmanager
def project_lock(output):
    """OS releases this advisory lock after a crash. Stale metadata cannot block retries."""
    import fcntl

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".polypangf.lock").open("a+") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise InputError(f"Another PolyPanGF process is writing {output}") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


class Runner:
    def __init__(self, cfg, resume=False, force=False, dry_run=False):
        self.cfg = cfg
        self.output = Path(cfg["output"])
        self.resume = resume
        self.force = force
        self.dry_run = dry_run

    def directory(self, stage):
        return self.output / stage

    def run(self, name, function, inputs, parameters):
        inputs = sorted({Path(path).resolve() for path in inputs})
        missing = [str(p) for p in inputs if not p.is_file()]
        if missing:
            raise InputError(f"{name}: missing input files {missing}")
        fingerprint = {
            "version": __version__,
            "parameters": parameters,
            "inputs": {str(path): sha256(path) for path in inputs},
            "code": {p.name: sha256(p) for p in Path(__file__).parent.glob("*.py")},
        }
        destination = self.directory(name)
        state_path = destination / "state.json"
        valid = False
        if state_path.exists():
            try:
                state = json.loads(state_path.read_text())
                valid = (
                    state.get("status") == "complete"
                    and bool(state.get("outputs"))
                    and state.get("fingerprint") == fingerprint
                    and all(
                        (destination / rel).is_file()
                        and sha256(destination / rel) == digest
                        for rel, digest in state.get("outputs", {}).items()
                    )
                )
            except (json.JSONDecodeError, OSError):
                valid = False
        if valid and self.resume and not self.force:
            logging.info("%s is current, reusing verified outputs", name)
            return destination
        if self.dry_run:
            print(
                json.dumps(
                    {
                        "stage": name,
                        "action": "reuse" if valid else "run",
                        "inputs": [str(p) for p in inputs],
                    }
                )
            )
            return destination
        if destination.exists() and not (self.force or self.resume):
            raise InputError(f"{destination} exists. Use --resume or --force.")
        work = self.output / ("." + name + "." + uuid.uuid4().hex)
        work.mkdir(parents=True)
        started = time.time()
        try:
            summary = function(work) or {}
            if any(
                not path.is_file() or sha256(path) != fingerprint["inputs"][str(path)]
                for path in inputs
            ):
                raise InputError(
                    f"{name}: an input changed during execution; retry after inputs are stable"
                )
            outputs = {
                str(p.relative_to(work)): sha256(p)
                for p in work.rglob("*")
                if p.is_file()
            }
            save_json(
                work / "state.json",
                {
                    "status": "complete",
                    "stage": name,
                    "fingerprint": fingerprint,
                    "outputs": outputs,
                    "summary": summary,
                    "elapsed_seconds": time.time() - started,
                },
            )
            backup = None
            if destination.exists():
                backup = self.output / ("." + name + ".previous." + uuid.uuid4().hex)
                destination.rename(backup)
            try:
                work.rename(destination)
            except BaseException:
                if backup:
                    backup.rename(destination)
                raise
            if backup:
                shutil.rmtree(backup)
            failure = self.output / "failures" / (name + ".json")
            failure.unlink(missing_ok=True)
            return destination
        except BaseException as exc:
            failed = self.output / "failures"
            failed.mkdir(exist_ok=True)
            retained = failed / (name + "." + uuid.uuid4().hex)
            if work.exists():
                work.rename(retained)
            save_json(
                failed / (name + ".json"),
                {"status": "failed", "error": str(exc), "work": str(retained)},
            )
            raise
