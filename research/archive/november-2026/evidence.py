"""Inspect or replay the November 2026 research evidence."""

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tempfile
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parent


class Evidence:
    def __init__(self, inventory, archive):
        self.archive = archive
        self.hashes = inventory["sha256"]
        self.suites = inventory["suites"]
        self.paths = {name: name for name in self.hashes}
        aliases = inventory["aliases"]
        if self.paths.keys() & aliases.keys() or set(aliases.values()) - self.paths.keys():
            raise ValueError("Invalid historical aliases")
        self.paths.update(aliases)
        for name in self.paths:
            path = PurePosixPath(name)
            if not path.parts or path.is_absolute() or ".." in path.parts or str(path) != name or "\\" in name:
                raise ValueError(f"Invalid historical path: {name}")
        if sorted(archive.namelist()) != sorted(self.hashes):
            raise ValueError("Archive and inventory disagree")

    def read(self, name):
        canonical = self.paths[name]
        data = self.archive.read(canonical)
        if hashlib.sha256(data).hexdigest() != self.hashes[canonical]:
            raise ValueError(f"Digest mismatch: {name}")
        return data

    def extract(self, destination):
        destination.mkdir()
        for name in sorted(self.paths):
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(self.read(name))

    def replay(self, name):
        case = self.suites[name]
        entry = case["entry"]
        package = case.get("package")
        with tempfile.TemporaryDirectory(prefix="hopf-evidence-") as temporary:
            work = Path(temporary)
            if package:
                prefix = package + "/"
                for path in self.paths:
                    if path.startswith(prefix):
                        target = work / path.removeprefix(prefix)
                        target.parent.mkdir(parents=True, exist_ok=True)
                        target.write_bytes(self.read(path))
                executable = entry.removeprefix(prefix)
            else:
                executable = "check.py"
                (work / executable).write_bytes(self.read(entry))
            environment = os.environ.copy()
            environment.pop("PYTHONHOME", None)
            environment.pop("PYTHONPATH", None)
            environment.update(PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1",
                               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
            command = [sys.executable, "-I", "-B", executable, *case.get("options", [])]
            print(f"Running {name}", flush=True)
            subprocess.run(command, cwd=work, env=environment, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)
    subcommands.add_parser("verify", help="check the inventory and every stored source digest")
    listing = subcommands.add_parser("list", help="list original paths beginning with PREFIX")
    listing.add_argument("prefix", nargs="?", default="")
    show = subcommands.add_parser("show", help="write one original record to standard output")
    show.add_argument("path")
    extract = subcommands.add_parser("extract", help="restore every historical path into a new directory")
    extract.add_argument("directory", type=Path)
    replay = subcommands.add_parser("replay", help="run selected suites, or all suites by default")
    replay.add_argument("suites", metavar="SUITE", nargs="*")
    args = parser.parse_args()
    inventory = json.loads((ROOT / "inventory.json").read_text(encoding="utf-8"))
    with ZipFile(ROOT / "sources.zip") as archive:
        evidence = Evidence(inventory, archive)
        if args.command == "verify":
            for path in evidence.hashes:
                evidence.read(path)
            print(f"Verified {len(evidence.hashes)} distinct records and {len(evidence.paths)} source paths.")
        elif args.command == "list":
            for path in sorted(evidence.paths):
                if path.startswith(args.prefix):
                    print(path)
        elif args.command == "show":
            if args.path not in evidence.paths:
                parser.error(f"Unknown source path: {args.path}")
            sys.stdout.buffer.write(evidence.read(args.path))
        elif args.command == "extract":
            evidence.extract(args.directory)
        elif args.command == "replay":
            selected = args.suites or list(evidence.suites)
            for name in selected:
                if name not in evidence.suites:
                    parser.error("Available suites: " + ", ".join(evidence.suites))
            for name in selected:
                evidence.replay(name)
            print(f"{len(selected)} suite entries passed.")


if __name__ == "__main__":
    main()
