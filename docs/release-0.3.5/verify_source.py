"""Verify frozen package inputs, statistical scope and optional archive/index."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/release-0.3.5"
parser = argparse.ArgumentParser()
parser.add_argument("--archive", type=Path)
parser.add_argument("--git-index", action="store_true")
args = parser.parse_args()


def normalized(name, data):
    # Include CSV fixtures and startup hooks as text; retain binary bytes.
    if b"\0" not in data:
        try:
            data.decode("utf-8")
        except UnicodeDecodeError:
            return data
        data = data.replace(b"\r\n", b"\n")
    return data


manifest = json.loads((OUT / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
assert manifest["version"] == "0.3.5" and manifest["status"] == "release"
for item in manifest["files"]:
    name = item["path"]
    assert hashlib.sha256(normalized(name, (ROOT / name).read_bytes())).hexdigest() == item["sha256"], name
    if args.git_index:
        data = subprocess.check_output(["git", "show", ":" + name], cwd=ROOT)
        assert hashlib.sha256(normalized(name, data)).hexdigest() == item["sha256"], "index: " + name
scope = json.loads((OUT / "STATISTICAL_SCOPE.json").read_text(encoding="utf-8"))
old_scope = json.loads((ROOT / "docs/release-0.3.4/final/STATISTICAL_SCOPE.json").read_text(encoding="utf-8"))
assert {k: v for k, v in scope.items() if k not in {"version", "status", "scope_unchanged_from"}} == {
    k: v for k, v in old_scope.items() if k != "version"}
validation = json.loads((OUT / "evidence/VALIDATION.json").read_text(encoding="utf-8"))
assert validation["check"] == "Status: OK"
assert validation["tests"]["fail"] == validation["tests"]["warn"] == 0
assert validation["public_interface_identical"] and validation["native_text_identical"]
if args.archive:
    assert args.archive.name == validation["source_archive"]
    assert hashlib.sha256(args.archive.read_bytes()).hexdigest() == validation["source_sha256"]
print(json.dumps(dict(version="0.3.5", frozen_inputs=len(manifest["files"]),
    scope_unchanged=True, git_index_checked=args.git_index, archive_checked=bool(args.archive))))
