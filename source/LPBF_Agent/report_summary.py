#!/usr/bin/env python3
"""Regenerate the concise human report from an existing assessment.json."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from lpbf_score.reporting import render_summary


# Backward-compatible import used by older notebooks/tests.
render = render_summary


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render assessment.json as a concise one-page summary."
    )
    parser.add_argument("results_dir", help="run_all.py output directory or assessment.json")
    parser.add_argument(
        "--print", dest="only_print", action="store_true", help="Print without writing a file."
    )
    args = parser.parse_args()

    root = Path(args.results_dir).expanduser()
    source = root if root.suffix.lower() == ".json" else root / "assessment.json"
    if not source.is_file():
        raise FileNotFoundError(f"Not found: {source}")
    text = render_summary(json.loads(source.read_text(encoding="utf-8")))
    print(text)
    if not args.only_print:
        target = source.parent / "00_READ_ME_FIRST.md"
        target.write_text(text, encoding="utf-8", newline="\n")
        print(f"[written: {target}]", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
