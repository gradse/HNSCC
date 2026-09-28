#!/usr/bin/env python3
"""Run one portable Python preprocessing step in a chosen data directory."""

import argparse
import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "scripts" / "python"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("step", help="Python filename under scripts/python")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT / "data" / "local",
        help="Directory containing local inputs and receiving generated outputs",
    )
    args = parser.parse_args()

    step = (SCRIPT_DIR / args.step).resolve()
    if step.parent != SCRIPT_DIR.resolve() or not step.is_file():
        parser.error(f"Unknown step: {args.step}")
    if not args.data_dir.is_dir():
        parser.error(f"Data directory does not exist: {args.data_dir}")

    os.chdir(args.data_dir.resolve())
    runpy.run_path(str(step), run_name="__main__")


if __name__ == "__main__":
    main()
