#!/usr/bin/env python3
"""Run the portable baseline and workflow experiments in dependency order."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


STAGES: list[tuple[str, list[str], bool]] = [
    ("baseline-polar", ["scripts/baselines/build_baseline_polar_mapping.py"], False),
    ("baseline-harmonic", ["scripts/baselines/build_baseline_harmonic_annulus.py"], False),
    ("baseline-conformal", ["scripts/baselines/build_baseline_conformal_annulus.py"], False),
    ("baseline-area", ["scripts/baselines/build_baseline_area_preserving_annulus.py"], True),
    ("step01", ["scripts/workflow/step01_build_center_polar_annulus.py"], False),
    ("step02", ["scripts/workflow/step02_build_center_topology.py"], False),
    ("step03", ["scripts/workflow/step03_build_line_md_scaffold.py"], False),
    ("step04-initial", ["scripts/workflow/step04_build_line_md_layer_column_partition.py"], False),
    (
        "step04-v2",
        [
            "scripts/workflow/step04_build_topology_aware_partition_v2.py",
            "--w-angle",
            "1.2",
            "--w-new",
            "4.5",
            "--w-angle-cap",
            "3.0",
            "--max-rank-shift",
            "2",
        ],
        True,
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=[name for name, _, _ in STAGES], help="Run one named stage only.")
    parser.add_argument("--skip-slow", action="store_true", help="Skip the area optimizer and topology-aware v2 stage.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    selected = [stage for stage in STAGES if args.only is None or stage[0] == args.only]
    for name, command, slow in selected:
        if args.skip_slow and slow:
            print(f"[skip] {name}")
            continue
        print(f"[run] {name}", flush=True)
        subprocess.run([sys.executable, *command], cwd=PROJECT_ROOT, check=True)
    print("All selected stages completed.")


if __name__ == "__main__":
    main()
