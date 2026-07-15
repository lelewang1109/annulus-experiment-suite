#!/usr/bin/env python3
"""Perform lightweight repository and result integrity checks."""

from __future__ import annotations

import ast
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    "README.md",
    "data/input/beijing_grid_cells.csv",
    "data/input/beijing_boundary_rays.csv",
    "results/baselines/polar/baseline_polar_mapping_metrics.csv",
    "results/baselines/harmonic/baseline_harmonic_annulus_metrics.csv",
    "results/baselines/conformal/baseline_conformal_annulus_metrics.csv",
    "results/baselines/area_preserving/baseline_area_preserving_annulus_metrics.csv",
    "results/workflow/step03/step03_topology_constraints.csv",
    "results/workflow/step04_topology_aware_v2/step04_topology_aware_metrics_v2.json",
]


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    for relative in REQUIRED:
        if not (ROOT / relative).is_file():
            fail(f"missing required file: {relative}")

    with (ROOT / "data/input/beijing_grid_cells.csv").open(encoding="utf-8") as handle:
        cells = list(csv.DictReader(handle))
    if len(cells) != 176 or len({row["cell_id"] for row in cells}) != 176:
        fail("canonical input must contain 176 unique cell_id values")

    with (ROOT / "results/workflow/step04_topology_aware_v2/step04_topology_aware_metrics_v2.json").open(encoding="utf-8") as handle:
        metrics = json.load(handle)
    optimized = metrics.get("optimized_effective_metrics", {})
    for key in ("adjacency_precision", "adjacency_recall", "adjacency_f1"):
        if key not in optimized:
            fail(f"optimized metrics missing key: {key}")

    expected_metrics = {
        "adjacency_precision": 0.595,
        "adjacency_recall": 0.7628205128205128,
        "adjacency_f1": 0.6685393258426966,
        "preserved_edge_count": 238,
        "lost_edge_count": 74,
        "new_edge_count": 162,
    }
    for key, expected in expected_metrics.items():
        actual = optimized.get(key)
        if actual is None or not math.isclose(float(actual), float(expected), rel_tol=1e-12, abs_tol=1e-12):
            fail(f"optimized metric {key}={actual!r}, expected {expected!r}")

    expected_parameters = {"w_new": 4.5, "w_angle": 1.2, "w_angle_cap": 3.0, "max_rank_shift": 2}
    parameters = metrics.get("parameters", {})
    for key, expected in expected_parameters.items():
        if parameters.get(key) != expected:
            fail(f"v2 parameter {key}={parameters.get(key)!r}, expected {expected!r}")
    for key in ("input_csv", "topology_csv", "output_dir"):
        if Path(str(parameters.get(key, ""))).is_absolute():
            fail(f"v2 parameter {key} must use a repository-relative path")

    with (ROOT / "results/workflow/step03/step03_topology_constraints.csv").open(encoding="utf-8") as handle:
        if sum(1 for _ in csv.DictReader(handle)) != 312:
            fail("topology constraints must contain 312 original edges")
    with (ROOT / "results/workflow/step04_topology_aware_v2/step04_topology_aware_assignment_v2.csv").open(encoding="utf-8") as handle:
        if sum(1 for _ in csv.DictReader(handle)) != 176:
            fail("topology-aware assignment must contain 176 cells")

    forbidden = []
    for path in ROOT.rglob("*"):
        if any(part in {".git", ".venv", "__pycache__"} for part in path.parts):
            continue
        if path.is_dir():
            continue
        if path.name == ".DS_Store":
            forbidden.append(str(path.relative_to(ROOT)))
        if path.suffix in {".py", ".md", ".toml", ".yaml", ".yml"}:
            text = path.read_text(encoding="utf-8", errors="ignore")
            if ("/" + "Users/") in text:
                forbidden.append(f"absolute path in {path.relative_to(ROOT)}")
    if forbidden:
        fail("repository contains non-portable artifacts: " + ", ".join(forbidden))

    for path in (ROOT / "scripts").rglob("*.py"):
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except SyntaxError as exc:
            fail(f"syntax error in {path.relative_to(ROOT)}: {exc}")

    oversized = [
        path.relative_to(ROOT)
        for path in ROOT.rglob("*")
        if path.is_file() and ".git" not in path.parts and path.stat().st_size >= 100_000_000
    ]
    if oversized:
        fail("files exceed GitHub's 100 MB limit: " + ", ".join(map(str, oversized)))

    print(f"OK: {len(cells)} cells; portable inputs and expected results are present.")


if __name__ == "__main__":
    main()
