#!/usr/bin/env python3
"""Run the two core analyses in the clean submission package."""

from __future__ import annotations

import argparse
from pathlib import Path

from scripts.host_network import main as run_host_network


def main() -> None:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-structure", action="store_true")
    parser.add_argument("--bootstraps", type=int, default=10_000)
    parser.add_argument("--permutations", type=int, default=5_000)
    args = parser.parse_args()
    host_dir = root / "results" / "host_network"
    print(run_host_network(bootstraps=args.bootstraps, permutations=args.permutations,
                           output_dir=host_dir))
    if not args.skip_structure:
        from scripts.structural_support import main as run_structural_support
        print(run_structural_support(host_dir=host_dir,
                                     output_dir=root / "results" / "structural_support"))


if __name__ == "__main__":
    main()
