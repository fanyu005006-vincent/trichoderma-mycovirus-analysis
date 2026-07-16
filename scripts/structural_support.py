#!/usr/bin/env python3
"""Add a secondary AlphaFold 3 and TM-align support layer to the host network."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

try:
    from tmtools import tm_align
    from tmtools.io import get_residue_data, get_structure
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Install tmtools from requirements.txt to run structural support.") from exc


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data" / "af3_structure_manifest.csv"
DEFAULT_STRUCTURES = ROOT / "data" / "af3_structures"
DEFAULT_HOST_NETWORK = ROOT / "results" / "host_network"
DEFAULT_OUTPUT = ROOT / "results" / "structural_support"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def number(value: str | None) -> float:
    return 0.0 if value in (None, "", "NA", "None") else float(value)


def structure_path(root: Path, row: dict[str, str]) -> Path:
    path = Path(row["structure_file"])
    return path if path.is_absolute() else root / path


def chain(path: Path) -> tuple[object, str]:
    structure = get_structure(str(path), format="mmcif")
    chains = []
    for value in structure.get_chains():
        coords, sequence = get_residue_data(value)
        if sequence:
            chains.append((len(sequence), coords, sequence))
    if not chains:
        raise ValueError(f"No protein chain found in {path}")
    _, coords, sequence = max(chains, key=lambda item: item[0])
    return coords, sequence


def compare(candidate: Path, reference: Path) -> dict[str, float]:
    coords_a, sequence_a = chain(candidate)
    coords_b, sequence_b = chain(reference)
    result = tm_align(coords_a, coords_b, sequence_a, sequence_b)
    aligned = [(a, b) for a, b in zip(result.seqxA, result.seqyA) if a != "-" and b != "-"]
    length = len(aligned)
    coverage = length / min(len(sequence_a), len(sequence_b)) if length else 0.0
    identity = sum(a == b for a, b in aligned) / length if length else 0.0
    tm_min = min(float(result.tm_norm_chain1), float(result.tm_norm_chain2))
    tm_mean = (float(result.tm_norm_chain1) + float(result.tm_norm_chain2)) / 2
    return {"candidate_length": len(sequence_a), "reference_length": len(sequence_b),
            "tm_score_min": tm_min, "tm_score_mean": tm_mean,
            "rmsd": float(result.rmsd), "aligned_length": length,
            "structure_coverage": coverage, "structure_identity": identity}


def support_score(metrics: dict[str, float], candidate: dict, reference: dict) -> float:
    """Fixed descriptive score; it is not a host-range probability."""
    tm = max(0.0, min(1.0, metrics["tm_score_min"]))
    coverage = max(0.0, min(1.0, metrics["structure_coverage"]))
    identity = max(0.0, min(1.0, metrics["structure_identity"]))
    confidence = math.sqrt(number(candidate.get("af3_ranking_score"))
                            * number(reference.get("af3_ranking_score")))
    return max(0.0, min(1.0, 0.50 * tm + 0.20 * coverage + 0.15 * identity + 0.15 * confidence))


def main(manifest_path: Path = DEFAULT_MANIFEST, structure_root: Path = DEFAULT_STRUCTURES,
         host_dir: Path = DEFAULT_HOST_NETWORK, output_dir: Path = DEFAULT_OUTPUT) -> dict:
    manifest = read_csv(manifest_path)
    candidates = [row for row in manifest if row.get("role") == "candidate"]
    references = [row for row in manifest if row.get("role") == "reference"]
    for row in manifest:
        path = structure_path(structure_root, row)
        if not path.exists():
            raise FileNotFoundError(path)

    network = read_csv(host_dir / "potential_host_network_edges.csv")
    accepted = read_csv(host_dir / "accepted_accession_evidence.csv")
    hosts_by_accession = {(row.get("query", ""), row.get("accession", "")): row.get("host_genus", "")
                          for row in accepted}
    pairs = []
    for candidate in candidates:
        for reference in references:
            if candidate.get("family") != reference.get("family"):
                continue
            profile = candidate.get("profile", "")
            accession = reference.get("accession", "")
            try:
                metrics = compare(structure_path(structure_root, candidate),
                                  structure_path(structure_root, reference))
            except Exception as exc:  # pragma: no cover
                print(f"[warn] skipped {profile} vs {accession}: {exc}")
                continue
            confidence = math.sqrt(number(candidate.get("af3_ranking_score"))
                                    * number(reference.get("af3_ranking_score")))
            host = hosts_by_accession.get((profile, accession), "")
            row = {"candidate_profile": profile, "candidate_labels": candidate.get("candidate_labels", ""),
                   "family": candidate.get("family", ""), "accession": accession,
                   "reference_virus": reference.get("source_virus", ""), "host_genus": host,
                   "candidate_structure": candidate.get("structure_file", ""),
                   "reference_structure": reference.get("structure_file", ""),
                   "af3_confidence_weight": confidence}
            row.update(metrics)
            row["af3_structure_support_score"] = support_score(metrics, candidate, reference)
            row["concordance_class"] = ("high_concordance" if metrics["tm_score_min"] >= 0.5
                                         and metrics["structure_coverage"] >= 0.5
                                         else "limited_concordance")
            pairs.append(row)

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "af3_structure_pair_scores.csv", pairs)
    pair_index = {(row["candidate_profile"], row["accession"]): row for row in pairs}
    annotated = []
    for edge in network:
        row = dict(edge)
        pair = pair_index.get((edge.get("query", ""), edge.get("accession", "")))
        row["af3_structure_status"] = "matched" if pair else "unavailable_for_edge"
        for key in ("tm_score_min", "tm_score_mean", "rmsd", "aligned_length", "structure_coverage",
                    "structure_identity", "af3_confidence_weight", "af3_structure_support_score"):
            row[f"af3_{key}"] = pair.get(key, "") if pair else ""
        annotated.append(row)
    write_csv(output_dir / "potential_host_network_edges_with_af3.csv", annotated)
    matched = sum(row["af3_structure_status"] == "matched" for row in annotated)
    summary = {"candidate_structure_count": len(candidates), "reference_structure_count": len(references),
               "structure_pair_count": len(pairs), "network_edge_count": len(network),
               "network_edges_with_structure_support": matched,
               "score": "0.50*TM-score + 0.20*coverage + 0.15*identity + 0.15*AF3 confidence",
               "interpretation": "secondary structural concordance, not a calibrated host-range probability"}
    (output_dir / "structural_support_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--structure-root", type=Path, default=DEFAULT_STRUCTURES)
    parser.add_argument("--host-dir", type=Path, default=DEFAULT_HOST_NETWORK)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(main(args.manifest, args.structure_root, args.host_dir, args.output_dir), indent=2))
