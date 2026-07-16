#!/usr/bin/env python3
"""Build the manuscript-facing potential host-network model.

The model uses saved BLASTp alignments plus NCBI ``source/host`` metadata.
Each observed query is kept as an independent profile.  It ranks up to ten
distinct host genera with ``bitscore * identity * query_coverage`` and reports
weighted bootstrap support and within-family permutation enrichment.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "host_prediction_inputs"
DEFAULT_METADATA = ROOT / "data" / "reference_host_metadata.csv"
DEFAULT_OUTPUT = ROOT / "results" / "host_network"

EVALUE = 1e-5
IDENTITY = 0.25
COVERAGE = 0.40
TOP_HOSTS = 10
BOOTSTRAPS = 10_000
PERMUTATIONS = 5_000
SEED = 20260715
IDENTICAL_PROFILE_SOURCES = {"CV2": "CV1"}

QUERY_RE = re.compile(r"^(?:CV|MIV|HV|MYV|BV)\d+(?:_RNA\d+)?$", re.I)
SEQ_RE = re.compile(r"^Sequence ID:\s*(\S+?)Length:\s*(\d+)")
RANGE_RE = re.compile(r"^Range\s+\d+:\s*(\d+)\s+to\s+(\d+)")
IDENTITY_RE = re.compile(r"(\d+)\s*/\s*(\d+)\((\d+)%\)")
EVALUE_RE = re.compile(r"bits\([^)]*\)\s+([^\s]+)")
BITS_RE = re.compile(r"^(\d+(?:\.\d+)?)\s+bits")
HOST_RE = re.compile(r"\[([^\]]+)\]\s*$")

NAMES = {
    "CV1": "TsCV1", "CV2": "TsCV2", "MIV1": "TsMiV1", "MIV2": "TsMiV2",
    "MIV3": "TsMiV3", "MIV5": "TsMiV5", "MIV7": "TsMiV7",
    "HV1": "TsHV1", "HV2": "TsHV2", "HV3": "TsHV3", "HV4": "TsHV4",
    "MYV2": "TsMyV2", "MYV3": "TsMyV3", "MYV5": "TsMyV5", "BV1": "TsBLV1",
}
FAMILIES = {
    "CV": "Chrysoviridae", "MIV": "Mitoviridae", "HV": "Hypoviridae",
    "MYV": "Mymonaviridae", "BV": "bunya-like lineage",
}
TICKS = {"Dermacentor", "Haemaphysalis", "Rhipicephalus", "Ixodes", "Amblyomma",
         "Hyalomma", "Argas", "Ornithodoros"}
OOMYCETES = {"Plasmopara", "Bremia", "Phytophthora"}
ARTHROPODS = {"Culex", "Diaphorina", "Apis"}
PLANT_HOSTS = {"Vitis", "Yellow silver pine"}
NON_SPECIFIC = {"plant", "grapevine", "unknown", "unidentified", "environment",
                "environmental", "soil", "water", "fungal", "unresolved"}


def family(query: str) -> str:
    upper = query.upper()
    return next((name for prefix, name in FAMILIES.items() if upper.startswith(prefix)), "unknown")


def category(genus: str) -> str:
    if genus in OOMYCETES:
        return "oomycete-associated"
    if genus in ARTHROPODS:
        return "animal/arthropod-associated"
    if genus in PLANT_HOSTS:
        return "plant-associated"
    return "fungal-associated"


def normalize_host(value: str) -> tuple[str, str]:
    text = " ".join(value.strip().split())
    lower = text.lower()
    if not text:
        return "", "source_host_missing"
    if text == "Vitis vinifera":
        return "Vitis", "reported_host_override"
    if text == "Yellow silver pine":
        return text, "reported_host_override"
    if (lower.startswith(("snownnecut", "unresolved", "mymonaviridae"))
            or lower in NON_SPECIFIC or lower.startswith(("plant ", "grapevine"))):
        return "", "nonspecific_host"
    first = re.split(r"\s+", text, maxsplit=1)[0].strip("[](),.;:")
    if not re.fullmatch(r"[A-Za-z][A-Za-z-]*", first or ""):
        return "", "unresolved_host_name"
    genus = first[0].upper() + first[1:]
    if genus in TICKS:
        return "", "tick_environment_host"
    return genus, "accepted"


def _query_context(lines: list[str], index: int) -> tuple[str, int]:
    query = lines[index].strip().upper()
    sequence: list[str] = []
    for line in lines[index + 1:]:
        value = line.strip()
        if value and "[" in value and "]" in value:
            break
        if value and not value.startswith("Sequence ID:"):
            sequence.append(re.sub(r"[^A-Za-z*]", "", value))
            if value.endswith("*"):
                break
    return query, len("".join(sequence).rstrip("*"))


def _description_index(lines: list[str], index: int) -> int | None:
    for position in range(index - 1, -1, -1):
        value = lines[position].strip()
        if value and "[" in value and "]" in value:
            return position
        if QUERY_RE.fullmatch(value) or value.startswith("Sequence ID:"):
            break
    return None


def parse_blast(path: Path) -> list[dict]:
    """Parse the saved BLASTp text without inferring host from virus names."""
    lines = path.read_text(errors="ignore").splitlines()
    boundaries = [i for i, line in enumerate(lines)
                  if QUERY_RE.fullmatch(line.strip()) or line.startswith("Sequence ID:")]
    query = ""
    query_length = 0
    rows: list[dict] = []
    for index, line in enumerate(lines):
        if QUERY_RE.fullmatch(line.strip()):
            query, query_length = _query_context(lines, index)
            continue
        match = SEQ_RE.match(line.strip())
        if not match or not query or not query_length:
            continue
        accession, hit_length = match.groups()
        description_index = _description_index(lines, index)
        if description_index is None:
            continue
        description = lines[description_index].strip()
        virus_match = HOST_RE.search(description)
        stats = None
        qstart = qend = None
        for lookahead in range(index + 1, min(index + 18, len(lines))):
            range_match = RANGE_RE.match(lines[lookahead].strip())
            if range_match:
                qstart, qend = map(int, range_match.groups())
            if "Identities" in lines[lookahead] and lookahead + 1 < len(lines):
                stats = lines[lookahead + 1].strip()
                break
        bits = BITS_RE.match(stats or "")
        evalue = EVALUE_RE.search(stats or "")
        identities = IDENTITY_RE.search(stats or "")
        if not (virus_match and bits and evalue and identities):
            continue
        identity_n = int(identities.group(1))
        identity_d = int(identities.group(2))
        identity = identity_n / identity_d if identity_d else 0.0
        coverage = ((qend - qstart + 1) / query_length if qstart is not None else
                    identity_d / query_length)
        next_boundary = next((b for b in boundaries if b > description_index), len(lines))
        rows.append({
            "query": query, "manuscript_virus": NAMES.get(query, query),
            "candidate_labels": NAMES.get(query, query), "family": family(query),
            "accession": accession, "reference_virus": virus_match.group(1).strip(),
            "query_length_aa": query_length, "reference_length_aa": int(hit_length),
            "alignment_length_aa": identity_d, "query_start": qstart, "query_end": qend,
            "bitscore": float(bits.group(1)), "evalue": float(evalue.group(1))
            if evalue.group(1).lower() not in {"0", "0.0"} else 0.0,
            "identity": identity, "identity_percent": identity * 100,
            "query_coverage": min(1.0, coverage), "query_coverage_percent": min(1.0, coverage) * 100,
            "raw_description": description, "source_file": path.name,
            "profile_provenance": "observed_input",
            "raw_block_start": description_index, "raw_block_end": next_boundary,
        })
    return rows


def add_identical_profiles(rows: list[dict]) -> list[dict]:
    """Create explicit independent profiles for manuscript-confirmed identical sequences."""
    observed = {row["query"] for row in rows}
    derived: list[dict] = []
    for target, source in IDENTICAL_PROFILE_SOURCES.items():
        if target in observed or source not in observed:
            continue
        for row in rows:
            if row["query"] != source:
                continue
            clone = dict(row)
            clone["query"] = target
            clone["manuscript_virus"] = NAMES.get(target, target)
            clone["candidate_labels"] = NAMES.get(target, target)
            clone["family"] = family(target)
            clone["profile_provenance"] = f"copied_from_{source}_identical_RdRp"
            derived.append(clone)
    return rows + derived


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["accession"]: row for row in csv.DictReader(handle)}


def annotate(rows: list[dict], metadata: dict[str, dict[str, str]]) -> None:
    for row in rows:
        record = metadata.get(row["accession"], {})
        row["source_host"] = record.get("source_host", "") or ""
        row["source_organism"] = record.get("source_organism", "")
        row["metadata_source"] = record.get("metadata_source", "")
        reported_host = record.get("reported_host_override", "") or ""
        row["reported_host"] = reported_host
        host_for_network = reported_host or row["source_host"]
        row["host_genus"], row["host_filter_reason"] = normalize_host(host_for_network)
        if reported_host:
            row["host_filter_reason"] = "reported_host_override"
        row["host_category"] = category(row["host_genus"]) if row["host_genus"] else "excluded"
        row["threshold_status"] = (row["evalue"] <= EVALUE and row["identity"] >= IDENTITY
                                    and row["query_coverage"] >= COVERAGE)


def strongest(rows: list[dict], key_fields: tuple[str, ...]) -> list[dict]:
    result: dict[tuple, dict] = {}
    for row in rows:
        key = tuple(row[field] for field in key_fields)
        if key not in result or row["raw_evidence_score"] > result[key]["raw_evidence_score"]:
            result[key] = row
    return list(result.values())


def score_and_select(rows: list[dict], top_n: int) -> tuple[list[dict], list[dict]]:
    accepted = [row for row in rows if row["threshold_status"] and row["host_genus"]]
    for row in accepted:
        row["raw_evidence_score"] = row["bitscore"] * row["identity"] * row["query_coverage"]
    accession_rows = strongest(accepted, ("query", "accession"))
    host_rows = strongest(accession_rows, ("query", "host_genus"))
    selected: list[dict] = []
    for query in sorted({row["query"] for row in host_rows}):
        query_rows = sorted((row for row in host_rows if row["query"] == query),
                            key=lambda r: (-r["raw_evidence_score"], r["host_genus"]))[:top_n]
        total = sum(row["raw_evidence_score"] for row in query_rows) or 1.0
        for rank, row in enumerate(query_rows, 1):
            row["evidence_weight"] = row["raw_evidence_score"] / total
            row["host_rank"] = rank
            selected.append(row)
    return accession_rows, selected


def bootstrap(rows: list[dict], repeats: int, rng: random.Random) -> list[dict]:
    output: list[dict] = []
    for query in sorted({row["query"] for row in rows}):
        query_rows = [row for row in rows if row["query"] == query]
        weights = [row["evidence_weight"] for row in query_rows]
        categories = sorted({row["host_category"] for row in query_rows})
        wins = Counter()
        for _ in range(repeats):
            sample = rng.choices(query_rows, weights=weights, k=len(query_rows))
            counts = Counter(row["host_category"] for row in sample)
            maximum = max(counts.values())
            tied = [name for name in categories if counts[name] == maximum]
            for name in tied:
                wins[name] += 1 / len(tied)
        for name in categories:
            output.append({"query": query, "candidate_labels": NAMES.get(query, query),
                           "family": family(query), "host_category": name,
                           "host_count": sum(row["host_category"] == name for row in query_rows),
                           "weighted_support_mass": sum(row["evidence_weight"] for row in query_rows
                                                        if row["host_category"] == name),
                           "bootstrap_support": wins[name] / repeats})
    return output


def bh(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    output = [1.0] * len(values)
    running = 1.0
    for rank in range(len(order), 0, -1):
        index = order[rank - 1]
        running = min(running, values[index] * len(values) / rank)
        output[index] = min(1.0, running)
    return output


def permutation(rows: list[dict], repeats: int, rng: random.Random) -> list[dict]:
    by_query = defaultdict(list)
    by_family = defaultdict(list)
    for row in rows:
        by_query[row["query"]].append(row)
        by_family[row["family"]].append(row)
    tests: list[dict] = []
    for query in sorted(by_query):
        family_name = by_query[query][0]["family"]
        categories = sorted({row["host_category"] for row in by_family[family_name]})
        observed = Counter(row["host_category"] for row in by_query[query])
        for name in categories:
            tests.append({"query": query, "candidate_labels": NAMES.get(query, query),
                          "family": family_name, "host_category": name,
                          "observed_count": observed.get(name, 0),
                          "query_n": len(by_query[query]), "family_background_n": len(by_family[family_name]),
                          "family_background_count": sum(row["host_category"] == name
                                                          for row in by_family[family_name])})
    exceed = [0] * len(tests)
    for _ in range(repeats):
        labels = {}
        for family_name, family_rows in by_family.items():
            values = [row["host_category"] for row in family_rows]
            rng.shuffle(values)
            labels.update({id(row): value for row, value in zip(family_rows, values)})
        for index, test in enumerate(tests):
            count = sum(labels[id(row)] == test["host_category"] for row in by_query[test["query"]])
            exceed[index] += count >= test["observed_count"]
    pvalues = [(1 + count) / (repeats + 1) for count in exceed]
    for test, pvalue, qvalue in zip(tests, pvalues, bh(pvalues)):
        test["permutation_p"] = pvalue
        test["fdr_q"] = qvalue
        test["observed_fraction"] = test["observed_count"] / test["query_n"]
        test["family_background_fraction"] = test["family_background_count"] / test["family_background_n"]
    return tests


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


def main(input_dir: Path = DEFAULT_INPUT, metadata_path: Path = DEFAULT_METADATA,
         output_dir: Path = DEFAULT_OUTPUT, top_n: int = TOP_HOSTS,
         bootstraps: int = BOOTSTRAPS, permutations: int = PERMUTATIONS,
         seed: int = SEED) -> dict:
    raw: list[dict] = []
    for path in sorted(input_dir.glob("*prediction analysis.txt")):
        raw.extend(parse_blast(path))
    if not raw:
        raise ValueError(f"No BLASTp records found in {input_dir}")
    raw = add_identical_profiles(raw)
    annotate(raw, read_metadata(metadata_path))
    accession_rows, top_rows = score_and_select(raw, top_n)
    rng = random.Random(seed)
    support = bootstrap(top_rows, bootstraps, rng)
    enrichment = permutation(top_rows, permutations, rng)
    output_dir.mkdir(parents=True, exist_ok=True)

    fields = ["query", "manuscript_virus", "candidate_labels", "family", "accession", "reference_virus",
              "source_host", "reported_host", "source_organism", "host_genus", "host_category", "host_filter_reason",
              "threshold_status", "evalue", "identity", "identity_percent", "query_coverage",
              "query_coverage_percent", "bitscore", "raw_evidence_score", "evidence_weight", "host_rank",
              "source_file", "profile_provenance", "metadata_source"]
    write_csv(output_dir / "potential_host_network_edges.csv",
              [{field: row.get(field, "") for field in fields} for row in top_rows])
    write_csv(output_dir / "accepted_accession_evidence.csv",
              [{field: row.get(field, "") for field in fields} for row in accession_rows])
    selected_keys = {(row["query"], row["accession"], row["host_genus"]) for row in top_rows}
    raw_fields = ["query", "manuscript_virus", "candidate_labels", "family", "accession",
                  "reference_virus", "source_host", "reported_host", "source_organism", "host_genus",
                  "host_category", "host_filter_reason", "threshold_status", "evalue",
                  "identity", "identity_percent", "query_coverage", "query_coverage_percent",
                  "bitscore", "raw_evidence_score", "source_file", "profile_provenance", "raw_description"]
    raw_export = []
    for row in raw:
        record = {field: row.get(field, "") for field in raw_fields}
        record["included_in_quantitative_network"] = (
            row.get("query"), row.get("accession"), row.get("host_genus")
        ) in selected_keys
        raw_export.append(record)
    write_csv(output_dir / "raw_blast_host_evidence.csv", raw_export)
    write_csv(output_dir / "host_category_support.csv", support)
    write_csv(output_dir / "host_category_enrichment.csv", enrichment)
    recurrence = []
    for genus in sorted({row["host_genus"] for row in top_rows}):
        genus_rows = [row for row in top_rows if row["host_genus"] == genus]
        recurrence.append({"host_genus": genus, "host_category": category(genus),
                           "profile_count": len({row["query"] for row in genus_rows}),
                           "edge_count": len(genus_rows),
                           "candidate_profiles": "; ".join(sorted({row["query"] for row in genus_rows}))})
    write_csv(output_dir / "host_taxon_recurrence.csv", recurrence)
    manifest = []
    for query in sorted({row["query"] for row in raw}):
        accepted = [row for row in accession_rows if row["query"] == query]
        selected = [row for row in top_rows if row["query"] == query]
        manifest.append({"query": query, "manuscript_virus": NAMES.get(query, query),
                         "family": family(query), "accepted_accession_count": len(accepted),
                         "top_host_count": len(selected),
                         "status": "included" if selected else "no valid source/host after filtering"})
    write_csv(output_dir / "candidate_manifest.csv", manifest)
    threshold_rows = [row for row in raw if row["threshold_status"]]
    excluded = Counter(row["host_filter_reason"] for row in threshold_rows if not row["host_genus"])
    summary = {
        "raw_alignment_count": len(raw), "threshold_hit_count": len(threshold_rows),
        "accepted_accession_count": len(accession_rows),
        "top_host_edge_count": len(top_rows), "unique_host_genera": len(recurrence),
        "query_profiles": sorted({row["query"] for row in raw}),
        "derived_profiles": {target: source for target, source in IDENTICAL_PROFILE_SOURCES.items()
                             if any(row["query"] == target and row["profile_provenance"].startswith("copied_from_")
                                    for row in raw)},
        "excluded_threshold_hit_reasons": dict(excluded),
        "thresholds": {"evalue": EVALUE, "identity": IDENTITY, "query_coverage": COVERAGE,
                        "top_distinct_host_genera_per_query": top_n},
        "statistics": {"weighted_bootstrap_replicates": bootstraps,
                        "within_family_permutations": permutations,
                        "permutation_p_value": "(1 + null_count_ge_observed) / (permutations + 1)",
                        "fdr": "Benjamini-Hochberg across all query-category tests", "seed": seed},
        "evidence_score": "bitscore * identity * query_coverage",
        "interpretation": "ranked host-association hypotheses, not infection probabilities or complete natural host ranges",
    }
    (output_dir / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--top-hosts", type=int, default=TOP_HOSTS)
    parser.add_argument("--bootstraps", type=int, default=BOOTSTRAPS)
    parser.add_argument("--permutations", type=int, default=PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    print(json.dumps(main(args.input_dir, args.metadata, args.output_dir, args.top_hosts,
                          args.bootstraps, args.permutations, args.seed), indent=2))
