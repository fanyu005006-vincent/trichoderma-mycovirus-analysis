# Proposed upload manifest

## Included

- `run_analysis.py`: one command entry point.
- `scripts/host_network.py`: BLASTp parsing, NCBI host filtering, weighted
  evidence ranking, bootstrap support and permutation enrichment.
- `scripts/structural_support.py`: AF3 and TM-align secondary support layer.
- `data/host_prediction_inputs/`: saved BLASTp inputs.
- `data/reference_host_metadata.csv`: accession to NCBI `source/host` mapping.
  It also records explicit literature based host overrides for Vitis vinifera
  and Yellow silver pine without changing the original source host field.
- `data/af3_structure_manifest.csv` and `data/af3_structures/`: supplied
  structural inputs.
- `README.md`, `CODE_AVAILABILITY.md`, `requirements.txt` and tests.

## Deliberately excluded

- The old `virome_workflow/` raw read and Snakemake tree.
- The old `filter_host_hits.py` top five host logic.
- The old full host network, figure redraw and Sankey scripts.
- `data/rdrp_candidates.fasta`, which is not the canonical input to the saved
  BLASTp host analysis and contains a different profile count.
- Local generated results, caches and machine specific files.
