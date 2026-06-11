# HIPC_cellannotationchallenge_Arpit_Vacc7_infection6studycode
This pipeline was used  for the HIPC cell-annotation challenge and run over the infection 6and vaccination 7 study cohorts.

## Annotation workflow

The required end-to-end annotation workflow is:

1. Run `scDownstream` on raw/quantified count `h5ad` inputs with `qc_only: true` (no filtering, no ambient correction, no doublet removal, no integration).
2. Generate first-pass labels from:
   - CellTypist (`Immune_All_Low`, `Immune_All_High`)
   - SingleR/celldex (`monaco_immune`, `dice`, `hpca`)
   and produce `reports/artifacts/concatenated.h5ad`.
3. Run AIFI CellTypist models (`aifi_annotation2.py`) to add `AIFI_L2` and `AIFI_L3`, producing `concatenated_AIFI2.h5ad`.
4. Run Azimuth (`panhumanpy ANNotate`) to add `azimuth_broad`, `azimuth_medium`, and `azimuth_fine`, producing `concatenated_AIFI2_ANN.h5ad`.
5. Export `obs` annotations to `concatenated_AIFI2_ANN_annotations.csv`.
6. Build the crosswalk from method labels to Cell Ontology nodes using `ct_aliases.py`, `ct_crosswalk_builder.py`, and the CT spreadsheet, producing `celltype_mapping_table.csv`.
7. Run majority voting per ontology tree level (`autoanno_consensus.py`) to compute consensus label, `agreementScore`, and confidence per level, then emit:
   - `suggestedCelltype` (deepest high-confidence level)
   - `cell_labels_<config>.csv`
   - `consensus_wide_<config>.csv`
