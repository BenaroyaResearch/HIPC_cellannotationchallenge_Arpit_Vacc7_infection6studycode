# HIPC Cell-Type Annotation Challenge — Consensus Pipeline

Multi-method ensemble cell-type annotation for the HIPC single-cell RNA-seq
annotation challenge, applied to two PBMC datasets:

- **vaccination_study_07** — 53,619 cells
- **infection_study_06** — 827,389 cells

Final submission uses the **fine6** 6-voter consensus configuration.
See `annotation.readme.docx` (in the submission ZIP) for the full methodology.

---

## Repository contents

| File | Purpose |
|---|---|
| `aifi_annotation2.py` | CellTypist pass with AIFI PBMC models — adds `AIFI_L2` / `AIFI_L3` to the scDownstream h5ad |
| `panhuman_azimuthannotatorrun.sh` | Batch Azimuth (panhumanpy) annotation across all study folders + CSV export |
| `ct_aliases.py` | Curated alias table: maps each tool's raw label vocabulary onto CT ontology nodes |
| `ct_crosswalk_builder.py` | Reads the CT ontology spreadsheet + alias table → per-method label → tree-level crosswalk |
| `autoanno_consensus.py` | Per-cell majority-vote consensus at each CT tree level; outputs final labels + confidence |
| `run_ct_consensus.sh` | Wrapper: runs crosswalk builder + consensus for one or more studies and configs |
| `params_vaccine7_annotation_full.yaml` | nf-core/scdownstream params for vaccination_study_07 |
| `params_infection6_annotation_full.yaml` | nf-core/scdownstream params for infection_study_06 |
| `samplesheet_vaccine7_annotation.csv` | scDownstream input samplesheet for vaccination_study_07 |
| `samplesheet_infection6_full.csv` | scDownstream input samplesheet for infection_study_06 (per-sample split) |
| `celldex_references.csv` | SingleR celldex reference list (monaco_immune, dice, hpca) |
| `CT_Ontology_Spreadsheet_20260526.xlsx` | HIPC Cell Ontology hierarchy (treeLevel1–6) |

---

## Pipeline overview

```
raw h5ad (per study)
       |
       v
[1] nf-core/scdownstream  (qc_only: true — annotation only, no filtering)
    CellTypist: Immune_All_Low, Immune_All_High
    SingleR:    monaco_immune, dice, hpca
       |
       v
[2] aifi_annotation2.py
    CellTypist AIFI PBMC models -> adds AIFI_L2, AIFI_L3
       |
       v
[3] panhuman_azimuthannotatorrun.sh
    panhumanpy annotate -> adds azimuth_broad, azimuth_medium, azimuth_fine
       |
       v
[4] run_ct_consensus.sh  (fine6 config)
    ct_crosswalk_builder.py  -> label x CT tree crosswalk
    autoanno_consensus.py    -> majority vote per tree level -> final labels
```

**Why no filtering?** Ambient correction, doublet detection, and integration are
all disabled so that every annotator sees the same set of cells — clean merge,
no reindexing needed.

---

## Voter configuration: fine6

| Voter | Tool | Reference |
|---|---|---|
| `celltypist_ImmLow` | CellTypist 1.6.3 | Immune_All_Low |
| `celltypist_AllenL2` | CellTypist 1.6.3 | AIFI PBMC L2 (2024-04-19) |
| `celltypist_AllenL3` | CellTypist 1.6.3 | AIFI PBMC L3 (2024-04-19) |
| `singleR_fine` | SingleR 2.12.0 | Monaco Immune (celldex 1.20.0) |
| `azimuth_l3` | panhumanpy 0.5.0 | AIFI pan-human M0.2 |
| `singleR_dice` | SingleR 2.12.0 | DICE (celldex 1.20.0) |

A cell is assigned the deepest CT tree level at which > 70% of voters agree
(`high` confidence). Cells with 50–70% agreement are assigned with `medium`
confidence. Cells with < 50% agreement fall back to their best-supported lineage
(e.g., Lymphoid Cell) — no cell is left unresolved.

Other available configs in `autoanno_consensus.py`:

| Config | Voters | Notes |
|---|---|---|
| `four` | 4 (base CellTypist + Monaco) | Original baseline |
| `all` | 10 (all annotators including coarse) | Broad coverage |
| `fine6` | 6 (fine-resolution only) | **Submission config** |
| `fine7` | 7 (fine6 + Seurat V5 PBMC 2023) | Experimental |

---

## Running the consensus

Requires the per-cell annotation CSVs produced by stages 1–3 and the CT
ontology spreadsheet. Update the `STUDIES` paths in `run_ct_consensus.sh`
to point at your CSVs, then:

```bash
# Run fine6 consensus for both datasets
BASE=$(pwd) ./run_ct_consensus.sh fine6 vaccine7
BASE=$(pwd) ./run_ct_consensus.sh fine6 infection6
```

Outputs go to `ct_consensus_runs/`:

| File | Contents |
|---|---|
| `crosswalk_fine6_<study>/celltype_mapping_table.csv` | Per-method label → CT tree-level crosswalk |
| `crosswalk_fine6_<study>/unmapped_labels_report.csv` | Labels that need curation in `ct_aliases.py` |
| `consensus_fine6_<study>/cell_labels_fine6.csv` | Final per-cell labels + confidence |
| `consensus_fine6_<study>/consensus_wide_fine6.csv` | Full wide table: all tree levels, scores, categories |

---

## Requirements

| Tool | Version |
|---|---|
| Nextflow | ≥ 26 |
| nf-core/scdownstream | v0.0.1dev-gbcb4c67 |
| CellTypist | 1.6.3 |
| SingleR | 2.12.0 |
| celldex | 1.20.0 |
| panhumanpy (Azimuth) | 0.5.0 |
| Python | ≥ 3.10 |
| pandas / numpy / anndata / openpyxl | 2.3 / 1.26 / 0.8 / 3.1 |
| R | 4.5.3 |

Python env: `pip install celltypist scanpy anndata pandas numpy openpyxl`

Azimuth env: `pixi` with `panhumanpy` — see `panhuman_azimuthannotatorrun.sh`
for `PIXI_MANIFEST` setup.

---

## Label curation

`ct_aliases.py` is the single file to extend when a new label vocabulary is
encountered. After adding aliases, re-run `run_ct_consensus.sh` to regenerate
the crosswalk — the unmapped labels report will confirm zero residual gaps.
