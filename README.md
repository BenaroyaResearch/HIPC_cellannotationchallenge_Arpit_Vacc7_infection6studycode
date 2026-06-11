# HIPC Cell-Type → Cell Ontology (CT) consensus mapping

Maps every unique annotation label in your per-cell CSV onto the HIPC **Cell
Ontology tree** (the 39-node `CT_Ontology_Spreadsheet`), then runs the autoAnno
majority-vote consensus. Two voter configurations: `four` (your original script)
and `all` (every available annotation column).

## Files

| File | What it does |
|------|--------------|
| `ct_aliases.py` | Hand-curated `method label → CT node` dictionary + normaliser + fuzzy fallback. **This is the file you extend.** |
| `ct_crosswalk_builder.py` | Reads the ontology, builds `treeLevel1..6` for each CT node, maps every unique label, emits the crosswalk + an unmapped-labels report. |
| `autoanno_consensus.py` | Your consensus pipeline, parameterised for both configs; reads the crosswalk. |
| `starter_all/`, `starter_four/` | Ready-made starter crosswalks built from the alias vocabulary (no CSV needed). Each contains `celltype_mapping_table.csv`, `unmapped_labels_report.csv`, `ct_tree_levels.csv`. |

## Voter configs

```
four : celltypist:Immune_All_Low (ImmLow), AIFI_L2 (AllenL2),
       AIFI_L3 (AllenL3), monaco_immune...labels (SingleR/fine)

all  : the four above + celltypist:Immune_All_High (ImmHigh),
       azimuth_broad/medium/fine (Azimuth l1/l2/l3),
       dice...labels (SingleR/dice), hpca...labels (SingleR/hpca)
```

## Workflow on your cluster

```bash
# 1. Build the crosswalk from the ACTUAL unique labels in your CSV.
#    (Run once per config. This also flags labels that need curation.)
python ct_crosswalk_builder.py \
    --ontology CT_Ontology_Spreadsheet_20260526.xlsx \
    --csv /nfs/.../concatenated_AIFI2_ANN_annotations.csv \
    --config all --outdir crosswalk_all

# 2. Open crosswalk_all/unmapped_labels_report.csv.
#    - rows with a fuzzyGuess: confirm or correct.
#    - rows with empty fuzzyGuess: add the label to the right CT node in ct_aliases.py.
#    Re-run step 1 until 'needs review' is acceptable.

# 3. Run the consensus.
python autoanno_consensus.py \
    --csv /nfs/.../concatenated_AIFI2_ANN_annotations.csv \
    --mapping crosswalk_all/celltype_mapping_table.csv \
    --config all --barcode-col auto --outdir consensus_all
```

Repeat with `--config four` for the 4-voter version.

## Outputs

- `cell_labels_<config>.csv` — `cellBarcode, suggestedCelltype, suggestedCelltypeLevel`
- `consensus_wide_<config>.csv` — per-level consensus label, agreement score, confidence category

## How the mapping works

Each method label is resolved to a CT node by (a) exact alias lookup in
`ct_aliases.ALIASES`, else (b) a conservative fuzzy match (token-set + difflib,
cutoff 0.84). The CT node is expanded to its full root→node path
(`Blood Cell > Leukocyte > … > node`) which becomes `treeLevel1..6`. The
consensus then majority-votes per tree level (level 1 dropped), assigns
`high` (>0.7), `medium` (≥0.5), or `low` confidence, and reports the deepest
high-confidence level as `suggestedCelltype`. Broad ancestors (Blood Cell,
Leukocyte, Lymphoid Cell, Myeloid Cell) are downgraded to `unresolved`, exactly
like autoAnno.

## Notes / caveats

- **The starter crosswalks cover common Monaco / Azimuth / CellTypist / AIFI /
  DICE / HPCA vocabularies.** Your real CSV will contain labels not yet aliased
  — that's what the unmapped report is for. Curating it is the one manual step.
- The join is normalisation-robust (case / punctuation / whitespace), so you
  don't need exact casing when adding aliases.
- `--barcode-col auto` tries `cell_barcode`, `cellBarcode`, `barcode`,
  `cell_id`, `index`, `Unnamed: 0`, then falls back to the row index. Pass an
  explicit column if your barcodes live elsewhere.
- Fuzzy matches are advisory. Anything below the cutoff is left unmapped rather
  than guessed, to avoid silently miscalling cells.
