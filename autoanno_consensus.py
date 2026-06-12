#!/usr/bin/env python3
"""
autoanno_consensus.py
=====================
Python port of the autoAnno.Rmd consensus pipeline, run directly on the
consolidated per-cell annotation CSV (concatenated_*_ANN_annotations.csv).

This version supports TWO voter configurations and reads the crosswalk that
ct_crosswalk_builder.py produces (celltype_mapping_table.csv with treeLevel1..6).

  --config four : the original 4 voters
        celltypist_ImmLow  <- celltypist:Immune_All_Low      (CellTypist/ImmLow)
        celltypist_AllenL2 <- AIFI_L2                         (CellTypist/AllenL2)
        celltypist_AllenL3 <- AIFI_L3                         (CellTypist/AllenL3)
        singleR_fine       <- monaco_immune.tar.labels...     (SingleR/fine)

  --config all  : the 4 above PLUS
        celltypist_ImmHigh <- celltypist:Immune_All_High      (CellTypist/ImmHigh)
        azimuth_l1/l2/l3   <- azimuth_broad/medium/fine       (Azimuth/l1,l2,l3)
        singleR_dice       <- dice.tar.labels_dice            (SingleR/dice)
        singleR_hpca       <- hpca.tar.labels_hpca            (SingleR/hpca)

  --config fine6 : four + fine-depth extras only (drops azimuth_l1/l2, ImmHigh, hpca)
        celltypist_ImmLow  <- celltypist:Immune_All_Low      (CellTypist/ImmLow)
        celltypist_AllenL2 <- AIFI_L2                         (CellTypist/AllenL2)
        celltypist_AllenL3 <- AIFI_L3                         (CellTypist/AllenL3)
        singleR_fine       <- monaco_immune.tar.labels...     (SingleR/fine)
        azimuth_l3         <- azimuth_fine                    (Azimuth/l3)
        singleR_dice       <- dice.tar.labels_dice            (SingleR/dice)

Reproduces the autoAnno chunks: mapToTree + bind, majority-vote consensus per
tree level, confidenceCategory, suggestedCelltype (deepest high-confidence
level), broad-ancestor downgrade, and CSV output.

USAGE
-----
  python autoanno_consensus.py \
      --csv /nfs/.../concatenated_AIFI2_ANN_annotations.csv \
      --mapping ./crosswalk_out/celltype_mapping_table.csv \
      --config all \
      --barcode-col auto \
      --outdir ./consensus_out
"""

import argparse
import os
import re
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Voter configs:  internal_voter_name -> (csv_column, method, methodLevel)
# The (method, methodLevel) tuple MUST match values in the crosswalk.
# ---------------------------------------------------------------------------
VOTERS_FOUR = {
    "celltypist_ImmLow":  ("celltypist:Immune_All_Low",              "CellTypist", "ImmLow"),
    "celltypist_AllenL2": ("AIFI_L2",                                "CellTypist", "AllenL2"),
    "celltypist_AllenL3": ("AIFI_L3",                                "CellTypist", "AllenL3"),
    "singleR_fine":       ("monaco_immune.tar.labels_monaco_immune", "SingleR",    "fine"),
}
VOTERS_ALL = {
    **VOTERS_FOUR,
    "celltypist_ImmHigh": ("celltypist:Immune_All_High", "CellTypist", "ImmHigh"),
    "azimuth_l1":         ("azimuth_broad",              "Azimuth",    "l1"),
    "azimuth_l2":         ("azimuth_medium",             "Azimuth",    "l2"),
    "azimuth_l3":         ("azimuth_fine",               "Azimuth",    "l3"),
    "singleR_dice":       ("dice.tar.labels_dice",       "SingleR",    "dice"),
    "singleR_hpca":       ("hpca.tar.labels_hpca",       "SingleR",    "hpca"),
}
VOTERS_FINE6 = {
    **VOTERS_FOUR,
    "azimuth_l3":   ("azimuth_fine",         "Azimuth", "l3"),
    "singleR_dice": ("dice.tar.labels_dice", "SingleR", "dice"),
}
VOTERS_FINE7 = {
    **VOTERS_FINE6,
    "sV5_pbmc":     ("sV5celltypeL2",        "Seurat",  "pbmc2023"),
}

BROAD_LABELS = {"Blood Cell", "Leukocyte", "Lymphoid Cell", "Myeloid Cell"}
# Nodes that are terminal at treeLevel2 (no children in the CT ontology).
# Deeper tree levels are meaningless for these — voters that don't know the
# cell type go silent at level 3+ and a single wrong voter creates spurious
# high-confidence calls. Accept the treeLevel2 majority vote directly.
TERMINAL_AT_L2 = {"Platelet", "HSC", "RBC", "Doublet"}
TREE_LEVELS = [2, 3, 4, 5, 6]                      # autoAnno drops level 1
TREE_COLS   = [f"treeLevel{i}" for i in range(1, 7)]


def _norm(s):
    """Match ct_aliases._norm so the crosswalk join is case/punctuation robust."""
    if s is None:
        return ""
    s = str(s).strip().lower()
    s = re.sub(r"[._/\\,\-+()]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def resolve_barcode(anno, barcode_col):
    """Return a Series of cell barcodes. 'auto' tries common names then index."""
    if barcode_col and barcode_col != "auto":
        if barcode_col not in anno.columns:
            raise ValueError(f"--barcode-col '{barcode_col}' not in CSV")
        return anno[barcode_col].astype(str).values
    for cand in ("cell_barcode", "cellBarcode", "barcode", "cell_id",
                 "CellID", "index", "Unnamed: 0"):
        if cand in anno.columns:
            print(f"[info] using '{cand}' as cell barcode")
            return anno[cand].astype(str).values
    print("[info] no barcode column found; using row index as barcode")
    return anno.index.astype(str).values


def map_to_tree(labels, mapping, method, method_level):
    """left_join one method's per-cell labels to the tree crosswalk.

    Join is done on a normalised label key so differences in case / punctuation
    / whitespace between the CSV and the crosswalk don't drop cells.
    """
    key = (mapping[(mapping["method"] == method) &
                   (mapping["methodLevel"] == method_level)]
           [["normLabel", *TREE_COLS]].drop_duplicates("normLabel"))
    left = pd.DataFrame({"normLabel": pd.Series(labels).map(_norm).values})
    return left.merge(key, on="normLabel", how="left").drop(columns="normLabel")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--mapping", required=True, help="celltype_mapping_table.csv")
    ap.add_argument("--config", choices=["four", "all", "fine6", "fine7"], default="all")
    ap.add_argument("--barcode-col", default="auto")
    ap.add_argument("--outdir", default="consensus_out")
    ap.add_argument("--fallback-to-lineage", action="store_true",
                    help="Assign unresolved cells to their best medium-confidence or "
                         "lineage-root label instead of 'unresolved'. Adds a "
                         "suggestedCelltypeConfidence column (high/medium/lineage_*).")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    voters = {"four": VOTERS_FOUR, "all": VOTERS_ALL,
              "fine6": VOTERS_FINE6, "fine7": VOTERS_FINE7}[args.config]

    anno = pd.read_csv(args.csv, dtype=str)
    mapping = pd.read_csv(args.mapping, dtype=str)
    for c in TREE_COLS:
        mapping[c] = mapping[c].fillna("").astype(str).str.strip()
    mapping["methodLabel"] = mapping["methodLabel"].astype(str).str.strip()
    if "normLabel" not in mapping.columns:        # tolerate hand-edited tables
        mapping["normLabel"] = mapping["methodLabel"].map(_norm)

    required = {"method", "methodLevel", "methodLabel", *TREE_COLS}
    missing = required - set(mapping.columns)
    if missing:
        raise ValueError(f"mapping table missing columns: {sorted(missing)}")

    barcodes = resolve_barcode(anno, args.barcode_col)
    annoMeta = pd.DataFrame({"cellBarcode": barcodes})

    # ---- build mapped tree-level blocks per available voter ----
    mapped_blocks = []
    used = []
    for vname, (csv_col, method, level) in voters.items():
        if csv_col not in anno.columns:
            print(f"[warn] column '{csv_col}' absent; voter '{vname}' skipped")
            continue
        used.append(vname)
        annoMeta[vname] = anno[csv_col].values
        block = map_to_tree(anno[csv_col], mapping, method, level)
        block.columns = [f"{vname}_{c}" for c in block.columns]
        mapped_blocks.append(block.reset_index(drop=True))
    print(f"[info] active voters ({len(used)}): {used}")

    annoMapped = pd.concat([annoMeta.reset_index(drop=True), *mapped_blocks], axis=1)

    # ---- naCheck (unmapped labels show up as NA at treeLevel4) ----
    print("\n=== naCheck: NA counts per *_treeLevel4 column (nonzero = unmapped) ===")
    for c in [c for c in annoMapped.columns if c.endswith("_treeLevel4")]:
        print(f"  {c}: {int(annoMapped[c].isna().sum())}")

    # ---- long + majority vote per (cell, level) ----
    value_cols = [c for c in annoMapped.columns if re.search(r"_treeLevel[2-6]$", c)]
    long = annoMapped.melt(id_vars="cellBarcode", value_vars=value_cols,
                           var_name="col", value_name="treeLabel")
    long[["methodLevel", "treeLevel"]] = long["col"].str.extract(r"(.+)_(treeLevel\d)")
    valid = long[long["treeLabel"].notna() & (long["treeLabel"] != "")].copy()

    counts = (valid.groupby(["cellBarcode", "treeLevel", "treeLabel"])
              .size().reset_index(name="n"))
    winners = counts.loc[counts.groupby(["cellBarcode", "treeLevel"])["n"].idxmax()]
    winners = winners.rename(columns={"treeLabel": "consensusLabel", "n": "nAgree"})
    totals = (counts.groupby(["cellBarcode", "treeLevel"])["n"].sum()
              .reset_index(name="nMethods"))
    cons = winners.merge(totals, on=["cellBarcode", "treeLevel"])
    cons["agreementScore"] = cons["nAgree"] / cons["nMethods"]

    cons["confidenceCategory"] = np.select(
        [cons["agreementScore"] > 0.7, cons["agreementScore"] >= 0.5],
        ["high", "medium"], default="low")

    print("\n=== cells per (treeLevel, confidenceCategory) ===")
    print(pd.crosstab(cons["treeLevel"], cons["confidenceCategory"]))

    # ---- pivot wide + suggestedCelltype ----
    wide = cons.pivot(index="cellBarcode", columns="treeLevel",
                      values=["consensusLabel", "agreementScore", "confidenceCategory"])
    wide.columns = [f"{val}_{lvl}" for val, lvl in wide.columns]
    wide = wide.reset_index()

    fallback = args.fallback_to_lineage

    def suggest(row):
        # Pre-check: terminal-at-level-2 nodes (Platelet, HSC, RBC).
        # Deeper tree levels are invalid for these — voters that don't recognise
        # the cell type go silent at level 3+ and a single wrong voter creates a
        # spurious high-confidence call.  Accept treeLevel2 if majority agree.
        l2_label = row.get("consensusLabel_treeLevel2", "")
        l2_cat   = row.get("confidenceCategory_treeLevel2", "")
        if l2_label in TERMINAL_AT_L2 and l2_cat in ("high", "medium"):
            return pd.Series([l2_label, "treeLevel2",
                              "high" if l2_cat == "high" else "medium"])

        # Step 1: deepest HIGH-confidence, non-broad label (primary call)
        for lvl in sorted(TREE_LEVELS, reverse=True):
            label = row.get(f"consensusLabel_treeLevel{lvl}")
            if (row.get(f"confidenceCategory_treeLevel{lvl}") == "high"
                    and label and label not in BROAD_LABELS):
                return pd.Series([label, f"treeLevel{lvl}", "high"])

        if not fallback:
            return pd.Series(["unresolved", "unresolved", "unresolved"])

        # Step 2: deepest MEDIUM-confidence, non-broad label
        # Correct at the lineage level; less confident at subtype resolution.
        for lvl in sorted(TREE_LEVELS, reverse=True):
            label = row.get(f"consensusLabel_treeLevel{lvl}")
            if (row.get(f"confidenceCategory_treeLevel{lvl}") == "medium"
                    and label and label not in BROAD_LABELS):
                return pd.Series([label, f"treeLevel{lvl}", "medium"])

        # Step 3: deepest HIGH or MEDIUM at any level including broad
        # (Lymphoid Cell / Myeloid Cell — always correct, just broad)
        for lvl in sorted(TREE_LEVELS, reverse=True):
            label = row.get(f"consensusLabel_treeLevel{lvl}")
            cat   = row.get(f"confidenceCategory_treeLevel{lvl}")
            if cat in ("high", "medium") and label:
                return pd.Series([label, f"treeLevel{lvl}", f"lineage_{cat}"])

        # Step 4: treeLevel2 raw consensus label — Blood Cell / Lymphoid / Myeloid
        # The majority vote at level 2 is always correct even at low confidence.
        label = row.get("consensusLabel_treeLevel2", "")
        if label:
            return pd.Series([label, "treeLevel2", "lineage_low"])

        return pd.Series(["unresolved", "unresolved", "unresolved"])

    wide[["suggestedCelltype", "suggestedCelltypeLevel",
          "suggestedCelltypeConfidence"]] = wide.apply(suggest, axis=1)

    if not fallback:
        # original behaviour: broad high-confidence calls → unresolved
        broad = wide["suggestedCelltype"].isin(BROAD_LABELS)
        wide.loc[broad, "suggestedCelltype"]           = "unresolved"
        wide.loc[broad, "suggestedCelltypeLevel"]      = "unresolved"
        wide.loc[broad, "suggestedCelltypeConfidence"] = "unresolved"

    print("\n=== suggestedCelltype counts ===")
    print(wide["suggestedCelltype"].value_counts(dropna=False))

    out_wide = os.path.join(args.outdir, f"consensus_wide_{args.config}.csv")
    out_labels = os.path.join(args.outdir, f"cell_labels_{args.config}.csv")
    wide.to_csv(out_wide, index=False)
    wide[["cellBarcode", "suggestedCelltype", "suggestedCelltypeLevel",
          "suggestedCelltypeConfidence"]].to_csv(out_labels, index=False)
    print(f"\nwrote {out_labels}  ({len(wide)} cells)")
    print(f"wrote {out_wide}")


if __name__ == "__main__":
    main()
