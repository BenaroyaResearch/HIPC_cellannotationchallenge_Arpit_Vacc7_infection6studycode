#!/usr/bin/env python3
"""
ct_crosswalk_builder.py
=======================
Build the `celltype_mapping_table.csv` crosswalk that autoanno_consensus.py
consumes -- but generated automatically for *all unique labels* in your
annotation CSV, mapped onto the HIPC Cell Ontology (CT) tree.

PIPELINE
--------
1. Load the CT ontology spreadsheet (Celltype, cl:ID, Parent Class,
   isTerminalNode, Definitions). Build the parent tree and compute, for every
   CT node, its full root->node ancestor path = treeLevel1..treeLevel6.
2. For each annotation method column, collect the UNIQUE labels (either from
   your CSV, or from a supplied list). Map each label -> a CT node using
   ct_aliases.exact_match, falling back to ct_aliases.fuzzy_match.
3. Emit:
     - celltype_mapping_table.csv   (method, methodLevel, methodLabel,
                                      treeLevel1..6, ctNode, clId, mapSource,
                                      mapScore)   <- feeds the consensus script
     - unmapped_labels_report.csv   (labels needing human curation, with the
                                      fuzzy best-guess + score)
     - ct_tree_levels.csv           (the flattened ontology, for reference)

USAGE
-----
  # run on your cluster where the CSV lives, all voters:
  python ct_crosswalk_builder.py \
      --ontology CT_Ontology_Spreadsheet_20260526.xlsx \
      --csv /nfs/.../concatenated_AIFI2_ANN_annotations.csv \
      --config all \
      --outdir ./crosswalk_out

  # or just the 4 voters in the original script:
  python ct_crosswalk_builder.py --ontology ... --csv ... --config four

  # no CSV handy? emit a starter crosswalk from the alias vocabularies only:
  python ct_crosswalk_builder.py --ontology ... --config all --starter

CONFIGS
-------
  four : celltypist:Immune_All_Low, AIFI_L2, AIFI_L3, monaco(fine)
  all  : the four above + Immune_All_High, azimuth broad/medium/fine,
         dice, hpca
"""

import argparse
import os
import sys
import pandas as pd

import ct_aliases as A

# ---------------------------------------------------------------------------
# Column  ->  (method, methodLevel) as used in the crosswalk + consensus script
# ---------------------------------------------------------------------------
COLUMN_METHODS_FOUR = {
    "celltypist:Immune_All_Low":              ("CellTypist", "ImmLow"),
    "AIFI_L2":                                ("CellTypist", "AllenL2"),
    "AIFI_L3":                                ("CellTypist", "AllenL3"),
    "monaco_immune.tar.labels_monaco_immune": ("SingleR",    "fine"),
}

COLUMN_METHODS_ALL = {
    **COLUMN_METHODS_FOUR,
    "celltypist:Immune_All_High":             ("CellTypist", "ImmHigh"),
    "azimuth_broad":                          ("Azimuth",    "l1"),
    "azimuth_medium":                         ("Azimuth",    "l2"),
    "azimuth_fine":                           ("Azimuth",    "l3"),
    "dice.tar.labels_dice":                   ("SingleR",    "dice"),
    "hpca.tar.labels_hpca":                   ("SingleR",    "hpca"),
}

TREE_COLS = [f"treeLevel{i}" for i in range(1, 7)]
NA_TOKENS = {"", "nan", "none", "na", "unknown", "unassigned", "unclassified",
             "unresolved", "doublet?", "low_quality", "lowquality",
             # non-lineage / non-cell-type tokens that must NOT be forced onto
             # the immune tree (booleans, tissue codes, lineage-free calls)
             "false", "true", "bm", "pb", "cb", "cycling cells",
             "cycling cell", "cells"}


# ---------------------------------------------------------------------------
# 1. Ontology -> tree levels
# ---------------------------------------------------------------------------
def load_ontology(path):
    """Return (df, parent_map, clid_map, levels_map, ct_names)."""
    if path.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(path, dtype=str)
    else:
        df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]

    # tolerate slight header variants
    def col(*cands):
        for c in cands:
            for actual in df.columns:
                if actual.lower().replace(" ", "") == c.lower().replace(" ", ""):
                    return actual
        raise KeyError(f"ontology missing column among {cands}")

    c_name = col("Celltype", "cell type", "celltype")
    c_id   = col("cl: ID", "cl:ID", "clID", "cl id", "id")
    c_par  = col("Parent Class", "parent", "parentclass")

    df = df[df[c_name].notna() & (df[c_name].astype(str).str.strip() != "")].copy()
    df[c_name] = df[c_name].astype(str).str.strip()
    df[c_par] = df[c_par].astype(str).str.strip()

    parent_map, clid_map = {}, {}
    for _, r in df.iterrows():
        nm = r[c_name]
        par = r[c_par]
        parent_map[nm] = None if par in ("None", "", "nan") else par
        clid_map[nm] = (None if str(r[c_id]).strip() in ("None", "", "nan")
                        else str(r[c_id]).strip())

    # compute root->node path for each node
    def path_to_root(nm):
        chain, seen = [], set()
        cur = nm
        while cur is not None and cur not in seen:
            seen.add(cur)
            chain.append(cur)
            cur = parent_map.get(cur)
        return list(reversed(chain))      # root first

    levels_map = {}
    for nm in parent_map:
        chain = path_to_root(nm)
        row = (chain + [""] * 6)[:6]      # pad/truncate to 6 levels
        if len(chain) > 6:
            sys.stderr.write(f"[warn] '{nm}' deeper than 6 levels; truncated\n")
        levels_map[nm] = row

    return df, parent_map, clid_map, levels_map, list(parent_map.keys())


# ---------------------------------------------------------------------------
# 2. Collect unique labels per column
# ---------------------------------------------------------------------------
def unique_labels_from_csv(csv_path, columns):
    """Stream the CSV in chunks; return {column: set(labels)} for present cols."""
    present = None
    out = {c: set() for c in columns}
    reader = pd.read_csv(csv_path, dtype=str, usecols=lambda c: True,
                         chunksize=200_000)
    for chunk in reader:
        if present is None:
            present = [c for c in columns if c in chunk.columns]
            missing = [c for c in columns if c not in chunk.columns]
            if missing:
                sys.stderr.write(f"[warn] columns not in CSV (skipped): {missing}\n")
        for c in present:
            out[c].update(chunk[c].dropna().astype(str).str.strip().unique())
    return {c: v for c, v in out.items() if c in (present or [])}


def labels_from_aliases(columns_methods):
    """Starter mode: use every alias string as a 'unique label' per column."""
    # map method/level -> the alias vocab is method-agnostic, so we just feed
    # all alias strings to every column. The mapping is identical regardless.
    all_labels = set()
    for ct, labs in A.ALIASES.items():
        all_labels.add(ct)
        all_labels.update(labs)
    return {col: set(all_labels) for col in columns_methods}


# ---------------------------------------------------------------------------
# 3. Map labels -> CT node -> tree levels
# ---------------------------------------------------------------------------
def build_rows(col_to_labels, columns_methods, levels_map, clid_map, ct_names,
               fuzzy_cutoff):
    rows, unmapped = [], []
    for col, labels in col_to_labels.items():
        method, level = columns_methods[col]
        for lab in sorted(labels):
            if A._norm(lab) in NA_TOKENS:
                continue
            ct = A.exact_match(lab)
            src, score = "alias", 1.0
            if ct is None:
                ct, score = A.fuzzy_match(lab, ct_names, cutoff=fuzzy_cutoff)
                src = "fuzzy" if ct else "unmapped"
            if ct is None:
                unmapped.append({"column": col, "method": method,
                                 "methodLevel": level, "methodLabel": lab,
                                 "fuzzyGuess": "", "fuzzyScore": score})
                continue
            lv = levels_map.get(ct, [""] * 6)
            row = {"method": method, "methodLevel": level, "methodLabel": lab,
                   "normLabel": A._norm(lab),
                   **{TREE_COLS[i]: lv[i] for i in range(6)},
                   "ctNode": ct, "clId": clid_map.get(ct, ""),
                   "mapSource": src, "mapScore": score}
            rows.append(row)
            if src == "fuzzy":
                unmapped.append({"column": col, "method": method,
                                 "methodLevel": level, "methodLabel": lab,
                                 "fuzzyGuess": ct, "fuzzyScore": score})
    cross = pd.DataFrame(rows).drop_duplicates(
        subset=["method", "methodLevel", "methodLabel"])
    unmapped_df = pd.DataFrame(unmapped)
    return cross, unmapped_df


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ontology", required=True, help="CT ontology .xlsx/.csv")
    ap.add_argument("--csv", help="annotation CSV to pull unique labels from")
    ap.add_argument("--config", choices=["four", "all"], default="all")
    ap.add_argument("--outdir", default="crosswalk_out")
    ap.add_argument("--starter", action="store_true",
                    help="ignore --csv; build crosswalk from alias vocab only")
    ap.add_argument("--fuzzy-cutoff", type=float, default=0.84)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    columns_methods = (COLUMN_METHODS_FOUR if args.config == "four"
                       else COLUMN_METHODS_ALL)

    df, parent_map, clid_map, levels_map, ct_names = load_ontology(args.ontology)

    # reference: flattened tree
    tl = pd.DataFrame([{"ctNode": nm, "clId": clid_map.get(nm, ""),
                        **{TREE_COLS[i]: levels_map[nm][i] for i in range(6)}}
                       for nm in ct_names])
    tl_path = os.path.join(args.outdir, "ct_tree_levels.csv")
    tl.to_csv(tl_path, index=False)

    # gather labels
    if args.starter or not args.csv:
        if not args.starter and not args.csv:
            sys.stderr.write("[info] no --csv given; running in --starter mode\n")
        col_to_labels = labels_from_aliases(columns_methods)
        mode = "starter (alias vocabulary)"
    else:
        col_to_labels = unique_labels_from_csv(args.csv, list(columns_methods))
        mode = f"CSV unique labels ({args.csv})"

    cross, unmapped = build_rows(col_to_labels, columns_methods, levels_map,
                                 clid_map, ct_names, args.fuzzy_cutoff)

    cw_path = os.path.join(args.outdir, "celltype_mapping_table.csv")
    um_path = os.path.join(args.outdir, "unmapped_labels_report.csv")
    cross.to_csv(cw_path, index=False)
    unmapped.to_csv(um_path, index=False)

    # ---- report ----
    print(f"config        : {args.config}")
    print(f"label source  : {mode}")
    print(f"ontology nodes : {len(ct_names)}")
    print(f"crosswalk rows : {len(cross)}  -> {cw_path}")
    if len(cross):
        print("  by mapSource :", cross["mapSource"].value_counts().to_dict())
        print("  by method    :",
              cross.groupby(['method', 'methodLevel']).size().to_dict())
    n_unmapped = int((unmapped["fuzzyGuess"] == "").sum()) if len(unmapped) else 0
    n_fuzzy = int((unmapped["fuzzyGuess"] != "").sum()) if len(unmapped) else 0
    print(f"needs review   : {len(unmapped)}  -> {um_path}")
    print(f"   fuzzy (confirm): {n_fuzzy}   |   unmapped (fill in): {n_unmapped}")
    print(f"tree reference : {tl_path}")
    print("\nNext: review unmapped_labels_report.csv, add fixes to ct_aliases.py,"
          " re-run, then feed celltype_mapping_table.csv to autoanno_consensus.py")


if __name__ == "__main__":
    main()
