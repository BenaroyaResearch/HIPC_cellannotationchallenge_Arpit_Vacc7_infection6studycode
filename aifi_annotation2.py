#!/usr/bin/env python3
"""
aifi_annotation2.py
===================
Run a second CellTypist pass using the AIFI PBMC reference models (L2 / L3)
on the scDownstream-annotated h5ad.  Adds AIFI_L2 and AIFI_L3 columns to .obs
and writes a new h5ad (*concatenated_AIFI2.h5ad*) beside the input.

Usage
-----
  python aifi_annotation2.py \
      --input  /path/to/concatenated.h5ad \
      --refs   /path/to/reference_models_dir

The reference_models_dir must contain:
  Immune_All_Low.pkl
  ref_pbmc_clean_celltypist_model_AIFI_L2_2024-04-19.pkl
  ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19.pkl
"""

import argparse
import gc
import os
import anndata as ad
import scanpy as sc
import celltypist

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Path to concatenated.h5ad")
    ap.add_argument("--refs",  required=True, help="Directory containing the .pkl model files")
    args = ap.parse_args()

    models = {
        "Immunelow": os.path.join(args.refs, "Immune_All_Low.pkl"),
        "AIFI_L2":   os.path.join(args.refs, "ref_pbmc_clean_celltypist_model_AIFI_L2_2024-04-19.pkl"),
        "AIFI_L3":   os.path.join(args.refs, "ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19.pkl"),
    }
    for name, path in models.items():
        if not os.path.exists(path):
            raise FileNotFoundError(f"Model not found: {path}  (key={name})")

    out = args.input.replace("concatenated.h5ad", "concatenated_AIFI2.h5ad")
    if out == args.input:
        raise ValueError("Input filename must contain 'concatenated.h5ad'")

    print("Loading h5ad...")
    a = ad.read_h5ad(args.input)
    print(f"Loaded: {a.shape}  X_max={float(a.X.max()):.2f}  layers={list(a.layers.keys())}")

    src = a.layers["counts"].copy() if "counts" in a.layers else a.X.copy()
    tmp = ad.AnnData(X=src, var=a.var.copy())
    if float(tmp.X.max()) > 30:
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)

    for name, path in models.items():
        print(f"Running {name}...")
        pred = celltypist.annotate(tmp, model=path, majority_voting=True)
        a.obs[name] = pred.predicted_labels["predicted_labels"].values
        del pred
        gc.collect()
        print(f"  {name} done")

    print("Saving...")
    a.write(out)
    print("Saved:", out)


if __name__ == "__main__":
    main()
