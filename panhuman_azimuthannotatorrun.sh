#!/usr/bin/env bash
set -o pipefail   # so annotate failures are caught even through the tee pipe

INPUT_ROOT="/homes/tedwards/projects/HIPC_annotation_challenge/data/inputData"
WORK_ROOT="$HOME/hipc_run"
PIXI_MANIFEST="$HOME/envs/azimuth/pixi.toml"

# -em / -umap SKIP embeddings + UMAP (per the help text).
# Leave empty if you want them generated.
ANNOTATE_FLAGS="-em -umap"

SKIP=("infection_study_06" "vaccination_study_07")

mkdir -p "$WORK_ROOT"

# --- write the obs -> CSV exporter once ---
EXPORTER="$WORK_ROOT/export_obs.py"
cat > "$EXPORTER" << 'PYEOF'
import sys
import anndata as ad

p = sys.argv[1]
obs = ad.read_h5ad(p, backed='r').obs          # backed: only obs loaded, cheap
keep = [c for c in obs.columns if '.scores_' not in c]
out = obs[keep].copy()
out.index.name = 'cell_barcode'
csv = p.replace('.h5ad', '_annotations.csv')
out.to_csv(csv)
print('wrote', csv, '|', out.shape)
print(out.columns.tolist())
PYEOF

pixi_run() { pixi run --manifest-path "$PIXI_MANIFEST" "$@"; }

for study_dir in "$INPUT_ROOT"/*/; do
    study=$(basename "$study_dir")

    for s in "${SKIP[@]}"; do
        if [[ "$study" == "$s" ]]; then
            echo ">>> skipping $study"
            continue 2
        fi
    done

    src=$(find "$study_dir" -maxdepth 1 -name '*_processed.h5ad' | head -n1)
    if [[ -z "$src" ]]; then
        echo ">>> no *_processed.h5ad in $study — skipping"
        continue
    fi

    dest_dir="$WORK_ROOT/$study"
    mkdir -p "$dest_dir"
    cp -n "$src" "$dest_dir/"
    dest="$dest_dir/$(basename "$src")"
    ann="${dest%.h5ad}_ANN.h5ad"

    echo ">>> [$(date +%H:%M:%S)] annotating $study"
    if ! pixi_run annotate $ANNOTATE_FLAGS "$dest" 2>&1 | tee "$dest_dir/annotate.log"; then
        echo ">>> FAILED annotate: $study (see log)"
        continue
    fi

    if [[ ! -f "$ann" ]]; then
        echo ">>> expected output missing: $ann — skipping CSV export"
        continue
    fi

    echo ">>> exporting CSV for $study"
    if pixi_run python "$EXPORTER" "$ann" 2>&1 | tee -a "$dest_dir/annotate.log"; then
        echo ">>> done: $study"
    else
        echo ">>> FAILED csv export: $study"
    fi
done
