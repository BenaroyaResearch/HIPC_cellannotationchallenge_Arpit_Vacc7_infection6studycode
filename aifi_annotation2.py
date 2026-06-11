import scanpy as sc
import anndata as ad
import celltypist
import gc

REF = '/homes/tedwards/projects/HIPC_annotation_challenge/data/inputData/references'
models = {
    'Immunelow': f'{REF}/Immune_All_Low.pkl',
    'AIFI_L2': f'{REF}/ref_pbmc_clean_celltypist_model_AIFI_L2_2024-04-19.pkl',
    'AIFI_L3': f'{REF}/ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19.pkl',
}

p = '/nfs/amishra/Arpit/hipc_challenge/infection6_annotation_full/reports/artifacts/concatenated.h5ad'
out = p.replace('concatenated.h5ad', 'concatenated_AIFI2.h5ad')

print('Loading h5ad...')
a = ad.read_h5ad(p)
print(f'Loaded: {a.shape} | X max: {float(a.X.max())} | layers: {list(a.layers.keys())}')

src = a.layers['counts'].copy() if 'counts' in a.layers else a.X.copy()
tmp = ad.AnnData(X=src, var=a.var.copy())
if float(tmp.X.max()) > 30:
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)

for name, path in models.items():
    print(f'Running {name}...')
    pred = celltypist.annotate(tmp, model=path, majority_voting=True)
    a.obs[name] = pred.predicted_labels['predicted_labels'].values
    del pred
    gc.collect()
    print(f'  {name} done')

print('Saving...')
a.write(out)
print('Saved:', out)
