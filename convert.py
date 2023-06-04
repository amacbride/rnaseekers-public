# data import for PDAC (Steele) data

import pandas as pd
from scipy.io import mmread
import scanpy as sc
from anndata import AnnData
import phate

# adjacent normal tissue
adjnorm = ['1', '2', '3']

# skip 14 and 16
tumor = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11A', '11B', '12', '13', '15'];

for a in adjnorm:
    # AdjNorm_TISSUE_*/filtered_feature_bc_matrix
    
    basedir = "AdjNorm_TISSUE_" + a
    dirname = basedir + "/" + "/filtered_feature_bc_matrix"

    print (basedir)
    # Load the data
    print("loading matrix...")
    matrix = mmread(dirname + '/' + 'matrix.mtx.gz').T.tocsc()  # Transpose the matrix, since Scanpy uses cells as rows
    print("loading cell barcodes...")
    barcodes = pd.read_csv(dirname + '/' + 'barcodes.tsv.gz', header=None, sep='\t', names=['cells'], index_col=0)
    print("loading features...")
    features = pd.read_csv(dirname + '/' + 'features.tsv.gz', header=None, sep='\t', names=['uniprot_id', 'gene_symbol', 'feature'], index_col=1)

    # Create AnnData object
    print("creating AnnData object...")
    adata = AnnData(X=matrix,
                    obs=barcodes,
                    var=features)

    # add some embeddings
    # Convert to float32
    adata.X = adata.X.astype('float32')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)  # Calculate the neighborhood graph
    sc.tl.umap(adata)  # Calculate the UMAP
    sc.tl.tsne(adata, n_pcs=40)

    phate_operator = phate.PHATE(n_components=3)
    adata_phate = phate_operator.fit_transform(adata)
    adata.obsm['X_phate'] = adata_phate

    # Save as .h5ad
    print("writing h5ad file...")
    adata.write(basedir + "/" + basedir + '_output.h5ad')

for t in tumor:
    # PDAC_TISSUE_*/filtered_feature_bc_matrix
    
    basedir = "PDAC_TISSUE_" + t
    dirname = basedir + "/" + "/filtered_feature_bc_matrix"

    print (basedir)
    # Load the data
    print("loading matrix...")
    matrix = mmread(dirname + '/' + 'matrix.mtx.gz').T.tocsc()  # Transpose the matrix, since Scanpy uses cells as rows
    print("loading cell barcodes...")
    barcodes = pd.read_csv(dirname + '/' + 'barcodes.tsv.gz', header=None, sep='\t', names=['cells'], index_col=0)
    print("loading features...")
    features = pd.read_csv(dirname + '/' + 'features.tsv.gz', header=None, sep='\t', names=['uniprot_id', 'gene_symbol', 'feature'], index_col=1)

    # Create AnnData object
    print("creating AnnData object...")
    adata = AnnData(X=matrix,
                    obs=barcodes,
                    var=features)

    # add some embeddings
    # Convert to float32
    adata.X = adata.X.astype('float32')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)  # Calculate the neighborhood graph
    sc.tl.umap(adata)  # Calculate the UMAP
    sc.tl.tsne(adata, n_pcs=40)

    phate_operator = phate.PHATE(n_components=3)
    adata_phate = phate_operator.fit_transform(adata)
    adata.obsm['X_phate'] = adata_phate

    # Save as .h5ad
    print("writing h5ad file...")
    adata.write(basedir + "/" + basedir + '_output.h5ad')

print("complete")
