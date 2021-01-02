import numpy as np
import pandas as pd
import scanpy as sc
from numpy.linalg import svd

# load the data from .MAT files
from scipy.io import loadmat

def loadmatfile(file_dir):
    
    f = loadmat(file_dir)
    names = f['data']['names'][0][0]
    name_list = [str(i[0]) for i in names[0]]
    intens_matrix = f['data']['intens'][0][0].T
    mzs = f['data']['mzs'][0][0][0]

    print('Loaded intensity matrix with shape {}'.format(intens_matrix.shape))
    intens_df = pd.DataFrame(intens_matrix)
    intens_df = intens_df.set_index([pd.Index(name_list)])
    #intens_df[intens_df==0]=1
    #intens_df = np.log(intens_df)
    intens_df.columns = np.round(mzs,2)
    
    return intens_df


def comp_lev(A,k,axis):
    
    U,D,V = svd(A, full_matrices=False)
    
    if axis==0:
        lev = np.sum(U[:,:k]**2,axis=1)
    if axis==1:
        lev = np.sum(V[:k,:]**2,axis=0)
    
    return lev


def cx_decomp(A,k,n_choose):
    
    lev = comp_lev(A,k,1)
    lev_rank_idx = lev.argsort()[::-1]
    
    C = np.dot(A[:,lev_rank_idx[:n_choose]], np.diag(1/lev[:n_choose]))
    
    X = np.dot(np.linalg.pinv(C), A)
    
    return C, X, lev


def cur_decomp(A,k,n_choose_col,n_choose_row):
    
    lev_col = comp_lev(A,k,1)
    lev_row = comp_lev(A,k,0)

    lev_col_idx = lev_col.argsort()[::-1]
    lev_row_idx = lev_row.argsort()[::-1]
    
    C = np.dot(A[:,lev_col_idx[:n_choose_col]], np.diag(1/lev_col_idx[:n_choose_col]))
    R = np.dot(np.diag(1/lev_row[:n_choose_row]), A[lev_row_idx[:n_choose_row],:])
    U = np.dot( np.dot(np.linalg.pinv(C), A), np.linalg.pinv(R) )
    
    return C, U, R


def process(adata, n_pcs, min_dist, resolution):
    #sc.pp.normalize_total(adata, target_sum=None)

    sc.pp.log1p(adata)
    print('pca..')
    sc.tl.pca(adata, svd_solver='arpack')

    print('getting neighbors..')
    sc.pp.neighbors(adata, n_neighbors=15, metric='cosine', n_pcs=n_pcs)
    print('tsne...')
    sc.tl.tsne(adata)  # , init_pos=adata.obsm['X_pca'])
    print('leiden...')
    sc.tl.leiden(adata, resolution=resolution)
