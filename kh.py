import multiprocessing as mp
from functools import partial
import numpy as np
from tqdm import tqdm
import anndata
from typing import Union
import scipy
import logging

def random_feats(X: np.ndarray,
                gamma: Union[int, float] = 1,
                frequency_seed: int = None):
    """Computes random Fourier frequency features: https://papers.nips.cc/paper/2007/hash/013a006f03dbc5392effeb8f18fda755-Abstract.html
    Parameters
    X: np.ndarray
        array of input data (dimensions = cells x features)
    gamma: Union([int, float]) (default = 1)
        scale for standard deviation of the normal distribution 
    frequency_seed: int (default = None):
        random state parameter     
    ----------
    Returns
    phi: np.ndarray
        random Fourier frequency features (dimensions = cells x 2000)
    ----------
    """
    scale = 1 / gamma

    if (frequency_seed is not None):
        np.random.seed(frequency_seed)
        W = np.random.normal(scale = scale, size = (X.shape[1], 1000))
    else:
        W = np.random.normal(scale = scale, size = (X.shape[1], 1000))

    XW = np.dot(X, W)
    sin_XW = np.sin(XW)
    cos_XW = np.cos(XW)
    phi = np.concatenate((cos_XW, sin_XW), axis=1)

    return phi

def kernel_herding(phi: np.ndarray,
                    num_subsamples: int = None):
    """Performs kernel herding subsampling: https://arxiv.org/abs/1203.3472
    Parameters
    phi: np.ndarray
        random Fourier frequency features (dimensions = cells x 2000)
    num_subsamples: int (default = None)
        number of cells to subsample  
    ----------
    Returns
    kh_indices: np.ndarray
        indices of subsampled cells
    ----------
    """
    w_t = np.mean(phi, axis=0)
    w_0 = w_t
    kh_indices = []
    while len(kh_indices) < num_subsamples:
        new_ind = np.argmax(np.dot(phi, w_t))
        w_t = w_t + w_0 - phi[new_ind]
        kh_indices.append(new_ind)
        kh_indices = list(set(kh_indices))

    return kh_indices

def _parse_input(adata: anndata.AnnData):
    """accesses and parses data from adata object
    Parameters
    adata: anndata.AnnData
        annotated data object where adata.X is the attribute for preprocessed data
    ----------
    Returns
    X: np.ndarray
        array of data (dimensions = cells x features)
    ----------
    """
    try:
        if isinstance(adata, anndata.AnnData):
            X = adata.X.copy()
        if isinstance(X, scipy.sparse.csr_matrix):
            X = np.asarray(X.todense())
    except NameError:
        pass
    
    return X

def kernel_herding_main(sample_set_ind,
                        X: np.ndarray = None,
                        gamma: Union[int, float] = 1,
                        frequency_seed: int = None,
                        num_subsamples: int = 500):
    """Performs kernel herding subsampling on a single sample-set using random features
    Parameters
    X: np.ndarray
        array of input data (dimensions = cells x features)
    gamma: Union([int, float]) (default = 1)
        scale for standard deviation of the normal distribution 
    frequency_seed: int (default = None):
        random state parameter     
    num_samples: int (default = None)
        number of cells to subsample 
    sample_set_ind: np.ndarray
        array containing the indices of the sample-set to subsample. if you'd like to use all cells within X, please pass in np.arange(0, X.shape[0])
    ----------
    Returns
    kh_indices: np.ndarray
        indices of subsampled cells within the sample-set
    ----------
    """
    X = X[sample_set_ind, :]
    phi = random_feats(X, gamma = gamma, frequency_seed = frequency_seed)
    kh_indices = kernel_herding(phi, num_subsamples)

    return kh_indices

def sketch(adata,
            sample_set_key: str = None,
            sample_set_inds = None,
            gamma: Union[int, float] = 1,
            frequency_seed: int = None,
            num_subsamples: int = 500,
            n_jobs: int = -1):
    """constructs a sketch using kernel herding and random Fourier frequency features
    Parameters
    adata: anndata.Anndata
        annotated data object (dimensions = cells x features)
    sample_set_key: str (default = None)
        string referring to the key within adata.obs that contains the sample-sets to subsample
            ~ if sample_set_key is None, will parse according to sample_set_inds
            ~ if sample_set_key is None and sample_set_inds is None, will use all cells as a single sample-set 
    sample_set_inds: list (default = None)
        list of arrays containig the indices of the sample-sets to subsample. (dimensions = len(sample_sets)) e.g. [np.array([]), np.array([]), ... , np.array([])]
            ~ if sample_set_key is None and sample_set_inds is None, will use all cells as a single sample-set 
    gamma: Union([int, float]) (default = 1)
        scale for standard deviation of the normal distribution within random Fourier frequency feature computation
    frequency_seed: int (default = None):
        random state parameter     
    num_samples: int (default = None)
        number of cells to subsample per sample-set
    n_jobs: int (default = -1)
        number of tasks
    ----------
    Returns
    kh_indices: np.ndarray
        list of indices of subsampled cells per sample-set e.g. [np.array(ind0_S0..indx_S0), np.array(ind0_S1..indx_S1), ... , np.array(ind0_SX..indx_SX)]
    adata_subsample: anndata.AnnData
        annotated data object containing subsampled data
    ----------
    """
    if n_jobs == -1:
        n_jobs = mp.cpu_count()
    elif n_jobs < -1:
        n_jobs = mp.cpu_count() + 1 + n_jobs

    if isinstance(adata, anndata.AnnData) and (sample_set_key is not None):
        sample_set_id = np.asarray(adata.obs[sample_set_key].cat.categories)
        sample_set_inds = [np.where(adata.obs[sample_set_key] == i)[0] for i in sample_set_id]
    elif sample_set_inds is None:
        sample_set_inds = [np.arange(0, adata.X.shape[0])]
        
    min_cell_size = min([len(i) for i in sample_set_inds])
    if num_subsamples > min_cell_size:
        logging.warning(f'Number of subsamples per sample-set {num_subsamples} is greater than the maximum number of cells in the smallest sample-set {min_cell_size}. \n Performing subsampling using {min_cell_size} cells per sample-set')
        num_subsamples = min_cell_size
        
    n_sample_sets = len(sample_set_inds)
    X = _parse_input(adata)

    p = mp.Pool(n_jobs)

    kh_indices = []
    for result in tqdm(p.imap(partial(kernel_herding_main, X = X, gamma = gamma, frequency_seed = frequency_seed, num_subsamples = num_subsamples), sample_set_inds),total = n_sample_sets, desc = 'performing subsampling'):
        kh_indices.append(result)

    adata_subsample = []
    for i in range(0, len(sample_set_inds)):
        sample_set = adata[sample_set_inds[i], :]
        subsampled_sample_set = sample_set[sample_set.obs.iloc[kh_indices[i]].index]
        adata_subsample.append(subsampled_sample_set)

    adata_subsample = anndata.concat(adata_subsample)

    return kh_indices, adata_subsample