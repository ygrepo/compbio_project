import numpy as np
import pandas as pd
import scipy.stats as stats
import warnings

#  This fucntion is from "https://github.com/broadinstitute/pyqtl/blob/master/qtl/norm.py"
def cpm(counts_df, lib_size=None, log=False, prior_count=0.25):
    if lib_size is None:
        lib_size = counts_df.sum(axis=0)
    if log:
        prior_count_scaled = lib_size/np.mean(lib_size) * prior_count
        lib_size <- lib_size + 2 * prior_count_scaled
    lib_size = 1e-6 * lib_size
    if log:
        return np.log2((counts_df + prior_count_scaled)/lib.size)
    else:
        return counts_df / lib_size

#  This fucntion is from "https://github.com/broadinstitute/pyqtl/blob/master/qtl/norm.py"
def calcNormfactors(counts_df, ref=None, logratio_trim=0.3,
                          sum_trim=0.05, acutoff=-1e10, verbose=False):
    # discard genes with all-zero counts
    Y = counts_df.values.copy()
    allzero = np.sum(Y>0,axis=1)==0
    if np.any(allzero):
        Y = Y[~allzero,:]

    # select reference sample
    if ref is None:  # reference sample index
        f75 = np.percentile(Y/np.sum(Y,axis=0), 75, axis=0)
        ref = np.argmin(np.abs(f75-np.mean(f75)))
        if verbose:
            print('Reference sample index: '+str(ref))

    N = np.sum(Y, axis=0)  # total reads in each library

    # with np.errstate(divide='ignore'):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # log fold change; Mg in [1]
        logR = np.log2((Y/N).T / (Y[:,ref]/N[ref])).T
        # average log relative expression; Ag in [1]
        absE = 0.5*(np.log2(Y/N).T + np.log2(Y[:,ref]/N[ref])).T
        v = (N-Y)/N/Y
        v = (v.T + v[:,ref]).T  # w in [1]

    ns = Y.shape[1]
    tmm = np.zeros(ns)
    for i in range(ns):
        fin = np.isfinite(logR[:,i]) & np.isfinite(absE[:,i]) & (absE[:,i] > acutoff)
        n = np.sum(fin)

        loL = np.floor(n*logratio_trim)+1
        hiL = n + 1 - loL
        loS = np.floor(n*sum_trim)+1
        hiS = n + 1 - loS
        rankR = stats.rankdata(logR[fin,i])
        rankE = stats.rankdata(absE[fin,i])
        keep = (rankR >= loL) & (rankR <= hiL) & (rankE >= loS) & (rankE <= hiS)
        # in [1], w erroneously defined as 1/v ?
        tmm[i] = 2**(np.nansum(logR[fin,i][keep]/v[fin,i][keep]) / np.nansum(1/v[fin,i][keep]))

    tmm = tmm / np.exp(np.mean(np.log(tmm)))
    return tmm


'''
Usage:

normalization = cpm(normalization)
tmm = calcNormfactors(normalization)
normalization = normalization / tmm
'''

