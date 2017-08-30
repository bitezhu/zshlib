import numpy as np
import sys

def cov(data):
    """
    Covariance matrix
    note: specifically for mean-centered data
    note: numpy's `cov` uses N-1 as normalization
    """
    X = data
    #return np.dot(X.T, X) / (X.shape[0] - 1) 
    return np.dot(X.T, X) / X.shape[0]


def pca(data, pc_count = 2,fraction=.95,force=0):
    """
    Principal component analysis using eigenvalues
    note: this mean-centers and auto-scales the data (in-place)
    """
    data = np.asarray(data,dtype=np.float64)
    data -= np.mean(data, 0)
    data /= np.std(data, 0)
    C = cov(data)
    E, Vector = np.linalg.eigh(C)
    sumvariance = np.cumsum(np.argsort(E)[::-1],dtype=float)
    sumvariance /= sumvariance[-1]
    if sumvariance[pc_count-1] < fraction:
        pc_count += 1
        sys.stderr.write("variance below 95%, will enable pc_count equals 3 automaticlly!!\nBut if force was not setted 1, pc_count will be unchanged.\n")
    if not force:
        pc_count -= 1 # in case you want two-dim plot anyway
    key = np.argsort(E)[::-1][:pc_count]
    U = np.dot(data,Vector)
    E, Vector = E[key], Vector[:, key]
    return U, E, Vector



