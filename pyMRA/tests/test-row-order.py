import scipy.optimize as opt
import gc
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pdb
import time
import sys
import scipy.linalg as lng

sys.path.append('..')

from MRA.MRATree import MRATree
import MRA.MRATools as mt




if __name__=='__main__':

    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S',level=logging.INFO)

    np.random.seed(12)

    
    frac_obs = 0.3

    dim_x = 9
    dim_y = 1
    M=1; J=2; r0=1
    critDepth = 5

    

    

    ### simulate data ###

    sig = 1.0
    me_scale=1e-6
    kappa = 1.0


    
    if dim_y==1:
        locs = mt.genLocations(dim_x)
    else:
        locs = mt.genLocations2d( Nx=dim_x, Ny=dim_y )

    Sig = sig*mt.ExpCovFun(locs, l=kappa)
    SigC = np.matrix(lng.cholesky(Sig))

    x_raw = np.matrix(np.random.normal(size=(locs.shape[0],1)))
    x = SigC.T * x_raw


    R = me_scale
    Rc = np.sqrt(R) if isinstance(me_scale, float) else np.linalg.cholesky(R)
    eps = Rc * np.matrix(np.random.normal(size=(locs.shape[0],1)))
    y = x + eps

    # introducing missing data
    if np.ndim(y)==1:
        square = np.where(np.logical_and(np.logical_and(locs[:,0]>0.2, locs[:,0]<0.6), np.logical_and(locs[:,1]>0.2, locs[:,1]<0.6)))[0]
        not_obs_inds = np.random.choice(square, size=int(0.9*len(square)), replace=False)
        obs_inds = np.setdiff1d(np.arange(len(locs)), not_obs_inds)
    else:
        obs_inds = np.random.choice(dim_x*dim_y, int(dim_x*dim_y*frac_obs), replace=False)
    obs_inds=np.sort(obs_inds)
    y_obs = np.empty(np.shape(y)); y_obs[:] = np.NAN; y_obs[obs_inds] = y[obs_inds]



       
    ### MRA ###

    cov = lambda _locs1, _locs2: mt.ExpCovFun(_locs1, _locs2, l=kappa)
    mraTree = MRATree(locs, M, J, r0, critDepth, cov, y_obs, R)
    BP = mraTree.getBasisFunctionsMatrix(distr="posterior")
    
    


    ### kriging ###

    Sig = mt.ExpCovFun(locs, l=kappa)       
    H = np.matrix(np.eye(len(locs)))
    obsBool = np.isfinite(y_obs).ravel()
    H = H[obsBool,:]
    varY = H * Sig * H.T + R
    SigP = np.linalg.inv(np.linalg.inv(Sig) + H.T * (1/R)*np.eye(len(obs_inds)) * H)


    order = mraTree.root.getOrderFromLeaves()

        
    BP = mraTree.getBasisFunctionsMatrix(distr="posterior")
    mt.dispMat(BP, title="B")
    mt.dispMat(BP*BP.T, title="BB'")
    mt.dispMat(SigP[np.ix_(order, order)], title="Sigma")
    mt.dispMat(BP*BP.T - SigP[np.ix_(order, order)])
