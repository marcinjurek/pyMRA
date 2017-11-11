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

from MRA.MRAGraph import MRAGraph
import MRA.MRATools as mt




if __name__=='__main__':

    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S',level=logging.INFO)
    logging.info("===Unit tests for MRA===")



    #################### #1. MRA in 1d exponential with 0 resolutions is the same as kriging ###########################
    dim_x = 2
    M=0; J=1; r0=dim_x
    frac_obs = 0.5

    sig = 1.0
    me_scale=1e-10
    kappa = 1.0

    locs = mt.genLocations(dim_x)
    Sig = sig*mt.ExpCovFun(locs, l=kappa)
    SigC = np.matrix(lng.cholesky(Sig))

    x_raw = np.matrix(np.random.normal(size=(locs.shape[0],1)))
    x = SigC.T * x_raw

    R = me_scale
    Rc = np.sqrt(R) if isinstance(me_scale, float) else np.linalg.cholesky(R)
    eps = Rc * np.matrix(np.random.normal(size=(locs.shape[0],1)))
    y = x + eps

    obs_inds = np.random.choice(dim_x, int(dim_x*frac_obs), replace=False)
    obs_inds=np.sort(obs_inds)
    y_obs = np.empty(np.shape(y)); y_obs[:] = np.NAN; y_obs[obs_inds] = y[obs_inds]
        
    
    ### MRA ###
    cov = lambda _locs1, _locs2: mt.ExpCovFun(_locs1, _locs2, l=kappa)
    MRATree = MRAGraph(locs, M, J, r0, M+1, cov, y_obs, R)

    xP, sdP = MRATree.predict()
    sdP = sdP.reshape((dim_x, 1), order='A')

    
    ### kriging ###
       
    H = np.matrix(np.eye(len(locs)))
    obsBool = np.isfinite(y_obs).ravel()
    H = H[obsBool,:]

    R = np.eye(len(obs_inds))*me_scale
        
    varY = H * (Sig + R) * H.T
    SigP = np.linalg.inv(np.linalg.inv(Sig) + H.T * np.linalg.inv(R) * H)
    sd = np.sqrt(np.diag(SigP))

    y_obs[np.where(np.isnan(y_obs))] = 0
    xk = Sig*H.T*np.linalg.inv(varY)*H*y_obs
    sdk = sd.reshape((dim_x, 1), order='A')

    order = np.floor(np.log10(np.mean(np.abs(x))))
    assert np.max(np.abs(xP - xk)) < (10**(order-3)), "test #1 (1d exponential, 0 resolutions) not passing"
    assert np.max(np.abs(sdP - sdk)) < 10**(order-3), "test #1 (1d exponential, 0 resolutions) not passing"

    #####################################################################################################################


    
    ############ #2. MRA in 1d exponential with 1 resolution, knots on the boundary the same as kriging #################

    
    dim_x = 3
    M=1; J=2; r0=1
    frac_obs = 0.5

    sig = 1.0
    me_scale=1e-10
    kappa = 1.0

    locs = mt.genLocations(dim_x)
    Sig = sig*mt.ExpCovFun(locs, l=kappa)
    SigC = np.matrix(lng.cholesky(Sig))

    x_raw = np.matrix(np.random.normal(size=(locs.shape[0],1)))
    x = SigC.T * x_raw

    R = me_scale
    Rc = np.sqrt(R) if isinstance(me_scale, float) else np.linalg.cholesky(R)
    eps = Rc * np.matrix(np.random.normal(size=(locs.shape[0],1)))
    y = x + eps

    obs_inds = np.random.choice(dim_x, int(dim_x*frac_obs), replace=False)
    obs_inds=np.sort(obs_inds)
    y_obs = np.empty(np.shape(y)); y_obs[:] = np.NAN; y_obs[obs_inds] = y[obs_inds]
        
    
    ### MRA ###
    cov = lambda _locs1, _locs2: mt.ExpCovFun(_locs1, _locs2, l=kappa)
    MRATree = MRAGraph(locs, M, J, r0, M+1, cov, y_obs, R)

    xP, sdP = MRATree.predict()
    sdP = sdP.reshape((dim_x, 1), order='A')
  

    ### kriging ###
       
    H = np.matrix(np.eye(len(locs)))
    obsBool = np.isfinite(y_obs).ravel()
    H = H[obsBool,:]

    R = np.eye(len(obs_inds))*me_scale
        
    varY = H * (Sig + R) * H.T
    SigP = np.linalg.inv(np.linalg.inv(Sig) + H.T * np.linalg.inv(R) * H)
    sd = np.sqrt(np.diag(SigP))

    y_obs[np.where(np.isnan(y_obs))] = 0
    xk = Sig*H.T*np.linalg.inv(varY)*H*y_obs
    sdk = sd.reshape((dim_x, 1), order='A')

    order = np.floor(np.log10(np.mean(np.abs(x))))
    assert np.max(np.abs(xP - xk)) < (10**(order-3)), "test #2 (1d exponential, 1 resolution) not passing"
    assert np.max(np.abs(sdP - sdk)) < 10**(order-3), "test #2 (1d exponential, 1 resolution) not passing"

    #####################################################################################################################


    
    
    ############################################# #3. MRA works in 2d ##################################################

    ### We don't plot the results here because we use the Matern covariance functions. It does not give us a screening
    ### effect so the MRA results are not exact. This is only to show that the code works properly (doesn't crash)
    ### when we deal with 2d domain.
    
    
    np.random.seed(12)
    
    frac_obs = 0.7

    dim_x = 10
    dim_y = 10
    M=2; J=3; r0=2
    critDepth = 6

    sig=1.0
    me_scale=1e-6
    kappa=0.5
    

    ### simulate data ###

    locs = mt.genLocations2d( Nx=dim_x, Ny=dim_y )

    Sig = sig*mt.Matern32(locs, l=kappa, sig=sig)
    SigC = np.matrix(lng.cholesky(Sig))

    x_raw = np.matrix(np.random.normal(size=(locs.shape[0],1)))
    x = SigC.T * x_raw
    del SigC

    R = me_scale
    Rc = np.sqrt(R) if isinstance(me_scale, float) else np.linalg.cholesky(R)
    eps = Rc * np.matrix(np.random.normal(size=(locs.shape[0],1)))
    y = x + eps

    # missing data

    obs_inds = np.random.choice(dim_x*dim_y, int(dim_x*dim_y*frac_obs), replace=False)
    obs_inds=np.sort(obs_inds)
    y_obs = np.empty(np.shape(y)); y_obs[:] = np.NAN; y_obs[obs_inds] = y[obs_inds]
    del Rc
       
    
    
    ### MRA ###

    cov = lambda _locs1, _locs2: mt.Matern32(_locs1, _locs2, l=kappa, sig=sig)
    MRATree = MRAGraph(locs, M, J, r0, critDepth, cov, y_obs, R)

    xP, sdP = MRATree.predict()
    sdP = sdP.reshape((dim_x, dim_y), order='A')
    


    ### kriging ###


    Sig = mt.Matern32(locs, l=kappa, sig=sig)
        
    H = np.matrix(np.eye(dim_x*dim_y))
    obsBool = np.isfinite(y_obs).ravel()
    H = H[obsBool,:]

    if np.ndim(R)==1:
        R = np.diag(R)
    elif isinstance(R, float):
        R = np.eye(len(obs_inds))*R
        
    varY = H * Sig * H.T + R
    SigP = np.linalg.inv(np.linalg.inv(Sig) + H.T * np.linalg.inv(R) * H)
    sd = np.sqrt(np.diag(SigP))

    y_obs[np.where(np.isnan(y_obs))] = 0
    xk = Sig*H.T*np.linalg.inv(varY)*H*y_obs
    sdk = sd.reshape((dim_x, dim_y), order='A')
            
    #####################################################################################################################

    
    logging.info("All tests passed")

    
