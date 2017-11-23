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

sys.path.append('../..')

from pyMRA.MRAGraph import MRAGraph
from pyMRA import MRATools as mt




if __name__=='__main__':

    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S',level=logging.INFO)

    #np.random.seed(12)

    test_small = True
    diagnose = True
    
    frac_obs = 0.4
    if test_small:
        dim_x = 10
        dim_y = 10
        M=2; J=4; r0=4
        critDepth = 5
    else:
        dim_x = 20
        dim_y = 20
        M=4; J=4; r0=4
        critDepth = 6
    
   

    ### simulate data ###

    sig = 1.0
    me_scale=1e-1
    kappa = 0.3
   
    logging.info('=== simulating data ===')
    logging.info('simulating locations')
    if dim_y==1:
        locs = mt.genLocations(dim_x)
    else:
        locs = mt.genLocations2d( Nx=dim_x, Ny=dim_y )

    logging.info('calculating covariance matrix')
    if test_small:
        #Sig = sig*mt.ExpCovFun(locs, l=kappa)
        Sig = sig*mt.Matern32(locs, l=kappa, sig=sig)
        SigC = np.matrix(lng.cholesky(Sig))
    else:
        np.save("SigC.npy", SigC)
        # SigC = np.load("SigC.npy")

    x_raw = np.matrix(np.random.normal(size=(locs.shape[0],1)))
    x = SigC.T * x_raw
    del SigC


    
    logging.info("simulating observations")
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
    del Rc
    gc.collect()

  
    
    
    ### MRA ###

    logging.info('=== starting MRA ===')
    start = time.time()
    #cov = lambda _locs1, _locs2: mt.ExpCovFun(_locs1, _locs2, l=kappa)
    cov = lambda _locs1, _locs2: mt.Matern32(_locs1, _locs2, l=kappa)
    MRATree = MRAGraph(locs, M, J, r0, critDepth, cov, y_obs, R)

    xP, sdP = MRATree.predict()
    sdP = sdP.reshape((dim_x, dim_y), order='A')
    #sdP = np.flipud(sdP.reshape((dim_x, dim_y), order='A'))
    MRATime = time.time()-start

    BP = MRATree.getBasisFunctionsMatrix(distr="posterior")
    
    logging.info('MRA predictions finished. It took {:.2}s'.format(MRATime))
    


   
    
    ### diagnostic plots ###

    if diagnose:        
        MRATree.drawBasisFunctions()
        MRATree.drawGridAndObs()
        MRATree.drawKnots()
