import logging
import numpy as np
import pdb
import sys
import scipy.linalg as lng

sys.path.append('../..')

from pyMRA.MRAGraph import MRAGraph
from pyMRA import MRATools as mt




if __name__=='__main__':

    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S',level=logging.INFO)


   
    frac_obs = 0.4

    dim_x = 100; dim_y = 1
    M=3; J=3; r0=2
    critDepth = M+1
    


    ### simulate data ###

    sig = 1.0; me_scale=1e-1; kappa = 0.3


    if dim_y==1:
        locs = mt.genLocations(dim_x)
    else:
        locs = mt.genLocations2d( Nx=dim_x, Ny=dim_y )


    Sig = sig*mt.Matern32(locs, l=kappa, sig=sig)
    SigC = np.matrix(lng.cholesky(Sig))

    x_raw = np.matrix(np.random.normal(size=(locs.shape[0],1)))
    x = SigC.T * x_raw

    R = me_scale
    Rc = np.sqrt(R) if isinstance(me_scale, float) else np.linalg.cholesky(R)
    eps = Rc * np.matrix(np.random.normal(size=(locs.shape[0],1)))
    y = x + eps

    # introducing missing data
    obs_inds = np.random.choice(dim_x*dim_y, int(dim_x*dim_y*frac_obs), replace=False)
    obs_inds=np.sort(obs_inds)
    y_obs = np.empty(np.shape(y)); y_obs[:] = np.NAN; y_obs[obs_inds] = y[obs_inds]
    
    
    
    ### MRA ###

    cov = lambda _locs1, _locs2: mt.Matern32(_locs1, _locs2, l=kappa)
    MRATree = MRAGraph(locs, M, J, r0, critDepth, cov, y_obs, R)
    xP, sdP = MRATree.predict()
   

    
    ### basis functions plots ###
        
    MRATree.drawBasisFunctions(distr="prior")
    MRATree.drawBasisFunctions(distr="posterior")
