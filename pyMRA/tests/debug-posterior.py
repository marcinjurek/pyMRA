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

    np.random.seed(11)

    diagnose = True
    krig = True
    compare=True
    
    frac_obs = 0.5

    dim_x = 2; dim_y = 1
    M=1; J=1; r0=1
    critDepth = M+1
    

    

    ### simulate data ###

    sig = 1.0
    me_scale=1e-4
    kappa = 0.3

    if dim_y==1:
        locs = mt.genLocations(dim_x)
    else:
        locs = mt.genLocations2d( Nx=dim_x, Ny=dim_y )

    Sig = sig*mt.ExpCovFun(locs, l=kappa)
    SigC = np.matrix(lng.cholesky(Sig))
        
    x_raw = np.matrix(np.random.normal(size=(locs.shape[0],1)))
    x = SigC.T * x_raw


    Rc = np.sqrt(me_scale)
    eps = Rc * np.matrix(np.random.normal(size=(locs.shape[0],1)))
    y = x + eps

    obs_inds = np.random.choice(dim_x*dim_y, int(dim_x*dim_y*frac_obs), replace=False)
    obs_inds=np.sort(obs_inds)
    y_obs = np.empty(np.shape(y)); y_obs[:] = np.NAN; y_obs[obs_inds] = y[obs_inds]






    
    
    
    
    ### MRA ###

    cov = lambda _locs1, _locs2: mt.ExpCovFun(_locs1, _locs2, l=kappa)
    MRATree = MRAGraph(locs, M, J, r0, critDepth, cov, y_obs, me_scale)

    xP, sdP = MRATree.predict()
    sdP = sdP.reshape((dim_x, dim_y), order='A')

    B = MRATree.getBasisFunctionsMatrix(distr="prior", timesKC=True)
    BP = MRATree.getBasisFunctionsMatrix(distr="posterior", timesKC=True)
    

    # some auxiliary matrices
    I = np.matrix(np.eye(dim_x*dim_y))

    H = np.matrix(np.eye(len(locs)))
    obsBool = np.isfinite(y_obs).ravel()
    H = H[obsBool,:]

    R = np.matrix(me_scale * np.eye(len(obs_inds)))



    # Sigma posterior reconstructed from B=MRD(Sigma)
    recSigP =  B * lng.inv(I + B.T * H.T * lng.inv(R) * H * B) * B.T

    





    ### kriging ###
      
    varY = H * Sig * H.T + R
    SigP = Sig - Sig * H.T * lng.inv(varY) * H * Sig
    

    sd = np.sqrt(np.diag(SigP))
    
    y_obs[np.where(np.isnan(y_obs))] = 0
    
    xk = Sig*H.T*np.linalg.inv(varY)*H*y_obs
    sdk = sd.reshape((dim_x, dim_y), order='A')
    #sdk = np.flipud(sd.reshape((dim_x, dim_y), order='A'))

  
    # illustrate what simple kriging is about
    # fig = plt.figure(figsize=plt.figaspect(0.2))
    # axMean = fig.add_subplot(111)
    # lineP = axMean.plot(locs, x, 'g-', linewidth=2, label='true state (X)')
    # lineO = axMean.plot(locs[obs_inds], y[obs_inds], 'ko', markersize='7', label='observations (Y)')
    
    # grid_ypos = np.min(np.vstack((np.min(x), y[obs_inds]))) - 0.1
    # axMean.plot(locs, np.ones(len(locs))*grid_ypos, marker='o', markersize='7', linestyle='None', label='grid')
    
    # axMean.legend(loc=(0.6, 0.2), fontsize='x-large')
    # plt.xticks(fontsize='x-large')
    # plt.yticks(fontsize='x-large')
    # plt.show()





   
    
    ### diagnostic plots ###

    if diagnose:
        #MRATree.drawBMatrix("prior")
        #MRATree.drawSparsityPat("prior")
        #MRATree.drawBMatrix("posterior")
        #MRATree.drawSparsityPat("posterior")
        
        #MRATree.drawBasisFunctions("prior")
        #MRATree.drawBasisFunctions("posterior")
        MRATree.drawGridAndObs()
        MRATree.drawKnots()



        
       

        
    
    ### compare results ###

    if krig and compare:
        if dim_y>1:
            mt.dispMat(xP.reshape((dim_x, dim_y)), cmap='RdYlBu', title="predicted x")
            mt.dispMat(x.reshape((dim_x, dim_y)), cmap='RdYlBu', title="true x")
            if test_small:
                mt.dispMat(xk.reshape((dim_x, dim_y), order='A'), cmap='RdYlBu', title='kriged x')
                mt.dispMat((xP-xk).reshape((dim_x, dim_y)), cmap='coolwarm', title="predicted x - true x")
            mt.dispMat(y_obs.reshape((dim_x, dim_y), order='A'), cmap='RdYlBu', title="observations")
            mt.dispMat(sdP, cmap='Reds', title="predicted sd")

        else:
            # compare MRA and kriging results
            fig = plt.figure(figsize=plt.figaspect(0.5))
            axMean = fig.add_subplot(121)
            lineP = axMean.plot(locs, xP, 'r-')
            lineO = axMean.plot(locs[obs_inds], y[obs_inds], 'ko', markersize='4')
            lineT = axMean.plot(locs, xk, linestyle='dashed', color="gray")
            axMean.set_title("mean")
               
            axSd = fig.add_subplot(122)
            #axSd.plot(locs, sdP, 'r-')
            #axSd.plot(locs, sd, linestyle='dashed', color="gray")
            axSd.plot(locs, sd-sdP.T.ravel(), 'r-')
            axSd.set_title("sd")
            plt.show()

            pdb.set_trace()
