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

#sys.path.append('../..')

from pyMRA.MRATree import MRATree
from pyMRA import MRATools as mt




if __name__=='__main__':

    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S',level=logging.INFO)

    np.random.seed(11)

    test_small = False
    diagnose = False
    krig = False
    compare=False
    find_params=False
    
    frac_obs = 0.4
    if test_small:
        dim_x = 20
        dim_y = 20
        M=5; J=2; r0=3
        critDepth = M+1
    else:
        dim_x = 100
        dim_y = 100
        M=5; J=4; r0=5
        critDepth = 0
    

    

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
        Sig = sig*mt.ExpCovFun(locs, l=kappa)
        #Sig = sig*mt.Matern32(locs, l=kappa, sig=sig)
        SigC = np.matrix(lng.cholesky(Sig))
    else:
        #Sig = sig*mt.ExpCovFun(locs, l=kappa)
        #SigC = np.matrix(lng.cholesky(Sig))
        #np.save("SigC.npy", SigC)
        SigC = np.load("SigC.npy")

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
    cov = lambda _locs1, _locs2: mt.ExpCovFun(_locs1, _locs2, l=kappa)
    #cov = lambda _locs1, _locs2: mt.Matern32(_locs1, _locs2, l=kappa)
    mraTree = MRATree(locs, M, J, r0, critDepth, cov, y_obs, R)
    logging.info( "leaf size: %f" % mraTree.avgLeafSize() )

    
    xP, sdP = mraTree.predict()
    sdP = sdP.reshape((dim_x, dim_y), order='A')
    #sdP = np.flipud(sdP.reshape((dim_x, dim_y), order='A'))
    MRATime = time.time()-start

    B = mraTree.getBasisFunctionsMatrix(distr="prior")

    logging.info('avg leaf size: %f' % mraTree.avgLeafSize())
    logging.info('max leaf size: %f' % mraTree.maxLeaf())
    logging.info('min leaf size: %f' % mraTree.minLeaf())
    logging.info('MRA predictions finished. It took {:.2}s'.format(MRATime))
    








    ### kriging ###

    if krig:
        logging.info('=== Starting ordinary kriging ===')

        start = time.time()
        Sig = mt.ExpCovFun(locs, l=kappa)
        
        H = np.matrix(np.eye(len(locs)))
        obsBool = np.isfinite(y_obs).ravel()
        H = H[obsBool,:]

        #if np.ndim(R)==1:
        #    R = np.diag(R)
        #elif isinstance(R, float):
        #    R = np.eye(len(obs_inds))*R

        varY = H * Sig * H.T + R

        SigP = np.linalg.inv(np.linalg.inv(Sig) + H.T * (1/R)*np.eye(len(obs_inds)) * H)
        sd = np.sqrt(np.diag(SigP))

        y_obs_k = np.zeros(np.shape(y_obs))
        y_obs_k[np.where(np.isnan(y_obs))] = 0

        xk = Sig*H.T*np.linalg.inv(varY)*H*y_obs_k
        sdk = sd.reshape((dim_x, dim_y), order='A')
        #sdk = np.flipud(sd.reshape((dim_x, dim_y), order='A'))
        regTime = time.time() - start

        logging.info('Kriging finished. It took {:.2}s'.format(regTime))
  
        #illustrate what simple kriging is about
        fig = plt.figure(figsize=plt.figaspect(0.2))
        axMean = fig.add_subplot(111)
        lineO = axMean.plot(locs[obs_inds], y[obs_inds], color='#deaa87', marker='o', linestyle='None', markersize='7', label='observations (Y)')
        lineP = axMean.plot(locs, x, color='#0000ff', linewidth=2, label='true state (X)')

        grid_ypos = np.min(np.vstack((np.min(x), y[obs_inds]))) - 0.1
        axMean.plot(locs, np.ones(len(locs))*grid_ypos, marker='o', color='#000000', markersize='7', linestyle='None', label='grid')
        
        axMean.legend(fontsize='x-large')
        plt.xticks(fontsize='x-large')
        plt.yticks(fontsize='x-large')
        plt.show()





   
    
    ### diagnostic plots ###

    if diagnose:
        mraTree.drawBMatrix("prior")
        mraTree.drawSparsityPat("prior")
        mraTree.drawBMatrix("posterior")
        mraTree.drawSparsityPat("posterior")
        
        mraTree.drawBasisFunctions("prior")
        mraTree.drawBasisFunctions("posterior")
        mraTree.drawGridAndObs()
        mraTree.drawKnots()



        


    ### parameter optimization ###

    def likelihood(kappa):
    #def likelihood(params):
        """
        Evaluates the likelihood of the model for a given set of parameters
        params=(kappa, sigma, R)
        """
        #par = dict(zip(('kappa', 'sigma', 'R'), np.abs(params)))      
        #logging.debug("kappa=%f, sig=%f, R=%f" % (par['kappa'], par['sigma'], par['R']))

        cov = lambda _locs1, _locs2: mt.ExpCovFun(_locs1, _locs2, l=kappa)
        #cov = lambda _locs1, _locs2: mt.Matern32(_locs1, _locs2, l=par['kappa'], sig=1.0)
        mraTree = MRATree(locs, M, J, r0, critDepth, cov, y_obs, me_scale)
        #mraTree = MRATree(locs, M, J, r0, critDepth, cov, y_obs, np.abs(par['R']))
        lik = mraTree.getLikelihood()
        return( lik )

    if find_params:
        logging.info("=== finding optimal parameter values ===")
        start = time.time()
        xmin = opt.minimize(likelihood, [kappa], method='nelder-mead', options={'xtol':1e-1, 'disp':False})
        #xmin = opt.minimize(likelihood, [kappa, sig, me_scale], method='nelder-mead', options={'xtol':1e-2, 'disp':False})
        opttime = time.time() - start

        
        
        logging.info("ML estimates found. It took %ds" % opttime)
        logging.info(str(xmin))
        #logging.info(dict(zip(('kappa'),xmin['x'])))
        #logging.info(dict(zip(('kappa', 'sigma', 'R'),xmin['x'])))
    


        
    
    ### compare results ###

    if krig and compare:
        if dim_y>1:


            my_map = mpl.cm.get_cmap('Spectral')
            my_map.set_bad(color='grey')

            vmin = np.min(np.vstack((y_obs_k, x)))
            vmax = np.max(np.vstack((y_obs_k, x)))
            
            mt.dispMat(xP.reshape((dim_x, dim_y)), cmap='RdYlBu', title="predicted x")          
            mt.dispMat(x.reshape((dim_x, dim_y)), cmap='Spectral', title="true state (X)", vmax=vmax, vmin=vmin, fontsize='xx-large')
            if test_small:
                mt.dispMat(xk.reshape((dim_x, dim_y), order='A'), cmap='RdYlBu', title='kriged x')
                mt.dispMat((xP-xk).reshape((dim_x, dim_y)), cmap='coolwarm', title="predicted x - true x")
            mt.dispMat(y_obs.reshape((dim_x, dim_y), order='A'), cmap=my_map, title="observations (Y)", vmax=vmax, vmin=vmin, fontsize='xx-large')
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
