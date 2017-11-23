import pdb
import scipy.optimize as opt
import logging
import numpy as np
import sys
import scipy.linalg as lng
from scipy.stats import multivariate_normal as mvn

sys.path.append('../..')

from pyMRA.MRAGraph import MRAGraph
import pyMRA.MRATools as mt



logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S',level=logging.DEBUG)

#np.random.seed(10)

###### set simulation parameters #####
frac_obs = 0.4 # the fraction of locations for which observations are available
dim_x = 5; dim_y = 5

sig = 1.0
me_scale=0.2
kappa = 0.1


### MRA parameters
M=0; J=3; r0=2
critDepth = M+1

    

##### simulate data #####

# generate the underlying process
locs = mt.genLocations(dim_x*dim_y)

Sig = sig*mt.Matern32(locs, l=kappa, sig=sig)
SigC = np.matrix(lng.cholesky(Sig))

x_raw = np.matrix(np.random.normal(size=(locs.shape[0],1)))
x = SigC.T * x_raw

# generate data
R = me_scale
Rc = np.sqrt(R) if isinstance(me_scale, float) else np.linalg.cholesky(R)
eps = Rc * np.matrix(np.random.normal(size=(locs.shape[0],1)))
y = x + eps


# introducing missing data
obs_inds = np.random.choice(dim_x*dim_y, int(dim_x*dim_y*frac_obs), replace=False)
obs_inds=np.sort(obs_inds)
y_obs = np.empty(np.shape(y));
y_obs[:] = np.NAN;
y_obs[obs_inds] = y[obs_inds]



##### parameter optimization #####

def MRALikelihood(kappa):
    
    cov = lambda _locs1, _locs2: mt.Matern32(_locs1, _locs2, l=kappa, sig=sig)
    MRATree = MRAGraph(locs, M, J, r0, critDepth, cov, y_obs, me_scale)
    lik = MRATree.getLikelihood()
    return( lik )


def TrueLikelihood(kappa):

    cov = lambda _locs1: mt.Matern32(_locs1, _locs1, l=kappa, sig=sig)
    obs = y_obs[obs_inds].ravel()
    obs_mat = np.matrix(obs).T
    obs_locs = locs[obs_inds,:]
    Sig = cov(obs_locs) + np.eye(len(obs))*me_scale

    sign, logdet = np.linalg.slogdet(Sig)
    full_hand_lik = -len(obs)*0.5*np.log(2*np.pi) - 0.5*logdet -0.5*(obs_mat.T*lng.inv(Sig)*obs_mat)
    hand_lik = logdet + obs_mat.T*lng.inv(Sig)*obs_mat
    
    model = mvn(mean=np.zeros(len(obs)), cov=Sig)
    numpy_lik = model.logpdf(obs)

    #logging.debug( "hand lik: %f" % hand_lik )
    #logging.debug( "numpy lik: %f" % numpy_lik )
    return( hand_lik )


logging.info("MRA likelihood: %f" % MRALikelihood(kappa))
logging.info("by-hand likelihood: %f" % TrueLikelihood(kappa))



#xMRA = opt.minimize(MRALikelihood, [kappa], method='nelder-mead', options={'xtol':1e-1, 'disp':False})
#xTrue = opt.minimize(TrueLikelihood, [kappa], method='nelder-mead', options={'xtol':1e-1, 'disp':False})

#logging.info("MRA estimate: %f" % xMRA.x[0])
#logging.info("True estimate: %f" % xTrue.x[0])

        
