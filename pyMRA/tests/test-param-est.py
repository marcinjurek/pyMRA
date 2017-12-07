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



logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S',level=logging.INFO)

np.random.seed(10)

###### set simulation parameters #####
frac_obs = 0.4 # the fraction of locations for which observations are available
dim_x = 30;  dim_y = 1

sig = 1.0
me_scale=1e-4
kappa = 0.3


### MRA parameters
M=2; J=3; r0=2
critDepth = M+1

    

##### simulate data #####

# generate the underlying process
locs = mt.genLocations(dim_x*dim_y)

Sig = mt.ExpCovFun(locs, l=kappa)
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
    
    cov = lambda _locs1, _locs2: mt.ExpCovFun(_locs1, _locs2, l=kappa)
    MRATree = MRAGraph(locs, M, J, r0, critDepth, cov, y_obs, me_scale)        
    lik = MRATree.getLikelihood()
    return( lik )


def TrueLikelihood(kappa):

    cov = lambda _locs1: mt.ExpCovFun(_locs1, _locs1, l=kappa)
    obs = y_obs[obs_inds].ravel()
    obs_mat = np.matrix(obs).T
    obs_locs = locs[obs_inds,:]
    varY = cov(obs_locs) + np.eye(len(obs))*me_scale

    
    sign, logdet = np.linalg.slogdet(varY)
    full_hand_lik = -len(obs)*0.5*np.log(2*np.pi) - 0.5*logdet -0.5*(obs_mat.T*lng.inv(varY)*obs_mat)
    quad_form = obs_mat.T*lng.inv(varY)*obs_mat
    hand_lik = logdet + quad_form
    
    model = mvn(mean=np.zeros(len(obs)), cov=varY)
    numpy_lik = model.logpdf(obs)

    return( logdet + quad_form )


logging.info("MRA likelihood: %f" % MRALikelihood(kappa))
logging.info("true likelihood: %f" % TrueLikelihood(kappa))



xMRA = opt.minimize(MRALikelihood, [kappa], method='nelder-mead', options={'xtol':1e-1, 'disp':False})
xTrue = opt.minimize(TrueLikelihood, [kappa], method='nelder-mead', options={'xtol':1e-1, 'disp':False})

logging.info("MRA estimate: %f" % xMRA.x[0])
logging.info("True estimate: %f" % xTrue.x[0])

        
