import pdb
import scipy.optimize as opt
import logging
import numpy as np
import sys
import scipy.linalg as lng
from scipy.stats import multivariate_normal as mvn

sys.path.append('../..')

from pyMRA.MRATree import MRATree
import pyMRA.MRATools as mt



logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S',level=logging.INFO)

np.random.seed(10)

###### set simulation parameters #####
frac_obs = 0.4 # the fraction of locations for which observations are available
dim_x = 30;  dim_y = 1

sig = 1.0
me_scale=0.2
kappa = 0.05


### MRA parameters
M=2; J=4; r0=2
critDepth = M+1



#filename = '/home/marcin/MRF/data/Exp_Theta0.1_X100_Y100.csv'
#filename = '/home/marcin/MRF/data/Exp_Theta0.1_X100_Y100_missing_all.csv'
#filename = '/home/marcin/MRF/data/sat_temps.csv'
filename = '/home/marcin/pyMRA/pyMRA/data/small/Exp_Theta0.1_X10_Y10.csv'
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S',level=logging.INFO)
    

    
with open(filename, 'r') as ipFile:
    csv = [line.strip().split(',') for line in ipFile.readlines()][1:]
    
N = len(csv)

locs=np.zeros((N, 2))
y_obs = np.zeros((N, 1))
y = np.zeros((N, 1))

    
for idx, line in enumerate(csv):

    locs[idx, 1] = float(line[0])
    locs[idx, 0] = float(line[1])
    y_obs[idx, 0] = float(line[2]) if line[2]!='NA' else np.NAN
    if len(line)>3:
        y[idx, 0] = float(line[3])

locs[:,0] = locs[:,0] - np.min(locs[:,0])
locs[:,1] = locs[:,1] - np.min(locs[:,1])
        
locs[:,0] = locs[:,0]/np.max(locs[:,0])
locs[:,1] = locs[:,1]/np.max(locs[:,1])

Nx = len(np.unique(locs[:,0])); Ny = len(np.unique(locs[:,1]))
    
obs_mean = np.nanmean(y_obs)
y_obs = y_obs - obs_mean
y = y - obs_mean
    
obs_inds = np.isfinite(y_obs).flatten()
R = me_scale
y_disp = y_obs.reshape((Nx, Ny))



##### parameter optimization #####

def MRALikelihood(params):

    kappa, sigma, me_scale = params
    cov = lambda _locs1, _locs2: sigma*mt.Matern32(_locs1, _locs2, l=np.abs(kappa))
    mraTree = MRATree(locs, M, J, r0, critDepth, cov, y_obs, me_scale)        
    lik = mraTree.getLikelihood()
    return( lik )


def TrueLikelihood(kappa):

    cov = lambda _locs1: mt.ExpCovFun(_locs1, _locs1, l=np.abs(kappa))
    obs = y_obs[obs_inds].ravel()
    obs_mat = np.matrix(obs).T
    obs_locs = locs[obs_inds,:]
    varY = cov(obs_locs) + np.eye(len(obs))*me_scale

    
    sign, logdet = np.linalg.slogdet(varY)
    const = len(obs)*np.log(2*np.pi)
    quad_form = obs_mat.T*lng.inv(varY)*obs_mat

    hand_lik = logdet + quad_form
    
    full_hand_lik = -0.5*( const + hand_lik )

    
    model = mvn(mean=np.zeros(len(obs)), cov=varY)
    numpy_lik = model.logpdf(obs)

    #print("full_hand_lik: %f" % full_hand_lik)
    #print("numpy_lik: %f" % numpy_lik)
    
    #return(numpy_lik)
    return( logdet + quad_form )



#logging.info("true likelihood: %f" % TrueLikelihood(kappa))
#logging.info("MRA likelihood: %f" % MRALikelihood(kappa))


xMRA = opt.minimize(MRALikelihood, [0.1, 0.8, 0.6], method='nelder-mead', options={'xtol':1e-3, 'disp':False})
#xTrue = opt.minimize(TrueLikelihood, [kappa], method='nelder-mead', options={'xtol':1e-1, 'disp':False})

logging.info("kappa: %f,\n sigma: %f,\n me_scale: %f" % tuple(xMRA.x))
#logging.info("True estimate: %f" % xTrue.x[0])

        
