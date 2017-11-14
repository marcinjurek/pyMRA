#import weakref
import os
import pdb
import numpy as np
from scipy.spatial.distance import squareform, pdist, cdist
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as sts
import numpy.linalg as lng




# Calculates the MSE between true state and predictions
# def MSE(xPred, Sigma, xTrue=0):
#     T_max = xPred.shape[1]
#     N = xPred.shape[0]
#     err = np.zeros(T_max)
#     for t in range(T_max):
#         err[t] = (lng.norm(xPred[:,t] - xTrue[:,t])**2)/N
#     return err






def MSE(xPred, Sigma, xTrue=0):

    N = xPred.shape[0]
    err = np.sqrt((lng.norm(xPred - xTrue)**2)/N)
    
    return err









def KLdiv(mu0, mu1, Sig0, Sig1):
    
    n = Sig0.shape[0]

    Sig0c = np.linalg.cholesky(Sig0)
    Sig1c = np.linalg.cholesky(Sig1)

    M = np.linalg.solve(Sig1c, Sig0c)
    traceterm = sum(np.diag(M.T*M)) - n
    
    logdetterm = 2*(np.log(np.linalg.det(Sig1c))) - 2*(np.log(np.linalg.det(Sig0c)))

    meandiff = mu1 - mu0
    temp = np.linalg.solve(Sig1c, meandiff)
    meanterm = np.linalg.norm(temp,2)**2

    kldiv = 0.5*( traceterm + logdetterm + meanterm)
    return kldiv
    



# Calculates the normal log score
def logscore(obs, muPred, SigPred):


    y = np.array(obs).ravel()
    obs_inds = np.logical_not(np.isnan(y))

    yObs = y[obs_inds]
    muObs = np.array(muPred).ravel()[obs_inds]
    SigObs = np.array(SigPred)[np.ix_(obs_inds, obs_inds)]

    # W, V = np.linalg.eigh(SigPred)
    
    # lW = np.log(W)

    # x = (V.T * (xPred-xTrue))
    # N = len(x)

    #return -N/2 * np.sum(np.log(W)) - x.T * np.diag(1/W) * x
    return sts.multivariate_normal.logpdf(yObs, mean=muObs, cov=SigObs)




# calculates a matrix with distances between points of v1 and v2
def dist(v1, v2=None):

    if not np.any(v2):
        v2=v1
        
    D = np.zeros([len(v1), len(v2)])
    for idx_k, k in enumerate(v1):
        for idx_l, l in enumerate(v2):
            D[idx_k, idx_l] = np.abs(k-l)
    return(D)



def filterNNZ(X):
    newX = X
    indNNZ = np.where(X!=0)
    newX[indNNZ] = 1
    return newX



# displays a matrix
def dispMat(mat, title="", cmap=None, fName=None):

    if cmap:
        fig = plt.matshow(mat, cmap=cmap)
    else:
        fig = plt.matshow(mat)
    plt.colorbar()
    if title:
        plt.title(title)
    if fName:
        plt.savefig(fName, dpi=400)
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
        

    plt.show()




    
    
# generates locations on a given domain
def genLocations(NGrid, lb=0, ub=1, random=False):

    if random:
        locs_all = np.random.uniform(lb, ub, NGrid)
    else:
        locs_all = np.linspace(lb, ub, num=NGrid+1)[1:]

    return(locs_all.reshape((NGrid, 1)))




def genLocations2d(Nx, lbx=0, ubx=1, Ny=0, lby=0, uby=1):

    if not Ny:
        Ny=Nx
    
    xx, yy = np.meshgrid(np.linspace(lbx, ubx, num=Nx), np.linspace(lby, uby, num=Ny))
    X = np.hstack((xx.flatten().reshape(Nx*Ny, 1), yy.flatten().reshape(Nx*Ny,1)))

    return X



def genClusters(N, k):

    NperK = int(N/k)
    points = np.empty((0,2))
    
    for cluster in range(k):
        pts = np.random.normal(loc=np.random.uniform(size=2), scale=np.random.uniform(low=0.1, high=0.2), size=(NperK,2))
        points = np.vstack((points, pts))

    for r in range(N- k*NperK):
        pts = np.random.uniform(size=2)
        points = np.vstack((points, pts))

    return points



def ExpCovFunCirc(locs, l=1):

    locs = locs if np.ndim(locs)==2 else np.reshape(locs, [len(locs), 1])

    xv, yv = np.meshgrid(locs, locs)
    m = np.minimum(xv, yv)
    M = np.maximum(xv, yv)
    dist = np.minimum(M - m, m + 1-M)

    covMat = np.exp(-dist/l)
    return(covMat)

    



# a covariance function

def ExpCovFun(locs, locs2=np.array([]), l=1, circular=False):

    locs = locs if np.ndim(locs)==2 else np.reshape(locs, [len(locs), 1])

    if circular:
        xv, yv = np.meshgrid(locs, locs)
        m = np.minimum(xv, yv)
        M = np.maximum(xv, yv)
        dist = np.minimum(M - m, m + 1-M)
    else:
        if len(locs2):
            dist = np.matrix(cdist(locs, locs2))
        else:
            dist = np.matrix(squareform(pdist(locs)))
    covMat = np.exp(-dist/l)


    return(covMat)





def Matern32(locs, locs2=np.array([]), l=1, sig=1):

    if len(locs2):
        dist = cdist(locs, locs2)
    else:
        dist = squareform(pdist(locs))

    covMat = sig*(1 + np.sqrt(3)*dist/l) * np.exp(-np.sqrt(3)*dist/l)

    return(np.matrix(covMat))




# a tapering function
def KanterCovFun(locs, radius=1):

    locs = locs if np.ndim(locs)==2 else np.reshape(locs, [len(locs), 1])
    
    D = squareform(pdist(locs))/radius
    piD2 = 2*np.pi*D
    R = (1 - D)*np.sin(piD2)/piD2 + 1/np.pi * (1-np.cos(piD2))/piD2
    R[D>1]=0
    R[D==0]=1
    return R




def simulate1D(locs, CFun, mean=None, seed=None, domainIsCircular=False):
    """
    Simulate from a 1-D random Gaussian field

    Parameters
    ----------
    locs : numpy.array
      1-D array with locations at which to simulate the process
    CFun : function
      covariance function to be used for simulation
    mean : function
      mean function of the process
    seed : int
      seed for numpy.random

    Returns
    -------
    y : numpy.array
      values of the process at the given locations
    """
    
    NGrid = len(locs)
    if seed:
        np.random.seed(seed)
    trueCovMat = CFun(locs, circular=True)
    covChol = lng.cholesky(trueCovMat)

    if not mean:
        mean = np.zeros([NGrid,1])
    else:
        assert(np.shape(mean)[0]==NGrid and np.shape(mean)[1]==1)
    
    innov = np.reshape(np.random.normal(size=NGrid), [NGrid, 1]) + mean
    y = np.dot(covChol, innov)
    return(y)
    
