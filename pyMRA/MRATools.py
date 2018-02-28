import os
import pdb
import numpy as np
from sklearn.gaussian_process.kernels import Matern as skMatern
from scipy.spatial.distance import squareform, pdist, cdist
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as sts
import numpy.linalg as lng





def get_layout(m,J,r):

    total = r*(J**m)
    if total<=2:
        tup = (1,2)
    elif total<=3:
        tup = (1,3)
    elif total<=4:
        tup = (1,4)
    elif total<=6:
        tup = (2,3)
    elif total<=8:
        tup = (2,4)
    elif total<=9:
        tup = (3,3)
    elif total<=12:
        tup = (3,4)
    elif total<=15:
        tup = (3,5)
    elif total==16:
        tup = (4,4)
    elif total<=18:
        tup = (3,6)
    elif total<=20:
        tup = (4,5)
    elif total<=24:
        tup = (4,6)
    elif total==25:
        tup = (5,5)
    elif total==28:
        tup = (4,7)
    elif total<=30:
        tup = (5,6)
    elif total<=35:
        tup = (5,7)
    elif total==36:
        tup = (6,6)
    else:
        raise ValueError("Too many functions to plot")
    return tup




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
def dispMat(mat, title="", cmap=None, fName=None, vmin=None, vmax=None, colorbar=True):

    fig = plt.matshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)
    if colorbar:
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
        if len(locs2):
            xv, yv = np.meshgrid(locs, locs2)
        else:
            xv, yv = np.meshgrid(locs, locs)
        m = np.minimum(xv, yv)
        M = np.maximum(xv, yv)
        dist = np.minimum(M - m, m + 1-M).T
    else:
        if len(locs2):
            dist = np.matrix(cdist(locs, locs2))
        else:
            dist = np.matrix(squareform(pdist(locs)))
    covMat = np.exp(-dist/l)


    return(covMat)




def Matern(locs, l=1, sig=1, nu=1.5):

    K = skMatern(nu=nu, length_scale=l)
    covMat = K(locs)*sig
    return( np.matrix(covMat) )







def Matern32(locs, locs2=np.array([]), l=1, sig=1):

    if len(locs2):
        dist = cdist(locs, locs2)
    else:
        dist = squareform(pdist(locs))

    covMat = sig*(1 + np.sqrt(3)*dist/l) * np.exp(-np.sqrt(3)*dist/l)

    return(np.matrix(covMat))




# a tapering function
def KanterCovFun(locs, locs2=np.array([]), radius=1, cir=False):

    # if cir:
    #     D = 1
    # else:
    #     locs = locs if np.ndim(locs)==2 else np.reshape(locs, [len(locs), 1])
    #     D = squareform(pdist(locs))/radius

    locs = locs if np.ndim(locs)==2 else np.reshape(locs, [len(locs), 1])

    if cir:
        if len(locs2):
            xv, yv = np.meshgrid(locs, locs2)
        else:
            xv, yv = np.meshgrid(locs, locs)
        m = np.minimum(xv, yv)
        M = np.maximum(xv, yv)
        D = np.minimum(M - m, m + 1-M).T
    else:
        if len(locs2):
            D = np.matrix(cdist(locs, locs2))
        else:
            D = np.matrix(squareform(pdist(locs)))
    D = D/radius
    
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
    




def simulateGRF(locs, CFun, mean=None, seed=None, domainIsCircular=False):
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

    Nx = len(np.unique(locs[:,0]))
    if locs.shape[1]==1:
        Ny = 1
    else:
        Ny = len(np.unique(locs[:,1]))
    
    if seed:
        np.random.seed(seed)
    if isinstance(CFun, np.ndarray):
        covChol = np.matrix(lng.cholesky(CFun))
    else:
        trueCovMat = CFun(locs)
        covChol = np.matrix(lng.cholesky(trueCovMat))



    if not mean:
        mean = np.zeros((Nx*Ny, 1))
    
    x_raw = np.matrix(np.random.normal(size=(Nx*Ny, 1)) + mean)
    y = covChol * x_raw

    if Ny>1:
        y = y.reshape((Ny, Nx))

    return(y)
