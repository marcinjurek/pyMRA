import os
import pdb
import scipy.sparse as sp
import numpy as np
from sklearn.gaussian_process.kernels import Matern as skMatern
from scipy.spatial.distance import squareform, pdist, cdist
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as sts
import numpy.linalg as lng
import logging


logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)


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


# scoring functions

def MSE(xPred, xTrue=0):

    N = xPred.shape[0]
    err = np.sqrt((lng.norm(xPred - xTrue)**2)/N)
    
    return err




def KLdivOld(mu0, mu1, Sig0, Sig1):
    
    n = Sig0.shape[0]

    Sig0c = np.linalg.cholesky(Sig0)
    Sig1c = np.linalg.cholesky(Sig1)

    M = np.linalg.solve(Sig1c, Sig0c)
    traceterm = sum(np.diag(M.T*M)) - n

    s, logdet1 = np.linalg.slogdet(Sig1c)
    s, logdet0 = np.linalg.slogdet(Sig0c)
    
    logdetterm = 2*(logdet1 - logdet0)

    meandiff = mu1 - mu0
    temp = np.linalg.solve(Sig1c, meandiff)
    meanterm = np.linalg.norm(temp,2)**2

    kldiv = 0.5*( traceterm + logdetterm + meanterm)
    return kldiv
    



def KLdiv(mu0, mu1, Sig0, Sig1):
    
    ndim = Sig0.shape[0]
    Sig1Inv = np.linalg.solve(Sig1, np.eye(ndim))
    traceterm = sum(np.diag(Sig1Inv * Sig0)) - ndim
    
    s, logdet1 = np.linalg.slogdet(Sig1)
    s, logdet0 = np.linalg.slogdet(Sig0)
    logdetterm = logdet1 - logdet0

    meandiff = mu1 - mu0
    sigMean = np.linalg.solve(Sig1, meandiff)
    meanterm = (meandiff.T * sigMean)[0,0]
    
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




def filterNNZ(X, tol=0.0):
    newX = np.zeros(X.shape)
    indNNZ = np.where(np.abs(X)>tol)
    newX[indNNZ] = 1
    return newX




# displays a matrix
def dispMat(mat, title="", cmap=None, fName=None, vmin=None, vmax=None, colorbar=True, pattern=False):

    if sp.issparse(mat):
        mat = mat.toarray()
    if pattern:
        mat = filterNNZ(mat)
        
    fig = plt.matshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
        
    if colorbar:
        plt.colorbar()
    if title:
        plt.title(title)
    if fName:
        plt.savefig(fName, dpi=400)
        

    plt.show()



    
    
# generates locations on a given domain
def genLocations(NGrid, lb=0, ub=1, random=False):

    if random:
        locs_all = np.random.uniform(lb, ub, NGrid)
    else:
        locs_all = np.linspace(lb, ub, num=NGrid+1)[1:]

    return(locs_all.reshape((NGrid, 1)))




def genLocations2d(Nx, lbx=0, ubx=1, Ny=0, lby=0, uby=1):
    """
    Generate locations on a 2d grid
    """
    
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







# calculates a matrix with distances between points of v1 and v2
def dist(locs, locs2=np.array([]), circular=False):

    locs = locs if np.ndim(locs)==2 else np.reshape(locs, [len(locs), 1])
    if circular:
        if len(locs2):
            xv, yv = np.meshgrid(locs, locs2)
        else:
            xv, yv = np.meshgrid(locs, locs)
        m = np.minimum(xv, yv)
        M = np.maximum(xv, yv)
        dist = np.matrix(np.minimum(M - m, m + 1-M).T)
    else:
        if len(locs2):
            dist = np.matrix(cdist(locs, locs2))
        else:
            dist = np.matrix(squareform(pdist(locs)))
    return dist







# covariance functions


def Iden(locs, locs2=np.array([]), l=1, circular=False):

    D = dist(locs, locs2, circular)
    inds = np.where(D==0)
    covMat = np.matrix(np.zeros(D.shape))
    covMat[inds] = 1
    return covMat


def ExpCovFun(locs, locs2=np.array([]), l=1, circular=False):

    D = dist(locs, locs2, circular)
    covMat = np.exp(-D/l)
    return(covMat)



def Matern(locs, l=1, sig=1, nu=1.5):

    K = skMatern(nu=nu, length_scale=l)
    covMat = K(locs)*sig
    return( np.matrix(covMat) )



def Matern52(locs, locs2=np.array([]), l=1, sig=1, circular=False):

    D = dist(locs, locs2, circular)
    covMat = sig*np.multiply(1 + np.sqrt(5)*D/l + (5/3)*np.square(D/l), np.exp(-np.sqrt(5)*D/l))
    return(np.matrix(covMat))

            

def Matern32(locs, locs2=np.array([]), l=1, sig=1, circular=False):
    
    D = dist(locs, locs2, circular)
    covMat = sig*np.multiply(1 + np.sqrt(3)*D/l, np.exp(-np.sqrt(3)*D/l))
    return(np.matrix(covMat))



def GaussianCovFun(locs, locs2=np.array([]), l=1, sig=1, circular=False):

    D = dist(locs, locs2, circular)
    covMat = sig*np.exp(-np.square(D)/(2*(l**2)))
    return(np.matrix(covMat))



def KanterCovFun(locs, locs2=np.array([]), radius=1.0, circular=False):

    if isinstance(radius, int):
        sort_x_locs = np.sort(np.unique(locs[:,0]))
        d = sort_x_locs[1] - sort_x_locs[0]

        if locs.shape[1]>1:
            ndim = len(np.unique(locs[:,1]))
        else:
            ndim = 1

        radius = determine_radius(radius, d, ndim=ndim)

    D = np.array(dist(locs, locs2, circular))
    D = D/radius
    piD2 = 2*np.pi*D
    R = (1 - D)*np.sin(piD2)/piD2 + 1/np.pi * (1-np.cos(piD2))/piD2
    R[D>1]=0
    R[D==0]=1
    return R




def determine_radius(k, h, ndim=2):
    """
    Determine the radius based on how many elements should be in the ensemble

    Parameters
    ----------
    k : int
      size of the ensemble
    h : float
      mesh diameter

    Returns
    -------
    radius : float
      tapering radius that ensures that the covariance matrix has approx. 
      k nonzero elements
    """

    if ndim==1:
        return( int(k/2)*h )
    

    
    if k==0:
        raise ValueError("Ensemble size must be stricly positive")
    s = np.floor(np.sqrt(k))

    if s % 2==0:
        sf = s-1
    else:
        sf = s

    if k==sf**2:
        return h*1.01*(sf-1)/2*np.sqrt(2)
        
    base = (sf-1)/2.0


    
    intervals = np.array([sf**2])
    while intervals[-1]<(sf+2)**2:
        if len(intervals)==1 or ((sf+2)**2 - intervals[-1]==4):
            intervals = np.append(intervals, intervals[-1] + 4)
        else:
            intervals = np.append(intervals, intervals[-1] + 8)

    ind = intervals.searchsorted(k)
    middle = (intervals[ind-1] + intervals[ind])/2.0

    if k<=middle:
        #print('k=%d, app=%d, ind=%d, base=%d' % (k, intervals[ind-1], ind-1, (sf-1)/2.0))
        app_ind = ind-1
    else:
        #print('k=%d, app=%d, ind=%d, base=%d' % (k, intervals[ind], ind, (sf-1)/2.0))
        app_ind = ind

    if app_ind==0:
        return h*base*np.sqrt(2) + h*0.01
    else:
        return h*np.sqrt((base+1)**2 + (app_ind-1)**2) + h*0.01






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
    




def simulateGRF(locs, CFun, mean=0, seed=None, domainIsCircular=False, reshape=True, CFunIsChol=False):
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
    if isinstance(CFun, np.ndarray) and CFunIsChol:
        covChol = CFun
    else:
        if isinstance(CFun, np.ndarray):
            covChol = np.matrix(lng.cholesky(CFun))
        else:
            trueCovMat = CFun(locs)
            covChol = np.matrix(lng.cholesky(trueCovMat))

    if np.ndim(mean)==1:
        mean = mean.reshape(len(mean),1)

    x_raw = np.matrix(np.random.normal(size=(len(locs), 1)))
    try:
        y = covChol * x_raw + mean
    except:
        pdb.set_trace()

    if Ny>1 and reshape:
        y = y.reshape((Ny, Nx))

    return(y)
