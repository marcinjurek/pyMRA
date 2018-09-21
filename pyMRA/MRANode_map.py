import gc
import numpy_indexed as npi
import multiprocessing as mp
import logging
import pdb
import scipy.linalg as lng
from numpy.linalg import slogdet
import numpy as np
from sklearn.cluster import KMeans
from scipy.spatial.distance import *
import scipy
from pyMRA import MRATools as mt
from pathos.multiprocessing import ProcessPool
#from memory_profiler import profile


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class Node(object):

    #@profile
    def __init__(self, parent, ID, locs, notKnots, levelsFromLeaves,
                 J, r, critDepth, cov, obs, R, pipe=None):

        self.parent = parent
        self.children = []
        self.ID = ID
        self.res = len(ID)-1
        self.N = len(locs)
        self.locs = locs
        self.inds = {}
        self.critDepth = critDepth
        self.leaf=(not bool(levelsFromLeaves))
        
        if not self.leaf and len(notKnots)>max(r,J):
            if len(notKnots)>1e2:
                self.knots, self.kInds = self._getKnotsInds(r, notKnots, random=True)
            else:
                self.knots, self.kInds = self._getKnotsInds(r, notKnots)
        else: 
            self.knots = notKnots
            self.leaf=True
            try:
                self.kInds = np.arange(self.N)[np.flatnonzero(npi.contains(notKnots, self.locs))]
            except:
                pdb.set_trace()

        self.calculatePrior(cov)
        
        if not self.leaf and len(notKnots)>max(r,J):
            
            newNotKnots = notKnots[~npi.contains(self.knots, notKnots),:]
            # if there is fewer "spare" locations than splits, then make as many splits as there are locations
            minJ = min(J, len(newNotKnots)) 
            if self.N>1e2:
                splits = self._getSplits()
            else:
                splits = self._getJSplits(minJ, newNotKnots)
            #splits = self._getSplitsRobust(newNotKnots)


            
            if self.res==self.critDepth:
                pipes=[]; procs=[]


            NCh = len(splits)
            data = {'chID' : [],
                    'chLocs' : [],
                    'chNotKnots' : [],
                    'lev' : [levelsFromLeaves-1]*NCh,
                    'J' : [J]*NCh,
                    'r' : [r]*NCh,
                    'critDepth' : [critDepth]*NCh,
                    'covCh' : [],
                    'chObs' : [],
                    'chR' : []}


            for j in range(NCh):
                chID = self.ID + str(j+1); data['chID'] += [chID]
                self.inds[chID] = splits[j]
                
                chLocs = locs[self.inds[chID]]; data['chLocs'] += [chLocs]

                if isinstance(cov, np.matrix):
                    C = cov[np.ix_(self.inds[chID], self.kInds)]
                    covCh = cov[np.ix_(self.inds[chID], self.inds[chID])] - C * self.k * C.T
                elif len(locs)<=1e4: #if the child covariance matrix will not be too large, calculate it directly
                    C = cov(chLocs, self.knots)
                    covCh = cov(chLocs, chLocs) - C * self.k * C.T
                else:
                    covCh = lambda _locs1, _locs2: cov(_locs1, _locs2) - cov(_locs1, self.knots) * self.k * cov(self.knots, _locs2)
                data['covCh'] += [covCh]
                chObs = obs[self.inds[chID]]; data['chObs'] += [chObs]
                chNotKnots = chLocs[npi.contains(newNotKnots, chLocs)]; data['chNotKnots'] += [chNotKnots]
                
                if isinstance(R, float):
                    chR = R
                else:
                    chR = R[np.ix_(self.inds[chID], self.inds[chID])]
                data['chR'] += [chR]

                if self.res != self.critDepth:
                    newChild = Node(self, chID, chLocs, chNotKnots, levelsFromLeaves-1, J, r, critDepth, covCh, chObs, chR)
                    self.children.append( newChild )
                   
        if self.res==self.critDepth:
            pp = ProcessPool(nodes=NCh)
            self.children = pp.map(Node, [self]*NCh, data['chID'], data['chLocs'], data['chNotKnots'], data['lev'], data['J'], data['r'], data['critDepth'], data['covCh'], data['chObs'], data['chR'])


        self.calculatePosterior(obs, R)
        
        while self.children:
            ch = self.children.pop()
            del ch
        gc.collect()

        
        

    def getLeaves(self):

        leaves = []
        if self.leaf:
            leaves.append(self)
        else:
            for ch in self.children:
                leaves.extend( ch.getLeaves() )

        return leaves

            


    
    def getOrderFromLeaves(self):

        leaves = self.getLeaves()

        
        def getLeafInds(leafID):
            path = list(leafID[(self.res+1):])[::-1]
            node = self
            inds = [np.arange(self.N)]
            
            while path:
                chId = node.children[int(path.pop())-1].ID
                inds.append( node.inds[chId] )
                node = node.children[int(chId[-1])-1]

            curInds = inds.pop()

            while inds:
                curInds = inds.pop()[curInds]

            return curInds

        leavesInds = [(leaf.ID, getLeafInds(leaf.ID)) for leaf in leaves]
        order = np.concatenate([tup[1] for tup in leavesInds])
        
        return order


            


    def getGrid(self):
        grid = self.knots
        for ch in self.children:
            chGrid = ch.getGrid()
            grid = np.vstack((grid, chGrid))
        return grid




            


    #@profile
    def _getKnotsInds(self, r, notKnots, random=False):

        if np.shape(self.locs)[1]==1:
            locs1d = self.locs.ravel()
            #knots = [np.percentile(locs1d, 100.0*i/(r+1), interpolation='nearest') for i in range(r+2)][1:-1]
            #knots = [np.percentile(notKnots1d, 100.0*i/(r+1), interpolation='nearest') for i in range(r+2)][1:-1]
            notKnots1d = notKnots.ravel()
            knots = [np.percentile(notKnots1d, 100.0*i/(r+1), interpolation='nearest') for i in range(r+2)][1:-1]
            kInds = np.arange(self.N)[np.flatnonzero(npi.contains(knots, locs1d))]
            

        else:
            if random and self.res>=0:
                inds = np.random.choice(np.arange(len(notKnots)), size=r, replace=False)
                knots = notKnots[inds]
            else:
                kmeans = KMeans(n_clusters=r, random_state=0).fit(notKnots)
                C = kmeans.cluster_centers_
                D = cdist(notKnots, C)
                knots = np.zeros((r, self.locs.shape[1]))
                for centr in range(r):
                    ind = np.argmin(D[:,centr])
                    knots[centr,:] = notKnots[ind]

            kInds = np.arange(self.N)[np.flatnonzero(npi.contains(knots, self.locs))]
        knots = self.locs[kInds] # this is to ensure that both kInds and knots are in the same order
        return knots, kInds
                                          

    
    
    

    #@profile
    def _getSplits(self):

        """
        If locations are 1D, splits into three segments with the same number of grid points.
        Otherwise splits into four segments, splitting in half along each axis.
        """

        if np.shape(self.locs)[1]==1:

            perc = np.percentile(self.locs, (33,66))
            
            locs_0 = np.where(self.locs[:,0]<perc[0])[0]
            locs_1 = np.where(np.logical_and(self.locs[:,0]>perc[0], self.locs[:,0]<perc[1]))[0]
            locs_2 = np.where(self.locs[:,0]>perc[1])[0]

            subdomains = [locs_0, locs_1, locs_2]

        else:
        
            means = np.mean(self.locs, axis=0)
        
            locs_00 = np.where(np.logical_and(self.locs[:,0]<=means[0], self.locs[:,1]<=means[1]))[0]
            locs_01 = np.where(np.logical_and(self.locs[:,0]<=means[0], self.locs[:,1]>means[1]))[0]
            locs_10 = np.where(np.logical_and(self.locs[:,0]>means[0], self.locs[:,1]<=means[1]))[0]
            locs_11 = np.where(np.logical_and(self.locs[:,0]>means[0], self.locs[:,1]>means[1]))[0]
            
            subdomains = [locs_00, locs_01, locs_10, locs_11]

        
        return subdomains






    def _getSplitsRobust(self, notKnots):

        """
        Works like _getSplits but accounts for the fact that there might be few grid points along
        a certain dimension. It first splits the points that are not the "middle" points (along x, y
        or both). Then it splits the "ties" into 4 parts. This is useful when we are dealing with, say
        a 3x3 grid and want to split it into four subregions. If we didn't use the robust method
        we might end up with one region that has 4 points, two regions with 2 points and one region with only
        a single point. 

        Useful only in 2D
        """
        med = np.median(notKnots, axis=0)

        all_ties = notKnots[np.where(np.logical_and(notKnots[:,0]==med[0], notKnots[:,1]==med[1]))]
        ties_inds = np.arange(self.N)[np.flatnonzero(npi.contains(all_ties, self.locs))]
        
        xy_ties = np.where(np.logical_and(notKnots[:,0]==med[0], notKnots[:,1]==med[1]))[0].ravel()
        x_ties = np.where(notKnots[:,0]==med[0])[0]
        y_ties = np.where(notKnots[:,1]==med[1])[0]
        
        pure_x_ties = np.setdiff1d(x_ties, xy_ties)
        pure_y_ties = np.setdiff1d(y_ties, xy_ties)
        array_xy_ties = np.array_split(xy_ties, 4)
                           
        locs_list = [np.array_split(xy_ties, 4), np.array_split(pure_x_ties, 4),  np.array_split(pure_y_ties, 4)[::-1]]
        locs = []
        for i in range(4):
            locs.append( np.hstack((locs_list[0][i], locs_list[1][i], locs_list[2][i])) )       
      
        locs[0] = np.setdiff1d(np.where(np.logical_and(self.locs[:,0] < med[0],  self.locs[:,1] < med[1]))[0], ties_inds)
        locs[1] = np.setdiff1d(np.where(np.logical_and(self.locs[:,0] < med[0],  self.locs[:,1] > med[1]))[0], ties_inds)
        locs[2] = np.setdiff1d(np.where(np.logical_and(self.locs[:,0] > med[0],  self.locs[:,1] < med[1]))[0], ties_inds)
        locs[3] = np.setdiff1d(np.where(np.logical_and(self.locs[:,0] > med[0],  self.locs[:,1] > med[1]))[0], ties_inds)       
            
        return locs


    

    def _getJSplits(self, J, notKnots):

        # if:
        #   * J=r+1, or if the number of knots is one less that the number
        #     of partitions,
        #   * we are in 1d,
        #   * there is enough grid points left,
        #   then we partition such that the knots are at the boundary;
        # else:
        #    do k-means etc.

        r = len(self.kInds)
        cond1 = J==r+1
        cond2 = self.locs.shape[1]==1
        cond3 = self.N>=(J+r)
        
        if cond1 and cond2 and cond3:
            subDomains = np.split(np.arange(self.N), self.kInds)
        else:
            subDomains = []

            # we need this to match the notKnots to their original indices
            allInds = np.arange(self.N)
            notKnotInds = np.arange(self.N)[np.flatnonzero(npi.contains(notKnots, self.locs))]

            #clustering
            nClusters = min(J, len(notKnots))
            kmeans = KMeans(n_clusters=nClusters, random_state=0).fit(notKnots)
            all_labels=kmeans.labels_

            #now assign knots from all previous resolutions to the nearby clusters
            centers = kmeans.cluster_centers_
            allKnotsInds = np.setdiff1d(np.arange(self.N), notKnotInds)
            allKnots = self.locs[ allKnotsInds, : ]
            D = cdist( allKnots, centers )
            knot_labels = np.argmin(D, axis=1)
            
            for j in range(J):
                indsInNotKnots = np.where(all_labels==j)[0]
                inds = notKnotInds[indsInNotKnots]
                # check if we should include a knot in this cluster
                
                kIndsInThisRegion = allKnotsInds[ np.where(knot_labels==j)[0] ]
                inds = np.hstack( (kIndsInThisRegion, inds) )
                inds = np.sort(inds)
                if len(inds):
                    subDomains.append(inds)

            if self.locs.shape[1]==1:
                subDomains = sorted(subDomains, key=lambda _arr: np.min(_arr))
                
        return subDomains





    def _getB_lk(self, l):
 
        node = self
        inds = np.arange(self.N)
        while node.res>l:
            chID = node.ID
            node = node.parent
            inds = node.inds[chID][inds]
 
        return node.B[inds,:]



    

    def _getRowOrder(self):

        
        if True:#not self.children:
            return np.arange(len(self.locs))
        
        allInds = np.array([])
        for ch in self.children:
            allInds = np.hstack((allInds, self.inds[ch.ID]))

        order = np.argsort(allInds)
        return order
    


        
    #@profile
    def calculatePrior(self, cov):
        #order = self._getRowOrder()
        logger.debug("Node %s: calculate prior" % self.ID)
        if isinstance(cov, np.matrix):
            self.B = cov[:,self.kInds]
        else:
            self.B = cov(self.locs, self.knots)
        self.kInv = self.B[self.kInds,:]
        try:
            self.k = np.linalg.inv(self.kInv)
        except:
            logger.critical("Problem with the knots!")
            pdb.set_trace()    
        self.kC = np.linalg.cholesky(self.k)

        #self.B = self.B[order,:]
        
        logger.debug("Node %s: finished calculating prior" % self.ID)






    #@profile
    def calculatePosterior(self, obs, R):
        logger.debug("Node %s: start calculating posterior" % self.ID)

        
        # calculate A and omega
        self.A = []
        omg = []

        if self.leaf:


            ################     proposed change to posterior calculations     ################
            obsInds = np.array(np.isfinite(obs)).ravel()
            
            H = np.matrix(np.eye(len(obs)))[obsInds,:]
            
            Rmat = np.matrix(R*np.eye(sum(obsInds)))
                
            HRinvH = (1/R)*H.T*H
            HRinvObs = H.T*(1/R)*obs[obsInds]

            B_lk = [ self._getB_lk(k) for k in range(self.res+1) ]
            
            for k in range(self.res+1):
                self.A.append( [] )
                omg.append( B_lk[k].T * HRinvObs )
                for l in range(self.res+1):
                    self.A[k].append( B_lk[k].T * HRinvH * B_lk[l] )
            
        else:           
                
            for k in range(self.res+1):    
                subList=[ch.omgTil[k] for ch in self.children]
                omg.append( sum(subList) )
                self.A.append([])
                for l in range(self.res+1):
                    subList = [ch.ATil[k][l] for ch in self.children]
                    self.A[k].append( sum(subList) )                    

   
        
        self.kTil = np.matrix(lng.inv(self.kInv + self.A[self.res][self.res]))
        self.kTilInv = lng.inv(self.kTil)

        
        
        #likelihood:
        if self.leaf:

            self.u = -omg[self.res].T * self.kTil * omg[self.res] + obs[obsInds].T * (1/R) * obs[obsInds]
            if np.max(np.shape(self.u))>1:
                pdb.set_trace()
            
            sgnTil, logdetTil = slogdet(self.kTilInv)
            sgn, logdet = slogdet(self.kInv)
            sgnR, logdetR = slogdet(Rmat)
            self.d = logdetTil - logdet + logdetR


        else:
            self.d = -np.log(np.linalg.det(self.kTil)) - np.log(np.linalg.det(self.kInv))
            self.u = - omg[self.res].T * self.kTil * omg[self.res]

            for ch in self.children:
                self.d += ch.d
                self.u += ch.u


        

        # calculate A-tilde and omega-tilde
        self.ATil = []
        self.omgTil = []
        for k in range(self.res):
            self.ATil.append( [] )
            self.omgTil.append( omg[k] - self.A[k][self.res] * self.kTil * omg[self.res] )
            for l in range(self.res+1):
                self.ATil[k].append( self.A[k][l] - self.A[k][self.res] * self.kTil * self.A[self.res][l] )



                
        # calculate B-tilde
        order = self._getRowOrder()
        self.BTil = []
        if self.leaf:
             self.BTil = [self._getB_lk(k) for k in range(self.res+1)]
        else:
            for k in range(self.res+1):
                self.BTil.append( np.matrix(np.zeros((self.N, len(self.kInds)))) )
                for ch in self.children:
                    chInds = self.inds[ch.ID]
                    self.BTil[k][chInds,:] = ch.BTil[k] - ch.BTil[ch.res] * ch.kTil * ch.A[ch.res][k]

      


        
                    

                
        W, V = np.linalg.eigh(self.kTil)
        negInd = np.where(W<0)[0]
        W[negInd] = -W[negInd]
        self.kTilC = np.matrix(V) * np.matrix(np.diag(np.sqrt(W)))

        # calculate the moments
        self.mean = self.BTil[self.res] * self.kTil * omg[self.res]
        self.var = np.linalg.norm(self.BTil[self.res] * self.kTilC, axis=1)**2 # we calculate only the diagonal to save memory


        
      
        # add contribution from children
        for ch in self.children:
            chInds = self.inds[ch.ID]
            self.mean[chInds,:] += ch.mean
            self.var[chInds] += ch.var

            
        logger.debug("Node %s: finished calculating posterior" % self.ID)





        
