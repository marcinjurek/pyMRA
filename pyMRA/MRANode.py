import gc
import numpy_indexed as npi
import multiprocessing as mp
import logging
import pdb
#import sys
import scipy.linalg as lng
import numpy as np
from sklearn.cluster import KMeans
from scipy.spatial.distance import *
import scipy
#sys.path.append('..')
#import MRA.MRATools as mt
from  pyMRA.MRA import MRATools as mt

class Node(object):



    #@profile
    def __init__(self, parent, ID, locs, notKnots, levelsFromLeaves, J, r, critDepth, cov, obs, R, pipe=None):

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
            self.knots, self.kInds = self._getKnotsInds(r, notKnots, random=False)
        else: 
            self.knots = notKnots
            try:
                self.kInds = np.arange(self.N)[np.flatnonzero(npi.contains(notKnots, self.locs))]
            except:
                logging.critical("CRITICAL: problem with knots")
                pdb.set_trace()
           

        self.calculatePrior(cov)
        
        if not self.leaf and len(notKnots)>max(r,J):
            
            newNotKnots = notKnots[~npi.contains(self.knots, notKnots),:]
            #minJ = min(J, len(newNotKnots)) # if there is fewer "spare" locations than splits, then make as many splits as there are locations
            splits = self._getJSplits(J)#, self.locs)
            #splits = self._getSplits()
            #splits = self._getSplitsRobust(newNotKnots)


            
            if self.res==self.critDepth:
                pipes=[]; procs=[]

            for j in range(len(splits)):
                chID = self.ID + str(j+1)
                self.inds[chID] = splits[j]
                
                chLocs = locs[self.inds[chID]]

                if isinstance(cov, np.matrix):
                    C = cov[np.ix_(self.inds[chID], self.kInds)]
                    covCh = cov[np.ix_(self.inds[chID], self.inds[chID])] - C * self.k * C.T
                elif len(locs)<=1e4: #if the child covariance matrix will not be too large, calculate it directly
                    C = cov(chLocs, self.knots)
                    covCh = cov(chLocs, chLocs) - C * self.k * C.T
                else:
                    covCh = lambda _locs1, _locs2: cov(_locs1, _locs2) - cov(_locs1, self.knots) * self.k * cov(self.knots, _locs2)

                chObs = obs[self.inds[chID]]
                chNotKnots = chLocs[npi.contains(newNotKnots, chLocs)]
                
                if isinstance(R, float):
                    chR = R
                else:
                    chR = R[np.ix_(self.inds[chID], self.inds[chID])]
                
                if self.res==self.critDepth:
                    parPipe, chPipe = mp.Pipe()
                    chP = mp.Process(target=Node, args=(self, chID, chLocs, chNotKnots, levelsFromLeaves-1, J, r, critDepth, covCh, chObs, chR, chPipe))
                    chP.start()
                    pipes.append(parPipe)
                    procs.append(chP)
                else:
                    newChild = Node(self, chID, chLocs, chNotKnots, levelsFromLeaves-1, J, r, critDepth, covCh, chObs, chR)
                    self.children.append( newChild )
                   

                    
        if self.res==self.critDepth:
            for idx, pipe in enumerate(pipes):
                self.children.append( pipe.recv() )
                procs[idx].join()

        self.calculatePosterior(obs, R)
        
        #while self.children:
        #    ch = self.children.pop()
        #    del ch
        #gc.collect()

        
        if self.res==(self.critDepth+1):
            pipe.send(self)

        

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
            # if leafID=='r12' and self.ID=='r':
            #     pdb.set_trace()
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

        # leavesLocs = sorted([(leaf.ID, leaf.locs) for leaf in leaves], key=lambda _t: _t[0])        
        # findIntersect = lambda _locs: np.arange(self.N)[np.in1d(self.locs, _locs)]
        # findIntersectTup = lambda _t: (_t[0], findIntersect([1]))
        # leavesLocsInds = list(map(findIntersectTup, leavesLocs))
        
        # if set(self.ID) <= set(['1', 'r']):
        #      print('===resolution %d===' % self.res)
        # print(self.ID)
        # print(self.locs[order])
        
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
            knots = [np.percentile(locs1d, 100.0*i/(r+1), interpolation='nearest') for i in range(r+2)][1:-1]
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


    

    


    def _getJSplits(self, J):
        
        subDomains = []
        kmeans = KMeans(n_clusters=J, random_state=0).fit(self.locs)
        all_labels=kmeans.labels_
        for j in range(J):
            inds = np.where(all_labels==j)[0]
            if len(inds):
                subDomains.append(inds)
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
        logging.debug("Node %s: calculate prior" % self.ID)
        if isinstance(cov, np.matrix):
            self.B = cov[:,self.kInds]
        else:
            self.B = cov(self.locs, self.knots)
        self.kInv = self.B[self.kInds,:]
        self.k = np.linalg.inv(self.kInv)
        self.kC = np.linalg.cholesky(self.k)

        #self.B = self.B[order,:]
        
        logging.debug("Node %s: finished calculating prior" % self.ID)






    #@profile
    def calculatePosterior(self, obs, R):
        logging.debug("Node %s: start calculating posterior" % self.ID)


        
        # calculate A and omega
        self.A = []
        omg = []            
        if self.leaf:

            finiteObsInds = np.isfinite(obs).ravel()
            finiteObs = obs[finiteObsInds,:]
            if isinstance(R, float):
                finiteR = np.matrix(R*np.eye(len(finiteObs)))
            else:
                finiteR = R[np.ix_(finiteObsInds, finiteObsInds)]
            invFinR = np.matrix(np.linalg.inv(finiteR))
            
            B_lf = [ self._getB_lk(k)[finiteObsInds,:] for k in range(self.res+1) ]
            
            for k in range(self.res+1):
                self.A.append( [] )
                if np.any(finiteObsInds):
                    omg.append( B_lf[k].T * invFinR * finiteObs )
                else:
                    omg.append( np.matrix(np.zeros((B_lf[k].shape[1],1))) )
                for l in range(self.res+1):
                    if np.any(finiteObsInds):
                        self.A[k].append( B_lf[k].T * invFinR * B_lf[l] )
                    else:
                        self.A[k].append( np.matrix(np.zeros((B_lf[k].shape[1],B_lf[l].shape[1]))) )
                    
        else:           
                
            for k in range(self.res+1):    
                subList=[ch.omgTil[k] for ch in self.children]
                omg.append( sum(subList) )
                self.A.append([])
                for l in range(self.res+1):
                    subList = [ch.ATil[k][l] for ch in self.children]
                    self.A[k].append( sum(subList) )                    

                    
        self.kTil = np.matrix(lng.inv(self.kInv + self.A[self.res][self.res]))
        W, V = np.linalg.eigh(self.kTil)
        negInd = np.where(W<0)[0]
        W[negInd] = -W[negInd]
        self.kTilC = np.matrix(V) * np.matrix(np.diag(np.sqrt(W)))


        if True:#likelihood:
            # calculate likelihood
            if self.leaf:
                if np.any(finiteObsInds):
                    self.u = float(finiteObs.T * np.linalg.inv(finiteR) * finiteObs)
                    sign, self.d = np.linalg.slogdet(finiteR)
                else:
                    self.u = 0
                    self.d = 0
            else:
                kTilInv = np.linalg.inv(self.kTil)
                assert np.abs(np.max(kTilInv * self.kTil - np.eye(len(self.knots))))<1e-10, "something went wrong with iversion in likelihood evaluation"
                self.d = np.log(np.linalg.det(kTilInv)) - np.log(np.linalg.det(self.kInv))
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
        

        # calculate the moments
        self.mean = self.BTil[self.res] * self.kTil * omg[self.res]
        self.var = np.linalg.norm(self.BTil[self.res] * self.kTilC, axis=1)**2 # we calculate only the diagonal to save memory

      
        # add contribution from children
        for ch in self.children:
            chInds = self.inds[ch.ID]
            self.mean[chInds,:] += ch.mean
            self.var[chInds] += ch.var

        logging.debug("Node %s: finished calculating posterior" % self.ID)





        
