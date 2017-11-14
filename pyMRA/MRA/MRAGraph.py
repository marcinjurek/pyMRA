from pyMRA.MRA.MRANode import Node
from functools import reduce
import scipy.linalg
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pdb

from pyMRA.MRA import MRATools as mt


class MRAGraph(object):


    def __init__(self, locs, M, J, r, critDepth, cov, obs, R):

        self.M = M
        self.J = J
        self.r = r
        self.locs = locs
        self.d = np.shape(locs)[1]
        if False:#isinstance(R, float):
            self.root = Node(None, 'r', locs, locs, M, J, r, critDepth, cov, obs, R*np.matrix(np.eye(len(locs))))
        else:
            self.root = Node(None, 'r', locs, locs, M, J, r, critDepth, cov, obs, R)
        self.obs_inds = np.where(np.logical_not(np.isnan(obs)))[0]
       
        #self.colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
        #self.colors = ['#377eb8','#4daf4a','#984ea3']
        #self.colors = ['#ffffff','#1f78b4','#b2df8a','#33a02c']
        self.colors = ['#000000','#1f78b4','#c71585','#33a02c']



        
        
    def getNodesBFS(self, groupByResolution=False):

        stack = [self.root]
        nodes = []
        while stack:
            n = stack.pop()
            
            if groupByResolution:
                if not nodes or nodes[-1][-1].res<n.res:
                    nodes.append([n])
                else:
                    nodes[-1].append(n)
            else:
                nodes.append(n)
                
            stack = n.children[::-1] + stack

        return nodes

    

    
    def getNodesDFS(self):
        
        stack = [self.root]
        nodes = []
        while stack:
            n = stack.pop()
            nodes.append(n)
            stack.extend(n.children[::-1])
            
        return nodes
    


    
    def drawKnots(self):

        fig = plt.figure()#figsize=plt.figaspect(0.2))
        nodes = self.getNodesBFS(groupByResolution=True)
        
        for m in range(self.M+1):

            if self.d==2:
                ax = fig.add_subplot(int(0.5*(self.M+1)), 2, m+1)
            else:
                ax = fig.add_subplot(self.M+1, 1, m+1)
                ax.set_ylim(-0.1, 2)

                
            # plot grid locations not used as knot before
            grid = [node.getGrid() for node in nodes[m]]

            for idx, region in enumerate(grid):
                col = self.colors[((idx+m) % len(self.colors))]
                if self.d==2:
                    ax.plot(region[:,0], region[:,1], marker='o', color=col, markersize='2', linestyle='None', label='grid')
                else:
                    ax.plot(region, np.zeros(len(region)), marker='o', color=col, markersize='5', linestyle='None', label='grid')
                    # if idx in [0, 1] and m>0:
                    #     dist = (self.locs[1]-self.locs[0])/2
                    #     ax.axvline(x=np.max(region)+dist, color='k', linestyle='dashed')

                        
            # plot knots
            knots = reduce(lambda a,b: np.vstack((a,b)), [node.locs[node.kInds] for node in nodes[m]])
            
            if self.d==2:
                ax.plot(knots[:,0], knots[:,1], marker='o', color='red', markersize='4', linestyle='None', label='knots')
            else:
                knots = ax.plot(knots[:,0], np.ones(len(knots)), marker='o', color='red', markersize='5', linestyle='None', label='knots')
                ax.get_yaxis().set_visible(False)
                ax.get_xaxis().set_visible(False)
            ax.set_title("resolution: %d" % m)

        plt.tight_layout()
        plt.show()      



    def drawBMatrix(self, distr="prior"):

        """
        Daws the B matrix. Its rows correspond to grid points and columns to
        hierarchically arranged basis functions
        """
        B = self.getBasisFunctionsMatrix(distr=distr)
        plt.matshow(B, cmap='Spectral')
        plt.colorbar()
        plt.tick_params(axis='both', which='both', bottom='off', top='off')
        plt.xlabel("%s basis functions" % distr, fontsize='x-large')
        plt.ylabel("grid points", fontsize='x-large')
        plt.show()





    def drawSparsityPat(self, distr="prior"):

        """
        draws the B matrix with ones whererver the corresponding element in the
        B matrix is nonzero and zeros otherwise.
        """
        if distr=="prior":
            B = self.getBasisFunctionsMatrix(distr="prior")
        elif distr=="posterior":
            B = self.getBasisFunctionsMatrix(distr="posterior")
        ind = np.where(np.abs(B)>1e-10)
        B1 = B; B1[ind] = 1
        mt.dispMat(B1, cmap='binary', title="%s sparsity pattern" % distr)




    def drawBasisFunctions(self):
        """
        plots basis functions
        ( so far for data in 1D only )
        """

        if np.shape(self.locs)[1]==1 and self.M>0:
            B = self.getBasisFunctionsMatrix(groupByResolution=True)
            fig = plt.figure()
            for m, Bm in enumerate(B[:-1]):
        
                ax = fig.add_subplot(self.M, 1, m+1)
                for col in range(Bm.shape[1]):
                    color = self.colors[col % len(self.colors)]
                    ax.plot(self.locs, Bm[:,col], color=color, linewidth=4)
                ax.tick_params(labelsize='x-large')
                if m<(len(B)-2):
                    ax.get_xaxis().set_visible(False)
                ax.set_title("resolution: %d" % m)#, fontsize='x-large')

            plt.tight_layout()
            plt.show()


        #if np.shape(self.locs)[2]==2 and self.M:
            
            


    def drawGridAndObs(self):

        """
        plot grid points and observation locations
        """
               
        if np.shape(self.locs)[1]==1:
            fig = plt.figure()
            ax = fig.add_subplot(111)

            ax.plot(self.locs, np.zeros(len(self.locs)), marker='o', color='black', markersize='3', linestyle='None', label="grid locations")
            ax.plot(self.locs[self.obs_inds], np.ones(len(self.obs_inds))*0.1, marker='o', color='red', markersize='4', linestyle='None', label="observations")
            ax.set_ylim([-.01, 0.2])
            ax.get_yaxis().set_visible(False)
            ax.legend(fontsize='large')
            plt.xticks(fontsize='large')
            plt.show()


            
        elif np.shape(self.locs)[1]==2 and self.M:

            if np.shape(self.locs)[1]>1:
                plt.scatter(self.locs[self.obs_inds,0], self.locs[self.obs_inds,1]);
            else:
                plt.scatter(self.locs[self.obs_inds,0], np.zeros(len(self.obs_inds)));
            plt.title("observation locations")
            plt.show()


            

        
    
    def getB_lk(self, callerID, k, l=None):
        """
        Queries node j_1, ... j_k for its covariance v_k with grid points in
        subregion j_1, ..., j_k, ..., j_l.

        Parameters
        ----------
        callerID : string
          ID of the node that asked for the information
        k : int
          resolution of the node on the caller path that contains the covariance
        l : int
          resolution of the node on the caller path with which the covariance
          should be calculated; l has to be between k and len(callerID)

        Returns
        -------
        Node:
          an instance of the Node class
        """

        kNode = self.getKNode(callerID, k)
        
        path = list(callerID[(k+1):] if not l else callerID[(k+1):(l+1)])[::-1]
        node = kNode
        inds = [np.arange(kNode.N)]
            
        while path:
            chId = node.children[int(path.pop())-1].ID
            inds.append( node.inds[chId] )
            node = node.children[int(chId[-1])-1]

        curInds = inds.pop()

        while inds:
            curInds = inds.pop()[curInds]
            
        return kNode.B[curInds,:]
    
        


    def getKNode(self, callerID, k):
        """
        Retrieves the k-th node on the path to the caller. For examples
        if caller ID is (j_1, j_2, j_3, j_4, j_5) and k=3 then the function
        returns (j_1, j_2, j_3). If k=0 then the root is returned. 
        This function should only be used by other nodes on the tree, not called
        from the outside.
        
        Parameters:
        -----------
        callerID : str
          ID of the node requesting the information
        k : int
          number of levels of the hierarchy down from the root

        Returns:
        --------
        node : Node
          the node with the requested ID
        """
        
        path = list(callerID[::-1])[:-1]
        node = self.root
        tempK = k
        
        while tempK>0:
            node = node.children[int(path.pop())-1]
            tempK = tempK-1
            
        return node




    
    def setPrior(self, xF, Sigma):

        self.root.calculatePrior(Sigma)




      

    def getBasisFunctionsMatrix(self, distr="prior", groupByResolution=False):

        nodes = self.getNodesBFS(groupByResolution=True)
        r = self.root
        
        if distr=="prior":
            B = r.B# * r.kC
        elif distr=="posterior":
            B = r.BTil[self.root.res]# * r.kTilC
        rootOrder = r.getOrderFromLeaves()
        B = B[rootOrder,:]
            
        if groupByResolution:
            B = [B]

            
        for m_nodes in nodes[1:]:

            # figure out the order in which the grid points in children will be sorted
            orders=[]
            for n in m_nodes:
                orders.append(n.getOrderFromLeaves())
            
            Jr = self.r*len(m_nodes)
            if distr=="prior":
                #block_list = [np.array(n.B*n.kC) for n in m_nodes]
                block_list = [np.array(n.B) for n in m_nodes]

            elif distr=="posterior":
                block_list = [np.array(n.BTil[n.res]) for n in m_nodes]
                #block_list = [np.array(n.BTil[n.res] * n.kTilC) for n in m_nodes]

            # reorder blocks
            for node_num in range(len(m_nodes)):
                block_list[node_num] = block_list[node_num][orders[node_num],:]

                
            Bm = scipy.linalg.block_diag(*block_list)
            if groupByResolution:
                B.append(Bm)
            else:
                B = np.hstack((B, Bm))
                
        return B



    

    def getLikelihood(self):

        return self.root.d + self.root.u


    
    

    def predict(self):

        xP = self.root.mean
        sdP = np.sqrt(self.root.var)
        return xP, sdP
