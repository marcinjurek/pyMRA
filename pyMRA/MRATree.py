from pyMRA.MRANode import Node
from functools import reduce
import scipy.linalg
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import pdb

from pyMRA import MRATools as mt


class MRATree(object):

    #@profile
    def __init__(self, locs, M, J, r, critDepth, cov, obs, R):

        if J>1:
            if len(locs) < r*(1-J**(M))/(1-J):
                raise ValueError("Not enough grid points for the M, J, r you specified.")
        else:
            if len(locs) < r*(M+1):
                raise ValueError("Not enough grid points for the M, J, r you specified.")
            
        self.M = M
        self.J = J
        self.r = r
        self.locs = locs
        self.d = np.shape(locs)[1]
        obsM = np.matrix(obs) # make sure observations are a matrix; otherwise many operations will not work
        self.root = Node(None, 'r', locs, locs, M, J, r, critDepth, cov, obsM, R)
        self.obs_inds = np.where(np.logical_not(np.isnan(obs)))[0]





    def getLikelihood(self):

        return self.root.d + self.root.u


    
    

    def predict(self):

        xP = self.root.mean
        sdP = np.sqrt(self.root.var)
        return xP, sdP


        

        
        
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
      
    

    def avgLeafSize(self):

        nodes = self.getNodesBFS(groupByResolution=True)
        leaves = nodes[-1]
        sizes = np.array([ len(leaf.kInds) for leaf in leaves])
        return( np.mean(sizes) )



    def minLeaf(self):
    
        nodes = self.getNodesBFS(groupByResolution=True)
        leaves = nodes[-1]
        sizes = np.array([ len(leaf.kInds) for leaf in leaves])
        return( np.min(sizes) )
        

    def maxLeaf(self):
        nodes = self.getNodesBFS(groupByResolution=True)
        leaves = nodes[-1]
        sizes = np.array([ len(leaf.kInds) for leaf in leaves])
        return( np.max(sizes) )
        
    
    
    def drawKnots(self):

        fig = plt.figure()#figsize=plt.figaspect(0.2))
        nodes = self.getNodesBFS(groupByResolution=True)
        colors = ['#a6cee3','#b2df8a','#fb9a99','#ff7f00','#6a3d9a','#b15928']

        
        for m in range(self.M+1):

            if self.d==2:
                ax = fig.add_subplot(int(0.5*(self.M))+1, 2, m+1)
            else:
                ax = fig.add_subplot(self.M+1, 1, m+1)
                ax.set_ylim(-0.1, 2)

                
            # plot grid locations not used as knot before
            grid = [node.getGrid() for node in nodes[m]]

            for idx, region in enumerate(grid):
                col = colors[((idx+m) % len(colors))]
                if self.d==2:
                    ax.plot(region[:,0], region[:,1], marker='s', color=col, markersize='6', linestyle='None', label='grid')
                else:
                    ax.plot(region, np.zeros(len(region)), marker='s', color=col, markersize='6', linestyle='None', label='grid')
                        
            # plot knots
            knots = reduce(lambda a,b: np.vstack((a,b)), [node.locs[node.kInds] for node in nodes[m]])
            
            if self.d==2:
                ax.plot(knots[:,0], knots[:,1], marker='s', color='red', markersize='6', linestyle='None', label='knots')
            else:
                knots = ax.plot(knots[:,0], np.ones(len(knots)), marker='s', color='red', markersize='6', linestyle='None', label='knots')
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
        Brr = B1
        #mt.dispMat(B1, cmap='binary', title="%s sparsity pattern" % distr, colorbar=False)
        mt.dispMat(Brr, cmap='binary', title="%s sparsity pattern" % distr, colorbar=False)
        #ind = np.where(np.abs(B.T*B)>1e-10)
        #BB1 = B.T * B; BB1[ind] = 1
        #Brr = BB1[:,::-1][::-1,:]
        #mt.dispMat(Brr, cmap='binary', colorbar=False)




    def drawBasisFunctions(self, distr="prior", order="root"):
        """
        plots basis functions
        """
        cmaps = [plt.cm.Blues, plt.cm.copper, plt.cm.Blues, plt.cm.copper]#plt.cm.Greys, plt.cm.YlOrBr, plt.cm.Purples, plt.cm.RdPu]

        ### 1D version
        if np.shape(self.locs)[1]==1 and self.M>0:
           
            B = self.getBasisFunctionsMatrix(groupByResolution=True, order=order, distr=distr)
            nodes = self.getNodesBFS(groupByResolution=True)
            fig = plt.figure(figsize=(8,6))
            
            for m in range(self.M+1):

                Bm = B[m]
                cmap = cmaps[m];
                ncol = Bm.shape[1]; offset=0.3*ncol
                
                ax = fig.add_subplot(self.M+1, 1, m+1)
                ax.set_ylim(top=1.1)
                ax.set_xlim(np.min(self.locs), np.max(self.locs))
                
                # draw partition lines
                for node_idx, node in enumerate(nodes[m]):
                    if np.min(node.locs)>np.min(self.locs):
                        ax.axvline(x=node.locs[0], linestyle='dashed', color='k', linewidth=1)

                # draw functions
                for col in range(Bm.shape[1]):
                    color = cmap(((offset+col)/(ncol+offset)))
                    ax.plot(self.locs, Bm[:,col], color=color, linewidth=2)
                    
                ax.tick_params(labelsize='large')
                ax.plot(self.locs, np.zeros(len(self.locs)), color='black', linewidth=2)
                
                if m<(len(B)-2):
                    ax.get_xaxis().set_visible(False)
                ax.set_title("resolution: %d" % m, fontsize='x-large')

            
            plt.tight_layout()
            plt.show()

        ### 2D version
        if np.shape(self.locs)[1]==2 and self.M:

            dim_x = len(np.unique(self.locs[:,0]))
            dim_y = len(np.unique(self.locs[:,1]))

            B = self.getBasisFunctionsMatrix(groupByResolution=True)

            for m, Bm in enumerate(B):

                if Bm.shape[1]>36:
                    continue
                nrows, ncols = mt.get_layout(m, self.J, self.r)
                
                fig = plt.figure()
                grid = AxesGrid(fig, 111,
                                nrows_ncols=(nrows, ncols),
                                axes_pad=0.1,
                                share_all=True,
                                label_mode="L",
                                cbar_location="right",
                                cbar_mode="single")

                for func, ax in zip(Bm.T, grid):
                    im = ax.imshow(np.array(func).reshape((dim_x, dim_y)), vmax=1, vmin=-0.1, cmap="coolwarm")
                    ax.get_xaxis().set_visible(False)
                    ax.get_yaxis().set_visible(False)

                plt.suptitle("resolution: %d" % m, fontsize="x-large")

                grid.cbar_axes[0].colorbar(im)

                plt.show()
                    


        

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




      

    def getBasisFunctionsMatrix(self, distr="prior", groupByResolution=False, order='root', timesKC=False):

        nodes = self.getNodesBFS(groupByResolution=True)
        r = self.root
        
        if distr=="prior":
            B = r.B
            if timesKC:
                B = B * r.kC
        elif distr=="posterior":
            B = r.BTil[self.root.res]
            if timesKC:
                B = B * r.kTilC

        if order=="leaves":
            leafOrder = r.getOrderFromLeaves()
            B = B[leafOrder,:]
            
        if groupByResolution:
            B = [B]

            
        for m_nodes in nodes[1:]:

            # figure out the order in which the grid points in children will be sorted
            if order=='leaves':
                orders=[]
                for n in m_nodes:
                    orders.append(n.getOrderFromLeaves())
            if order=='root':
                smallest_locs = np.array([np.min(n.locs[:,0]) for n in m_nodes])
                orders = np.argsort(smallest_locs)

                    
            Jr = self.r*len(m_nodes)
            if distr=="prior":
                if timesKC:
                    block_list = [np.array(n.B*n.kC) for n in m_nodes]
                else:
                    block_list = [np.array(n.B) for n in m_nodes]

            elif distr=="posterior":
                if timesKC:
                    block_list = [np.array(n.BTil[n.res] * n.kTilC) for n in m_nodes]
                else:
                    block_list = [np.array(n.BTil[n.res]) for n in m_nodes]


            # reorder blocks
            if order=="leaves":
                for node_num in range(len(m_nodes)):
                    block_list[node_num] = block_list[node_num][orders[node_num],:]
            if order=='root':
                block_list = [block_list[i] for i in orders]

                    
                
            Bm = scipy.linalg.block_diag(*block_list)
            if groupByResolution:
                B.append(Bm)
            else:
                B = np.hstack((B, Bm))
                
        return B
