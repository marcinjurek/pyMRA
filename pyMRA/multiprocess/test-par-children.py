import numpy as np
import sys
import pdb
import multiprocessing as mp
from mpi4py import MPI


class Node(object):

    
    def __init__(self, ID, NCh, levelsFromLeaves):

        self.ID = ID
        self.children = []
        self.depth = len(self.ID)-1

        dim = 300
        self.mat = np.random.normal(size=dim*dim).reshape((dim, dim))
        self.inv = np.linalg.inv(self.mat)

        
        if levelsFromLeaves:

            if self.depth==2:
                self.pool = mp.Pool(processes=NCh)
                self.children = self.pool.starmap(Node, [(self.ID + str(j+1), NCh, levelsFromLeaves-1) for j in range(NCh)])

            else:
                for j in range(NCh):
                    chID = self.ID + str(j+1)
                    ch = Node(chID, NCh, levelsFromLeaves-1)
                    self.children.append(ch)
        


        

def times(x, y):

    return x*y

                    



if __name__=='__main__':

    M = 4
    J = 2

    root = Node('r', J, M)

    pdb.set_trace()
    sys.exit()
    
    pool = mp.Pool(processes=4)
    A = list([(i, i+1) for i in range(10)])
    results = pool.starmap(times, A)
    
    print(results)
