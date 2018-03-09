import logging
import numpy as np
import sys
sys.path.append('..')
import MRATools as mt



if __name__=='__main__':



    dim_x = 100

    locs = np.matrix(np.linspace(0, 1, dim_x).reshape((dim_x, 1)))

    # specify the bandwith
    b = 10
    h = locs[1] - locs[0]
    radius = float(b*h)
    print("radius=%f" % radius)

    
    Sig = mt.KanterCovFun(locs, radius=radius, circular=False)
    mt.dispMat(Sig)
    nnz = np.count_nonzero(Sig[0,:])

