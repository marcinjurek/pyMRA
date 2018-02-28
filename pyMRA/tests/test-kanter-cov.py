import numpy as np
import sys
sys.path.append('..')
import MRATools as mt



if __name__=='__main__':



    dim_x = 100

    locs = np.linspace(0, 1, 100)

    # specify the bandwith
    b = 10
    h = locs[1] - locs[0]
    radius = b*h
    
    Sig = mt.KanterCovFun(locs, radius=b*h, cir=True)
    mt.dispMat(Sig)
    nnz = np.count_nonzero(Sig[0,:])

