import pdb
import logging
import numpy as np
import sys
sys.path.append('..')
import MRATools as mt






def determine_radius(k, h):

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
        print('k=%d, app=%d, ind=%d, base=%d' % (k, intervals[ind-1], ind-1, (sf-1)/2.0))
        app_ind = ind-1
    else:
        print('k=%d, app=%d, ind=%d, base=%d' % (k, intervals[ind], ind, (sf-1)/2.0))
        app_ind = ind

    if app_ind==0:
        return h*base*np.sqrt(2) + h*0.01
    else:
        return h*np.sqrt((base+1)**2 + (app_ind-1)**2) + h*0.01
        




if __name__=='__main__':



    dim_x = 10
    dim_y = 10
    
    locs = mt.genLocations2d( Nx=dim_x, Ny=dim_y )

    for Nens in range(36,50):
        radius = determine_radius(Nens, 1.0/(dim_x-1))
        print(radius)
        Sig = mt.KanterCovFun(locs, radius=radius, circular=False)
        onerow = mt.filterNNZ(Sig[55,:].reshape((dim_x, dim_y)))
        mt.dispMat(onerow, title="Nens=%d" % Nens)


