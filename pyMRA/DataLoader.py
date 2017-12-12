import pyMRA
import os
import numpy as np

def load_data(size="small", include_truth=False):

    if size not in ["small", "large"]:
        raise ValueError("size has to be 'small' or 'large'")

    packagedir = pyMRA.__path__[0]
    dirname = os.path.join(os.path.dirname(packagedir), 'pyMRA', 'data', size)
    y_obs = np.load(os.path.join(dirname, 'y_obs.npy'))
    locs = np.load(os.path.join(dirname, 'locs.npy'))

    if include_truth:
        y = np.load(os.path.join(dirname, 'y.npy'))
        return( y, locs, y_obs )

    return( y, locs )
