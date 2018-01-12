### pyMRA
Multi-resolution approximation for Gaussian processes

This package implements the multiresolution approximation algorithm described in [1]. It is based on the recursive application of the predictive process approximation [2]. The main functionalities include parameter estimation and prediction. See technical note for more details.

### Installation

Using pip is the easiest way to install the package:

```
pip install pyMRA
```
It is recommended to use Python 3.5.2.



### Examples

The first example shows how to load a sample data set supplied with the package and make predictions at unobserved locations.
```python
from pyMRA.MRATree import MRATree
import pyMRA.MRATools as mt
import pyMRA.DataLoader as dl

M=4; J=4; r0=4
me_scale=1e-4
critDepth = 0

y, locs, y_obs = dl.load_data("large", True)
    
Nx = y.shape[0]; Ny = y.shape[1]
y_obs = y_obs.reshape((Nx*Ny,1))
   
cov = lambda _locs1, _locs2: mt.ExpCovFun(_locs1, _locs2, l=2)
mraTree = MRATree(locs, M, J, r0, critDepth, cov, y_obs, me_scale)      
yP, sdP = mraTree.predict()
    
sdP = sdP.reshape((Nx, Ny))
yP = yP.reshape((Nx, Ny))
      
### compare results
mt.dispMat(yP, cmap="Spectral", title="prediction")
y = y.reshape((Nx, Ny), order='A')
mt.dispMat(y_obs.reshape((Nx, Ny)), cmap="Spectral", title="observations")
```




Now we demonstrate how to perform maximum likelihood estimation to obtain point estimates of unknown parameters
```python
import scipy.optimize as opt
import logging
import numpy as np
import scipy.linalg as lng
import pyMRA.MRATools as mt
from pyMRA.MRATree import MRATree

logging.basicConfig(format='%(asctime)s %(message)s', \
			datefmt='%H:%M:%S',level=logging.INFO)

###### set simulation parameters #####
frac_obs = 0.4 # the fraction of locations for which observations are available
dim_x = 100; dim_y = 1
sig = 1.0; me_scale=1e-1; kappa = 0.3

### MRA parameters
M=3; J=3; r0=2
critDepth = M+1

##### simulate data #####
# generate the underlying process
locs = mt.genLocations(dim_x)
Sig = sig*mt.Matern32(locs, l=kappa, sig=sig)
SigC = np.matrix(lng.cholesky(Sig))
x_raw = np.matrix(np.random.normal(size=(locs.shape[0],1)))
x = SigC.T * x_raw

# generate data
R = me_scale
Rc = np.sqrt(R) if isinstance(me_scale, float) else np.linalg.cholesky(R)
eps = Rc * np.matrix(np.random.normal(size=(locs.shape[0],1)))
y = x + eps

# introducing missing data
obs_inds = np.random.choice(dim_x*dim_y, int(dim_x*dim_y*frac_obs), replace=False)
obs_inds=np.sort(obs_inds)
y_obs = np.empty(np.shape(y)); y_obs[:] = np.NAN; y_obs[obs_inds] = y[obs_inds]

##### parameter optimization #####
def likelihood(kappa):
    cov = lambda _locs1, _locs2: mt.Matern32(_locs1, _locs2, l=kappa, sig=sig)
    mraTree = MRATree(locs, M, J, r0, critDepth, cov, y_obs, me_scale)
    lik = mraTree.getLikelihood()
    return( lik )

xmin = opt.minimize(likelihood, [kappa], method='nelder-mead', \
                    options={'xtol':1e-2, 'disp':False})
logging.info(str(xmin))
```


Finally we show how to produce diagnostic plots
```python
import logging
import numpy as np
import scipy.linalg as lng

from pyMRA.MRATree import MRATree
import pyMRA.MRATools as mt

logging.basicConfig(format='%(asctime)s %(message)s', \
                    datefmt='%H:%M:%S',level=logging.INFO)

frac_obs = 0.4

dim_x = 100; dim_y = 1
M=3; J=3; r0=2
critDepth=M+1

### simulate data ###

sig = 1.0; me_scale=1e-1; kappa = 0.3

locs = mt.genLocations(dim_x)
Sig = sig*mt.ExpCovFun(locs, l=kappa)
SigC = np.matrix(lng.cholesky(Sig))
x_raw = np.matrix(np.random.normal(size=(locs.shape[0],1)))
x = SigC.T * x_raw
eps = np.sqrt(me_scale) * np.matrix(np.random.normal(size=(locs.shape[0],1)))
y = x + eps

    
# introducing missing data
obs_inds = np.random.choice(dim_x*dim_y, int(dim_x*dim_y*frac_obs), replace=False)
obs_inds=np.sort(obs_inds)
y_obs = np.empty(np.shape(y)); y_obs[:] = np.NAN; y_obs[obs_inds] = y[obs_inds]

    
### MRA ###
cov = lambda _locs1, _locs2: mt.ExpCovFun(_locs1, _locs2, l=kappa)
mraTree = MRATree(locs, M, J, r0, critDepth, cov, y_obs, me_scale)
xP, sdP = mraTree.predict()
sdP = sdP.reshape((dim_x, dim_y), order='A')

    
### diagnostic plots ###
mraTree.drawBMatrix("prior")
mraTree.drawSparsityPat("prior")
mraTree.drawBMatrix("posterior")
mraTree.drawSparsityPat("posterior")
mraTree.drawGridAndObs()
mraTree.drawKnots()
```






### Acknowledgements

The package was mainly developed during the SiParCS 2017 program at the National Center For Atmospheric Research. The author gratefully acknowledges funding from NCAR and helpful comments offered by NCAR employees and interns. Special thanks to Dr. Dorit Hammerling and Dr. Matthias Katzfuss.


### References
[1] Katzfuss, M. (2017). A multi-resolution approximation for massive spatial datasets. Journal of the American Statistical Association, 112(517), 201-214.

[2] Banerjee, S., Gelfand, A. E., Finley, A. O., & Sang, H. (2008). Gaussian predictive process models for large spatial data sets. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(4), 825-848.