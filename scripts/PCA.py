#!/usr/bin/env python
#
# Code module for calculating the PCA matrices of the
# SalishSeaCast surface nitrate and temperature records.
#
# required for the analyses presented in:
# 
# B. Moore-Maley and S. E. Allen: Wind-driven upwelling and
# surface nutrient delivery in a semi-enclosed coastal sea,
# Ocean Sci., 2021.
#
# $ cd scripts
# $ python3 PCA.py /path/to/files

import numpy as np
import xarray as xr
import sys
from datetime import datetime, timedelta
from scipy import signal, fft
from tqdm import tqdm

import tools


def regrid(xflat, yflat, valuesflat):
    """Regrid a flattened array with the landpoints removed according to the
    corresponding xflat, yflat coordinate arrays. This function assumes a
    2D input shape for valuesflat of [space, mode].
    """
    
    shape = (max(yflat)+1, max(xflat)+1, valuesflat.shape[1])
    valuesgridded = np.zeros(shape)
    for y, x, row in zip(yflat, xflat, valuesflat):
        valuesgridded[y, x, :] = row
    
    return valuesgridded


def varimax(A, maxiter=40, tol=1e-5):
    """Calculate the varimax rotation matrix H from the n x p PC loadings matrix A. H is determined
    iteratively from the Lagrange multiplier optimization of the varimax criterion.
    
    Adapted from:
    
    Horst, P. (1965) Factor Analysis of Data Matrices. Holt, Rinehart and Winston. New York, USA.
    Chapter 18: Analytical Rotations
        - Section 18.4: Simultaneous Factor Varimax Solution, Equations 18.4.1-10, pp. 428-429
        - Section 18.7: Mathematical Proofs, Equations 18.7.32-54, pp. 437-438
    
    The algorithm described in Section 18.4 has been reformulated to use SVD based on equivalent
    definitions for the rotation matrix described in Section 18.7. The eigenvalue matrix is used
    to evaluate convergence.
    
    This version of the varimax algorithm is functionally identical to those found in Sci-kit learn,
    Matlab, R, and presumably others.
    """
    
    # Initialization
    n, p = A.shape
    H = np.eye(p)
    d = 0
    
    # Iteration
    for i in tqdm(range(maxiter), desc='Calculating rotation matrix'):
        d_old = d
        B = A.dot(H)  # -------------------------------------------------------- 18.4.5
        beta = B * B * B - B.dot(np.diag(np.diag(B.T.dot(B)))) / n  # ---------- 18.4.6
        P, Delta, Q_T = np.linalg.svd(A.T.dot(beta))  # ------------------------ 18.7.42
        H = P.dot(Q_T)  # ------------------------------------------------------ 18.7.45
        d = sum(Delta)
        
        # Convergence
        if d_old != 0 and d/d_old < 1 + tol: break

    return H


def calc_PCA(z):
    """Calculate EOF matrices of n x p data matrix z using SVD
    and optional varimax rotation
    """
    
    # Calculate orthogonal PCA matrices
    A_prime, sqrtL, E_T = np.linalg.svd(z, full_matrices=False)
    A = A_prime.dot(np.diag(sqrtL))
    A2 = A * A
    var = A2.sum(axis=0) / A2.sum()
    E = E_T.T

    # Get varimax rotation matrix
    R = varimax(A)

    # Rotate matrices
    B = A.dot(R)
    B2 = B * B
    var_rot = B2.sum(axis=0) / B2.sum()
    U = E.dot(R)

    # Sort rotated matrices
    isort = var_rot.argsort()[::-1]

    # Return xarray-compatible netCDF dict
    PCA = {'A': A, 'E': E, 'var': var, 'B': B[:, isort], 'U': U[:, isort], 'var_rot': var_rot[isort]}
    
    return PCA


def build_PCA_files(results_path, subsample=5, cutoff=1235):
    """Call the principal component analysis and varimax rotation
    functions and build them to netCDF output
    """

    # Load aggregated results file
    slc = slice(None, None, subsample)
    with xr.open_dataset(results_path + 'EOF_paper_model_fields.nc') as ds:
        data = {var: ds[var].values for var in ('temperature', 'nitrate')}
        coords = {var: ds[var].values[slc] for var in ('x', 'y')}
        coords['time'] = tools.formattime(ds.time.values)
        landmask = ds.landmask.values

    # Calculate seasonal indices
    waterpoints = tools.openwaterpoints(landmask)
    nitrate = np.median(tools.flatten(data['nitrate'], waterpoints), axis=1)
    isegments, iseason = tools.calc_seasonal_indices(coords['time'], nitrate)
    coords['time'] = coords['time'][iseason]

    # Build flattened, subsampled coordinate arrays
    landmask = landmask[slc, slc]
    y, x = [range(dim) for dim in landmask.shape]
    maskflat = landmask.ravel().astype(bool)
    xflat, yflat = [var.ravel()[maskflat] for var in np.meshgrid(x, y)]

    # Calculate EOFs
    for var in ['temperature', 'nitrate']:

        # Subsample and flatten
        raw = tools.flatten(data[var][:, slc, slc], landmask)

        # Subtract lowpass filter and extract productive season
        z = np.vstack([col - tools.lowpass(col, cutoff) for col in raw.T]).T[iseason, :]

        # Subtract mean and calculate PCA
        PCA = calc_PCA(z - z.mean(axis=0)[None, :])

        # Build PCA results as xarray Dataset and save to netCDF
        variables = {
            'landmask': (['y', 'x'], landmask),
            'median': (['y', 'x'], np.median(data[var][iseason, slc, slc], axis=0)),
            'A': (['time', 'mode'], PCA['A']),
            'B': (['time', 'mode'], PCA['B']),
            'E': (['y', 'x', 'mode'], regrid(xflat, yflat, PCA['E'])),
            'U': (['y', 'x', 'mode'], regrid(xflat, yflat, PCA['U'])),
            'var': ('mode', PCA['var']),
            'var_rot': ('mode', PCA['var_rot']),
        }
        xr.Dataset(variables, coords).to_netcdf(results_path + var + '_PCA.nc')


if __name__ == "__main__":
    build_PCA_files(sys.argv[1])