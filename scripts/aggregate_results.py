#!/usr/bin/env python
#
# Code module for aggregating SalishSeaCast and HRDPS results
# from the SalishSeaCast ERDDAP server
# https://salishsea.eos.ubc.ca/erddap/griddap
#
# required for the analyses presented in:
# 
# B. Moore-Maley and S. E. Allen: Wind-driven upwelling and
# surface nutrient delivery in a semi-enclosed coastal sea,
# Ocean Sci., 2021.

import numpy as np
import xarray as xr
import sys
from datetime import datetime
from calendar import monthrange
from tqdm import tqdm


def build_HRDPS_mask(
    subdomain=[102, 154, 120, 186],
    erddap_url='https://salishsea.eos.ubc.ca/erddap/griddap/',
    grid_file='ubcSSaAtmosphereGridV1',
    mask_file='ubcSSn3DMeshMaskV17-02',
    gridref_path='../grid/grid_from_lat_lon_mask999.nc',
):
    '''Build HRDPS mask over the specified subdomain using the 
    `grid_from_lat_lon_mask999.nc` lookup reference file and return
    the flattened boolean mask array `wmask`.
    '''

    # Load netCDF files and variables
    slc = {'gridY': slice(*subdomain[2:]), 'gridX': slice(*subdomain[:2])}
    grid = xr.open_dataset(erddap_url + grid_file).isel(slc)
    mask = xr.open_dataset(erddap_url + mask_file)
    gridref = xr.open_dataset(gridref_path)
    lons, lats = [grid[var].values for var in ('longitude', 'latitude')]
    tmask = mask.tmask[0, 0, ...].values
    
    # Build mask from lonlat comparison to NEMO
    wmask = []
    print('Loading HRDPS mask ...')
    for lon, lat in zip(tqdm(lons.ravel() - 360), lats.ravel()):
        j, i = [gridref[var].sel(lats=lat, lons=lon, method='nearest').item() for var in ('jj', 'ii')]
        if (j > 0) & (i > 0):
            wmask.append(tmask[j, i])
        else:
            wmask.append(0)

    # Reshape and clip Gulf Islands and inlets
    wmask = np.array(wmask).reshape(lons.shape)
    wmask[:10, :49] = 0
    wmask[:20, :40] = 0
    wmask[50:, 26:] = 0
    wmask[28:, 46:] = 0

    # Flatten to logical
    wmask = wmask.ravel().astype(bool)
    
    return wmask


def load_HRDPS(
    rotation=55.5, subdomain=[102, 154, 120, 186],
    erddap_url='https://salishsea.eos.ubc.ca/erddap/griddap/',
    HRDPS_file='ubcSSaSurfaceAtmosphereFieldsV1',
):
    '''Load HRDPS velocities and calculate along-axis velocity and
    wind stress over the open-water points determined by `wmask`.
    Return results as the dictionary `data`. The `subdomain` is
    chosen to limit the results coverage to a box around the SoG.
    '''
    
    # Wind processing parameters
    rho_a = 1.225     # kg/m3
    theta = np.deg2rad(rotation)

    # Cd coefficients from Hellerman and Rosenstein (1983) JPO (neglect T)
    coeff = [0.934e-3, 0.788e-4, -0.616e-6]

    # Lazy load and slice HRDPS netCDF file
    slc = {'gridY': slice(*subdomain[2:]), 'gridX': slice(*subdomain[:2])}
    HRDPS = xr.open_dataset(erddap_url + HRDPS_file).isel(slc)
    
    # Build mask
    wmask = build_HRDPS_mask(subdomain=subdomain)

    # Load wind velocities by year
    data = {'v_along': [], 'tau_along': []}
    print('Loading HRDPS velocities ...')
    for year in tqdm(range(2015, 2020)):
        tslc = slice(*[datetime(year, 1, 1), datetime(year, 12, 31, 23, 59)])
        u, v = [HRDPS[var].sel(time=tslc).values.reshape(-1, len(wmask))[:, wmask] for var in ('u_wind', 'v_wind')]
        wspd = np.sqrt(u**2 + v**2)
        C_d = coeff[0] + coeff[1] * wspd + coeff[2] * wspd**2
        v_along = v * np.cos(theta) - u * np.sin(theta)
        tau_along = rho_a * C_d * v_along * wspd
        data['v_along'].append(np.median(v_along, axis=1))
        data['tau_along'].append(np.median(tau_along, axis=1))
    
    # Concatenate results and add time
    for var in data:
        data[var] = np.hstack(data[var])
    tslc = slice(*[datetime(2015, 1, 1), datetime(2019, 12, 31, 23, 59)])
    data['time'] = HRDPS.time.sel(time=tslc).values.astype('datetime64[s]').astype(datetime)
    
    return data


def load_NEMO(
    subdomain=[110, 370, 300, 850],
    erddap_url='https://salishsea.eos.ubc.ca/erddap/griddap/',
    tracers_file='ubcSSg3DTracerFields1hV19-05',
    biology_file='ubcSSg3DBiologyFields1hV19-05',
    grid_file='ubcSSnBathymetryV17-02',
    mask_file='ubcSSn3DMeshMaskV17-02',
):
    '''Load NEMO surface temperature and nitrate fields over
    `subdomain`. Return results as the dictionary `data`.
    '''

    # Lazy load and slice NEMO netCDF files
    slc = {'gridY': slice(*subdomain[2:]), 'gridX': slice(*subdomain[:2])}
    tracers = xr.open_dataset(erddap_url + tracers_file).isel(depth=0, **slc)
    biology = xr.open_dataset(erddap_url + biology_file).isel(depth=0, **slc)
    grid = xr.open_dataset(erddap_url + grid_file).isel(slc)
    mask = xr.open_dataset(erddap_url + mask_file).isel(slc)
    
    # Load results by month
    data = {'temperature': [], 'nitrate': []}
    for year in range(2015, 2020):
        print(f'Loading NEMO surface tracers {year} ...')
        for month in tqdm(range(1, 13)):
            day = monthrange(year, month)[1]
            tslc = slice(*[datetime(year, month, 1), datetime(year, month, day, 23, 59)])
            data['temperature'].append(tracers.temperature.sel(time=tslc).values)
            data['nitrate'].append(biology.nitrate.sel(time=tslc).values)
    
    # Concatenate results and add time, lon, lat, landmask
    for var in data:
        data[var] = np.concatenate(data[var])
    tslc = slice(*[datetime(2015, 1, 1), datetime(2019, 12, 31, 23, 59)])
    data['time'] = tracers.time.sel(time=tslc).values.astype('datetime64[s]').astype(datetime)
    for var in ['longitude', 'latitude']:
        data[var] = grid[var].values
    data['landmask'] = mask.tmask[0, 0, ...].values
    
    return data


def build_results_file(savepath, subdomain=[110, 370, 300, 850]):
    '''Call the HRDPS and NEMO load functions and build results
    netCDF file for SoG upwelling EOF paper.
    '''
    
    # Load HRDPS results (~5-10 min)
    HRDPS = load_HRDPS()
    
    # Load NEMO results (~5-10 hours)
    NEMO = load_NEMO(subdomain=subdomain)

    # Build xarray dataset and export to netCDF
    coords = {'time': NEMO['time'], 'y': range(*subdomain[2:]), 'x': range(*subdomain[:2])}
    variables = {
        'longitude'  : (['y', 'x'], NEMO['longitude']),
        'latitude'   : (['y', 'x'], NEMO['latitude']),
        'landmask'   : (['y', 'x'], NEMO['landmask']),
        'v_along'    : ('time', HRDPS['v_along']),
        'tau_along'  : ('time', HRDPS['tau_along']),
        'temperature': (['time', 'y', 'x'], NEMO['temperature']),
        'nitrate'    : (['time', 'y', 'x'], NEMO['nitrate']),
    }
    xr.Dataset(variables, coords).to_netcdf(savepath + 'EOF_paper_model_fields.nc')


if __name__ == "__main__":
    build_results_file(*[str(arg) for arg in sys.argv[1:]])