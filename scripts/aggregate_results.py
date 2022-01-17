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
#
# This routine can be run to aggregate the files with
#
# $ cd scripts
# $ python3 aggregate_results.py /path/to/files

import numpy as np
import xarray as xr
import sys
import yaml
from datetime import datetime, timedelta, timezone
from calendar import monthrange
from tqdm import tqdm

import tools


def build_HRDPS_mask(
    subdomain=[102, 154, 120, 186],
    erddap_url='https://salishsea.eos.ubc.ca/erddap/griddap/',
    grid_file='ubcSSaAtmosphereGridV1',
    mask_file='ubcSSn3DMeshMaskV17-02',
    gridref_path='../grid/grid_from_lat_lon_mask999.nc',
):
    """Build HRDPS mask over the specified subdomain using the 
    `grid_from_lat_lon_mask999.nc` lookup reference file and return
    the flattened boolean mask array `wmask`.
    """

    # Load netCDF files and variables
    slc = {'gridY': slice(*subdomain[2:]), 'gridX': slice(*subdomain[:2])}
    grid = xr.open_dataset(erddap_url + grid_file).isel(slc)
    mask = xr.open_dataset(erddap_url + mask_file)
    gridref = xr.open_dataset(gridref_path)
    lons, lats = [grid[var].values for var in ('longitude', 'latitude')]
    tmask = mask.tmask[0, 0, ...].values
    
    # Build mask from lonlat comparison to NEMO
    wmask = []
    for lon, lat in zip(tqdm(lons.ravel() - 360, desc='Loading HRDPS mask'), lats.ravel()):
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
    
    return wmask


def load_HRDPS(
    attrs_ref, rotation=55.5, subdomain=[102, 154, 120, 186],
    erddap_url='https://salishsea.eos.ubc.ca/erddap/griddap/',
    HRDPS_file='ubcSSaSurfaceAtmosphereFieldsV1',
):
    """Load HRDPS velocities and calculate along-axis velocity and
    wind stress over the open-water points determined by `wmask`.
    Return results as the dictionary `data`. The `subdomain` is
    chosen to limit the results coverage to a box around the SoG.
    """
    
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

    # Initialize variables dict and lists
    variables = {'v_along': [], 'tau_along': []}
    
    # Load wind velocities by year
    for year in tqdm(range(2015, 2020), desc='Loading HRDPS velocities'):
        tslc = slice(datetime(year, 1, 1), datetime(year, 12, 31, 23, 59))
        u, v = [tools.flatten(HRDPS[var].sel(time=tslc).values, wmask) for var in ('u_wind', 'v_wind')]
        wspd = np.sqrt(u**2 + v**2)
        C_d = coeff[0] + coeff[1] * wspd + coeff[2] * wspd**2
        v_along = v * np.cos(theta) - u * np.sin(theta)
        tau_along = rho_a * C_d * v_along * wspd
        variables['v_along'].append(np.median(v_along, axis=1))
        variables['tau_along'].append(np.median(tau_along, axis=1))
    
    # Concatenate results and add attributes
    for var in ('v_along', 'tau_along'):
        variables[var] = ('time', np.hstack(variables[var]), attrs_ref[var])
    
    return variables


def load_NEMO(
    variables, attrs_ref,
    percentiles=[5, 25, 50, 75, 95], subdomain=[110, 370, 300, 850],
    erddap_url='https://salishsea.eos.ubc.ca/erddap/griddap/',
    tracers_file='ubcSSg3DTracerFields1hV19-05',
    biology_file='ubcSSg3DBiologyFields1hV19-05',
    grid_file='ubcSSnBathymetryV17-02',
    mask_file='ubcSSn3DMeshMaskV17-02',
):
    """Load NEMO surface temperature and nitrate fields over
    `subdomain`. Return results as the dictionary `data`.
    """

    # Lazy load and slice NEMO netCDF files
    slc = {'gridY': slice(*subdomain[2:]), 'gridX': slice(*subdomain[:2])}
    tracers = xr.open_dataset(erddap_url + tracers_file).isel(slc)
    biology = xr.open_dataset(erddap_url + biology_file).isel(slc)
    grid = xr.open_dataset(erddap_url + grid_file).isel(slc)
    mask = xr.open_dataset(erddap_url + mask_file).isel(slc)
    
    # Required variables
    tslc = slice(datetime(2015, 1, 1), datetime(2019, 12, 31, 23, 59))
    time = tools.formattime(tracers.time.sel(time=tslc).values)
    landmask = mask.tmask[0, 0, ...].values
    
    # Initialize lists
    coords = {}
    for var in ('temperature', 'nitrate', 'thermocline', 'nitracline'):
        variables[var] = []
    
    # Load results by month (~5-10 hours)
    for year in range(2015, 2020):
        for month in tqdm(range(1, 13), desc=f'Loading surface tracers {year}'):
            day = monthrange(year, month)[1]
            tslc = slice(datetime(year, month, 1), datetime(year, month, day, 23, 59))
            for ds, var in zip([tracers, biology], ['temperature', 'nitrate']):
                variables[var].append(ds[var].isel(depth=0).sel(time=tslc).values)

    # Concatenate results
    for var in ('temperature', 'nitrate'):
        variables[var] = np.concatenate(variables[var])
    
    # Calculate "productive season" date cutoff indices
    waterpoints = tools.openwaterpoints(landmask)
    values = np.median(tools.flatten(variables['nitrate'], waterpoints), axis=1)
    isegments, iseason = tools.calc_seasonal_indices(time, values)
    
    # Loop through productive seasons, loading one day at a time (~5-10 hours)
    for year, segment in zip(range(2015, 2020), isegments):
        season = [time[i].replace(hour=0, minute=0) for i in segment]
        for day in tqdm(range(np.diff(season)[0].days+1), desc=f'Loading profiles {year}'):
            start = season[0] + timedelta(days=day)
            tslc = slice(start, start + timedelta(days=1))
            for ds, key, var in zip([tracers, biology], ['temperature', 'nitrate'], ['thermocline', 'nitracline']):
                values = tools.flatten(ds[key].sel(time=tslc).values, waterpoints)
                variables[var].append([np.percentile(values, p, axis=2) for p in percentiles])

    # Calculate percentiles
    for var in ('thermocline', 'nitracline'):
        variables[var] = np.vstack([np.percentile(row, p, axis=0) for row, p in zip(np.hstack(variables[var]), percentiles)])
    
    # Prepare variables and coordinates for xarray netCDF format
    # -- time, depth, percentiles ----
    attrs = tracers.time.attrs.copy()
    attrs.update({
        'actual_range': tools.formattime(time[[0, -1]], target='int64'),
        'comment': attrs['comment'] + '\n' + attrs_ref['time']['comment'],
    })
    coords['time'] = ('time', time, attrs)
    attrs = tracers.depth.attrs.copy()
    attrs.pop('_ChunkSizes')
    coords['depth'] = ('depth', tracers.depth.values.copy(), attrs)
    coords['percentile'] = ('percentile', percentiles, attrs_ref['percentile'])
    
    # ------ x, y --------------------
    for var in ['x', 'y']:
        key = f'grid{var.capitalize()}'
        coord = grid[key].values
        attrs = grid[key].attrs.copy()
        attrs['actual_range'] = coord[[0, -1]]
        coords[var] = (var, coord, attrs)
    
    # ------ lon, lat ----------------
    for var in ['longitude', 'latitude']:
        attrs = grid[var].attrs.copy()
        attrs['_ChunkSizes'] = np.array(landmask.shape)
        variables[var] = (['y', 'x'], grid[var].values, attrs)
    
    # ------ landmask ----------------
    attrs = mask.tmask.attrs.copy()
    attrs.update({
        '_ChunkSizes': np.array(landmask.shape),
        'comment': 'surface mask',
    })
    variables['landmask'] = (['y', 'x'], landmask, attrs)
    
    # ------ surface tracers ---------
    for ds, var in zip([tracers, biology], ['temperature', 'nitrate']):
        attrs = ds[var].attrs.copy()
        attrs.update({
            '_ChunkSizes': np.array((1,) + landmask.shape),
            'comment': 'surface fields',
        })
        variables[var] = (['time', 'y', 'x'], variables[var], attrs)
    
    # ------ tracer profiles ---------
    for ds, key, var in zip([tracers, biology], ['temperature', 'nitrate'], ['thermocline', 'nitracline']):
        attrs = ds[key].attrs.copy()
        attrs.pop('_ChunkSizes')
        attrs['comment'] = attrs_ref['profiles']['comment']
        variables[var] = (['percentile', 'depth'], variables[var], attrs)
    
    return variables, coords


def build_results_file(savepath):
    """Call the HRDPS and NEMO load functions and build results
    netCDF file for SoG upwelling EOF paper.
    """
    
    # Load attributes from YAML
    with open('attrs.yaml', 'r') as f:
        attrs_ref = yaml.safe_load(f)
    
    # Load HRDPS results (~5-10 min)
    variables = load_HRDPS(attrs_ref['variable_attributes'])
    
    # Load NEMO results (~12-20 hours)
    variables, coords = load_NEMO(variables, attrs_ref['variable_attributes'])
    
    # Set up time attributes
    datefmt = '%Y-%m-%dT%H:%M:%SZ'
    timenow = datetime.now(tz=timezone.utc)
    timestr, timestamp = [timenow.strftime(fmt) for fmt in (datefmt, '%Y-%b-%d %H:%M:%S UTC')]
    start, end = [t.strftime(datefmt) for t in tools.formattime(coords['time'][2]['actual_range'])]
    
    # Define global attributes
    attrs = attrs_ref['dataset_attributes']
    attrs.update({
        'history': (
            timestr + ' multiple files at https://salishsea.eos.ubc.ca/erddap/griddap/\n' +
            timestr + ' https://github.com/SalishSeaCast/SoG_upwelling_EOF_paper/blob/master/scripts/aggregate_results.py'
        ),
        'time_coverage_end': end,
        'time_coverage_start': start,
        'timeStamp': timestamp,
    })
    encoding = {var: {'zlib': True} for var in ('v_along', 'tau_along', 'temperature', 'nitrate')}
    
    # Create xarray dataset and save to netCDF
    xr.Dataset(variables, coords, attrs).to_netcdf(savepath + 'EOF_paper_model_fields.nc', encoding=encoding)


if __name__ == "__main__":
    build_results_file(sys.argv[1])