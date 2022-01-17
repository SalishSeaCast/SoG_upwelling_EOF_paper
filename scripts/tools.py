#!/usr/bin/env python
#
# Tools module required for the analyses presented in:
# 
# B. Moore-Maley and S. E. Allen: Wind-driven upwelling and
# surface nutrient delivery in a semi-enclosed coastal sea,
# Ocean Sci., 2021.
#
# Recommended use is to add this module path to your
# PYTHONPATH environment variable
#
# $ export PYTHONPATH=$PYTHONPATH:/path/to/SoG_upwelling_EOF_paper/scripts

import numpy as np
from scipy import signal
from datetime import datetime


def formattime(time, target=datetime):
    """Format time to target format
    """
    
    return time.astype('datetime64[s]').astype(target)


def lowpass(data, cutoff, window_type='blackman'):
    """Apply a Finite Impulse Response (FIR) lowpass filter according
    to the window_type and cutoff using a convolution algorithm
    """
    
    window = signal.get_window(window_type, cutoff)
    filtered = np.convolve(data, window / sum(window), mode='same')
    
    return filtered


def flatten(values, mask):
    """Flatten the spatial dimensions of the input array
    and removed masked points
    """
    
    # Flatten and mask the input array
    args = [-1] if values.ndim < 4 else values.shape[:2]
    maskflat = mask.ravel().astype(bool)
    valuesflat = values.reshape(*args, len(maskflat))[..., maskflat]
    
    return valuesflat


def get_seasonal_indices(time=None, values=None, threshold=2, lag=5):
    """Get seasonal indices of a SalishSeaCast time series based
    on a threshold value and lag time. Return the hardcoded indices
    if the inputs are None.
    """
    
    # Get hardcoded seasonal bounds
    if (time is None) | (values is None):
        seasonbounds = [
            [ 2419,  5657], # 2015
            [11107, 14741], # 2016
            [20327, 23494], # 2017
            [28533, 32219], # 2018
            [37227, 41053], # 2019
        ]
    
    # Get seasonal bounds from values threshold
    else:
        seasonbounds = []
        for year in range(2015, 2020):
            index = np.array([datetime(year, 1, 1) <= t < datetime(year+1, 1, 1) for t in time])
            bounds = np.ma.where(np.ma.masked_where(index==0, values) < threshold)[0][[0, -1]]
            bounds = [bound + lag * hour for bound, hour in zip(bounds, [24, -24])]
            seasonbounds.append(bounds)
    
    # Calculate index array from bounds
    seasonindex = np.hstack([np.arange(*bounds) for bounds in seasonbounds])
    
    return seasonbounds, seasonindex


def openwaterpoints(landmask):
    """Define open water mask from NEMO landmask by removing
    specified islands, inlets and passages
    """

    # Initiate waterpoints from landmask and remove specified areas
    waterpoints = np.copy(landmask)
    waterpoints[:120, :] = 0     # Southern end
    waterpoints[410:, :] = 0     # Northern end
    waterpoints[:170, :130] = 0  # Gulf Islands
    waterpoints[:, 190:] = 0     # Fraser River/Howe Sound
    waterpoints[280:, 130:] = 0  # Jervis Inlet
    
    return waterpoints