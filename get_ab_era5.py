#! /usr/bin/env python

#==============================================================================
#   get_ab_era5.py
"""
Python module.

This module returns a/b arrays for units of ps in Pa for ERA5.

a/b's stored in order levels 1 to 137; 1: top, 137: surface.
"""
#==============================================================================

import numpy
n = numpy


def get_ab_era5():
#------------
    """
    See file header.
    """
    #The following two lines may need to be changed to where the files are
    full_filename = 'zdata/era5_137lev_ab_full.dat'
    half_filename = 'zdata/era5_137lev_ab_half.dat'

    full = n.loadtxt(full_filename)
    half = n.loadtxt(half_filename)

    hyam = full[:,0]
    hybm = full[:,1]

    hyai = half[:,0]
    hybi = half[:,1]

    p0 = 1.0

    return p0, hyam, hybm, hyai, hybi
