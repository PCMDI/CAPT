#! /usr/bin/env python

#==============================================================================
#   get_abera5.py
"""
Python module.

This module returns a/b arrays for units of ps in Pa.

a/b's stored in order levels 1 to 137; 1: top, 137: surface.
"""
#==============================================================================

import numpy
n = numpy


def get_abera5():
#------------
    """
    See file header.
    """

    full_filename = '/global/homes/h/hyma/cssef/utils/zdata/era5_137lev_ab_full.dat'
    half_filename = '/global/homes/h/hyma/cssef/utils/zdata/era5_137lev_ab_half.dat'

    full = n.loadtxt(full_filename)
    half = n.loadtxt(half_filename)

    hyam = full[:,0]
    hybm = full[:,1]

    hyai = half[:,0]
    hybi = half[:,1]

    p0 = 1.0

    return p0, hyam, hybm, hyai, hybi
