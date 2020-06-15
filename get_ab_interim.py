#! /usr/bin/env python

#==============================================================================
#   get_ab_interim.py
"""
Python module.

This module returns a/b arrays for units of ps in Pa for ERA Interim.

a/b's stored in order levels 1 to 60; 1: top, 60: surface.
"""
#==============================================================================

import numpy
n = numpy


def get_ab_interim():
#------------
    """
    See file header.
    """

    full_filename = 'erainterim_60lev_ab_full.dat'
    half_filename = 'erainterim_60lev_ab_half.dat'

    full = n.loadtxt(full_filename)
    half = n.loadtxt(half_filename)

    hyam = full[:,0]
    hybm = full[:,1]

    hyai = half[:,0]
    hybi = half[:,1]

    p0 = 1.0

    return p0, hyam, hybm, hyai, hybi
