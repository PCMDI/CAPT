#! /usr/bin/env python

#==============================================================================
#   make_ab.py
"""
Python script.

This script is to generate the a, b coefficients for ERA5 vertical levels

Interactive usage:
    make_ab_ERA5.py
"""
#==============================================================================

import numpy
n = numpy

PSFC = 101325.0  # Units of Pa.


def _make_ab_ERA5():
#--------------
    """
    See file header.
    """

    #abh = n.loadtxt('/global/homes/h/hyma/cssef/utils/zdata/ERA5_L137.txt', skiprows=1)
    abh = n.loadtxt('/global/homes/h/hyma/cssef/utils/zdata/ERA5_L137.txt')
    ah = abh[:,0]
    bh = abh[:,1]

    nhalf_levs = ah.shape[0]
    nfull_levs = nhalf_levs - 1

    af = n.zeros((nfull_levs,), n.float32)
    bf = n.zeros((nfull_levs,), n.float32)

    for ii in range(nfull_levs):
        af[ii] = (ah[ii] + (ah[ii+1]) / 2.0)
        bf[ii] = (bh[ii] + (bh[ii+1]) / 2.0)

    #pfc = af + (bf * PSFC)
    #phc = ah + (bh * PSFC)

    #print n.sqrt(n.average((pf - (pfc / 100.0))**2)), \
    #      n.sqrt(n.average((ph - (phc / 100.0))**2))

    fh = open('era5_137lev_ab_half.dat', 'wr')

    for ii in range(nhalf_levs):
        fh.write('%f %f \n' % (ah[ii], bh[ii]))

    fh.close()

    ff = open('era5_137lev_ab_full.dat', 'wr')

    for ii in range(nfull_levs):
        ff.write('%f %f \n' % (af[ii], bf[ii]))

    ff.close()

    return


if __name__ == '__main__':
#-------------------------
    _make_ab_ERA5()
