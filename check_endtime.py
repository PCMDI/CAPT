#! /usr/bin/env python

#==============================================================================
#   check_endtime.py
"""
Python script.

This script checks on the last sst time point.

Interactive usage:
    check_endtime.py <filename>

Arguments:
    filename  # File name to check.
"""
#==============================================================================

import sys

import cdms2
import cdtime


def _check_endtime(filename):
#----------------------------
    """
    See file header.
    """

    fin = cdms2.openDataset(filename)

    sst = fin.variables['sst']

    t = sst.getTime()

    print ('Final time: ', cdtime.reltime(t[-1], t.units).tocomp())

    fin.close()

    return


if __name__ == '__main__':
#-------------------------
    filename = sys.argv[1]

    _check_endtime(filename)
