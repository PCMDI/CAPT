#! /usr/bin/env python

#==============================================================================
#   make_scrip.py
"""
Python script.

This script is to generate a SCRIP grid descriptor file.

Interactive usage:
    make_scrip.py <src_dir src_file src_var src_grid_title>

Arguments:
    src_dir         # Source directory.
    src_file        # Source file.
    src_var         # Source variable.
    src_grid_title  # Source grid title.
"""
#==============================================================================

import sys

import cdms2


def _make_scrip(src_dir, src_file, src_var, src_grid_title):
#-------------------------------------------------------------------
    """
    See file header.
    """

    scrip_filename = src_grid_title + '_scrip.nc'

    fin = cdms2.openDataset(src_dir + '/' + src_file)
    src_grid = fin.variables[src_var].getGrid()

    cdms2.writeScripGrid(scrip_filename, src_grid, src_grid_title)

    return


if __name__ == '__main__':
#-------------------------
    src_dir  = sys.argv[1]
    src_file = sys.argv[2]
    src_var  = sys.argv[3]
    src_grid_title = sys.argv[4]

    _make_scrip(src_dir, src_file, src_var, src_grid_title)
