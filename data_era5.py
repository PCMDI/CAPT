#! /usr/bin/env python

#==============================================================================
#   data_era5.py
"""
Python module.

This module creates a dataera5 class to retrieve data from ERA5 grib
files.
"""
#==============================================================================

import numpy
n = numpy

import cdms2

import get_ab_era5


class dataera5:
    """
    TBD.
    """

    def __init__(self, state_xmlfile_dir, state_xmlfile_name):
    #---------------------------------------------------------
        """
        Initializes dataera5 class.

        state_files_dir  # TBD.
        state_file_name  #
        """

        self.state_xmlfile_dir  = state_xmlfile_dir
        self.state_xmlfile_name = state_xmlfile_name

        # These are to conform with ERA5 names.
        #---------------------------------------
        self.var_dict = {'PHIS': 'PHIS', 'PS': 'PS', 'Q': 'Q',
                         'V': 'V', 'U': 'U', 'T': 'T'}

        fin = \
            cdms2.open(self.state_xmlfile_dir + '/' + self.state_xmlfile_name)

        zsfc = fin.variables['PHIS']

        # Get coordinates.
        #-----------------
        lats = zsfc.getLatitude()[:]
        lons = zsfc.getLongitude()[:]

        self.lons = lons
        self.lats = lats

        # Create a CDMS grid object, so can use later for regridding.
        #------------------------------------------------------------
        self.grid = cdms2.createGenericGrid(self.lats, self.lons)

        # Get hybrid coordinate information.
        #-----------------------------------
        p0, hyam, hybm, hyai, hybi = get_ab_era5.get_ab_era5()

        self.p0 = p0

        self.hyam = hyam
        self.hybm = hybm
        self.hyai = hyai
        self.hybi = hybi

        self.fin = fin

    def get_const_data(self, var):
    #-----------------------------
        """
        TBD.

        Only constant is surface geopotential.

        var  # TBD.
        """

        const_data = None

        if var != 'PHIS':
            mystr = '\nONLY DOES PHIS - EXITING!\n\n'
            raise SystemExit(mystr)

        phis = self.fin.variables[self.var_dict[var]][0]

        const_data = n.squeeze(phis)

        return const_data

    def get_state_data(self, var, time, regrid_func=None):
    #-----------------------------------------------------
        """
        Returns ERA5 state data.

        State_data assumed to be (time, level, lat, lon).

        var          # State variable to get.
        time         # cdtime component time.
        regrid_func  # CDMS regrid function.
        """

        state_data = None

        if var[:8] != 'pressure':
            state_data = self.fin(self.var_dict[var], time=time, squeeze=1)

            if var == 'PS':
                state_data = n.exp(state_data).filled()

            if regrid_func:
                state_data = regrid_func(state_data.filled())

        else:
            ps = self.fin('PS', time=time, squeeze=1)
            ps = n.exp(ps)

            if regrid_func:
                ps = regrid_func(ps.filled())

            if var == 'pressure_interface':
                edges = (self.hyai[:,n.newaxis] * self.p0) + \
                        n.outer(self.hybi, ps)
                edges.shape = (self.hyai.shape[0], ps.shape[0], ps.shape[1])

            else:
                edges = (self.hyam[:,n.newaxis] * self.p0) + \
                        n.outer(self.hybm, ps)
                edges.shape = (self.hyam.shape[0], ps.shape[0], ps.shape[1])

            state_data = edges

        return state_data
