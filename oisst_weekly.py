#! /usr/bin/env python

#==============================================================================
#   oisst_weekly.py
"""
Python script.

This script is to combine SST and sea ice into one file for model runs.

This file can be used by both E3SM and CESM. 

Interactive usage:
    oisst_weekly.py <weekly_sst_dir weekly_sst_filename>

Arguments:
    weekly_sst_dir       # 
    weekly_sst_filename  #
"""
#==============================================================================

import sys

import numpy
n = numpy

import cdms2
from cdms2 import Cdunif
import cdtime
from regrid2 import Regridder


def _oisst_weekly(weekly_sst_dir, weekly_sst_filename):
#---------------------------------------------------------
    """
    See file header.
    """

    time_span_tag = weekly_sst_filename.split('.')[2]
    time_units = 'days since ' + time_span_tag.split('-')[0] + '-1-1 00:00:00'

    year_start = int(time_span_tag.split('-')[0])

    weekly_ice_filename = weekly_sst_filename.replace('sst', 'icec')

    # Resolution can be varied if desired by changing the arguments to
    # createUniformGrid below.
    #
    # For a "1x1" resolution use:
    #     targ_grid = cdms2.createUniformGrid(-90., 181, 1., 0., 360, 1.)
    # This will create a grid with:
    #     latitude  starting at -90 & going north 181 points,
    #         with an increment of 1 degree;
    #     longitude starting at  0E & going east  360 points,
    #         with an increment of 1 degree.
    # Be sure to change the out_filename to reflect the new resolution.
    #--------------------------------------------------------------------
    targ_grid = cdms2.createUniformGrid(-90., 181, 1., 0., 360, 1.)

    out_filename = 'sst_weekly_cdcunits_1x1_' + time_span_tag + '.nc'

    fweekly_sst = cdms2.open(weekly_sst_dir + '/' + weekly_sst_filename)
    fweekly_ice = cdms2.open(weekly_sst_dir + '/' + weekly_ice_filename)

    input_grid = fweekly_sst.variables['sst'].getGrid()

    rg_in2targ = Regridder(input_grid, targ_grid)

    # Create file and variables for output.
    #--------------------------------------
    #fout = NetCDF.NetCDFFile(out_filename, 'w')
    fout = Cdunif.CdunifFile(out_filename, 'w')

    lons = targ_grid.getLongitude()[:]
    lats = targ_grid.getLatitude()[:]

    fout.createDimension('lon', len(lons))
    fout.createDimension('lat', len(lats))
    fout.createDimension('time', None)

    sst_cpl = fout.createVariable('sst', 'f', ('time', 'lat', 'lon'))
    sst_cpl.long_name = 'sea surface temperature'
    sst_cpl.units = 'degrees_C'

    ifrac = fout.createVariable('ifrac', 'f', ('time', 'lat', 'lon'))
    ifrac.long_name = 'ice fraction'
    ifrac.units = 'fraction'

    lat = fout.createVariable('lat', 'd', ('lat',))
    lat.long_name = 'latitude of grid cell center'
    lat.units = 'degrees_north'

    lon = fout.createVariable('lon', 'd', ('lon',))
    lon.long_name = 'longitude of grid cell center'
    lon.units = 'degrees_east'

    time = fout.createVariable('time', 'd', ('time',))
    time.long_name = 'time'
    time.units = time_units
    time.calendar = 'noleap'

    date = fout.createVariable('date', 'i', ('time',))
    date.long_name = 'calendar date (YYYYMMDD)'

    datesec = fout.createVariable('datesec', 'i', ('time',))
    datesec.long_name = 'seconds elapsed on calendar date'
    datesec.units = 'seconds'

    # Coordinate data.
    #-----------------
    lat[:] = lats
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'

    lon[:] = lons
    lon.long_name = 'longitude'
    lon.units = 'degrees_east'

    sst_w     = fweekly_sst.variables['sst']
    ice_cov_w = fweekly_ice.variables['icec']

    ntimes = sst_w.shape[0]

    intime        = sst_w.getTime()
    intime_units  = intime.units
    intime_bounds = fweekly_sst.variables['time_bnds'][:]

    # Time loop.
    #-----------
    time_idx_out = -1

    for time_idx in range(ntimes):
        # Data NOT centered on time in file?
        #-----------------------------------
        mid_intime = (intime_bounds[time_idx,0] + \
                      intime_bounds[time_idx,1]) / 2.0

        rtime = cdtime.reltime(mid_intime, intime_units)
        ctime = rtime.tocomp()

        new_reltime = ctime.torel(time_units, cdtime.NoLeapCalendar)
        new_ctime   = new_reltime.tocomp()

        year = ctime.year

        if year < year_start:
            #=======
            continue
            #=======

        month  = ctime.month
        day    = ctime.day
        hour   = ctime.hour
        minute = ctime.minute
        second = ctime.second

        # Change time units.
        #-------------------
        print ('time_idx_out, ctime, new_ctime: ', \
              time_idx_out, ctime, new_ctime)

        time[time_idx_out] = new_reltime.value
        print ('time[time_idx_out]: ', time[time_idx_out])

        date[time_idx_out]    = (year * 10000) + (month * 100) + day
        datesec[time_idx_out] = (hour * 60 * 60) + (minute * 60) + second

        out_sst = rg_in2targ(sst_w[time_idx])
        print ('out_sst min,max,mean: ', \
              out_sst.min(), out_sst.max(), out_sst.mean())

        out_ice = rg_in2targ(ice_cov_w[time_idx]) / 100.0
        print ('out_ice min,max,mean: ', \
              out_ice.min(), out_ice.max(), out_ice.mean())

        # Set ice to zero where missing - over land.
        #-------------------------------------------
        out_ice = n.where(n.greater(out_ice, 1.0e4), 0.0, out_ice)
        print ('out_ice min,max,mean: ', \
              out_ice.min(), out_ice.max(), out_ice.mean())

        sst_cpl[time_idx_out,:,:] = out_sst
        ifrac[time_idx_out,:,:]   = out_ice

        time_idx_out = time_idx_out + 1

        fout.sync()

    fout.close()

    return


if __name__ == '__main__':
#-------------------------
    weekly_sst_dir      = sys.argv[1]
    weekly_sst_filename = sys.argv[2]

    _oisst_weekly(weekly_sst_dir, weekly_sst_filename)
