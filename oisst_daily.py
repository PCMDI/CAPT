#! /usr/bin/env python

#==============================================================================
#   oisst_daily.py
"""
Python script.

This script is to combine SST and sea ice into one file for model runs.

This file can be used by both E3SM and CESM. 

Interactive usage:
    oisst_daily.py <daily_sst_dir daily_sst_filename>
                   [optional arguments]

Required arguments:
    daily_sst_dir       # TBD.
    daily_sst_filename  #

Optional arguments:
    targ_grid_res       # Target grid resolution; choices are 'HALFXHALF',
                          'ONEXONE', or 'TWOXTWO' ['ONEXONE'].
"""
#==============================================================================

import sys

import numpy
n = numpy

import cdms2
from cdms2 import Cdunif
import cdtime
from regrid2 import Regridder

import fill_msg_grid


TARG_GRID_RES_DICT = \
    {'HALFXHALF': {'label': '_0.5x0.5_',
                   'args':  (-90.0, 361, 0.5, 0.0, 720, 0.5)},
     'ONEXONE':   {'label': '_1x1_',
                   'args':  (-90.0, 181, 1.0, 0.0, 360, 1.0)},
     'TWOXTWO':   {'label': '_2x2_',
                   'args':  (-90.0,  91, 2.0, 0.0, 180, 2.0)}}

TARG_GRID_RES_DEF = 'ONEXONE'  # ['HALFXHALF'|'ONEXONE'|'TWOXTWO']


def fill_msg(data, epsx=0.01, nscan=100):
#----------------------------------------
    """
    TBD.

    data   # TBD.
    epsx   #
    nscan  #
    """

    gtype = 1
    guess = 1
    ier   = 0
    mscan = 10000
    relc  = 0.6
    xmsg  = -999.0

    sstf = data.filled(xmsg)
    sstf = n.array(sstf, n.float64)

    dummy = fill_msg_grid.poisxy1(sstf.T, xmsg, guess, gtype, nscan, epsx,
                                  relc, mscan, ier)

    return sstf


def _oisst_daily(daily_sst_dir, daily_sst_filename, targ_grid_res):
#---------------------------------------------------------------------
    """
    See file header.
    """

    print ('targ_grid_res: ', targ_grid_res)

    time_span_tag = daily_sst_filename.split('.')[2]
    time_units = 'days since ' + time_span_tag.split('-')[0] + '-1-1 00:00:00'

    year_start = int(time_span_tag.split('-')[0])

    daily_ice_filename = daily_sst_filename.replace('sst', 'icec')

    # Create target grid.
    # For a ONEXONE target grid resolution with arguments:
    #     (-90., 181, 1., 0., 360, 1.)
    # A grid will be created with:
    #     latitude  starting at -90 & going north 181 points,
    #         with an increment of 1 degree;
    #     longitude starting at  0E & going east  360 points,
    #         with an increment of 1 degree.
    # The out_filename will reflect the designated resolution.
    #---------------------------------------------------------
    args = TARG_GRID_RES_DICT[targ_grid_res]['args']
    targ_grid = cdms2.createUniformGrid(args[0], args[1], args[2], args[3],
                                        args[4], args[5])

    label = TARG_GRID_RES_DICT[targ_grid_res]['label']
    out_filename = 'sst_daily_cdcunits' + label + time_span_tag + '.nc'

    fdaily_sst = cdms2.open(daily_sst_dir + '/' + daily_sst_filename)
    fdaily_ice = cdms2.open(daily_sst_dir + '/' + daily_ice_filename)

    input_grid = fdaily_sst.variables['sst'].getGrid()

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

    sst_w = fdaily_sst.variables['sst']

    ntimes = sst_w.shape[0]

    intime       = sst_w.getTime()
    intime_units = intime.units
    intimes      = intime[:]

    # Time loop.
    #-----------
    time_idx_out = -1

    for time_idx in range(ntimes - 1):
        # Data is centered on time in file.
        #----------------------------------
        mid_intime = intimes[time_idx]

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

        data = fdaily_sst('sst', time=slice(time_idx, (time_idx + 1)), raw=1,
                          squeeze=1)

        data_f = fill_msg(data, nscan=200)
        data_f = n.array(data_f, n.float32)
        print ('data_f min,max,mean: ', \
                data_f.min(), data_f.max(), data_f.mean() )

        data_f = rg_in2targ(data_f).filled()
        out_sst = data_f
        print ('out_sst min,max,mean: ', \
                out_sst.min(), out_sst.max(), out_sst.mean() )

        data = fdaily_ice('icec', time=slice(time_idx, (time_idx + 1)), raw=1,
                          squeeze=1)
        data_f = data * 1.0
        print ('data_f min,max,mean: ', \
                data_f.min(), data_f.max(), data_f.mean() )

        # Set ice to zero where missing - over land.
        #-------------------------------------------
        data_f = rg_in2targ(data_f).filled(0.0)
        out_ice = data_f
        print ('out_ice min,max,mean: ', \
                out_ice.min(), out_ice.max(), out_ice.mean() )

        sst_cpl[time_idx_out,:,:] = out_sst
        ifrac[time_idx_out,:,:]   = out_ice

        time_idx_out = time_idx_out + 1

        fout.sync()

    fout.close()

    return


if __name__ == '__main__':
#-------------------------
    targ_grid_res = TARG_GRID_RES_DEF

    daily_sst_dir      = sys.argv[1]
    daily_sst_filename = sys.argv[2]

    if len(sys.argv) > 3:
        targ_grid_res = sys.argv[3]

    _oisst_daily(daily_sst_dir, daily_sst_filename, targ_grid_res)
