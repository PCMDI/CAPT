#! /usr/bin/env python

#==============================================================================
#   get_interim.py
"""
Python script.

This script is to download ERA Interim model level data. 
Please refer this website for more information on how to download the data.

https://confluence.ecmwf.int/display/WEBAPI/Python+ERA-interim+examples#PythonERAinterimexamples-ERA-InterimModelLevels

Interactive usage:
    get_interim.py <year_start month_start day_start
                    year_end   month_end   day_end>

Arguments:
    year_start   # 
    month_start  #
    day_start    #
    year_end     #
    month_end    #
    day_end      #
"""
#==============================================================================

import datetime
import os
import sys
import time

from ecmwf import ECMWFDataServer


def _get_interim(year_start, month_start, day_start, year_end, month_end,
                 day_end):
#------------------------------------------------------------------------
    """
    See file header.
    """

    login_name = os.getlogin()

    ## You have to change the following manually to your own credential for the server.
    server = \
        ECMWFDataServer('http://data-portal.ecmwf.int/data/d/dataserver/',
                        '538c43fc1d36a64fc183d2a71301343f',
                        '%s@llnl.gov' % login_name)

    dt     = datetime.datetime(year_start, month_start, day_start, 0)
    dt_end = datetime.datetime(year_end, month_end, day_end, 18)

    delt_6hr = datetime.timedelta(hours=6)

    while dt <= dt_end:
        ymd  = '%04d%02d%02d' % (dt.year, dt.month, dt.day)
        hour = '%02d' % dt.hour
        ymdh = ymd + hour

        ###data_file = 'UVTQPS_' + ymdh + '.grib'
        ###data_file = 'UVTQPS-' + ymdh + '.grib'
        data_file = 'UVTQPS-' + ymd + '_' + hour + '.grib'
        print (data_file)

        data_dict = \
            {'stream':   'oper',
             'step':     0,
             'dataset':  'interim_full_daily',
             'time':     hour,
             'levtype':  'ml',
             'date':     ymd,
             'type':     'an',
             'class':    'ei',
             'param':    '129.128/152.128/130.128/133.128/131.128/132.128',
             'levelist': 'all',
             'target':   data_file}

        try:
            server.retrieve(data_dict)

        except:
            print ('Try again.')

            time.sleep(5 * 60)

            server.retrieve(data_dict)

        dt = dt + delt_6hr

    return


if __name__ == '__main__':
#-------------------------
    year_start  = int(sys.argv[1])
    month_start = int(sys.argv[2])
    day_start   = int(sys.argv[3])
    year_end    = int(sys.argv[4])
    month_end   = int(sys.argv[5])
    day_end     = int(sys.argv[6])

    _get_interim(year_start, month_start, day_start, year_end, month_end,
                 day_end)
