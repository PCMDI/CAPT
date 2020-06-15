#! /usr/bin/env python

#==============================================================================
#   get_ERA5.py
"""
Python script.

2020/05/13
This script is to download ERA5 reanalysis from ECMWF server.

One has to register an account first over https://cds.climate.copernicus.eu/#!/home
Then follow the instruction on how to setup the download on your local server.

The current setup is to download model level U, V, T, Q, VOR, DIV, and ln(ps)

Interactive usage:
    get_ERA5.py <year_start month_start day_start
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
#One has to change the following line for your own setup
sys.path.insert(1, '/global/homes/h/hyma/ECMWF')
import time
import cdsapi

def _get_ERA5(year_start, month_start, day_start, hr_start, year_end, month_end, day_end, hr_end):

  dt     = datetime.datetime(year_start, month_start, day_start, hr_start)
  dt_end = datetime.datetime(year_end, month_end, day_end, hr_end)
  #delt_1d = datetime.timedelta(hours=24)
  #dt_end = datetime.datetime(year_end, month_end, day_end, 23)
  delt_1hr = datetime.timedelta(hours=1)
  print dt,dt_end
  print delt_1hr
  print ' '

  while dt <= dt_end :
    ymd  = '%04d%02d%02d' % (dt.year, dt.month, dt.day)
    filedate  = '%04d-%02d-%02d' % (dt.year, dt.month, dt.day)
    hour = '%02d' % dt.hour
    #ymdh = ymd + hour

    data_file = 'UVTQPS-' + ymd + '_' + hour + '.grib'
    print data_file
    print filedate,ymd
    print ' '

    data_dict = \
     { 'class':'ea',
       'expver':'1',
       'stream':'oper',
       'type':'an',
       'param':'129.128/152.128/130.128/133.128/131.128/132.128/138.128/155.128',
       'levtype':'ml',
       'levelist':'all',
       'date':filedate,
       'time':hour }

    try:
      c = cdsapi.Client()
      c.retrieve('reanalysis-era5-complete',data_dict,data_file)

    except:
      print 'Try again.'
      time.sleep(1 * 60)
      c = cdsapi.Client()
      c.retrieve('reanalysis-era5-complete',data_dict,data_file)

    dt = dt + delt_1hr

  return

if __name__ == '__main__':
#-------------------------
    year_start  = int(sys.argv[1])
    month_start = int(sys.argv[2])
    day_start   = int(sys.argv[3])
    hr_start    = int(sys.argv[4])
    year_end    = int(sys.argv[5])
    month_end   = int(sys.argv[6])
    day_end     = int(sys.argv[7])
    hr_end      = int(sys.argv[8])

    _get_ERA5(year_start, month_start, day_start, hr_start, year_end, month_end, day_end, hr_end)

