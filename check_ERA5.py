#! /usr/bin/env python

#==============================================================================
#   check_ERA5.py
"""
Python script.

2020/5/13
This script is to check if all the necessary ERA5 files are in place.

Interactive usage:
    check_ERA5.py <year_start month_start year_end month_end data_dir> [optional arguments]

Required arguments:
    year_start   # 
    month_start  #
    year_end     #
    month_end    #
    data_dir     #

Optional arguments:
    cdo_path     # This must be set if cdo is not in your executable
                   search path [''].
"""
#==============================================================================

import calendar
import os
import subprocess
import sys


CDO_PATH_DEF = ''

MONTH_END_DATE = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def _check_ERA5(year_start, month_start, year_end, month_end,
                  data_dir, cdo_path):
#--------------------------------------------------------------
    """
    See file header.
    """

    for year in range(year_start, (year_end + 1)):
        year_str = '%04d' % year

        is_leap_yr = calendar.isleap(year)

        for month in range(month_start, (month_end + 1)):
            month_str = '%02d' % month

            #dirname = data_dir + '/' + year_str + '/' + month_str
            dirname = data_dir + '/' 
            os.chdir(dirname)

            ndays = MONTH_END_DATE[month-1]

            if is_leap_yr and (month == 2):
                print 'Leap year.'
                ndays = 29

            day_start = 1
            day_end = ndays

            missing_files = []

            for day in range(day_start, (day_end + 1)):
                day_str = '%02d' % day
                ymd = year_str + month_str + day_str

                for hour in [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]:
                    hour_str = '%02d' % hour
                    in_filename = 'UVTQPS-' + ymd + '_' + hour_str + '.grib'

                    if not os.path.isfile(in_filename):
                        missing_files.append(in_filename)

            if len(missing_files) > 0:
                print 'File(s) missing:'
                for mf in missing_files:
                    print '    ', mf

                mystr = '\nONE OR MORE FILES MISSING - EXITING!\n\n'
                raise SystemExit(mystr)

            else:
                print 'All files in place.'

    return


if __name__ == '__main__':
#-------------------------
    cdo_path = CDO_PATH_DEF

    year_start  = int(sys.argv[1])
    month_start = int(sys.argv[2])
    year_end    = int(sys.argv[3])
    month_end   = int(sys.argv[4])
    data_dir    = sys.argv[5]

    if len(sys.argv) > 6:
        cdo_path = sys.argv[6]

        if not cdo_path.endswith('/'):
            cdo_path = cdo_path + '/'

    _check_ERA5(year_start, month_start, year_end, month_end,
                  data_dir, cdo_path)
