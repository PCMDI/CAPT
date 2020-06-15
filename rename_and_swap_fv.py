#! /usr/bin/env python

#==============================================================================
#   rename_and_swap_fv.py
"""
Python script.

This script TBD.

This script should be run in the directory where the nudging files are
found.

Interactive usage:
    rename_and_swap_fv.py <case_id nudge_xml_file targ_dir
                           year month day_start day_end hours>

Arguments:
    case_id         # TBD.
    nudge_xml_file  #
    targ_dir        #
    year            #
    month           #
    day_start       #
    day_end         #
    hours           #
"""
#==============================================================================

import subprocess
import sys

import numpy as n

import cdms2
import cdtime


VAR_LIST = ['US', 'VS', 'T', 'Q']


def _make_dtg_tag(time):
#-----------------------
    """
    TBD.

    time  # TBD.
    """

    seconds = time.hour * 60 * 60

    ymdh = '%d-%02d-%02d-%05d' % (time.year, time.month, time.day, seconds)

    return(ymdh)


def _swap(src_filename, targ_filename, time, targ_type='initial'):
#-----------------------------------------------------------------
    """
    TBD.

    src_filename   # TBD.
    targ_filename  #
    time           #
    targ_type      #
    """

    cdms2.axis.longitude_aliases.append('slon')
    cdms2.axis.latitude_aliases.append('slat')

    print 'File names: ', src_filename, targ_filename

    fsrc  = cdms2.open(src_filename)
    ftarg = cdms2.open(targ_filename, 'r+')

    if targ_type == 'initial':
        targ_data_order = ('tyzx')
    elif targ_type == 'history':
        targ_data_order = ('tzyx')

    for myvar in VAR_LIST:
        vsrc  = fsrc(myvar, time=time, order=targ_data_order)
        vtarg = ftarg.variables[myvar]

        print 'time: ', time
        print 'src_filename:  ', src_filename
        print 'targ_filename: ', targ_filename

        try:
            vtarg[0] = vsrc[0]
        except TypeError:
            vtarg[0] = vsrc[0].astype(n.float)

    vsrc  = fsrc('PS', time=time)
    vtarg = ftarg.variables['PS']

    try:
        vtarg[0] = vsrc[0]
    except TypeError:
        vtarg[0] = vsrc[0].astype(n.float)

    fsrc.close()
    ftarg.close()

    return


def _rename_and_swap_fv(case_id, nudge_xml_file, targ_dir,
                        year, month, day_start, day_end, hours):
#---------------------------------------------------------------
    """
    See file header.
    """

    for day in range(day_start, (day_end + 1)):
        for hour in hours:
            time_suff = cdtime.comptime(year, month, day, hour)
            time_pref = time_suff.sub(6, cdtime.Hours)

            tag_pref = _make_dtg_tag(time_pref)
            tag_suff = _make_dtg_tag(time_suff)

            # Copy new atm, land files to correct names for running.
            # If '.clm2.i.' does not work, try '.clm2.r.'.
            #-------------------------------------------------------
            for model in ['.cam.i.', '.clm2.i.']:
                src_file  = case_id + model + tag_suff + '.nc'
                targ_file = 'fa_' + tag_pref + model + tag_suff + '.nc'

                cmd = 'cp ' + src_file + ' ' + targ_file
                print cmd

                ret = subprocess.Popen(cmd, shell=True)
                ret.wait()

            # Swap in the new fields in the atm file.
            #----------------------------------------
            model = '.cam.i.'
            targ_file = 'fa_' + tag_pref + model + tag_suff + '.nc'

            _swap(nudge_xml_file, targ_file, time_suff, targ_type='history')

            # Move the renamed, swapped files to the proper directory.
            # If '.clm2.i.' does not work, try '.clm2.r.'.
            #---------------------------------------------------------
            for model in ['.cam.i.', '.clm2.i.']:
                targ_file = 'fa_' + tag_pref + model + tag_suff + '.nc'

                cmd = 'mv ' + targ_file + ' ' + targ_dir
                print cmd

                ret = subprocess.Popen(cmd, shell=True)
                ret.wait()

    return


if __name__ == '__main__':
#-------------------------
    case_id        = sys.argv[1]
    nudge_xml_file = sys.argv[2]
    targ_dir       = sys.argv[3]
    year           = int(sys.argv[4])
    month          = int(sys.argv[5])
    day_start      = int(sys.argv[6])
    day_end        = int(sys.argv[7])

    _rename_and_swap_fv(case_id, nudge_xml_file, targ_dir,
                        year, month, day_start, day_end, hours=[0,6,12,18])
