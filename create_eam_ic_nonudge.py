#! /usr/bin/env python

#==============================================================================
#   create_eam_ic_nonudge.py
"""
Python script.

This script is to replace state variables and surface pressure (U, V, T, Q, PS) 
in a EAM initial condition file with reanalysis data.

The Renalaysis has to be processed to the model grid (e.g., the nudging files)

You need to have a initial condition file from any model runs (e.g., .i file).
One way to obtain this without doing a nudging run is to perform a AMIP 
run and save the ".i" file at the desired time.

This script should be run in the directory where the renalysis data is.

Interactive usage:
    create_eam_ic_nonudge.py <case_id nudge_xml_file targ_dir
                           year month day hour>
                           year month day_start day_end hours>

Arguments:
    case_id         # input file name
    nudge_xml_file  # reanalysis file (nc or xml)
    targ_dir        # move new file to the target directory (have to create first)
    year            # reanalysis year
    month           # reanalysis month
    day             # reanalysis day
    hour            # reanalysis hour
"""
#==============================================================================

import subprocess
import sys

import numpy as n

import cdms2
import cdtime


VAR_LIST = ['U', 'V', 'T', 'Q']


def _make_dtg_tag(time):
#-----------------------
    """
    time set up 
    """

    seconds = time.hour * 60 * 60

    ymdh = '%d-%02d-%02d-%05d' % (time.year, time.month, time.day, seconds)

    return(ymdh)


def _swap(src_filename, targ_filename, time):
#--------------------------------------------
    """
    src_filename   # source file name
    targ_filename  # target file name
    time           # time
    """

    print ('File names: ', src_filename, targ_filename)

    fsrc  = cdms2.open(src_filename)
    ftarg = cdms2.open(targ_filename, 'r+')

    for myvar in VAR_LIST:
        vsrc  = fsrc(myvar, time=time)
        vtarg = ftarg.variables[myvar]

        print ('time: ', time)
        print ('src_filename:  ', src_filename)
        print ('targ_filename: ', targ_filename)

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


def _create_eam_ic_nonudge(case_id, nudge_xml_file, targ_dir,
                           year, month, day, hour):
#---------------------------------------------------------------
    """
    See file header.
    """

    time_suff = cdtime.comptime(year, month, day, hour)

    tag_suff = _make_dtg_tag(time_suff)

    # Copy new atm file to correct name for running.
    #-------------------------------------------------------
    src_file  = case_id + '.nc' 
    targ_file = 'fa_' + case_id + '_' + tag_suff + '.nc'

    cmd = 'cp ' + src_file + ' ' + targ_file
    print (cmd)

    ret = subprocess.Popen(cmd, shell=True)
    ret.wait()

    # Swap in the new fields in the atm file.
    #----------------------------------------

    _swap(nudge_xml_file, targ_file, time_suff)

    # Move the renamed, swapped files to the proper directory.
    #---------------------------------------------------------
    #cmd = 'mv ' + targ_file + ' ' + targ_dir
    #print (cmd)

    #ret = subprocess.Popen(cmd, shell=True)
    #ret.wait()

    return


if __name__ == '__main__':
#-------------------------
    case_id        = sys.argv[1]
    nudge_xml_file = sys.argv[2]
    targ_dir       = sys.argv[3]
    year           = int(sys.argv[4])
    month          = int(sys.argv[5])
    day            = int(sys.argv[6])
    hour           = int(sys.argv[7])

    _create_eam_ic_nonudge(case_id, nudge_xml_file, targ_dir,
                           year, month, day, hour)
