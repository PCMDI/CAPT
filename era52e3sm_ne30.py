#! /usr/bin/env python

#==============================================================================
#   era52e3sm_ne30.py
"""
Python script.

This script interpolates ERA5 data to a E3SM ne30 grid. 
The output file(s) can then be used as nudging file(s) for nudging runs.
This should also be usalbe for CESM2 ne30(SE) grid.

Interactive usage:
    era52e3sm_ne30.py  <year month day hour
                        ntimes
                        state_xmlfile_dir state_xmlfile_name
                        map_file_dir      map_file_name
                        hybrid_file       topo_file>
                       [optional arguments]

Required arguments:
    year                # Start year.
    month               # Start month.
    day                 # Start day.
    hour                # Start hour.
    ntimes              # Number of 1 hour time intervals.
    state_xmlfile_dir   # xml file directory
    state_xmlfile_name  # xml file name
    map_file_dir        # mapping file directory
    map_file_name       # mapping file name
    hybrid_file         # hybrid file name 
    topo_file           # topography file name

Other required environment/scripts:
    1. load CDAT
    2. vertical_interpolation_code (vertical_interpolation_code.f90)   
    3. data_era5.py
    4. get_ab_era5.py

Optional arguments:
    -ft=$ | --file-tag=$  # outfile = 'era5' + <file_tag> + '_' + \
                            <start> + '_' + <end> + '_' [FILE_TAG_DEF].

Notes:
    * If the first character of hybrid_file or topo_file is not a '/',
      INFILE_PREF will be prepended to them.

    * 20200108, to run this script on NERSC Cori, one has to do the following
      to load the gfortran library for the vertical_interpolation_code.so 
      to work. The idea is that the .so was compiled with gfortran so the 
      gfortran libraries have to be loaded while this script is running. The 
      default fortran is ifort so it has to be unloaded first.

      1) put the vertical_interpolation_code.so in the same directory with this 
         script
      2) module unload PrgEnv-intel/6.0.5
      3) module load PrgEnv-gnu
"""
#==============================================================================

import sys

import numpy
n = numpy

import cdms2
import cdtime
import regrid2

#path of the vertical_interpolation_code.so 
sys.path.insert(1, '/global/homes/h/hyma/CAPT/utils/zfortran')
import vertical_interpolation_code
vic = vertical_interpolation_code

# Access to era5 raw reanalysis data is through dataera5 class.
#--------------------------------------------------------------------
import data_era5
di = data_era5


INFILE_PREF = '/p/lscratchc/cuqrusr/ccsm3data/inputdata/'
OUTFILE_PREF = './'

FILE_TAG_DEF = '_ne30'

#DELTAT = 6
DELTAT = 1
GEE = 9.81
SURF_GEOPOTA = 'PHIS'
SURF_HT_IS_GEOPOT = True
SURF_PRES_SCALE = 1.0     # gds ps is in cb. cb*1000. = Pa ECMWF is Pa.
SURF_PRESA = 'PS'

FAC_DICT    = {'T': 1.0, 'U': 1.0, 'V': 1.0, 'Q': 1.0}
OFFSET_DICT = {'T': 0.0, 'U': 0.0, 'V': 0.0, 'Q': 0.0}

# Mappings of names & units from ECMWF to CCSM on staggered & unstaggered.
#-------------------------------------------------------------------------
#VAR_DICT = {'T': 'T', 'U': 'U', 'V': 'V', 'Q': 'QV'}
VAR_DICT = {'T': 'T', 'U': 'U', 'V': 'V', 'Q': 'Q'}

VAR_LIST = ['T', 'Q', 'U', 'V']


def _get_hybrid_vars(hybrid_file):
#---------------------------------
    """
    Returns various E3SM ne30 horizontal coordinate variables.

    hybrid_file  # TBD.
    """

    try:
        hfin = cdms2.open(hybrid_file)

    except:
        mystr = 'hybrid_file NOT FOUND IN _era52e3sm_ne30()!:\n    %s\n' % \
                hybrid_file
        raise Exception('\n  %s' % mystr)

    ntrk = hfin.variables['ntrk'].getValue()
    ntrm = hfin.variables['ntrm'].getValue()
    ntrn = hfin.variables['ntrn'].getValue()

    return hfin, ntrk, ntrm, ntrn


def _get_hybridv_vars(hybridv_file):
#-----------------------------------
    """
    Returns various E3SM ne30 vertical hybrid coordinate variables.

    hybridv_file  # 
    """

    try:
        # cdms2.open() does not permit access to non-variable dimension
        # parameters (e.g., hyai, hyam, hybi, hybm), so use cdms2.Cdunif.
        #----------------------------------------------------------------
        hvfin = cdms2.Cdunif.CdunifFile(hybridv_file)

    except:
        mystr = 'hybridv_file NOT FOUND IN _era52e3sm_ne30()!:\n    %s\n' % \
                hybridv_file
        raise Exception('\n  %s' % mystr)

    ilevs = hvfin.variables['ilev'][:]
    levs = hvfin.variables['lev'][:]

    # p = p0*hyai + hybi*ps (E3SM ne30 model's version of era5 ak/bk).
    #------------------------------------------------------------------
    hyai = hvfin.variables['hyai'][:]
    hyam = hvfin.variables['hyam'][:]
    hybi = hvfin.variables['hybi'][:]
    hybm = hvfin.variables['hybm'][:]

    p0 = hvfin.variables['P0'].getValue()

    if type(p0) == type(levs):
        p0 = p0[0]

    if p0 > 1.0e20:
        p0 = 100000.0

    hvfin.close()

    return ilevs, levs, hyai, hyam, hybi, hybm, p0


def _get_topo_vars(topo_file):
#-----------------------------
    """
    Returns various topography variables.

    topo_file  # TBD.
    """

    # Get terrain.
    #-------------
    try:
        tfin = cdms2.open(topo_file)

    except:
        mystr = 'topo_file NOT FOUND IN _era52e3sm_ne30()!:\n    %s\n' % \
                topo_file
        raise Exception('\n  %s' % mystr)

    lats = tfin('lat', raw=1, squeeze=1)
    lons = tfin('lon', raw=1, squeeze=1)

    ncol = len(lons)

    # Topo heights.
    #--------------
    phisfcm_in = tfin('PHIS', squeeze=1)
    print ('phisfcm_in type: ', type(phisfcm_in))
    print ('phisfcm_in shape: ', phisfcm_in.shape)

    return lats, lons, ncol, phisfcm_in, tfin


def _setup_outfile(cbaset, cendt, file_tag, hyai, hyam, hybi, hybm, ilevs,
                   lats, levs, lons, ncol, ntrk, ntrm, ntrn, p0):
#-------------------------------------------------------------------------
    """
    Sets up E3SM ne30 output file (same format as E3SM .i files);
    returns various CDMS variables.

    cbaset  # Starting time of data in E3SM ne30 output file.
    cendt   # Ending   time of data in E3SM ne30 output file.
    """

    start = '%04d%02d%02d%02d' % \
            (cbaset.year, cbaset.month, cbaset.day, cbaset.hour)

    end = '%04d%02d%02d%02d' % \
          (cendt.year, cendt.month, cendt.day, cendt.hour)

    outfile_name = 'era5' + file_tag + '_' + start + '_' + end + '_'

    for var in VAR_LIST:
        outfile_name = outfile_name + var

    outfile_name = OUTFILE_PREF + outfile_name + '.nc'

    print ('outfile_name: ', outfile_name)

    fout = cdms2.Cdunif.CdunifFile(outfile_name, 'w')

    fout.createDimension('ncol', ncol)
    fout.createDimension('lev', len(levs))
    fout.createDimension('ilev', len(ilevs))
    fout.createDimension('time', None)

    dim_tuple = ('time', 'lev', 'ncol')

    if 'T' in VAR_LIST:
        t = fout.createVariable('T', 'f', dim_tuple)
        t.long_name = 'temperature'
        t.units = 'K'

    if 'Q' in VAR_LIST:
        q = fout.createVariable('Q', 'f', dim_tuple)
        q.long_name = 'specific humidity'
        q.units = 'kg/kg'

    if 'U' in VAR_LIST:
        u = fout.createVariable('U', 'f', dim_tuple)
        u.long_name = 'zonal wind component'
        u.units = 'm/s'

    if 'V' in VAR_LIST:
        v = fout.createVariable('V', 'f', dim_tuple)
        v.long_name = 'meridional wind component'
        v.units = 'm/s'

    # Create surface pressure grids.
    #-------------------------------
    ps = fout.createVariable('PS', 'f', ('time', 'ncol'))
    ps.long_name = 'surface pressure'
    ps.units = 'Pa'

    phis = fout.createVariable('PHIS', 'f', ('time', 'ncol'))
    phis.long_name = 'surface geopotential'
    phis.units = 'm2/s2'

    lon = fout.createVariable('lon', 'd', ('ncol',))
    lon[:] = lons
    lon.long_name = 'longitude'
    lon.units = 'degrees_east'

    lat = fout.createVariable('lat', 'd', ('ncol',))
    lat[:] = lats
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'

    lev = fout.createVariable('lev', 'd', ('lev',))
    lev[:] = levs
    lev.long_name = 'hybrid level at layer midpoints (100*(A+B))'
    lev.units = 'hybrid_sigma_pressure'
    lev.positive = 'down'
    lev.A_var = 'hyam'
    lev.B_var = 'hybm'
    lev.P0_var = 'P0'
    lev.PS_var = 'PS'
    lev.edges = 'ilev'

    ilev = fout.createVariable('ilev', 'd', ('ilev',))
    ilev[:] = ilevs
    ilev.long_name = 'hybrid level at layer interfaces (1000*(A+B))'
    ilev.units = 'hybrid_sigma_pressure'
    ilev.positive = 'down'
    ilev.A_var = 'hyai'
    ilev.B_var = 'hybi'
    ilev.P0_var = 'P0'
    ilev.PS_var = 'PS'

    hyamv = fout.createVariable('hyam', 'd', ('lev',))
    hyamv[:] = hyam
    hyamv.long_name = 'hybrid A coefficient at layer midpoints'

    hybmv = fout.createVariable('hybm', 'd', ('lev',))
    hybmv[:] = hybm
    hybmv.long_name = 'hybrid B coefficient at layer midpoints'

    hyaiv = fout.createVariable('hyai', 'd', ('ilev',))
    hyaiv[:] = hyai
    hyaiv.long_name = 'hybrid A coefficient at layer interfaces'

    hybiv = fout.createVariable('hybi', 'd', ('ilev',))
    hybiv[:] = hybi
    hybiv.long_name = 'hybrid B coefficient at layer interfaces'

    p0v = fout.createVariable('P0', 'f', ())
    p0v.assignValue(p0)
    p0v.long_name = 'reference pressure'
    p0v.units = 'Pa'

    ntrmv = fout.createVariable('ntrm', 'i', ())
    ntrmv.assignValue(ntrm)
    ntrmv.long_name = 'spectral truncation parameter M'

    ntrnv = fout.createVariable('ntrn', 'i', ())
    ntrnv.assignValue(ntrn)
    ntrnv.long_name = 'spectral truncation parameter N'

    ntrkv = fout.createVariable('ntrk', 'i', ())
    ntrkv.assignValue(ntrk)
    ntrkv.long_name = 'spectral truncation parameter K'

    time = fout.createVariable('time', 'f', ('time',))
    time.long_name = 'simulation time'

    new_units = 'days since ' + \
                str(cbaset.year) + '-' + \
                str(cbaset.month) + '-' + \
                str(cbaset.day) + ' ' + \
                str(cbaset.hour) + ':00:00'

    time.units = new_units
    time.calendar = 'gregorian'

    ndbase = fout.createVariable('ndbase', 'i', ())
    ndbase.assignValue((year * 10000) + (month * 100) + day)
    ndbase.long_name = 'base day for this case'

    nsbase = fout.createVariable('nsbase', 'i', ())
    nsbase.assignValue(hour * 60 * 60)
    nsbase.long_name = 'seconds to complete base day'

    nbdate = fout.createVariable('nbdate', 'i', ())
    nbdate.assignValue((year * 10000) + (month * 100) + day)
    nbdate.long_name = 'base date (YYYYMMDD)'

    nbsec = fout.createVariable('nbsec', 'i', ())
    nbsec.assignValue(hour * 60 * 60)
    nbsec.long_name = 'seconds to complete base date'

    date = fout.createVariable('date', 'i', ('time',))
    date.long_name = 'current date (YYYYMMDD)'

    datesec = fout.createVariable('datesec', 'i', ('time',))
    datesec.long_name = 'seconds to complete current date'
    datesec.units = 'seconds'

    return date, datesec, fout, lev, phis, ps, q, t, time, u, v


def _get_map(map_file_dir, map_file_name):
#-----------------------------------------
    """
    Returns fv2se_map.

    map_file_dir   # TBD.
    map_file_name  #
    """

    map_file = map_file_dir + '/' + map_file_name

    try:
        fv2se_map = cdms2.open(map_file)

    except:
        mystr = 'map_file NOT FOUND IN _era52e3sm_ne30()!:\n    %s\n' % \
                map_file
        raise Exception('\n  %s' % mystr)

    return fv2se_map


def _get_phisfcr_in(ddi, rgfa2m):
#--------------------------------
    """
    Returns era5 geopotential on E3SM ne30 grid.

    ddi     # era5 surface geopotential.
    rgfa2m  # era5 to E3SM SE regrid function.
    """

    phisfca = ddi.get_const_data(SURF_GEOPOTA)
    print ('phisfca type: ', type(phisfca))
    print ('phisfca shape: ', phisfca.shape)

    if not SURF_HT_IS_GEOPOT:
        phisfca = phisfca * GEE

    phisfcr_in = rgfa2m(cdms2.asVariable(phisfca))
    print ('phisfcr_in type: ', type(phisfcr_in))
    print ('phisfcr_in shape: ', phisfcr_in.shape)
    print (phisfcr_in.iscontiguous())

    return phisfcr_in


def _do_one_timestep(ctime, ddi, hyam, hybm, levs, ncol, p0, phis, phisfcm_in,
                     phisfcr_in, plat, plevm, plevr, plevrp1, plon, ps, q,
                     rgfa2m, t, tidx, u, v):
#------------------------------------------------------------------------------
    """
    Does one timestep of putting era5 raw reanalysis state variables
    onto E3SM ne30 grid.

    ctime  # cdtime component time.

    Suffix a=era5 raw reanalysis,
           m=E3SM ne30 model,
           r=era5 raw reanalysis regridded to E3SM ne30 model
             (but still era5).
    """

    psfca = ddi.get_state_data(SURF_PRESA, ctime)

    psfcr = rgfa2m(cdms2.asVariable(psfca))
    psfcr = psfcr * SURF_PRES_SCALE
    print ('psfcr shape: ', psfcr.shape)
    print ('psfcr min/max: ', psfcr.min(), psfcr.max())

    ta = ddi.get_state_data(VAR_DICT['T'], ctime)
    tr = rgfa2m(cdms2.asVariable(ta))
    print ('tr shape: ', tr.shape)

    pressr = ddi.get_state_data('pressure_mid', ctime)
    pressr = rgfa2m(cdms2.asVariable(pressr))

    pressir = ddi.get_state_data('pressure_interface', ctime)
    pressir = rgfa2m(cdms2.asVariable(pressir))

    # In order to use the same Fortran, we need to dummy in two dimensions,
    # (lev, ncol, 1) analogous to (lev, lat, lon), that the Fortran is
    # expecting.
    #----------------------------------------------------------------------
    ncol_lat = ncol
    ncol_lon = 1

    # Fortran wants numpy arrays.
    #----------------------------
    phisfcm = phisfcm_in.filled()
    phisfcr = phisfcr_in.filled()

    psfcr = psfcr.filled()
    tr = tr.filled()

    nlevs = tr.shape[0]
    nlevsp1 = nlevs + 1

    phisfcm = phisfcm.reshape((ncol_lat, ncol_lon))
    phisfcr = phisfcr.reshape((ncol_lat, ncol_lon))

    psfcr = psfcr.reshape((ncol_lat, ncol_lon))
    tr = tr.reshape((nlevs, ncol_lat, ncol_lon))

    pressr = pressr.reshape((nlevs, ncol_lat, ncol_lon))
    pressr = pressr.filled()

    pressir = pressir.reshape((nlevsp1, ncol_lat, ncol_lon))
    pressir = pressir.filled()

    print ('plevr, plevrp1, plat, plon: ', plevr, plevrp1, plat, plon)
    print ('tr pressr pressir shapes: ', \
          tr.shape, pressr.shape, pressir.shape)
    print ('phisfcr, phisfcm, psfcr shapes: ', \
          phisfcr.shape, phisfcm.shape, psfcr.shape)

    print ('Executing vic.')

    # Adjust lower level temperature and pressue using hydrostatic approximation
    # to account for differences between ERA5 and E3SM topography.
    #
    # In Fortran, the arrays are defined:
    #   real*8 press_m(plon,plat,plev)   ! analysis pressures
    #   real*8 press_i(plon,plat,plevp1) ! analysis pressures (interfaces)
    # Since the lat, lon are arbitrary, we just want to get the vertical in
    # the correct place - transpose since lev in last index.
    #----------------------------------------------------------------------
    psfcm = vic.psadj(plevr, plevrp1, plat, plon,
                      n.transpose(tr[:,:,:], (2, 1, 0)),
                      n.transpose(pressr[:,:,:], (2, 1, 0)),
                      n.transpose(pressir[:,:,:], (2, 1, 0)),
                      n.transpose(phisfcr),
                      n.transpose(phisfcm),
                      n.transpose(psfcr))

    print ('psfcm shape: ', psfcm.shape)

    psfcm = n.transpose(psfcm)
    print ('psfcm shape: ', psfcm.shape)

    # Write new surface pressure to file.
    #------------------------------------
    phis[tidx,:] = (phisfcm.astype(n.float32)).reshape((ncol,))

    ps[tidx,:] = (psfcm.astype(n.float32)).reshape((ncol,))

    # Now compute vertical pressures using new surface pressure.
    #-----------------------------------------------------------
    new_shape = (len(hyam), psfcm.shape[0], psfcm.shape[1])

    pressm = (hyam[:,n.newaxis] * p0) + n.outer(hybm, psfcm)
    pressm.shape = new_shape
    print ('pressm shape: ', pressm.shape)

    # Need to flip 3D arrays to new.
    #-------------------------------
    print ('pressr, pressir, pressm shapes: ', \
          pressr.shape, pressir.shape, pressm.shape)

    pressir = n.transpose(pressir[:,:,:], (0, 2, 1))
    pressr = n.transpose(pressr[:,:,:], (0, 2, 1))

    pressm = n.transpose(pressm, (0, 2, 1))

    print ('pressr, pressir, pressm shapes: ', \
          pressr.shape, pressir.shape, pressm.shape)

    print ('phisfcr, psfcr shapes: ', phisfcr.shape, psfcr.shape)

    phisfcr = n.transpose(phisfcr)
    psfcr = n.transpose(psfcr)

    print ('phisfcr, psfcr shapes: ', phisfcr.shape, psfcr.shape)

    for var in VAR_LIST:
        ta = ddi.get_state_data(VAR_DICT[var], ctime)
        tr = rgfa2m(cdms2.asVariable(ta))

        nlevsa = tr.shape[0]

        fac = FAC_DICT[var]
        offset = OFFSET_DICT[var]

        tr = tr.reshape((nlevsa, ncol, 1)).filled()
        tr = n.transpose(tr[:,:,:], (0, 2, 1))

        if var == 'T':
            print ('**T**')
            print ('tr, pressr, pressir shapes: ', \
                  tr.shape, pressr.shape, pressir.shape)
            print ('pressm, phisfcr, psfcr shapes: ', \
                  pressm.shape, phisfcr.shape, psfcr.shape)

            tint = vic.vert_quad_opt1(plevr, plevrp1, plevm, plat, plon, tr,
                                      pressr, pressir, pressm, phisfcr, psfcr,
                                      0)

        if var == 'Q':
            print ('**Q**')
            print ('tr, pressr, pressm shapes: ', \
                  tr.shape, pressr.shape, pressm.shape)

            tint = vic.vert_int_opt1(plat, plon, plevr, plevm,
                                     pressr, pressm, tr, 1)

        if (var == 'U') or (var == 'V'):
            print ('**U/V**')
            print ('tr, pressr, pressm shapes: ', \
                  tr.shape, pressr.shape, pressm.shape)

            tint = vic.vert_int_opt2(plat, plon, plevr, plevm,
                                     pressr, pressm, tr, 0)

        print ('tint shape: ', tint.shape)

        tint = (fac * n.transpose(tint, (0, 2, 1))) + offset

        print ('tint shape: ', tint.shape)

        if var == 'T':
            tint = (tint.astype(n.float32)).reshape((len(levs), ncol))
            t[tidx,:,:] = tint[:,:]

        if var == 'Q':
            # Ensure q is >= 0.0.
            #--------------------
            tint = n.where(n.less(tint, 0.0), 0.0, tint)
            tint = (tint.astype(n.float32)).reshape((len(levs), ncol))

            q[tidx,:,:] = tint[:,:]

        if var == 'U':
            tint = (tint.astype(n.float32)).reshape((len(levs), ncol))
            u[tidx,:,:] = tint[:,:]

        if var == 'V':
            tint = (tint.astype(n.float32)).reshape((len(levs), ncol))
            v[tidx,:,:] = tint[:,:]

    return


def _era52e3sm_ne30(year, month, day, hour, ntimes, state_xmlfile_dir,
                     state_xmlfile_name, map_file_dir, map_file_name,
                     hybrid_file, topo_file, file_tag):
#----------------------------------------------------------------------
    """
    See file header.
    """

    # Base time is set from first time slice selected.
    #-------------------------------------------------
    cbaset = cdtime.comptime(year, month, day, hour)
    cendt = cbaset.add(((ntimes - 1) * DELTAT), cdtime.Hours)

    ddi = di.dataera5(state_xmlfile_dir=state_xmlfile_dir,
                      state_xmlfile_name=state_xmlfile_name)

    lats, lons, ncol, phisfcm_in, tfin = \
        _get_topo_vars(topo_file)

    plat = ncol
    plon = 1

    hfin, ntrk, ntrm, ntrn = \
        _get_hybrid_vars(hybrid_file)

    hybridv_file = hybrid_file

    ilevs, levs, hyai, hyam, hybi, hybm, p0 = \
        _get_hybridv_vars(hybridv_file)

    date, datesec, fout, lev, phis, ps, q, t, time, u, v = \
        _setup_outfile(cbaset, cendt, file_tag, hyai, hyam, hybi, hybm, ilevs,
                       lats, levs, lons, ncol, ntrk, ntrm, ntrn, p0)

    # Now access ERA5 data to obtain grid/level parameters of the analysis.
    #-------------------------------------------------------------------------
    ta = ddi.get_state_data('T', cbaset)
    print ('ta shape: ', ta.shape)

    # Model levels.
    #--------------
    plevm = len(lev)

    # Analysis levels.
    #-----------------
    plevr = ta.shape[0]
    plevrp1 = plevr + 1

    # Set up regridder functions for all needed tranformations.
    #----------------------------------------------------------
    fv2se_map = _get_map(map_file_dir, map_file_name)

    rgfa2m = regrid2.readRegridder(fv2se_map)

    fv2se_map.close()

    phisfcr_in = _get_phisfcr_in(ddi, rgfa2m)

    # Time loop.
    #-----------
    for tidx in range(ntimes):
        ctime = cbaset.add((tidx * DELTAT), cdtime.Hours)
        print ('In time loop, tidx, ctime: ', tidx, ctime)

        time[tidx] = (tidx * DELTAT) / 24.0
        print ('time[tidx], DELTAT: ', time[tidx], DELTAT)

        date[tidx] = (ctime.year * 10000) + (ctime.month * 100) + ctime.day
        datesec[tidx] = (ctime.hour * 60 * 60) + (ctime.minute * 60) + \
                        ctime.second
        print ('date[tidx], datesec[tidx]: ', date[tidx], datesec[tidx])

        _do_one_timestep(ctime, ddi, hyam, hybm, levs, ncol, p0, phis,
                         phisfcm_in, phisfcr_in, plat, plevm, plevr, plevrp1,
                         plon, ps, q, rgfa2m, t, tidx, u, v)

        fout.sync()

    fout.close()
    hfin.close()
    tfin.close()

    return


if __name__ == '__main__':
#-------------------------
    args = sys.argv[1:]

    if len(args) < 11:
        print (__doc__)
        print ('\nFirst 11 arguments required - EXITING!\n')
        raise SystemExit(1)

    else:
        year = int(args[0])
        month = int(args[1])
        day = int(args[2])
        hour = int(args[3])
        ntimes = int(args[4])
        state_xmlfile_dir = args[5]
        state_xmlfile_name = args[6]
        map_file_dir = args[7]
        map_file_name = args[8]
        hybrid_file = args[9]
        topo_file = args[10]

        args = args[11:]

        if not hybrid_file.startswith('/'):
            hybrid_file = INFILE_PREF + hybrid_file

        if not topo_file.startswith('/'):
            topo_file = INFILE_PREF + topo_file

        file_tag = FILE_TAG_DEF

        while len(args) >= 1:
            arg = args[0]
            del args[0]

            if arg.startswith('-ft=') or \
                 arg.startswith('--file-tag='):
                file_tag = arg[(arg.find('=')+1):]

            else:
                print (__doc__)
                print ('Invalid option: %s\n' % arg)
                raise SystemExit(1)

        print (year, month, day, hour)
        print (ntimes)
        print (state_xmlfile_dir)
        print (state_xmlfile_name)
        print (map_file_dir)
        print (map_file_name)
        print (hybrid_file)
        print (topo_file)
        print (file_tag)

        _era52e3sm_ne30(year, month, day, hour, ntimes, state_xmlfile_dir,
                         state_xmlfile_name, map_file_dir, map_file_name,
                         hybrid_file, topo_file, file_tag)
