#! /usr/bin/env python

#==============================================================================
#   interim2cesm_fv.py
"""
Python script.

This script interpolates Interim data to a CESM Finite Volume (FV) grid.

Interactive usage:
    interim2cesm_fv.py <year month day hour
                        ntimes
                        state_xmlfile_dir state_xmlfile_name
                        hybrid_file       topo_file>
                       [optional arguments]

Required arguments:
    year                # Start year.
    month               # Start month.
    day                 # Start day.
    hour                # Start hour.
    ntimes              # Number of 6 hour time intervals.
    state_xmlfile_dir   # TBD.
    state_xmlfile_name  #
    hybrid_file         #
    topo_file           #

Optional arguments:
    -ft=$ | --file-tag=$  # outfile = 'interim' + <file_tag> + '_' + \
                            <start> + '_' + <end> + '_' [FILE_TAG_DEF].

Notes:
    * If the first character of hybrid_file or topo_file is not a '/',
      INFILE_PREF will be prepended to them.
"""
#==============================================================================

import sys

import numpy
n = numpy

import cdms2
import cdtime
import regrid2

import vertical_interpolation_code
vic = vertical_interpolation_code

# Access to Interim raw reanalysis data is through dataInterim class.
#--------------------------------------------------------------------
import data_interim
di = data_interim

import filter_field
ff = filter_field


INFILE_PREF = '/p/lscratchc/cuqrusr/ccsm3data/inputdata/'
OUTFILE_PREF = './'

FILE_TAG_DEF = '_fv'

DELTAT = 6
FIL_IDX = 1
FIL_TYPE = 'gerrity'
GEE = 9.81
SURF_GEOPOTA = 'PHIS'
SURF_HT_IS_GEOPOT = True
SURF_PRES_SCALE = 1.0     # gds ps is in cb. cb*1000. = Pa ECMWF is Pa.
SURF_PRESA = 'PS'

FAC_DICT    = {'T': 1.0, 'U': 1.0, 'V': 1.0, 'Q': 1.0, 'US': 1.0, 'VS': 1.0}
OFFSET_DICT = {'T': 0.0, 'U': 0.0, 'V': 0.0, 'Q': 0.0, 'US': 0.0, 'VS': 0.0}

# Mappings of names & units from GDAS to CCSM on staggered & unstaggered.
#------------------------------------------------------------------------
VAR_DICT = {'T': 'T', 'U': 'U', 'V': 'V', 'Q': 'QV', 'US': 'U', 'VS': 'V'}

VAR_LIST = ['T', 'Q', 'US', 'VS']


def _get_hybrid_vars(hybrid_file):
#---------------------------------
    """
    Returns various CESM FV horizontal coordinate variables.

    hybrid_file  # TBD.
    """

    # Need to add the lat/lon names of the staggered grid to cdms,
    # so it can manipulate the grid.
    #-------------------------------------------------------------
    cdms2.axis.longitude_aliases.append('slon')
    cdms2.axis.latitude_aliases.append('slat')

    try:
        hfin = cdms2.open(hybrid_file)

    except:
        mystr = 'hybrid_file NOT FOUND IN _interim2cesm_fv()!:\n    %s\n' % \
                hybrid_file
        raise Exception('\n  %s' % mystr)

    gwts = hfin.variables['gw'][:]

    # Staggered lats and lons from staggered winds;
    # this only works if aliases are set above.
    #----------------------------------------------
    model_grid_us = hfin.variables['US'].getGrid()
    model_grid_vs = hfin.variables['VS'].getGrid()

    ntrk = hfin.variables['ntrk'].getValue()
    ntrm = hfin.variables['ntrm'].getValue()
    ntrn = hfin.variables['ntrn'].getValue()

    slats = model_grid_us.getLatitude()[:]
    slons = model_grid_vs.getLongitude()[:]

    w_stagwts = hfin.variables['w_stag'][:]

    return gwts, hfin, model_grid_us, model_grid_vs, ntrk, ntrm, ntrn, slats, \
           slons, w_stagwts


def _get_hybridv_vars(hybridv_file):
#-----------------------------------
    """
    Returns various CESM FV vertical hybrid coordinate variables.

    hybridv_file  # TBD.
    """

    try:
        # cdms2.open() does not permit access to non-variable dimension
        # parameters (e.g., hyai, hyam, hybi, hybm), so use cdms2.Cdunif.
        #----------------------------------------------------------------
        hvfin = cdms2.Cdunif.CdunifFile(hybridv_file)

    except:
        mystr = 'hybridv_file NOT FOUND IN _interim2cesm_fv()!:\n    %s\n' % \
                hybridv_file
        raise Exception('\n  %s' % mystr)

    ilevs = hvfin.variables['ilev'][:]
    levs = hvfin.variables['lev'][:]

    # p = p0*hyai + hybi*ps (CESM FV model's version of Interim ak/bk).
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
        mystr = 'topo_file NOT FOUND IN _interim2cesm_fv()!:\n    %s\n' % \
                topo_file
        raise Exception('\n  %s' % mystr)

    # Grid on which the topography is defined.
    #-----------------------------------------
    model_grid = tfin.variables['PHIS'].getGrid()

    lats = model_grid.getLatitude()[:]
    lons = model_grid.getLongitude()[:]

    # Topo heights.
    #--------------
    phisfcm_in = tfin.variables['PHIS'][:]
    print 'phisfcm_in type: ', type(phisfcm_in)
    print 'phisfcm_in shape: ', phisfcm_in.shape

    return lats, lons, model_grid, phisfcm_in, tfin


def _setup_outfile(cbaset, cendt, file_tag, gwts, hyai, hyam, hybi, hybm,
                   ilevs, lats, levs, lons, ntrk, ntrm, ntrn, p0, slats,
                   slons, w_stagwts):
#------------------------------------------------------------------------
    """
    Sets up CESM FV output file (same format as CESM .i files);
    returns various CDMS variables.

    cbaset  # Starting time of data in CESM FV output file.
    cendt   # Ending   time of data in CESM FV output file.
    """

    start = '%04d%02d%02d%02d' % \
            (cbaset.year, cbaset.month, cbaset.day, cbaset.hour)

    end = '%04d%02d%02d%02d' % \
          (cendt.year, cendt.month, cendt.day, cendt.hour)

    outfile_name = 'interim' + file_tag + '_' + start + '_' + end + '_'

    for var in VAR_LIST:
        outfile_name = outfile_name + var

    outfile_name = OUTFILE_PREF + outfile_name + '.nc'

    print 'outfile_name: ', outfile_name

    fout = cdms2.Cdunif.CdunifFile(outfile_name, 'w')

    fout.createDimension('lon', len(lons))
    fout.createDimension('lat', len(lats))
    fout.createDimension('slon', len(slons))
    fout.createDimension('slat', len(slats))

    fout.createDimension('lev', len(levs))
    fout.createDimension('ilev', len(ilevs))
    fout.createDimension('time', None)

    dim_tuple = ('time', 'lev', 'lat', 'lon')
    dim_tuple_us = ('time', 'lev', 'slat', 'lon')
    dim_tuple_vs = ('time', 'lev', 'lat', 'slon')

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

    if 'US' in VAR_LIST:
        us = fout.createVariable('US', 'f', dim_tuple_us)
        us.long_name = 'zonal wind, staggered'
        us.units = 'm/s'

        u = fout.createVariable('U', 'f', dim_tuple)
        u.long_name = 'zonal wind component'
        u.units = 'm/s'

    if 'V' in VAR_LIST:
        v = fout.createVariable('V', 'f', dim_tuple)
        v.long_name = 'meridional wind component'
        v.units = 'm/s'

    if 'VS' in VAR_LIST:
        vs = fout.createVariable('VS', 'f', dim_tuple_vs)
        vs.long_name = 'meridional wind, staggered'
        vs.units = 'm/s'

        v = fout.createVariable('V', 'f', dim_tuple)
        v.long_name = 'meridional wind component'
        v.units = 'm/s'

    # Create two staggered surface pressure grids - for convenience.
    #---------------------------------------------------------------
    ps = fout.createVariable('PS', 'f', ('time', 'lat', 'lon'))
    ps.long_name = 'surface pressure'
    ps.units = 'Pa'

    psus = fout.createVariable('PSUS', 'f', ('time', 'slat', 'lon'))
    psus.long_name = 'surface pressure on US grid'
    psus.units = 'Pa'

    psvs = fout.createVariable('PSVS', 'f', ('time', 'lat', 'slon'))
    psvs.long_name = 'surface pressure on VS grid'
    psvs.units = 'Pa'

    phis = fout.createVariable('PHIS', 'f', ('time', 'lat', 'lon'))
    phis.long_name = 'surface geopotential'
    phis.units = 'm2/s2'

    lon = fout.createVariable('lon', 'f', ('lon',))
    lon[:] = n.array(lons, n.float32)
    lon.long_name = 'longitude'
    lon.units = 'degrees_east'

    lat = fout.createVariable('lat', 'f', ('lat',))
    lat[:] = n.array(lats, n.float32)
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'

    slon = fout.createVariable('slon', 'f', ('slon',))
    slon[:] = n.array(slons, n.float32)
    slon.long_name = 'staggered longitude'
    slon.units = 'degrees_east'

    slat = fout.createVariable('slat', 'f', ('slat',))
    slat[:] = n.array(slats, n.float32)
    slat.long_name = 'staggered latitude'
    slat.units = 'degrees_north'

    lev = fout.createVariable('lev', 'f', ('lev',))
    lev[:] = n.array(levs, n.float32)
    lev.long_name = 'hybrid level at layer midpoints (100*(A+B))'
    lev.units = 'hybrid_sigma_pressure'
    lev.positive = 'down'
    lev.A_var = 'hyam'
    lev.B_var = 'hybm'
    lev.P0_var = 'P0'
    lev.PS_var = 'PS'
    lev.edges = 'ilev'

    ilev = fout.createVariable('ilev', 'f', ('ilev',))
    ilev[:] = n.array(ilevs, n.float32)
    ilev.long_name = 'hybrid level at layer interfaces (1000*(A+B))'
    ilev.units = 'hybrid_sigma_pressure'
    ilev.positive = 'down'
    ilev.A_var = 'hyai'
    ilev.B_var = 'hybi'
    ilev.P0_var = 'P0'
    ilev.PS_var = 'PS'

    hyamv = fout.createVariable('hyam', 'f', ('lev',))
    hyamv[:] = n.array(hyam, n.float32)
    hyamv.long_name = 'hybrid A coefficient at layer midpoints'

    hybmv = fout.createVariable('hybm', 'f', ('lev',))
    hybmv[:] = n.array(hybm, n.float32)
    hybmv.long_name = 'hybrid B coefficient at layer midpoints'

    hyaiv = fout.createVariable('hyai', 'f', ('ilev',))
    hyaiv[:] = n.array(hyai, n.float32)
    hyaiv.long_name = 'hybrid A coefficient at layer interfaces'

    hybiv = fout.createVariable('hybi', 'f', ('ilev',))
    hybiv[:] = n.array(hybi, n.float32)
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

    gw = fout.createVariable('gw', 'f', ('lat',))
    gw[:] = n.array(gwts.filled(), n.float32)
    gw.long_name = 'gauss weights'

    w_stag = fout.createVariable('w_stag', 'f', ('slat',))
    w_stag[:] = n.array(w_stagwts.filled(), n.float32)
    w_stag.long_name = 'staggered latitude weights'

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

    return date, datesec, fout, lev, phis, ps, psus, psvs, q, t, time, u, us, \
           v, vs


def _get_phisfcr_in(ddi, rgfa2m):
#--------------------------------
    """
    Returns Interim geopotential on CESM FV grid.

    ddi     # Interim surface geopotential.
    rgfa2m  # Interim to CESM FV regrid function.
    """

    phisfca = ddi.get_const_data(SURF_GEOPOTA)
    print 'phisfca type: ', type(phisfca)
    print 'phisfca shape: ', phisfca.shape

    if not SURF_HT_IS_GEOPOT:
        phisfca = phisfca * GEE

    print 'phisfca type: ', type(phisfca)
    print 'phisfca shape: ', phisfca.shape
    print 'rgfa2m: ', rgfa2m

    phisfcr_in = rgfa2m(phisfca.filled())
    print 'phisfcr_in type: ', type(phisfcr_in)
    print 'phisfcr_in shape: ', phisfcr_in.shape
    print phisfcr_in.iscontiguous()

    return phisfcr_in


def _do_one_timestep(ctime, ddi, hyam, hybm, p0, phis, phisfcm_in, phisfcr_in,
                     plat, plevm, plevr, plevrp1, plon, ps, pslat, pslon, psus,
                     psvs, q, rgfa2m, rgfa2mus, rgfa2mvs, rgfm2us, rgfm2vs,
                     rgfus2m, rgfvs2m, t, tidx, u, us, v, vs):
#------------------------------------------------------------------------------
    """
    Does one timestep of putting Interim raw reanalysis state variables
    onto CESM FV grid.

    ctime  # cdtime component time.

    Suffix a=Interim raw reanalysis,
           m=CESM FV model,
           r=Interim raw reanalysis regridded to CESM FV model
             (but still Interim).
    """

    psfca = ddi.get_state_data(SURF_PRESA, ctime)

    psfcr = rgfa2m(psfca)
    psfcr = ff.filter_field(psfcr.filled(), FIL_TYPE, FIL_IDX)
    psfcr = psfcr * SURF_PRES_SCALE
    print 'psfcr shape: ', psfcr.shape
    print 'psfcr min/max: ', psfcr.min(), psfcr.max()

    ta = ddi.get_state_data(VAR_DICT['T'], ctime)
    tr = rgfa2m(ta).filled()
    print 'tr shape: ', tr.shape

    print 'plevr: ', plevr

    for ii in range(plevr):
        tr[ii,:,:] = ff.filter_field(tr[ii,:,:], FIL_TYPE, FIL_IDX)

    print 'tr[0] min/max: ', tr[0].min(), tr[0].max()

    pressr = ddi.get_state_data('pressure_mid', ctime,
                                regrid_func=rgfa2m)

    pressir = ddi.get_state_data('pressure_interface', ctime,
                                 regrid_func=rgfa2m)

    print 'plevr, plevrp1, plat, plon: ', plevr, plevrp1, plat, plon
    print 'tr pressr pressir shapes: ', \
          tr.shape, pressr.shape, pressir.shape
    print 'phisfcr_in, phisfcm_in, psfcr shapes: ', \
          phisfcr_in.shape, phisfcm_in.shape, psfcr.shape

    print 'Executing vic.'

    # Adjust lower level temperature and pressue to account for
    # differences between Interim and CESM topography.
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
                      n.transpose(phisfcr_in),
                      n.transpose(phisfcm_in),
                      n.transpose(psfcr))

    print 'psfcm shape: ', psfcm.shape

    psfcm = n.transpose(psfcm)
    print 'psfcm shape: ', psfcm.shape

    psfcm = n.log(psfcm)
    psfcm = ff.filter_field(psfcm, FIL_TYPE, FIL_IDX)

    psfcm = n.exp(psfcm)

    # Write new surface pressure to file.
    #------------------------------------
    phis[tidx,:,:] = n.array(phisfcm_in, n.float32)

    ps[tidx,:,:] = n.array(psfcm, n.float32)

    # Now compute vertical pressures using new surface pressure.
    #-----------------------------------------------------------
    new_shape = (len(hyam), psfcm.shape[0], psfcm.shape[1])

    pressm = (hyam[:,n.newaxis] * p0) + n.outer(hybm, psfcm)
    pressm.shape = new_shape
    print 'pressm shape: ', pressm.shape

    # At this point we have surface pressure defined for the model &
    # analysis; we choose to interpolate these final products to the
    # staggered grids - it is only a small shift so further smoothing is
    # not needed.
    #-------------------------------------------------------------------
    psfcmus = rgfm2us(psfcm)
    psfcmvs = rgfm2vs(psfcm)

    psus[tidx,:,:] = psfcmus
    psvs[tidx,:,:] = psfcmvs

    new_shape = (len(hyam), psfcmus.shape[0], psfcmus.shape[1])

    pressmus = (hyam[:,n.newaxis] * p0) + n.outer(hybm, psfcmus)
    pressmus.shape = new_shape

    new_shape = (len(hyam), psfcmvs.shape[0], psfcmvs.shape[1])

    pressmvs = (hyam[:,n.newaxis] * p0) + n.outer(hybm, psfcmvs)
    pressmvs.shape = new_shape

    psfcrus = rgfm2us(psfcr)
    psfcrvs = rgfm2vs(psfcr)

    pressrus = ddi.get_state_data('pressure_mid', ctime,
                                  regrid_func=rgfa2mus)
    pressrvs = ddi.get_state_data('pressure_mid', ctime,
                                  regrid_func=rgfa2mvs)

    # Need to flip the 3D arrays to new.
    #-----------------------------------
    print 'pressr, pressir, pressm shapes: ', \
          pressr.shape, pressir.shape, pressm.shape

    pressir = n.transpose(pressir[:,:,:], (0, 2, 1))
    pressr = n.transpose(pressr[:,:,:], (0, 2, 1))

    pressm = n.transpose(pressm, (0, 2, 1))

    print 'pressr, pressir, pressm shapes: ', \
          pressr.shape, pressir.shape, pressm.shape

    print 'phisfcr_in, psfcr shapes: ', phisfcr_in.shape, psfcr.shape

    phisfcr = n.transpose(phisfcr_in)
    psfcr = n.transpose(psfcr)

    print 'phisfcr, psfcr shapes: ', phisfcr.shape, psfcr.shape

    pressrus = n.transpose(pressrus[:,:,:], (0, 2, 1))
    pressrvs = n.transpose(pressrvs[:,:,:], (0, 2, 1))

    pressmus = n.transpose(pressmus, (0, 2, 1))
    pressmvs = n.transpose(pressmvs, (0, 2, 1))

    psfcrus = n.transpose(psfcrus)
    psfcrvs = n.transpose(psfcrvs)

    for var in ['T', 'Q', 'US', 'VS']:
        ta = ddi.get_state_data(VAR_DICT[var], ctime)

        if var == 'US':
            tr = rgfa2mus(ta).filled()
        elif var == 'VS':
            tr = rgfa2mvs(ta).filled()
        else:
            tr = rgfa2m(ta).filled()

        for ii in range(plevr):
            tr[ii,:,:] = ff.filter_field(tr[ii,:,:], FIL_TYPE, FIL_IDX)

        fac = FAC_DICT[var]
        offset = OFFSET_DICT[var]

        tr = n.transpose(tr[:,:,:], (0, 2, 1))

        if var == 'T':
            print '**T**'
            print 'tr, pressr, pressir shapes: ', \
                  tr.shape, pressr.shape, pressir.shape
            print 'pressm, phisfcr, psfcr shapes: ', \
                  pressm.shape, phisfcr.shape, psfcr.shape

            tint = vic.vert_quad_opt1(plevr, plevrp1, plevm, plat, plon, tr,
                                      pressr, pressir, pressm, phisfcr, psfcr,
                                      0)

        if var == 'Q':
            print '**Q**'
            print 'tr, pressr, pressm shapes: ', \
                  tr.shape, pressr.shape, pressm.shape

            tint = vic.vert_int_opt1(plat, plon, plevr, plevm,
                                     pressr, pressm, tr, 1)

        if (var == 'U') or (var == 'V'):
            print '**U/V**'
            print 'tr, pressr, pressm shapes: ', \
                  tr.shape, pressr.shape, pressm.shape

            tint = vic.vert_int_opt2(plat, plon, plevr, plevm,
                                     pressr, pressm, tr, 0)

        if var == 'US':
            print '**US**'
            print 'tr, pressrus, pressmus shapes: ', \
                  tr.shape, pressrus.shape, pressmus.shape

            tint = vic.vert_int_opt2(pslat, plon, plevr, plevm,
                                     pressrus, pressmus, tr, 0)

        if var == 'VS':
            print '**VS**'
            print 'tr, pressrvs, pressmvs shapes: ', \
                  tr.shape, pressrvs.shape, pressmvs.shape

            tint = vic.vert_int_opt2(plat, pslon, plevr, plevm,
                                     pressrvs, pressmvs, tr, 0)

        print 'tint shape: ', tint.shape

        tint = (fac * n.transpose(tint, (0, 2, 1))) + offset

        print 'tint shape: ', tint.shape

        if var == 'T':
            t[tidx,:,:,:] = n.array(tint, n.float32)

        if var == 'Q':
            # Ensure q is >= 0.0.
            #--------------------
            tint = n.where(n.less(tint, 0.0), 0.0, tint)

            q[tidx,:,:,:] = n.array(tint, n.float32)

        if var == 'U':
            u[tidx,:,:,:] = n.array(tint, n.float32)

        if var == 'V':
            v[tidx,:,:,:] = n.array(tint, n.float32)

        if var == 'US':
            us[tidx,:,:,:] = n.array(tint, n.float32)
            uu = rgfus2m(tint)
            u[tidx,:,:,:] = n.array(uu, n.float32)

        if var == 'VS':
            vs[tidx,:,:,:] = n.array(tint, n.float32)
            vv = rgfvs2m(tint)
            v[tidx,:,:,:] = n.array(vv, n.float32)

    return


def _interim2cesm_fv(year, month, day, hour, ntimes, state_xmlfile_dir,
                     state_xmlfile_name, hybrid_file, topo_file, file_tag):
#--------------------------------------------------------------------------
    """
    See file header.
    """

    # Base time is set from first time slice selected.
    #-------------------------------------------------
    cbaset = cdtime.comptime(year, month, day, hour)
    cendt = cbaset.add(((ntimes - 1) * DELTAT), cdtime.Hours)

    ddi = di.dataInterim(state_xmlfile_dir=state_xmlfile_dir,
                         state_xmlfile_name=state_xmlfile_name)

    lats, lons, model_grid, phisfcm_in, tfin = \
        _get_topo_vars(topo_file)

    plat = len(lats)
    plon = len(lons)

    gwts, hfin, model_grid_us, model_grid_vs, ntrk, ntrm, ntrn, slats, slons, \
    w_stagwts = \
        _get_hybrid_vars(hybrid_file)

    pslat = len(slats)
    pslon = len(slons)

    hybridv_file = hybrid_file

    ilevs, levs, hyai, hyam, hybi, hybm, p0 = \
        _get_hybridv_vars(hybridv_file)

    date, datesec, fout, lev, phis, ps, psus, psvs, q, t, time, u, us, v, \
    vs = \
        _setup_outfile(cbaset, cendt, file_tag, gwts, hyai, hyam, hybi, hybm, \
                       ilevs, lats, levs, lons, ntrk, ntrm, ntrn, p0, slats, \
                       slons, w_stagwts)

    # Now access Interim data to obtain grid/level parameters of the analysis.
    #-------------------------------------------------------------------------
    ta = ddi.get_state_data('T', cbaset)
    print 'ta shape: ', ta.shape

    gdas_grid = ddi.grid

    # Model levels.
    #--------------
    plevm = len(lev)

    # Analysis levels.
    #-----------------
    plevr = ta.shape[0]
    plevrp1 = plevr + 1

    # Set up regridder functions for all needed tranformations.
    #----------------------------------------------------------
    rgfa2m = regrid2.Regridder(gdas_grid, model_grid)

    rgfa2mus = regrid2.Regridder(gdas_grid, model_grid_us)
    rgfa2mvs = regrid2.Regridder(gdas_grid, model_grid_vs)

    rgfm2us = regrid2.Regridder(model_grid, model_grid_us)
    rgfm2vs = regrid2.Regridder(model_grid, model_grid_vs)

    rgfus2m = regrid2.Regridder(model_grid_us, model_grid)
    rgfvs2m = regrid2.Regridder(model_grid_vs, model_grid)

    phisfcm_in = phisfcm_in.filled()
    phisfcm_in = n.squeeze(phisfcm_in)

    phisfcr_in = _get_phisfcr_in(ddi, rgfa2m)
    phisfcr_in = ff.filter_field(phisfcr_in.filled(), FIL_TYPE, FIL_IDX)

    # Time loop.
    #-----------
    for tidx in range(ntimes):
        ctime = cbaset.add((tidx * DELTAT), cdtime.Hours)
        print 'In time loop, tidx, ctime: ', tidx, ctime

        time[tidx] = (tidx * DELTAT) / 24.0
        print 'time[tidx], DELTAT: ', time[tidx], DELTAT

        date[tidx] = (ctime.year * 10000) + (ctime.month * 100) + ctime.day
        datesec[tidx] = (ctime.hour * 60 * 60) + (ctime.minute * 60) + \
                        ctime.second
        print 'date[tidx], datesec[tidx]: ', date[tidx], datesec[tidx]

        _do_one_timestep(ctime, ddi, hyam, hybm, p0, phis, phisfcm_in,
                         phisfcr_in, plat, plevm, plevr, plevrp1, plon, ps,
                         pslat, pslon, psus, psvs, q, rgfa2m, rgfa2mus,
                         rgfa2mvs, rgfm2us, rgfm2vs, rgfus2m, rgfvs2m, t,
                         tidx, u, us, v, vs)

        fout.sync()

    fout.close()
    hfin.close()
    tfin.close()

    return


if __name__ == '__main__':
#-------------------------
    args = sys.argv[1:]

    if len(args) < 9:
        print __doc__
        print '\nFirst 9 arguments required - EXITING!\n'
        raise SystemExit(1)

    else:
        year = int(args[0])
        month = int(args[1])
        day = int(args[2])
        hour = int(args[3])
        ntimes = int(args[4])
        state_xmlfile_dir = args[5]
        state_xmlfile_name = args[6]
        hybrid_file = args[7]
        topo_file = args[8]

        args = args[9:]

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
                print __doc__
                print 'Invalid option: %s\n' % arg
                raise SystemExit(1)

        print year, month, day, hour
        print ntimes
        print state_xmlfile_dir
        print state_xmlfile_name
        print hybrid_file
        print topo_file
        print file_tag

        _interim2cesm_fv(year, month, day, hour, ntimes, state_xmlfile_dir,
                         state_xmlfile_name, hybrid_file, topo_file, file_tag)
