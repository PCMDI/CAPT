#! /usr/bin/env python

#==============================================================================
#   ERA5toN480nc.py
"""
Python script.

This script is to convert ERA5 (grib2 format) from spectral grid to 
N480 full Gaussian grid (netcdf format) 

Interactive usage:
    ERA5toN480nc.py <year_start month_start day_start hr_start 
                     year_end month_end day_end hr_end
                     data_dir>
                    [optional arguments]

Required variables: 
    Input: Vorticity, Divergence, Temperature, Specific Humidity, 
           Surface Presurre (lnps), Surface Geopotential
    Output: Zonal Wind, Meridiornal Wind, Temperature, Specific Humidity, 
           Surface Presurre (lnps), Surface Geopotential

Note that this script requires both cdo and CDAT loaded.

Required arguments:
    year_start   # 
    month_start  #
    hr_start     #
    year_end     #
    month_end    #
    hr_end       #
    data_dir     #

Optional arguments:
    cdo_path     # This must be set if cdo is not in your executable
                   search path [''].
"""
#==============================================================================
import datetime
import calendar
import os
import subprocess
import sys
import cdms2,MV2,genutil,cdutil,numpy,cdtime,MV2

CDO_PATH_DEF = ''

def _ERA5_2_N480(year_start, month_start, day_start, hr_start, year_end, month_end, day_end, hr_end,
                  data_dir, cdo_path):
#--------------------------------------------------------------

  dt     = datetime.datetime(year_start, month_start, day_start, hr_start)
  dt_end = datetime.datetime(year_end, month_end, day_end, hr_end)
  delt_1hr = datetime.timedelta(hours=1)
  print(dt)
  print(dt_end)
  print(delt_1hr)
  print(' ')

  while dt <= dt_end :
    ymd  = '%04d%02d%02d' % (dt.year, dt.month, dt.day)
    filedate  = '%04d-%02d-%02d' % (dt.year, dt.month, dt.day)
    hour_str = '%02d' % dt.hour

    in_filename = 'UVTQPS-' + ymd + '_' + hour_str + '.grib'

    print(in_filename)
    print(filedate)
    print(ymd)

    dirname = data_dir + '/' 
    os.chdir(dirname)

    # cdo sinfon *.grib    #check grib variable information
    # select spectral fields of sfc z(129), ln sfc p(152) 
    #--------------------------------
    sel2d_call = cdo_path + 'cdo selvar,z,lnsp' + ' ' + \
                 in_filename + ' ' + 'ml_2d.grib'
    print(sel2d_call)
    sel2d = subprocess.Popen(sel2d_call, shell=True)
    sel2d.wait()

    # select spectral fields of T(130)
    #--------------------------------
    #selSp_call = cdo_path + 'cdo selcode,129,152,130' + ' ' + \
    selt_call = cdo_path + 'cdo selvar,t' + ' ' + \
                in_filename + ' ' + 'ml_t.grib'
    print(selt_call)
    selt = subprocess.Popen(selt_call, shell=True)
    selt.wait()

    # Convert spectral to full gg
    #--------------------------------
    sp2gp_call = cdo_path + 'cdo sp2gp' + ' ' + \
                 'ml_2d.grib' + ' ' + 'ml_2d_fgg.grib'
    print(sp2gp_call)
    spp = subprocess.Popen(sp2gp_call, shell=True)
    spp.wait()

    sp2gp_call = cdo_path + 'cdo sp2gp' + ' ' + \
                 'ml_t.grib' + ' ' + 'ml_t_fgg.grib'
    print(sp2gp_call)
    spp = subprocess.Popen(sp2gp_call, shell=True)
    spp.wait()

    # select div(155), vor(138)  in spectral
    #--------------------------------
    divVor_call = cdo_path + 'cdo selvar,d,vo' + ' ' + \
                  in_filename + ' ' + 'ml_divvo.grib'
    print(divVor_call)
    divv = subprocess.Popen(divVor_call, shell=True)
    divv.wait()

    # convert div, vor spectral to u,v full gaussian
    #--------------------------------
    dv2uv_call = cdo_path + 'cdo dv2uv' + ' ' + \
                 'ml_divvo.grib' + ' ' + 'ml_uv_fgg.grib'
    print(dv2uv_call)
    dvuv = subprocess.Popen(dv2uv_call, shell=True)
    dvuv.wait()

    # select reduced gg fields Q (133)
    #--------------------------------
    selRgg_call = cdo_path + 'cdo selvar,q' + ' ' + \
                  in_filename + ' ' + 'ml_q_rgg.grib'
    print(selRgg_call)
    slgg = subprocess.Popen(selRgg_call, shell=True)
    slgg.wait()

    # convert reduced gg to full gg
    #----------------------------
    #fgg_call = cdo_path + 'cdo -R copy' + ' ' + \
    fgg_call = cdo_path + 'cdo setgridtype,regular' + ' ' + \
               'ml_q_rgg.grib' + ' ' + 'ml_q_fggN320.grib'
    print(fgg_call)
    fgg = subprocess.Popen(fgg_call, shell=True)
    fgg.wait()

    # Remap Q onto N480 grid (Q is on a different grid).
    #------------------------------------------------
    remap_call = cdo_path + 'cdo remapbil,n480' + ' ' + \
                 'ml_q_fggN320.grib' + ' ' + 'ml_q_fgg.grib'
    print(remap_call)
    rmp = subprocess.Popen(remap_call, shell=True)
    rmp.wait()

    # Combine all the files
    #------------------------------------------------
    cat2_call = 'cat ml_uv_fgg.grib ml_t_fgg.grib ml_2d_fgg.grib > ' + \
                'UVTQPS-' + ymd + '_' + hour_str + '_N480_UVTn2D.grib'
    print(cat2_call)
    ct2 = subprocess.Popen(cat2_call, shell=True)
    ct2.wait()

    cp1_call = 'cp ml_q_fgg.grib' + ' ' + \
               'UVTQPS-' + ymd + '_'  + hour_str + '_N480_Q.grib'
    print(cp1_call)
    cp1 = subprocess.Popen(cp1_call, shell=True)
    cp1.wait()

    # Convert grib2 to grib1 and make grads ctl for hourly files.
    #----------------------------------
    # Use cdo gradsdes to geneate ctls files
    # This only works for grib1 format. Have to convert grib2 to grib1
    # In the processes, cdo -f grb copy would mess up the variable names so have to separate the variables into two files
    grib2to1_call = cdo_path + 'cdo -f grb copy' + ' ' + \
            'UVTQPS-' + ymd + '_' + hour_str + '_N480_UVTn2D.grib' + ' ' + \
            'UVTQPS-' + ymd + '_' + hour_str + '_N480_UVTn2D.grib1'
    print(grib2to1_call)
    grb21 = subprocess.Popen(grib2to1_call, shell=True)
    grb21.wait()

    grib2to1_call = cdo_path + 'cdo -f grb copy' + ' ' + \
            'UVTQPS-' + ymd + '_'  + hour_str + '_N480_Q.grib' + ' ' + \
            'UVTQPS-' + ymd + '_'  + hour_str + '_N480_Q.grib1'
    print(grib2to1_call)
    grb21 = subprocess.Popen(grib2to1_call, shell=True)
    grb21.wait()

    gradsdes1_call = cdo_path + 'cdo gradsdes' + ' ' + \
                     'UVTQPS-' + ymd + '_' + hour_str + '_N480_UVTn2D.grib1'
    print(gradsdes1_call)
    grd1 = subprocess.Popen(gradsdes1_call, shell=True)
    grd1.wait()

    gradsdes2_call = cdo_path + 'cdo gradsdes' + ' ' + \
                     'UVTQPS-' + ymd + '_'  + hour_str + '_N480_Q.grib1'
    print(gradsdes2_call)
    grd2 = subprocess.Popen(gradsdes2_call, shell=True)
    grd2.wait()

    # Generate NC file
    #----------------------------------

    infile1='UVTQPS-' + ymd + '_' + hour_str + '_N480_UVTn2D.grib1.ctl'
    infile2='UVTQPS-' + ymd + '_'  + hour_str + '_N480_Q.grib1.ctl'
    out_filename = 'UVTQPS-' + ymd + '_' + hour_str + '.nc'

    print(infile1)
    print(infile2)
    print(out_filename)

    fid1 = cdms2.open(infile1)
    fid2 = cdms2.open(infile2)

    var_u = fid1('var131')
    var_v = fid1('var132')
    var_t = fid1('var0')
    var_lnsp = fid1('var25')
    var_z = fid1('var4')
    var_q = fid2('var0')

    lev=var_u.getLevel()
    lat=var_u.getLatitude()
    lon=var_u.getLongitude()

    cdms2.setNetcdfShuffleFlag(0)
    cdms2.setNetcdfDeflateFlag(0)
    cdms2.setNetcdfDeflateLevelFlag(0)

    c = cdms2.open(out_filename,'w')

    time=cdms2.createAxis([0.])
    time.id = 'time'
    time.units = 'hours since ' + filedate + ' ' + hour_str + ':00:00'
    time.long_name = 'time'
    time.axis = 'T'
    time.calendar = 'gregorian'

    print(time.units)

    var1=cdms2.createVariable(var_u,axes=(time,lev,lat,lon),typecode='f',id='U')
    var1.long_name='zonal wind component'
    var1.units='m/s'

    var2=cdms2.createVariable(var_v,axes=(time,lev,lat,lon),typecode='f',id='V')
    var2.long_name='meridional wind component'
    var2.units='m/s'

    var3=cdms2.createVariable(var_t,axes=(time,lev,lat,lon),typecode='f',id='T')
    var3.long_name='temperature'
    var3.units='K'

    var4=cdms2.createVariable(var_q,axes=(time,lev,lat,lon),typecode='f',id='Q')
    var4.long_name='specific humidity'
    var4.units='kg/kg'

    var5=cdms2.createVariable(var_lnsp,axes=(time,lat,lon),typecode='f',id='PS')
    var5.long_name='surface pressure'
    var5.units='Pa'

    var6=cdms2.createVariable(var_z,axes=(time,lat,lon),typecode='f',id='PHIS')
    var6.long_name='surface geopotential'
    var6.units='m2/s2'

    c.write(var1)
    c.write(var2)
    c.write(var3)
    c.write(var4)
    c.write(var5)
    c.write(var6)

    fid1.close()
    fid2.close()
    c.close()

    ######### remove all temporary files and continue onto next day
    remove_call1 = 'rm -f' + ' ' + 'ml_2d.grib' + ' ' + 'ml_t.grib' + ' ' + \
                      'ml_2d_fgg.grib' + ' ' + 'ml_t_fgg.grib' + ' ' + \
                      'ml_divvo.grib' + ' ' + 'ml_uv_fgg.grib' + ' ' + \
                      'ml_q_rgg.grib' + ' ' + 'ml_q_fggN320.grib' + ' ' + 'ml_q_fgg.grib'
    print(remove_call1)
    rmc1 = subprocess.Popen(remove_call1, shell=True)
    rmc1.wait()

    remove_call2 = 'rm -f' + ' ' + '*_N480_*grib*'
    print(remove_call2)
    rmc2 = subprocess.Popen(remove_call2, shell=True)
    rmc2.wait()
    ##################################
    dt = dt + delt_1hr

  return

if __name__ == '__main__':
#-------------------------
    cdo_path = CDO_PATH_DEF

    year_start  = int(sys.argv[1])
    month_start = int(sys.argv[2])
    day_start   = int(sys.argv[3])
    hr_start    = int(sys.argv[4])
    year_end    = int(sys.argv[5])
    month_end   = int(sys.argv[6])
    day_end     = int(sys.argv[7])
    hr_end      = int(sys.argv[8])
    data_dir    = sys.argv[9]

    if len(sys.argv) > 10:
        cdo_path = sys.argv[10]

        if not cdo_path.endswith('/'):
            cdo_path = cdo_path + '/'

    _ERA5_2_N480(year_start, month_start, day_start, hr_start, year_end, month_end, day_end, hr_end,
                 data_dir, cdo_path)
