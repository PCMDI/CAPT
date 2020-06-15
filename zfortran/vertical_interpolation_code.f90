c    Here are the routines that do various vertical interpolations.  Any
c  comments referring to input "analyses" assume the input fields have 
c  already been horizontally interpolated to the model's horizontal grid.
c
c   Note that the routine "psadj" (adjusts the newly interpolated surface
c  pressure) must be called *first* before any vertical interpolation is
c  done.
c
c
c    This is what each routine does:
c
c     psadj         : adjusts surface pressure for differences in surface height
c                     between model and analysis.
c
c     tsadj         : adjusts surface temperature for differences in surface
c                     height between model and analysis.
c
c     vert_quad_opt1: Quadratic vertical interpolation designed for Temperature
c                     interpolation.  We set the parameter "loglin" to 0 to
c                     interpolate in ln(P).
c
c     vert_int_opt1 : Linear interpolation designed for moisture fields like
c                     q, cloud water, cloud ice, cloud frac, etc.  We set the
c                     parameter "loglin" to 1 to interpolate in P.
c
c     vert_int_opt2 : Combination of Linear and Quadratic interpolation designed
c                     for vertical interpolation of U/V.  We set the
c                     parameter "loglin" to 0 to interpolate in ln(P).
c
c
c
c---------
c
c
C NCLFORTSTART
      subroutine tsadj(plat    ,plon    ,phis_old,phis_new,ts      )
C
C-----------------------------------------------------------------------
C
C Adjust Ts based on difference between old and new phis.
C
C-----------------------------------------------------------------------
C
      implicit none
C
C-----------------------------------------------------------------------
C
C     INPUTS
C
      integer plat !  latitude dimension
      integer plon !  longitude dimension
C
      real*8 phis_old(plon,plat) ! analysis phis (e.g., ECMWF)
      real*8 phis_new(plon,plat) ! model phis
C
C     INPUT/OUTPUT
C
      real*8 ts      (plon,plat) ! Surface Temp
C
C NCLEND
C
C---------------------------Local workspace-----------------------------
C
      real*8 dtdz
      real*8 gravit
      real*8 del_z
      integer i, j, k             ! Indices
C
      dtdz    = -0.0065           ! -6.5 deg/km
      gravit  = 9.80616           ! acceleration of gravity ~ m/s^2

      do j = 1,plat
        do i = 1,plon

          del_z = ( phis_new(i,j) - phis_old(i,j) )/gravit
          ts(i,j) = ts(i,j) + dtdz*del_z

        end do
      end do

      return
      end

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C NCLFORTSTART
      subroutine vert_quad_opt1(plevi , plevip1 , plev , plat  ,plon,
     $                 t_old ,pressi_m ,pressi_i ,presso_m, phis_old,
     $                 ps_old  , t_new , loglin )
C
C-----------------------------------------------------------------------
C
C Quadratic interpolation (designed for Temperature interpolation)
C
C                         (if "loglin" == 1: in P
C                          if "loglin" /= 1: in ln(P) )
C
C  Above input top           :  quadratic using top, levels 1 and 2 (top
C                               defined as 1.e-10 Pa for now).  Top
C                               value set to value at level 1
C  Between levels 1 and "bot":  quadratic interp using 3 closest levels
C  Between levels "bot"
C                 and surface:  linear interpolation using "Tbot"and
C                               "Tsurf".
C  Below surface             :  You don't wanna know (see doc)...
C
C
C-----------------------------------------------------------------------
C
      implicit none
C
C-----------------------------------------------------------------------
C
C     INPUTS
C
      integer plevi     ! vertical dimension of input analysis fields
      integer plevip1   ! "plevi+1 (vert dimension of input interfaces)
      integer plev      ! vertical dimension of model fields
      integer plat      ! latitude dimension
      integer plon      ! longitude dimension
C
      real*8 t_old   (plevi  ,plon,plat)   ! analysis tempertatures
      real*8 pressi_m (plevi  ,plon,plat)  ! analysis pressures at all levels
      real*8 pressi_i (plevip1,plon,plat)  ! analysis interface pressures
      real*8 presso_m (plev   ,plon,plat)  ! model pressures (based on adjusted PS)
      real*8 phis_old(plon,     plat)      ! analysis phis
      real*8 phis_new(plon,     plat)      ! model phis
      real*8 ps_old  (plon,     plat)      ! analysis surface pressure
      real*8 ps_new  (plon,     plat)      ! "adjusted" model surface pressure
      integer loglin                       ! interpolation flag
C
C     OUTPUT
C
      real*8 t_new   (plev   ,plon,plat)   ! Interpolated Temperatures
C
C
C NCLEND
C
C---------------------------Local workspace-----------------------------
C
      real*8 tsurf
      real*8 t0
      real*8 t_ref1
      real*8 t_ref2
      real*8 t_ref3
      real*8 t_ref3_top
      real*8 t_ref3_bot
      real*8 z_ref_top
      real*8 z_ref_bot
      real*8 tbot
      real*8 pbot
      real*8 psurf
      real*8 dtdz
      real*8 lapse
      real*8 boltz
      real*8 avogad
      real*8 mwdair
      real*8 rgas
      real*8 rdair
      real*8 gravit
      real*8 x
      real*8 threshold
      real*8 tmp
      real*8 z
      real*8 z_min
      real*8 z_incr
      real*8 hkk
      real*8 x1
      real*8 x2
      real*8 x3
      real*8 p1
      real*8 p2
      real*8 p3
      real*8 px
      real*8 pt                       ! top input pressure (ghost)
      real*8 xt                       ! top input value (linear extrap)
      real*8 tmp1
      real*8 tmp2
      real*8 tmp3
      real*8 tmpt
      real*8 beta
      integer i, j, k, kk, kkp1, kkp2 ! Indices
      integer k_bot
C
      dtdz    = -0.0065           ! -6.5 deg/km
      gravit  = 9.80616           ! acceleration of gravity ~ m/s^2
      boltz   = 1.38065e-23       ! boltzmann's constant ~ J/k/molecule
      avogad  = 6.02214e26        ! avogadro's number ~ molecules/kmole
      mwdair  = 28.966            ! molecular weight dry air ~ kg/kmole

      rgas    = avogad*boltz      ! universal gas constant ~ J/k/kmole
      rdair   = rgas/mwdair       ! constant for dry air   ~ J/k/kg

      t_ref1    = 290.5
      t_ref2    = 255.0
      t_ref3    = 298.0
      z_ref_bot = 2000.
      z_ref_top = 2500.

      threshold = 0.001
      pt = 1.e-10
      if(loglin .ne. 1) pt = log(pt)

      do j = 1,plat
        do i = 1,plon
C
C Tbot and Pbot are determined from the first model level that is at
C least 150m above the surface
C  ... will be used later
C
            z_min = 150.
            z     = 0.

            do k = plevi,1,-1
              k_bot  = k
              hkk    = 0.5*( pressi_i(k+1,i,j) - pressi_i(k,i,j) )/
     $                                                   pressi_m(k,i,j)
              z_incr = (rdair/gravit)*t_old(k,i,j)*hkk
              z      = z + z_incr
              if(z .gt. z_min) go to 10
              z      = z + z_incr
            end do

            write(6,*) 'Error:  could not find model level above ',z_min
            call abort

   10       continue
            lapse = -dtdz

            tbot  = t_old  (k_bot,i,j)
            pbot  = pressi_m(k_bot,i,j)
            tmp   = lapse*(rdair/gravit)*(ps_old(i,j)/pbot - 1.)
            tsurf = tbot*(1. + tmp)

          kk   = 1
          kkp1 = kk + 1
          kkp2 = kk + 2
C
C Find bracketting input levels if possible
C
          do k = 1,plev
   20       continue
            if(pressi_m(kkp1,i,j) .le. presso_m(k,i,j) ) then
              beta = presso_m(k,i,j) - pressi_m(kkp1,i,j)
              beta = beta/( pressi_m(kkp2,i,j) - pressi_m(kkp1,i,j) )
              if(beta .ge. 0.5) then
                kk   = kk + 1
                kkp1 = kk + 1
                kkp2 = kk + 2
                if(kkp2 .ge. plevi+1) then
                  kk   = plevi - 2
                  kkp1 = kk + 1
                  kkp2 = kk + 2
                  go to 30
                endif
                go to 20
              endif
            endif

   30       continue

            x1    = t_old   (kk   ,i,j)
            x2    = t_old   (kkp1 ,i,j)
            x3    = t_old   (kkp2 ,i,j)
            p1    = pressi_m(kk   ,i,j)
            p2    = pressi_m(kkp1 ,i,j)
            p3    = pressi_m(kkp2 ,i,j)
            pbot  = pressi_m(k_bot,i,j)
            psurf = ps_old  (      i,j)

            px    = presso_m(k    ,i,j)
C
C Convert to log(P) for log(P) interpolation
C
            if(loglin .ne. 1) then
              p1    = log(p1)
              p2    = log(p2)
              p3    = log(p3)
              pbot  = log(pbot)
              psurf = log(psurf)
              px    = log(px)
            endif
C
C If above 1st analysis level: quadratic interp
C
            if    (px .lt. p1 ) then
              xt         = x1

              tmpt       = ( (px-p1)*(px-p2) )/( (pt-p1)*(pt-p2) )
              tmp1       = ( (px-pt)*(px-p2) )/( (p1-pt)*(p1-p2) )
              tmp2       = ( (px-pt)*(px-p1) )/( (p2-pt)*(p2-p1) )

              t_new(k,i,j) = xt*tmpt + x1*tmp1 + x2*tmp2
C
C Elseif between "pbot" and analysis surface: linear interp
C
            elseif(px .ge. pbot .and. px .le. psurf ) then
              t_new(k,i,j) = ( tbot*(psurf - px) + tsurf*(px - pbot) )
     $                                                   /(psurf - pbot)
C
C Elseif below analysis surface: special case - see documentation
C
            elseif(px .gt. psurf ) then
              t_ref3_bot =      tsurf - z_ref_bot*dtdz
              t_ref3_top = min( tsurf - z_ref_top*dtdz, t_ref3 )
              z = phis_old(i,j)/gravit

              if(z .ge. z_ref_bot) then
                if(z .ge. z_ref_top) then
                  t0 = min( tsurf - z*dtdz, t_ref3 )
                else
                  t0 = ( t_ref3_bot*(z_ref_top - z        )   +
     $                   t_ref3_top*(z         - z_ref_bot) ) /
     $                              (z_ref_top - z_ref_bot)
                endif
                lapse = max ( (t0 - tsurf)/z , 0.)
              else
                lapse = -dtdz
              endif

              x   = lapse*rdair/gravit*log( presso_m(k    ,i,j)/
     $                                      ps_old  (      i,j) )
              t_new(k,i,j) = tsurf*(1. + x + x*x/2. + x*x*x/6.)
C
C Should never happen
C
            elseif(px .ge. p3 ) then
              write(6,*) 'Error:  extrapolation below input levels'
              write(6,*) 'not allowed'
              call abort
C
C Else between 1st analysis level and "pbot":  quadratic interp
C
            else
              tmp1       = ( (px-p2)*(px-p3) )/( (p1-p2)*(p1-p3) )
              tmp2       = ( (px-p1)*(px-p3) )/( (p2-p1)*(p2-p3) )
              tmp3       = ( (px-p1)*(px-p2) )/( (p3-p1)*(p3-p2) )

              t_new(k,i,j) = x1*tmp1 + x2*tmp2 + x3*tmp3

            endif
          end do
        end do
      end do

      return
      end

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C NCLFORTSTART
      subroutine vert_int_opt2(plat    ,plon    ,plevi   ,plevo   ,
     $                         pressi  ,presso  ,xxi     ,xxo     ,
     $                         loglin  )
C
C-----------------------------------------------------------------------
C
C Designed for vertical interpolation of U/V
C
C Linear and Quadratic interpolation (if "loglin" == 1: in P
C                                     if "loglin" /= 1: in ln(P) )
C
C  Above input top        :  quadratic using top, levels 1 and 2 (top
C                            defined as 1.e-10 Pa for now).  Top value
C                            determined from linear extrapolation from
C                            levels 1 and 2
C  Between levels 1 and 2 :  quadratic interp using levels 1,2, & 3
C  Between levels 2 and K :  linear interpolation using adjacent levels
C  Below level K          :  set equal to level K
C
C-----------------------------------------------------------------------
C
      implicit none
C
C-----------------------------------------------------------------------
C
C     INPUTS
C
      integer plat  ! latitude dimension
      integer plon  ! longitude dimension
      integer plevi ! vertical dimension of analysis fields
      integer plevo ! vertical dimension of model fields
C
      real*8 pressi  (plevi,plon,plat) ! analysis pressures
      real*8 presso  (plevo,plon,plat) ! model pressures (based on adjusted PS)
      real*8 xxi     (plevi,plon,plat) ! input analysis field
      integer loglin                   ! interpolation flag
C
C     OUTPUTS
C
      real*8 xxo     (plevo,plon,plat) ! model field
C
C NCLEND
C
C---------------------------Local workspace-----------------------------
C
      real*8 x1
      real*8 x2
      real*8 x3
      real*8 p1
      real*8 p2
      real*8 p3
      real*8 px
      real*8 pt                       ! top input pressure (ghost)
      real*8 xt                       ! top input value (linear extrap)
      real*8 tmp1
      real*8 tmp2
      real*8 tmp3
      real*8 tmpt
      integer i, j, k, kk, kkp1, kkp2 ! Indices
C
C-----------------------------------------------------------------------
C
      pt = 1.e-10
      if(loglin .ne. 1) pt = log(pt)
C
      do j = 1,plat
        do i = 1,plon

          kk   = 1
          kkp1 = kk + 1
          kkp2 = kk + 2
C
C Find bracketting analysis pressure levels
C
          do k = 1,plevo
   10       continue
            if(pressi(kkp1,i,j) .le. presso(k,i,j) ) then
              kk   = kk + 1
              kkp1 = kk + 1
              kkp2 = kk + 2
              if(kkp1 .eq. plevi+1) then
                kk   = plevi - 1
                kkp1 = kk + 1
                kkp2 = kk + 2
                go to 20
              endif
              go to 10
            endif

   20       continue

            x1 = xxi   (kk  ,i,j)
            x2 = xxi   (kkp1,i,j)
            p1 = pressi(kk  ,i,j)
            p2 = pressi(kkp1,i,j)
            px = presso(k   ,i,j)

            if(kkp2 .le. plevi) then
              x3 = xxi   (kkp2,i,j)
              p3 = pressi(kkp2,i,j)
            else
              x3 = 1.e+36
              p3 = 1.e+36
            endif
C
C Convert to log(P) for log(P) interpolation
C
            if(loglin .ne. 1) then
              p1 = log(p1)
              p2 = log(p2)
              p3 = log(p3)
              px = log(px)
            endif
C
C If above 1st analysis level:  quadratic interp
C
            if    (px .lt. p1 ) then
              xt         = ( x1*(p2 - pt) - x2*(p1 - pt) )/(p2 - p1)

              tmpt       = ( (px-p1)*(px-p2) )/( (pt-p1)*(pt-p2) )
              tmp1       = ( (px-pt)*(px-p2) )/( (p1-pt)*(p1-p2) )
              tmp2       = ( (px-pt)*(px-p1) )/( (p2-pt)*(p2-p1) )

              xxo(k,i,j) = xt*tmpt + x1*tmp1 + x2*tmp2
C
C Elseif below bottom analysis level:  output = bottome analysis field value
C
            elseif(px .ge. p2 ) then
              xxo(k,i,j) = x2
C
C Elseif between 1st and 2nd analysis levels:  quadratic interp
C
            elseif(kk .eq. 1 ) then
              tmp1       = ( (px-p2)*(px-p3) )/( (p1-p2)*(p1-p3) )
              tmp2       = ( (px-p1)*(px-p3) )/( (p2-p1)*(p2-p3) )
              tmp3       = ( (px-p1)*(px-p2) )/( (p3-p1)*(p3-p2) )

              xxo(k,i,j) = x1*tmp1 + x2*tmp2 + x3*tmp3
C
C Else, Linear interpolation
C
            else
              xxo(k,i,j) = ( x1*(p2 - px) + x2*(px - p1) )/(p2 - p1)

            endif
          end do
        end do
      end do

      return
      end

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C NCLFORTSTART
      subroutine vert_int_opt1(plat    ,plon    ,plevi   ,plevo   ,
     $                         pressi  ,presso  ,xxi     ,xxo     ,
     $                         loglin  )
C
C-----------------------------------------------------------------------
C
C Designed for moisture fields like q, cloud water, cloud ice, cloud frac, etc.
C
C Linearly interpolate (if "loglin" == 1: in P
C                       if "loglin" /= 1: in ln(P) )
C
C  Above input top        :  set equal to level 1
C  Between levels 1 and K :  linear interpolation using adjacent levels
C  Below level K          :  set equal to level K
C
C-----------------------------------------------------------------------
C
      implicit none
C
C-----------------------------------------------------------------------
C
C     INPUTS
C
      integer plat   ! latitude dimension
      integer plon   ! longitude dimension
      integer plevi  ! vertical dimension of analysis fields
      integer plevo  ! vertical dimension of model fields
C
      real*8 pressi  (plevi,plon,plat)  ! analysis pressures
      real*8 presso  (plevo,plon,plat)  ! model pressures (based on adjusted PS)
      real*8 xxi     (plevi,plon,plat)  ! analysis field
      integer loglin                    ! interpolation flag
C
C     OUTPUTS
C
      real*8 xxo     (plevo,plon,plat)  ! model field
C
C NCLEND
C
C---------------------------Local workspace-----------------------------
C
      real*8 p1
      real*8 p2
      real*8 px
      integer i, j, k, kk, kkp1         ! Indices
C
C-----------------------------------------------------------------------
C
      do j = 1,plat
        do i = 1,plon
C
C Find bracketting analysis pressure levels
C
          kk   = 1
          kkp1 = kk + 1

          do k = 1,plevo
   10       continue
            if(pressi(kkp1,i,j) .le. presso(k,i,j) ) then
              kk   = kk + 1
              kkp1 = kk + 1
              if(kkp1 .eq. plevi+1) then
                kk   = plevi - 1
                kkp1 = kk + 1
                go to 20
              endif
              go to 10
            endif

   20       continue
C
C If above 1st analysis level:  output = top analysis field value
C
            if    (presso(k,i,j) .lt. pressi(kk,i,j) ) then
              xxo(k,i,j) = xxi(kk  ,i,j)
C
C If below bottom analysis level:  output = bottom analysis field value
C
            elseif(presso(k,i,j) .ge. pressi(kkp1,i,j) ) then
              xxo(k,i,j) = xxi(kkp1,i,j)
C
C Else, Linear interpolation
C
            else
              p1 = pressi(kk  ,i,j)
              p2 = pressi(kkp1,i,j)
              px = presso(k   ,i,j)
              if(loglin .ne. 1) then
                p1 = log(p1)
                p2 = log(p2)
                px = log(px)
              endif
              xxo(k,i,j) = xxi(kk  ,i,j)*(p2 - px) + 
     $                     xxi(kkp1,i,j)*(px - p1)
              xxo(k,i,j) = xxo(k   ,i,j)/(p2 - p1)
            endif
          end do
        end do
      end do

      return
      end

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C NCLFORTSTART
      subroutine psadj(plev    ,plevp1  ,plat    ,plon    ,t       ,
     $                 press_m ,press_i ,phis_old,phis_new,ps_old  ,
     $                 ps_new  )
C
C-----------------------------------------------------------------------
C
C Adjust Ps based on difference between "analysis" phis and model phis.
C Also uses T and P arrays
C
C-----------------------------------------------------------------------
C
      implicit none
C
C-----------------------------------------------------------------------
C
C     INPUTS
C
      integer plev      ! vertical dimension
      integer plevp1    ! "plev+1"
      integer plat      ! latitude dimension
      integer plon      ! longitude dimension
C
      real*8 t       (plon,plat,plev)    ! analysis Temperatures
      real*8 press_m (plon,plat,plev)    ! analysis pressures
      real*8 press_i (plon,plat,plevp1)  ! analysis pressures (interfaces)
      real*8 phis_old(plon,     plat)    ! analysis phis
      real*8 phis_new(plon,     plat)    ! model phis
      real*8 ps_old  (plon,     plat)    ! analysis Ps (horizontally
C                                        ! interpolated to model grid)
C
C     OUTPUTS
C
      real*8 ps_new  (plon,     plat)    ! adjusted model Ps
C
C NCLEND
C
C---------------------------Local workspace-----------------------------
C
      real*8 tsurf
      real*8 t_ref1
      real*8 t_ref2
      real*8 t0
      real*8 tbot
      real*8 pbot
      real*8 dtdz
      real*8 lapse
      real*8 boltz
      real*8 avogad
      real*8 mwdair
      real*8 rgas
      real*8 rdair
      real*8 gravit
      real*8 del_phis
      real*8 x
      real*8 threshold
      real*8 tmp
      real*8 z
      real*8 z_min
      real*8 z_incr
      real*8 hkk
      integer i, j, k, kk         ! Indices
C
      dtdz    = -0.0065           ! -6.5 deg/km
      gravit  = 9.80616           ! acceleration of gravity ~ m/s^2
      boltz   = 1.38065e-23       ! boltzmann's constant ~ J/k/molecule
      avogad  = 6.02214e26        ! avogadro's number ~ molecules/kmole
      mwdair  = 28.966            ! molecular weight dry air ~ kg/kmole

      rgas    = avogad*boltz      ! universal gas constant ~ J/k/kmole
      rdair   = rgas/mwdair       ! constant for dry air   ~ J/k/kg

      t_ref1    = 290.5
      t_ref2    = 255.0
      threshold = 0.001

      do j = 1,plat
        do i = 1,plon

          del_phis = phis_old(i,j) - phis_new(i,j)
C
C If difference between analysis and model phis is negligible,
C then set model Ps = analysis
C
          if(abs(del_phis) .le. threshold) then
            ps_new(i,j) = ps_old(i,j)
C
C Else, go nuts...

          else
C
C Tbot and Pbot are determined from the first model level that is at
C least 150m above the surface
C
            z_min = 150.
            z     = 0.

            do k = plev,1,-1
              kk     = k
              hkk    = 0.5*( press_i(i,j,k+1) - press_i(i,j,k) )/
     $                                                    press_m(i,j,k)
              z_incr = (rdair/gravit)*t(i,j,k)*hkk
              z      = z + z_incr
              if(z .gt. z_min) go to 10
              z      = z + z_incr
            end do

            write(6,*) 'Error:  could not find model level above ',z_min
            call abort

   10       continue
            lapse = -dtdz
            k     = kk
C
C Define Tbot & Pbot
C
            tbot  = t      (i,j,k)
            pbot  = press_m(i,j,k)
            tmp   = lapse*(rdair/gravit)*(ps_old(i,j)/pbot - 1.)
            tsurf = tbot*(1. + tmp)
            t0    = tsurf + lapse*phis_old(i,j)/gravit
C
C See documentation for explanation of following code
C
            if     (t0 .gt. t_ref1 .and. tsurf .le. t_ref1) then
              lapse = (t_ref1 - tsurf)*gravit/phis_old(i,j)
            elseif (t0 .gt. t_ref1 .and. tsurf .gt. t_ref1) then
              lapse = 0.
              tsurf = (t_ref1 + tsurf)*0.5
            endif

            if(tsurf .lt. t_ref2) then
              lapse = -dtdz
              tsurf = (t_ref2 + tsurf)*0.5
            end if              

            x   = lapse*del_phis/(gravit*tsurf)
            tmp = 1. - x/2. + x**2./3.
            tmp = del_phis/(rdair*tsurf)*tmp
            ps_new(i,j) = ps_old(i,j)*exp(tmp)

          endif
        end do
      end do

      return
      end
