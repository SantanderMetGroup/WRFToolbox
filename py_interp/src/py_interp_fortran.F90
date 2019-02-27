MODULE routines
IMPLICIT NONE
CONTAINS
!  Subroutines extracted from p_interp.F90 Authors:
!
!  November 2007 - Cindy Bruyere
!  December 2009 - Lluis Fita (Santander Meteorolgy Group, Univ. Cantabria)
!  January 2012 - J. Fernandez (Santander Meteorolgy Group, Univ. Cantabria) -> jf
!
SUBROUTINE interp (data_out, data_in, pres_field, interp_levels, psfc, ter, tk, qv, ix, iy, iz, it, &
                     num_metgrid_levels, LINLOG, extrapolate, GEOPT, MISSING)
    IMPLICIT NONE
    ! Input fields:
    !
    ! ix, iy, iz, it     NetCDF dimensions
    ! num_metgrid_levels Number pf pressure levels.
    ! LINLOG             Interpolation option (linear or log)
    ! data_in            Data of the input variable
    ! pres_field, tk, qv Pressure, temperature and specific humidity fields.
    ! psfc               Surface pressure
    ! ter                Terrain height
    ! interp_levels      Pressure levels to interpolate
    ! GEOPT              Boolean: True if we are interpolating the geopotential
    ! extrapolate        Boolean: Extrapolate below ground
    !
    INTEGER, INTENT(IN)                              :: ix, iy, iz, it    
    INTEGER, INTENT(IN)                              :: num_metgrid_levels, LINLOG
    REAL(8), DIMENSION(ix, iy, iz, it), INTENT(IN)      :: data_in, pres_field, tk, qv 

    REAL(8), DIMENSION(ix, iy, it), INTENT(IN)          :: psfc
    REAL(8), DIMENSION(ix, iy), INTENT(IN)              :: ter
    REAL(8), DIMENSION(num_metgrid_levels), INTENT(IN)  :: interp_levels
    REAL(8)                                             :: ptarget, dp, dpmin, expon, pbot, zbot, tbotextrap, tvbotextrap
    LOGICAL, INTENT(IN)                                 :: GEOPT
    INTEGER, INTENT(IN)                                 :: extrapolate
    REAL(8), INTENT(IN)                                 :: MISSING
    
    REAL(8), DIMENSION(ix, iy, num_metgrid_levels, it), INTENT(OUT)  :: data_out
    INTEGER                                          :: i, j, itt, k, kk, kin, kupper
    REAL(8), DIMENSION(num_metgrid_levels)              :: data_out1D
    REAL(8), DIMENSION(iz)                              :: data_in1D, pres_field1D

!     REAL, DIMENSION(ix, iy, num_metgrid_levels, it)  :: N
!     REAL                                             :: sumA, sumN, AVE_geopt


!    N = 1.0

    expon=287.04*.0065/9.81


    do itt = 1, it
        do j = 1, iy
            do i = 1, ix
               data_in1D(:)    = data_in(i,j,:,itt)
               pres_field1D(:) = pres_field(i,j,:,itt)
               CALL int1D (data_out1D, data_in1D, pres_field1D, interp_levels, iz, num_metgrid_levels, LINLOG, MISSING)
               data_out(i,j,:,itt) = data_out1D(:)
            end do
        end do
    end do


    ! Fill in missing values
    IF ( extrapolate == 0 ) RETURN       !! no extrapolation - we are out of here

    ! First find where about 400 hPa is located
    kk = 0
    find_kk : do k = 1, num_metgrid_levels
        kk = k
        if ( interp_levels(k) <= 40000. ) exit find_kk
    end do find_kk


    IF ( GEOPT ) THEN     !! geopt is treated different below ground

    do itt = 1, it
       do k = 1, kk
          do j = 1, iy
              do i = 1, ix
                 IF ( data_out(i,j,k,itt) == MISSING .AND. interp_levels(k) < psfc(i,j,itt) ) THEN

!                We are below the first model level, but above the ground 

                    data_out(i,j,k,itt) = ((interp_levels(k) - pres_field(i,j,1,itt))*ter(i,j)*9.81 +  &
                                           (psfc(i,j,itt) - interp_levels(k))*data_in(i,j,1,itt) ) /   &
                                          (psfc(i,j,itt) - pres_field(i,j,1,itt))

                 ELSEIF ( data_out(i,j,k,itt) == MISSING ) THEN

!                We are below both the ground and the lowest data level.

!                First, find the model level that is closest to a "target" pressure
!                level, where the "target" pressure is delta-p less that the local
!                value of a horizontally smoothed surface pressure field.  We use
!                delta-p = 150 hPa here. A standard lapse rate temperature profile
!                passing through the temperature at this model level will be used
!                to define the temperature profile below ground.  This is similar
!                to the Benjamin and Miller (1990) method, except that for
!                simplicity, they used 700 hPa everywhere for the "target" pressure.
!                Code similar to what is implemented in RIP4

                    ptarget = (psfc(i,j,itt)*.01) - 150.
                    dpmin=1.e4
                    kupper = 0
                    loop_kIN : do kin=iz,1,-1
                       kupper = kin
                       dp=abs( (pres_field(i,j,kin,itt)*.01) - ptarget )
                       if (dp.gt.dpmin) exit loop_kIN
                       dpmin=min(dpmin,dp)
                    enddo loop_kIN

                    pbot=max(pres_field(i,j,1,itt),psfc(i,j,itt))
                    zbot=min(data_in(i,j,1,itt)/9.81,ter(i,j))

                    tbotextrap=tk(i,j,kupper,itt)*(pbot/pres_field(i,j,kupper,itt))**expon
                    tvbotextrap=virtual(tbotextrap,qv(i,j,1,itt))

                    data_out(i,j,k,itt) = (zbot+tvbotextrap/.0065*(1.-(interp_levels(k)/pbot)**expon))*9.81
               
                 ENDIF
              enddo
          enddo
       enddo
    enddo


    !!! Code for filling missing data with an average - we don't want to do this
    !!do itt = 1, it
       !!loop_levels : do k = 1, num_metgrid_levels
          !!sumA = SUM(data_out(:,:,k,itt), MASK = data_out(:,:,k,itt) /= MISSING)
          !!sumN = SUM(N(:,:,k,itt), MASK = data_out(:,:,k,itt) /= MISSING)
          !!IF ( sumN == 0. ) CYCLE loop_levels
          !!AVE_geopt = sumA/sumN
          !!WHERE ( data_out(:,:,k,itt) == MISSING )
             !!data_out(:,:,k,itt) = AVE_geopt
          !!END WHERE
       !!end do loop_levels
    !!end do

    END IF

    !!! All other fields and geopt at higher levels come here
    do itt = 1, it
    do j = 1, iy
    do i = 1, ix
      do k = 1, kk
         if ( data_out(i,j,k,itt) == MISSING ) data_out(i,j,k,itt) = data_in(i,j,1,itt)
      end do
      do k = kk+1, num_metgrid_levels
         if ( data_out(i,j,k,itt) == MISSING ) data_out(i,j,k,itt) = data_in(i,j,iz,itt)
      end do
    end do
    end do
    end do

END SUBROUTINE interp
 
SUBROUTINE int1D(xxout, xxin, ppin, ppout, npin, npout, LINLOG, MISSING)
    IMPLICIT NONE
!
! Modified from int2p - NCL code
! routine to interpolate from one set of pressure levels
! .   to another set  using linear or ln(p) interpolation
!
! NCL: xout = int2p (pin,xin,pout,linlog)
! This code was originally written for a specific purpose.
! .   Several features were added for incorporation into NCL's
! .   function suite including linear extrapolation.
!
! nomenclature:
!
! .   ppin   - input pressure levels. The pin can be
! .            be in ascending or descending order
! .   xxin   - data at corresponding input pressure levels
! .   npin   - number of input pressure levels >= 2
! .   ppout  - output pressure levels (input by user)
! .            same (ascending or descending) order as pin
! .   xxout  - data at corresponding output pressure levels
! .   npout  - number of output pressure levels
! .   linlog - if abs(linlog)=1 use linear interp in pressure
! .            if abs(linlog)=2 linear interp in ln(pressure)
! .   missing- missing data code. 
!                                                ! input types
    INTEGER, INTENT(IN)     :: npin,npout,linlog
    real(8), INTENT(IN)     :: ppin(npin),xxin(npin),ppout(npout)
    real(8), INTENT(IN)     :: MISSING
    !                                                ! output
    real(8),     INTENT(OUT)  :: xxout(npout)
    INTEGER                :: np,nl,nlmax
    INTEGER                :: nlsave,nlstrt
    real(8)                :: slope,pa,pb,pc

    ! automatic arrays
    real(8)      :: p(npin),x(npin)
    real(8)    :: pout(npout),xout(npout)


    xxout = MISSING
    pout  = ppout
    p     = ppin
    x     = xxin
    nlmax = npin

    ! exact p-level matches
    nlstrt = 1
    nlsave = 1
    do np = 1,npout
        xout(np) = MISSING
        do nl = nlstrt,nlmax
            if (pout(np).eq.p(nl)) then
                xout(np) = x(nl)
                nlsave = nl + 1
                go to 10
            end if
        end do
    10     nlstrt = nlsave
    end do

    if (LINLOG.eq.1) then
        do np = 1,npout
            do nl = 1,nlmax - 1
                if (pout(np).lt.p(nl) .and. pout(np).gt.p(nl+1)) then
                  slope = (x(nl)-x(nl+1))/ (p(nl)-p(nl+1))
                  xout(np) = x(nl+1) + slope* (pout(np)-p(nl+1))
                end if
            end do
        end do
    elseif (LINLOG.eq.2) then
        do np = 1,npout
            do nl = 1,nlmax - 1
            if (pout(np).lt.p(nl) .and. pout(np).gt.p(nl+1)) then
              pa = log(p(nl))
              pb = log(pout(np))
            ! special case: in case someone inadvertently enter p=0.
              if (p(nl+1).gt.0.d0) then
                  pc = log(p(nl+1))
              else
                  pc = log(1.d-4)
              end if

              slope = (x(nl)-x(nl+1))/ (pa-pc)
              xout(np) = x(nl+1) + slope* (pb-pc)
            end if
            end do
        end do
    end if

! place results in the return array;
    xxout = xout

END SUBROUTINE int1D
!------------------------------------------------------------------------------
SUBROUTINE compute_mslp(data_out, pres_field, psfc, ter, tk, qv, ix, iy, iz, it)
! New subroutine to compute mean_sealevel pressure values 

    INTEGER, INTENT(IN)                                           :: ix, iy, iz, it
    REAL(8), DIMENSION(ix, iy, it), INTENT(OUT)                   :: data_out
    REAL(8), DIMENSION(ix, iy, iz, it), INTENT(IN)                :: pres_field, tk, qv
    REAL(8), DIMENSION(ix, iy, it), INTENT(IN)                    :: psfc
    REAL(8), DIMENSION(ix, iy), INTENT(IN)                        :: ter
    INTEGER                                                       :: i, j, itt, kin, kupper
    REAL(8)                                                       :: ptarget, dpmin, dp, tkint, tbotextrap
    REAL(8)                                                       :: tvbotextrap, pbot, expon

!    N = 1.0
    expon = 287.04*.0065/9.81

!    data_out = 0.
!    We are below both the ground and the lowest data level.

!    First, find the model level that is closest to a "target" pressure
!    level, where the "target" pressure is delta-p less that the local
!    value of a horizontally smoothed surface pressure field.  We use
!    delta-p = 150 hPa here. A standard lapse rate temperature profile
!    passing through the temperature at this model level will be used
!    to define the temperature profile below ground.  This is similar
!    to the Benjamin and Miller (1990) method, using  
!    700 hPa everywhere for the "target" pressure.
    do itt = 1, it
        do j = 1, iy
            do i = 1, ix
                ptarget = (psfc(i,j,itt)*.01) - 150.
                dpmin = 1.e4
                kupper = 0
                loop_kIN : do kin = iz,1,-1
                kupper = kin
                dp=abs( (pres_field(i, j, kin, itt)*.01) - ptarget )
                if (dp.gt.dpmin) exit loop_kIN
                dpmin=min(dpmin, dp)
            enddo loop_kIN
            ptarget=ptarget*100.
            !          
            if (pres_field(i, j, kupper + 1, itt) - ptarget .ge. 0) then
                kupper = kupper + 1
            endif

!           kupper = 8
!           ptarget = pres_field(i, j, kupper, itt)
!
!           García-Díez 2012-06
!           The reference level temperature is computed by
!           linear interpolation, so there is no jump when the selected eta level changes.
!
            tkint = (tk(i, j, kupper, itt)*abs(ptarget - pres_field(i, j, kupper, itt)) + &
            tk(i, j, kupper + 1, itt)*abs(ptarget - pres_field(i, j, kupper + & 
            1, itt)))/abs(pres_field(i, j, kupper, itt) - pres_field(i, j, kupper + 1, itt))
!           tkint = tk(i, j, kupper, itt)

            pbot = pres_field(i,j,1,itt)
            
            !         tbotextrap = tkint*(psfc(i, j, itt)/ptarget)**expon
            tkint = tk(i, j, kupper, itt)
            tbotextrap = tkint*(pbot/ptarget)**expon
            tvbotextrap = virtual(tbotextrap, qv(i,j,1,itt))
            data_out(i, j, itt) = pbot*((tvbotextrap + 0.0065*ter(i, j))/tvbotextrap)**(1/expon)

!         IF (i==INT(ix/2) .AND. j==INT(iy/2) ) THEN
!         IF (i==54 .AND. j==134) THEN
!           IF (ter(i,j) > 2300.) THEN
!             PRINT *, 'pkupper - 1:', pres_field(i,j,kupper - 1,itt), 'pkupper:', pres_field(i,j,kupper,itt), 'pkupper + 1:', pres_field(i,j,kupper + 1,itt), 'pkupper + 2:', pres_field(i,j,kupper + 2,itt)
!             PRINT *, 'tk kupper:', tk(i,j,kupper,itt), 'tk kupper + 1:', tk(i,j,kupper + 1,itt), 'tkint:', tkint
!             PRINT *, 'qv kupper:', qv(i,j,kupper,itt), 'qv kupper + 1:', qv(i,j,kupper + 1,itt), 'qvint:', qvint
!             PRINT *,itt,' ptarget',ptarget,'kupper:',kupper
!             PRINT *,'tk:',tk(i,j,kupper,itt),'psfc:',psfc(i,j,itt)
!             PRINT *,'tbot:',tbotextrap,'tvbot:',tvbotextrap,'ter:',ter(i,j)
!             PRINT *,'qv:',qv(i,j,kupper,itt),'mslp:',data_out(i,j,itt,1)
!           ENDIF
            enddo ! i
        enddo ! j
    enddo ! itt

END SUBROUTINE compute_mslp

!---------------------------------------------------------------------
SUBROUTINE calcslptwo(slp, PP, P_s, PHI_s, T_L, nz, ns, ew, nt)
! Josipa's subroutine
IMPLICIT NONE

INTEGER,                INTENT(IN)  :: nz, ns, ew, nt !(input) dimensions: vertical, north-south, east-west
real, DIMENSION(ew, ns, nt),   INTENT(OUT) :: slp        !(output) sea level pressure
real, DIMENSION(ew, ns, nz, nt), INTENT(IN)  :: PP         !(input) 3D pressure
real, DIMENSION(ew, ns, nt),   INTENT(IN)  :: P_s        !(input) pressure at surface
real, DIMENSION(ew, ns),   INTENT(IN)  :: PHI_s      !(input) 2D geopotential of the surface
real, DIMENSION(ew, ns, nt),   INTENT(IN)  :: T_L        !(input) temperature at lowest level


real, DIMENSION(ew, ns, nt)         :: P_L        !(calculated) pressure at lowest level
real                                :: T_surf     !(calculated) surface temperature
real                                :: gamma_mod  !(calculated) modified lapse rate that is actually used for calculations
real                                :: x          !(calculated) expansion coefficient of (eq. 9)
real                                :: T_0        !(calculated) auxiliary variable

real, PARAMETER                     :: Rd    = 287.04 ![J kg-1 K-1] (const.) dry air constant
real, PARAMETER                     :: g     = 9.81   ![m s-2] (const.) acceleration due to gravity
real, PARAMETER                     :: gamma = 0.0065 ![K m-1] (const.) lapse rate at const. 0.0065 K/m, also denoted as (dT/dz)_st
integer                             :: j,k,t !loop parameters

!!PRINT *,'this is calcslptwo !!!'

P_L = PP(:,:,1,:)  !extract lowest layer

DO t = 1, nt
  DO j = 1 , ew
    DO k = 1 , ns

      !always assume none of the IF conditions trigger, then (according to step 5)
      gamma_mod = gamma

      !(0) if abs(PHI_s)<0.001 ("sea" grid cells) then set slp to surface pressure
      IF (ABS(PHI_s(j,k)).lt.0.001) THEN

        slp(j,k,t) = P_s(j,k,t) !in this case we are done with this grid cell

      ELSE !else apply the following algorithm

        !always assume none of the following IF conditions trigger
        !then we use the constant camma (according to step 5)
        gamma_mod = gamma

        !(1) compute T_surf according to eq (1)
        !T_surf = T_L + gamma*(Rd/g)*(P_s/P_L-1.0)*T_L
        T_surf = T_L(j,k,t) + gamma * (Rd/g) * ( P_s(j,k,t)/P_L(j,k,t) - 1.0 ) * T_L(j,k,t)

        !(2) compute T_0=T_surf+gamma*PHI_s/g
        T_0 = T_surf + gamma * PHI_s(j,k) / g

        !(3) to avoid extrapolation of too low pressures over high and warm surfaces:
        ! if (T_0 > 290.5):
        !    if (T_surf <= 290.5): gamma_mod=(290.5-T_surf)*g/PHI_s (eq. 7)
        !    else: gamma_mod=0.0, T_surf=0.5*(290.5+T_surf)
        IF (T_0 .gt. 290.5) THEN
          IF (T_surf .le. 290.5) THEN
            gamma_mod = ( 290.5 - T_surf ) * g / PHI_s(j,k)
          ELSE
            gamma_mod = 0.0
            T_surf = 0.5 * ( 290.5 + T_surf )
          END IF
        END IF

        !(4) to avoid extrapolation of too high pressures over cold surfaces:
        ! if T_surf < 255: gamma_mod=gamma, T_surf=0.5*(255+T_surf)
        IF (T_surf .lt. 255.0) THEN
          gamma_mod = gamma
          T_surf = 0.5 * ( 255.0 + T_surf )
        END IF

        !(5) in other cases set gamma_mod=gamma
        !this was already done in the beginning of the loop!

        !(6) compute mean sea level pressure (eq. 8) using the above determined parameters
        !x=gamma_mod*PHI_s/(g*T_surf)
        !slp=P_s*exp(PHI_s/(R_d*T_surf)*(1-x/2.+x*x/3.))
        x = gamma_mod * PHI_s(j,k) / ( g * T_surf )
        slp(j,k,t) = P_s(j,k,t) * EXP( PHI_s(j,k) / (Rd * T_surf ) * (1.0 - x/2. + x*x/3.) )

        !done.

      END IF

    END DO
  END DO
END DO
END SUBROUTINE calcslptwo

SUBROUTINE clt_sundqvist(dx, dy, dz, dt, cldfra, totcloudfr)
!  Subroutine to compute total cloud cover in base 1. BY LLUIS

    IMPLICIT NONE
    INTEGER,INTENT(IN)                                     :: dx, dy, dz, dt
    REAL(8), DIMENSION(dx,dy,dz,dt), INTENT(IN)           :: cldfra
    REAL(8), DIMENSION(dx,dy,dt),    INTENT(OUT)          :: totcloudfr
    ! Local
    INTEGER                                                 :: k
    REAL(8), DIMENSION(dx,dy,dz,dt)                        :: cldfram1, cldmax

    !!!!!!!!!!!!!! Variables
    ! dx, dy, dz, dt: dimensions of fields
    ! cldfra: cloud fraction at each level
    ! totcfr: total cloud fraction
    !
    ! rcode = nf_inq_varid(ncid, cldfraname, idcldfra)
    ! rcode = nf_get_var_real(ncid, idcldfra, cldfra)

    totcloudfr = 1.
    cldfram1(:,:,2:dz,:) = cldfra(:,:,1:dz-1,:)
    cldfram1(:,:,1,:) = 0.
    cldmax = MAX(cldfram1,cldfra)

    WHERE (cldfram1 == 1.) cldfram1 = 0.99

    vertical_levels: DO k=1, dz 
        totcloudfr=totcloudfr*((1. - cldmax(:,:,k,:))/(1. - cldfram1(:,:,k,:)))
!        PRINT *, "totcloudfr(dx/, dy/2, 2):", totcloudfr(50, 50, 2)
!        PRINT *, "cldfram1(dx/, dy/2, 2):", cldfram1(50,50,k,2)
!        PRINT *, "cldmax(dx/, dy/2, 2):", cldmax(50,50,k,2)
    END DO vertical_levels
    totcloudfr = 1. - totcloudfr

END SUBROUTINE clt_sundqvist

!---------------------------------------------------------------------

SUBROUTINE clt_maxrand(tot_cldfra, cldfra, ix, iy, iz, it)
!
!  Subroutine to compute total cloud cover in base 1 using maximum-random overlapping
!
    IMPLICIT NONE
    INTEGER,INTENT(IN)                                     :: ix, iy, iz, it
    REAL(8), DIMENSION(ix,iy,iz,it), INTENT(IN)           :: cldfra
    REAL(8), DIMENSION(ix,iy,it),    INTENT(OUT)          :: tot_cldfra
    ! Local
    REAL(8)                                                 :: tot_col_cldfra
    INTEGER                                                 :: i, j, t
    !
    ! ix, iy, iz, it: dimensions of fields
    ! cldfra: cloud fraction at each model level
    ! tot_cldfra: total cloud fraction (ix, iy, it) array
    ! tot_col_cldfra: column total cloud fraction real
    ! minpres, maxpres, optional bounds for pressure levels in hPa
    ! pres_field, optional pressure field
    !
    do t = 1, it
        do j = 1, iy
            do i = 1, ix
                call clt_maxrand_column(tot_col_cldfra, cldfra(i,j,:,t), iz)
                tot_cldfra(i,j,t) = tot_col_cldfra
            end do ! i
        end do ! j
    end do ! t
    
END SUBROUTINE clt_maxrand

!---------------------------------------------------------------------

SUBROUTINE clt_maxrand_levels(tot_cldfra, cldfra, pres_field, maxpres, minpres, ix, iy, iz, it)
!
! Subroutine to compute total cloud cover in base 1 using maximum-random overlapping.
! This one is able to compute high, low and medium clouds. Its sepparated from the total
! because optional arguments are not working and then it is not efficient.
!
    IMPLICIT NONE
    INTEGER,INTENT(IN)                                     :: ix, iy, iz, it
    REAL(8), DIMENSION(ix,iy,iz,it), INTENT(IN)           :: cldfra
    REAL(8), DIMENSION(ix,iy,iz,it), INTENT(IN)           :: pres_field
    REAL(8), INTENT(IN)                                    :: minpres, maxpres
    REAL(8), DIMENSION(ix,iy,it),    INTENT(OUT)          :: tot_cldfra
    ! Local
    REAL(8), DIMENSION(ix,iy,iz,it)                        :: masked_cldfra
    REAL(8)                                                 :: tot_col_cldfra
    INTEGER                                                 :: i, j, t
    !
    ! ix, iy, iz, it: dimensions of fields
    ! cldfra: cloud fraction at each model level
    ! tot_cldfra: total cloud fraction (ix, iy, it) array
    ! tot_col_cldfra: column total cloud fraction real
    ! minpres, maxpres, bounds for pressure levels in hPa
    ! pres_field, pressure field
    !
    masked_cldfra = cldfra
    where (pres_field <= minpres*100.) masked_cldfra(:,:,:,:) = 0.
    where (pres_field >= maxpres*100.) masked_cldfra(:,:,:,:) = 0.

    do t = 1, it
        do j = 1, iy
            do i = 1, ix
                call clt_maxrand_column(tot_col_cldfra, masked_cldfra(i,j,:,t), iz)
                tot_cldfra(i,j,t) = tot_col_cldfra
            end do ! i
        end do ! j
    end do ! t
    
END SUBROUTINE clt_maxrand_levels

!---------------------------------------------------------------------

SUBROUTINE clt_maxrand_column(tot_col_cldfra, col_cldfra, iz)
!
! Subroutine for performing maximum-random overlapping in one column.
!
    INTEGER,INTENT(IN)                                     :: iz
    REAL(8), DIMENSION(iz), INTENT(IN)                    :: col_cldfra
    REAL(8), INTENT(OUT)                                   :: tot_col_cldfra
    INTEGER                                                 :: k, nseq
    REAL(8),  ALLOCATABLE, DIMENSION(:)                  :: seq_cldfra
    LOGICAL                                                 :: in_cloudy_section
    !
    ! iz (in): Number of vertical levels.
    ! col_cldfra (in): Cloud fraction in model levels.
    ! tot_col_cldfra (out) : Column total cloud fraction real
    ! k: Vertical level counter
    ! nseq: Cloudy section number (levels with cldfra > 1 sepparated by cloud free levels)
    ! seq_cldfra: Allocatable array to store each sections cloud fraction
    ! in_cloudy_section: Flag that tells us if we are inside a cloudy section or not.
    !
    ! First, we need to find the number of cloudy sections.
    !
    nseq = 0
    in_cloudy_section = .FALSE.
    do k = 1, iz
        if (col_cldfra(k) > 0.) then
            !
            ! Check if we are already in a cloudy section
            !
            if (in_cloudy_section) then
                !
                ! We are inside a section, don't count it
                !
                cycle
            else
                !
                ! We found a new cloudy section
                !
                in_cloudy_section = .TRUE.
                nseq = nseq + 1
            endif
        else
            !
            ! We are outside a cloudy section.
            !
            in_cloudy_section = .FALSE.
        endif
    end do ! k
    !debug
    !print *, "Number of cloudy sections found:", nseq
    !
    ! Allocate the vector with cloud fractions of each section
    !
    allocate(seq_cldfra(nseq))
    !
    ! Then loop again over the vertical axis storing the cloud fraction
    ! of each section, computed assuming maximum overlapping.
    !
    nseq = 0
    in_cloudy_section = .FALSE.
    do k = 1, iz
        !
        ! Check if we are already in a cloudy section
        !
        if (col_cldfra(k) > 0.) then
            if (in_cloudy_section) then
                !
                ! We are inside a section, check if this level cltfra is
                ! larger than previous and store it.
                !
                seq_cldfra(nseq) = max(seq_cldfra(nseq), col_cldfra(k))
            else
                !
                ! We found a new cloudy section. Move the counter and 
                ! save its cldfra.
                !
                in_cloudy_section = .TRUE.
                nseq = nseq + 1
                seq_cldfra(nseq) = col_cldfra(k)
            endif
        else
            !
            ! We are outside a cloudy section.
            !
            in_cloudy_section = .FALSE.
        endif
    end do ! k

    !print *, "Cloud section cldfra:", seq_cldfra
    !print *, "All levels cloud fraction:" , col_cldfra
    !
    ! Compute the total cloud cover assuming random overlapping
    ! between cloudy sections.
    !
    tot_col_cldfra = 1. - product(1. - seq_cldfra)
END SUBROUTINE clt_maxrand_column
!------------------------------------------------------------------------------
SUBROUTINE pbl_height(nx, ny, nz, nt, ght, pt, ter, pblh)
  !
  !   Subroutine to compute the PBL height using the 1.5-theta-increase
  !   method by Nielsen-Gammon et al. 2008, JAMC 47:27-43
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN)                    :: nx, ny, nz, nt
  REAL,DIMENSION(nx,ny,nz,nt),INTENT(IN) :: ght ! geopotential
  REAL,DIMENSION(nx,ny,nz,nt),INTENT(IN) :: pt ! potential temperature
  REAL, DIMENSION(nx,ny), INTENT(IN)     :: ter ! terrain height
  REAL, DIMENSION(nx,ny,nt), INTENT(OUT) :: pblh
  INTEGER                                :: ix,iy,iz,it, rcode
  INTEGER, DIMENSION(nx,ny,nt)           :: minlevel  ! vertical level where minimal pot. temp. is reached
  REAL, DIMENSION(nx,ny,nt)              :: threshold ! minimum pot. temp. in the column plus a delta of 1.5 K


  minlevel = minloc(pt, 3)
  threshold = minval(pt, 3) + 1.5
  DO it=1,nt
    DO ix=1,nx
      DO iy=1,ny
        DO iz=minlevel(ix,iy,it),nz
          if (pt(ix,iy,iz,it).ge.threshold(ix,iy,it)) EXIT ! iz contains the level
        END DO
        ! Interpolate between levels iz-1 and iz
        pblh(ix,iy,it) = ght(ix,iy,iz,it)+ &
                            (threshold(ix,iy,it)-pt(ix,iy,iz,it)) * &
                            (ght(ix,iy,iz,it) - ght(ix,iy,iz-1,it)) / &
                            (pt(ix,iy,iz,it) - pt(ix,iy,iz-1,it))
        ! Convert to meters and remove terrain height
        pblh(ix,iy,it) = pblh(ix,iy,it)/9.81 - ter(ix,iy)
      END DO
    END DO
  END DO
  ! jf: Uncomment the following to dump into the PBLHGT variable the level
  !     where the minimum theta occurs (only for debug).
  !pblh = minlevel
END SUBROUTINE pbl_height

SUBROUTINE cape_3d(nx, ny, nz, nt, pres_field, tk, qv, cape)
  !
  !   Subroutine to compute the Convective Available Potential Energy
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN)                    :: nx, ny, nz, nt
  REAL,DIMENSION(nx,ny,nz,nt),INTENT(IN) :: pres_field ! pressure (Pa)
  REAL,DIMENSION(nx,ny,nz,nt),INTENT(IN) :: tk ! temperature (K)
  REAL,DIMENSION(nx,ny,nz,nt),INTENT(IN) :: qv ! specific humidity (kg kg-1)
  INTEGER                                :: ix,iy,iz,it, num_threads
  REAL                                   :: cape_column, cin, zlcl, zlfc 
  REAL                                   :: zel, psource, tsource, qvsource
  REAL,DIMENSION(nx,ny,nt), INTENT(OUT)  :: cape ! CAPE (J kg-1)
  
  DO it=1,nt
    DO iy=1,ny
      DO ix=1,nx
        call getcape(1, nz , pres_field(ix, iy, :, it), tk(ix, iy, :, it), &
               qv(ix, iy, :, it), cape_column , cin,  zlcl, zlfc, zel, psource, &
               tsource , qvsource )
        cape(ix, iy, it) = cape_column
      END DO
    END DO
  END DO
END SUBROUTINE cape_3d

subroutine getcape( source, nk , p_in , t_in , q_in, cape , cin,   &
                        zlcl, zlfc, zel , psource , tsource , qvsource )
    implicit none

    integer, intent(in) :: source,nk
    real, dimension(nk), intent(in) :: p_in,t_in,q_in
    real, intent(out) :: cape,cin,psource,tsource,qvsource

!-----------------------------------------------------------------------
!
!  getcape - a fortran90 subroutine to calculate Convective Available
!            Potential Energy (CAPE) from a sounding.
!
!   *** Modified version to calculate Lifted Condensation Level (LCL)  ***
!   *** Level of Free Convection (LFC), and Equilibrium Level (EL)     ***
!
!  Version 1.04                           Last modified:  8 October 2010
!
!  Author:  George H. Bryan
!           Mesoscale and Microscale Meteorology Division
!           National Center for Atmospheric Research
!           Boulder, Colorado, USA
!           gbryan@ucar.edu
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!
!  References:  Bolton (1980, MWR, p. 1046) (constants and definitions)
!               Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
!
!-----------------------------------------------------------------------
!
!  Input:     nk - number of levels in the sounding (integer)
!
!           p_in - one-dimensional array of pressure (mb) (real)
!
!           t_in - one-dimensional array of temperature (C) (real)
!
!           q_in - one-dimensional array of water vapor mixing ratio (kg/kg) (real)
!
!  Output:  cape - Convective Available Potential Energy (J/kg) (real)
!
!            cin - Convective Inhibition (J/kg) (real)
!
!-----------------------------------------------------------------------
!  User options:

    real, parameter :: pinc = 1000.0  ! Pressure increment (Pa) 
                                      ! Original was 100 but 1000 looks enough for WRF
                                      ! and it is way faster.
                                      ! (smaller number yields more accurate
                                      !  results,larger number makes code 
                                      !  go faster)

!!!    integer, parameter :: source = 2    ! Source parcel:
!!!                                        ! 1 = surface
!!!                                        ! 2 = most unstable (max theta-e)
!!!                                        ! 3 = mixed-layer (specify ml_depth)

    real, parameter :: ml_depth =  500.0  ! depth (m) of mixed layer 
                                          ! for source=3

    integer, parameter :: adiabat = 1   ! Formulation of moist adiabat:
                                        ! 1 = pseudoadiabatic, liquid only
                                        ! 2 = reversible, liquid only
                                        ! 3 = pseudoadiabatic, with ice
                                        ! 4 = reversible, with ice

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!            No need to modify anything below here:
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    logical :: doit,ice,cloud,not_converged
    integer :: k,kmax,n,nloop,i,orec
    real, dimension(nk) :: p,t,td,pi,q,th,thv,z,pt,pb,pc,pn,ptv,ptd,pqv,pql

    real :: the,maxthe,parea,narea,lfc
    real :: th1,p1,t1,qv1,ql1,qi1,b1,pi1,thv1,qt,dp,dz,ps,frac
    real :: th2,p2,t2,qv2,ql2,qi2,b2,pi2,thv2
    real :: thlast,fliq,fice,tbar,qvbar,qlbar,qibar,lhv,lhs,lhf,rm,cpm
    real*8 :: avgth,avgqv
!   real :: getqvl,getqvi,getthx,gettd
    real :: ee,zlcl,zlfc,zel,plcl,plfc,pel,zh

!-----------------------------------------------------------------------

    real, parameter :: g     = 9.81
    real, parameter :: p00   = 100000.0
    real, parameter :: cp    = 1005.7
    real, parameter :: rd    = 287.04
    real, parameter :: rv    = 461.5
    real, parameter :: xlv   = 2501000.0
    real, parameter :: xls   = 2836017.0
    real, parameter :: t0    = 273.15
    real, parameter :: cpv   = 1875.0
    real, parameter :: cpl   = 4190.0
    real, parameter :: cpi   = 2118.636
    real, parameter :: lv1   = xlv+(cpl-cpv)*t0
    real, parameter :: lv2   = cpl-cpv
    real, parameter :: ls1   = xls+(cpi-cpv)*t0
    real, parameter :: ls2   = cpi-cpv

    real, parameter :: rp00  = 1.0/p00
    real, parameter :: eps   = rd/rv
    real, parameter :: reps  = rv/rd
    real, parameter :: rddcp = rd/cp
    real, parameter :: cpdrd = cp/rd
    real, parameter :: cpdg  = cp/g

    real, parameter :: converge = 0.0002

    integer, parameter :: debug_level =   0

!-----------------------------------------------------------------------

!---- convert p,t,td to mks units; get pi,q,th,thv ----!

    do k=1,nk
        p(k) = p_in(k)
        t(k) = t_in(k)
        q(k) = q_in(k)
!!!       td(k) = 273.15+td_in(k)
       ee = alog((q(k)/eps)*p(k)/100.0/(1.0+(q(k)/eps)))
       td(k) = 273.15+(243.5*ee-440.8)/(19.48-ee)
       pi(k) = (p(k)*rp00)**rddcp
!!!        q(k) = getqvl(p(k),td(k))
       th(k) = t(k)/pi(k)
      thv(k) = th(k)*(1.0+reps*q(k))/(1.0+q(k))
!!!      print *,k,t(k)-273.15,td(k)-273.15
    enddo

!---- get height using the hydrostatic equation ----!

    z(1) = 0.0
    do k=2,nk
      dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))
      z(k) = z(k-1) + dz
    enddo

!----------------------------------------------------------------

!---- find source parcel ----!

  IF(source.eq.1)THEN
    ! use surface parcel
    kmax = 1

  ELSEIF(source.eq.2)THEN
    ! use most unstable parcel (max theta-e)

    IF(p(1).lt.50000.0)THEN
      ! first report is above 500 mb ... just use the first level reported
      kmax = 1
      maxthe = getthx(p(1),t(1),td(1),q(1))
    ELSE
      ! find max thetae below 500 mb
      maxthe = 0.0
      do k=1,nk
        if(p(k).ge.50000.0)then
          the = getthx(p(k),t(k),td(k),q(k))
          if( the.gt.maxthe )then
            maxthe = the
            kmax = k
          endif
        endif
      enddo
    ENDIF
    if(debug_level.ge.100) print *,'  kmax,maxthe = ',kmax,maxthe

  ELSEIF(source.eq.3)THEN
    ! use mixed layer

    IF( (z(2)-z(1)).gt.ml_depth )THEN
      ! the second level is above the mixed-layer depth:  just use the
      ! lowest level

      avgth = th(1)
      avgqv = q(1)
      kmax = 1

    ELSEIF( z(nk).lt.ml_depth )THEN
      ! the top-most level is within the mixed layer:  just use the
      ! upper-most level

      avgth = th(nk)
      avgqv = q(nk)
      kmax = nk

    ELSE
      ! calculate the mixed-layer properties:

      avgth = 0.0
      avgqv = 0.0
      k = 2
      if(debug_level.ge.100) print *,'  ml_depth = ',ml_depth
      if(debug_level.ge.100) print *,'  k,z,th,q:'
      if(debug_level.ge.100) print *,1,z(1),th(1),q(1)

      do while( (z(k).le.ml_depth) .and. (k.le.nk) )

        if(debug_level.ge.100) print *,k,z(k),th(k),q(k)

        avgth = avgth + 0.5*(z(k)-z(k-1))*(th(k)+th(k-1))
        avgqv = avgqv + 0.5*(z(k)-z(k-1))*(q(k)+q(k-1))

        k = k + 1

      enddo

      th2 = th(k-1)+(th(k)-th(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))
      qv2 =  q(k-1)+( q(k)- q(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))

      if(debug_level.ge.100) print *,999,ml_depth,th2,qv2

      avgth = avgth + 0.5*(ml_depth-z(k-1))*(th2+th(k-1))
      avgqv = avgqv + 0.5*(ml_depth-z(k-1))*(qv2+q(k-1))

      if(debug_level.ge.100) print *,k,z(k),th(k),q(k)

      avgth = avgth/ml_depth
      avgqv = avgqv/ml_depth

      kmax = 1

    ENDIF

    if(debug_level.ge.100) print *,avgth,avgqv

  ELSE

    print *
    print *,'  Unknown value for source'
    print *
    print *,'  source = ',source
    print *
    stop

  ENDIF

!---- define parcel properties at initial location ----!
    narea = 0.0

  if( (source.eq.1).or.(source.eq.2) )then
    k    = kmax
    th2  = th(kmax)
    pi2  = pi(kmax)
    p2   = p(kmax)
    t2   = t(kmax)
    thv2 = thv(kmax)
    qv2  = q(kmax)
    b2   = 0.0
  elseif( source.eq.3 )then
    k    = kmax
    th2  = avgth
    qv2  = avgqv
    thv2 = th2*(1.0+reps*qv2)/(1.0+qv2)
    pi2  = pi(kmax)
    p2   = p(kmax)
    t2   = th2*pi2
    b2   = g*( thv2-thv(kmax) )/thv(kmax)
  endif

    psource = p2
    tsource = t2
   qvsource = qv2

    ql2 = 0.0
    qi2 = 0.0
    qt  = qv2

    cape = 0.0
    cin  = 0.0
    lfc  = 0.0

    doit = .true.
    cloud = .false.
    if(adiabat.eq.1.or.adiabat.eq.2)then
      ice = .false.
    else
      ice = .true.
    endif

      the = getthx(p2,t2,t2,qv2)
      if(debug_level.ge.100) print *,'  the = ',the

      pt(k) = t2
      if( cloud )then
        ptd(k) = t2
      else
        ptd(k) = gettd(p2,t2,qv2)
      endif
      ptv(k) = t2*(1.0+reps*qv2)/(1.0+qv2)
      pb(k) = 0.0
      pqv(k) = qv2
      pql(k) = 0.0

      zlcl = -1.0
      zlfc = -1.0
      zel  = -1.0

!---- begin ascent of parcel ----!

      if(debug_level.ge.100)then
        print *,'  Start loop:'
        print *,'  p2,th2,qv2 = ',p2,th2,qv2
      endif

    do while( doit .and. (k.lt.nk) )

        k = k+1
       b1 =  b2

       dp = p(k-1)-p(k)

      if( dp.lt.pinc )then
        nloop = 1
      else
        nloop = 1 + int( dp/pinc )
        dp = dp/float(nloop)
      endif

      do n=1,nloop

         p1 =  p2
         t1 =  t2
        pi1 = pi2
        th1 = th2
        qv1 = qv2
        ql1 = ql2
        qi1 = qi2
        thv1 = thv2

        p2 = p2 - dp
        pi2 = (p2*rp00)**rddcp

        thlast = th1
        i = 0
        not_converged = .true.

        do while( not_converged )
          i = i + 1
          t2 = thlast*pi2
          if(ice)then
            fliq = max(min((t2-233.15)/(273.15-233.15),1.0),0.0)
            fice = 1.0-fliq
          else
            fliq = 1.0
            fice = 0.0
          endif
          qv2 = min( qt , fliq*getqvl(p2,t2) + fice*getqvi(p2,t2) )
          qi2 = max( fice*(qt-qv2) , 0.0 )
          ql2 = max( qt-qv2-qi2 , 0.0 )

          tbar  = 0.5*(t1+t2)
          qvbar = 0.5*(qv1+qv2)
          qlbar = 0.5*(ql1+ql2)
          qibar = 0.5*(qi1+qi2)

          lhv = lv1-lv2*tbar
          lhs = ls1-ls2*tbar
          lhf = lhs-lhv

          rm=rd+rv*qvbar
          cpm=cp+cpv*qvbar+cpl*qlbar+cpi*qibar
          th2=th1*exp(  lhv*(ql2-ql1)/(cpm*tbar)     &
                       +lhs*(qi2-qi1)/(cpm*tbar)     &
                       +(rm/cpm-rd/cp)*alog(p2/p1) )

          if(i.gt.90) print *,i,th2,thlast,th2-thlast
          if(i.gt.100)then
            print *
            print *,'  Error:  lack of convergence'
            print *
            print *,'  ... stopping iteration '
            print *
            stop 1001
          endif
          if( abs(th2-thlast).gt.converge )then
            thlast=thlast+0.3*(th2-thlast)
          else
            not_converged = .false.
          endif
        enddo

        ! Latest pressure increment is complete.  Calculate some
        ! important stuff:

        if( ql2.ge.1.0e-10 ) cloud = .true.
        if( cloud .and. zlcl.lt.0.0 )then
           zlcl = z(k-1)+(z(k)-z(k-1))*float(n)/float(nloop)
           plcl = p(k-1)+(p(k)-p(k-1))*float(n)/float(nloop)
        endif

        IF(adiabat.eq.1.or.adiabat.eq.3)THEN
          ! pseudoadiabat
          qt  = qv2
          ql2 = 0.0
          qi2 = 0.0
        ELSEIF(adiabat.le.0.or.adiabat.ge.5)THEN
          print *
          print *,'  Undefined adiabat'
          print *
          stop 10000
        ENDIF

      enddo

      thv2 = th2*(1.0+reps*qv2)/(1.0+qv2+ql2+qi2)
        b2 = g*( thv2-thv(k) )/thv(k)
        dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))

      if( zlcl.gt.0.0 .and. zlfc.lt.0.0 .and. b2.gt.0.0 )then
        if( b1.gt.0.0 )then
          zlfc = zlcl
          plfc = plcl
        else
          zlfc = z(k-1)+(z(k)-z(k-1))*(0.0-b1)/(b2-b1)
          plfc = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        endif
      endif

      if( zlfc.gt.0.0 .and. zel.lt.0.0 .and. b2.lt.0.0 )then
        zel = z(k-1)+(z(k)-z(k-1))*(0.0-b1)/(b2-b1)
        pel = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
      endif

      the = getthx(p2,t2,t2,qv2)

      pt(k) = t2
      if( cloud )then
        ptd(k) = t2
      else
        ptd(k) = gettd(p2,t2,qv2)
      endif
      ptv(k) = t2*(1.0+reps*qv2)/(1.0+qv2)
      pb(k) = b2
      pqv(k) = qv2
      pql(k) = ql2

      ! Get contributions to CAPE and CIN:

      if( (b2.ge.0.0) .and. (b1.lt.0.0) )then
        ! first trip into positive area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b2/(b2-b1)
        parea =  0.5*b2*dz*frac
        narea = narea-0.5*b1*dz*(1.0-frac)
        if(debug_level.ge.200)then
          print *,'      b1,b2 = ',b1,b2
          print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
          print *,'      frac = ',frac
          print *,'      parea = ',parea
          print *,'      narea = ',narea
        endif
        cin  = cin  + narea
        narea = 0.0
      elseif( (b2.lt.0.0) .and. (b1.gt.0.0) )then
        ! first trip into neg area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b1/(b1-b2)
        parea =  0.5*b1*dz*frac
        narea = -0.5*b2*dz*(1.0-frac)
        if(debug_level.ge.200)then
          print *,'      b1,b2 = ',b1,b2
          print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
          print *,'      frac = ',frac
          print *,'      parea = ',parea
          print *,'      narea = ',narea
        endif
      elseif( b2.lt.0.0 )then
        ! still collecting negative buoyancy
        parea =  0.0
        narea = narea-0.5*dz*(b1+b2)
      else
        ! still collecting positive buoyancy
        parea =  0.5*dz*(b1+b2)
        narea =  0.0
      endif

      cape = cape + max(0.0,parea)
      pc(k) = cape

      if(debug_level.ge.200)then
        write(6,102) p2,b1,b2,cape,cin,cloud
102     format(5(f13.4),2x,l1)
      endif

      if( (p(k).le.10000.0).and.(b2.lt.0.0) )then
        ! stop if b < 0 and p < 100 mb
        doit = .false.
      endif

    enddo

!!!    print *,'  zlcl,zlfc,zel = ',zlcl,zlfc,zel
!!!    print *,'  plcl,plfc,pel = ',plcl,plfc,pel

!---- All done ----!

    return
    end subroutine getcape

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

subroutine get_tlapser(ni, nj, nt, t2, topo, lmask, halo, tlapse)
  implicit none
  INTEGER, INTENT(IN)                    :: ni, nj, nt
  REAL,DIMENSION(ni,nj,nt),INTENT(IN)    :: t2 ! 2m temperature (C)
  REAL,DIMENSION(ni,nj),INTENT(IN)       :: topo ! topography
  LOGICAL, DIMENSION(ni, nj), INTENT(IN) :: lmask ! land sea mask
  INTEGER                                :: i, j, t, halo, imax, imin, jmax, jmin, tile_width, tile_size
  REAL ( kind = 8 )                                  :: a, b ! regression coefficients
  REAL( kind = 8 ),DIMENSION(:, :), ALLOCATABLE      :: t2_tile, topo_tile
  LOGICAL, DIMENSION(:, :), ALLOCATABLE              :: lmask_tile
 ! REAL( kind = 8 ),DIMENSION(:), ALLOCATABLE         :: t2_tilef, topo_tilef
  REAL( kind = 8 ),DIMENSION(ni,nj,nt),INTENT(OUT)   :: tlapse

  tile_width = 2*halo + 1
  tile_size = tile_width*tile_width
  ALLOCATE(t2_tile(tile_width, tile_width))
  ALLOCATE(topo_tile(tile_width, tile_width))
  ALLOCATE(lmask_tile(tile_width, tile_width))
  !ALLOCATE(t2_tilef(tile_size))
 ! ALLOCATE(topo_tilef(tile_size))

  DO t=1,nt
    DO j=1,nj
      DO i=1,ni
        !WRITE(*, *) t, j, i
        imin = max(i - halo, 1)
        jmin = max(j - halo, 1)
        imax = min(i + halo, ni)
        jmax = min(j + halo, nj)
        t2_tile = t2(imin:imax, jmin:jmax, t)
        topo_tile = topo(imin:imax, jmin:jmax)
        lmask_tile = lmask(imin:imax, jmin:jmax)

        !t2_tilef = PACK(t2_tile, mask=lmask_tile)
        !topo_tilef = PACK(topo_tile, mask=lmask_tile)
        CALL llsq(COUNT(lmask_tile), PACK(topo_tile, mask=lmask_tile), PACK(t2_tile, mask=lmask_tile), a, b )
        tlapse(i, j, t) = a

      END DO
    END DO
  END DO

  DEALLOCATE(t2_tile)
  DEALLOCATE(topo_tile)
  DEALLOCATE(lmask_tile)

end subroutine get_tlapser

subroutine llsq ( n, x, y, a, b )

!*****************************************************************************80
!
!! LLSQ solves a linear least squares problem matching a line to data.
!
!  Discussion:
!
!    A formula for a line of the form Y = A * X + B is sought, which
!    will minimize the root-mean-square error to N data points ( X(I), Y(I) );
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the data points.
!
!    Output, real ( kind = 8 ) A, B, the slope and Y-intercept of the
!    least-squares approximant to the data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bot
  real ( kind = 8 ) top
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xbar
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ybar
!
!  Special case.
!
  if ( n == 1 ) then
    a = 0.0D+00
    b = y(1)
    return
  end if
!
!  Average X and Y.
!
  xbar = sum ( x(1:n) ) / real ( n, kind = 8 )
  ybar = sum ( y(1:n) ) / real ( n, kind = 8 )
!
!  Compute Beta.
!
  top = dot_product ( x(1:n) - xbar, y(1:n) - ybar )
  bot = dot_product ( x(1:n) - xbar, x(1:n) - xbar )

  a = top / bot

  b = ybar - a * xbar

  return
end subroutine llsq

real function getqvl(p,t)
    implicit none
    
    real, intent(in) :: p,t
    real             :: es
    
    real, parameter :: eps = 287.04/461.5
    
    es = 611.2*exp(17.67*(t-273.15)/(t-29.65))
    getqvl = eps*es/(p-es)
    
    return
end function getqvl

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

real function getqvi(p,t)
    implicit none
    
    real :: p,t,es
    
    real, parameter :: eps = 287.04/461.5
    
    es = 611.2*exp(21.8745584*(t-273.15)/(t-7.66))
    getqvi = eps*es/(p-es)
    
    return
end function getqvi

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

real function getthx(p,t,td,q)
    implicit none
    
    real :: p,t,td,q
    real :: tlcl
    
    if( (td-t).ge.-0.1 )then
      tlcl = t
    else
      tlcl = 56.0 + ( (td-56.0)**(-1) + 0.00125*alog(t/td) )**(-1)
    endif
    
    getthx=t*( (100000.0/p)**(0.2854*(1.0-0.28*q)) )   &
            *exp( ((3376.0/tlcl)-2.54)*q*(1.0+0.81*q) )
    
    return
end function getthx

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

real function gettd(p,t,q)
    implicit none
    
    real :: p,t,q
    
    real :: el
    real, parameter :: eps = 287.04/461.5
    
    el = alog((q/eps)*p/100.0/(1.0+(q/eps)))
    gettd = 273.15+(243.5*el-440.8)/(19.48-el)
    
    return
end function gettd
  
!
!------------------------------------------------------------------------------
FUNCTION virtual (tmp,rmix)
    IMPLICIT NONE
!      This function returns virtual temperature in K, given temperature
!      in K and mixing ratio in kg/kg.

    real(8), intent(IN)                  :: tmp, rmix
    real(8)                              :: virtual

    virtual=tmp*(0.622+rmix)/(0.622*(1.+rmix))

END FUNCTION virtual

end MODULE
