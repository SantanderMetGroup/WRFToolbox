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
FUNCTION virtual (tmp,rmix)
    IMPLICIT NONE
!      This function returns virtual temperature in K, given temperature
!      in K and mixing ratio in kg/kg.

    real(8), intent(IN)                  :: tmp, rmix
    real(8)                              :: virtual

    virtual=tmp*(0.622+rmix)/(0.622*(1.+rmix))

END FUNCTION virtual

end MODULE
