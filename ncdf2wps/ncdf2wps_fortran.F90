MODULE routines
IMPLICIT NONE
CONTAINS
!   This subroutine in fortran 90 is taken from http://www2.mmm.ucar.edu/wrf/OnLineTutorial/Basics/IM_files/sample.f90

SUBROUTINE writeWPS (OFILE, SLAB, NX, NY, HDATE, XFCST, STARTLOC, FIELD, UNITS, DESC, MAP_SOURCE, XLVL, IPROJ, STARTLAT, &
                     STARTLON, DELTALAT, DELTALON, DX, DY, XLONC, TRUELAT1, TRUELAT2, NLATS)
	IMPLICIT NONE
	!
	! HDATE Date of the record, written as

    character(len=*), intent(in)   :: OFILE
	integer, intent(in)            :: NX
	integer, intent(in)            :: NY
	real, dimension(NX, NY), intent(in) :: SLAB
	character(len=24), intent(in)  :: HDATE
	real, intent(in)               :: XFCST
	character(len=8) , intent(in)  :: STARTLOC
	character(len=9) , intent(in)  :: FIELD
	character(len=25), intent(in)  :: UNITS
	character(len=46), intent(in)  :: DESC
	character(len=32), intent(in)  :: MAP_SOURCE
	real, intent(in)               :: XLVL
	integer, intent(in)            :: IPROJ
	real, intent(in)               :: STARTLAT
	real, intent(in)               :: STARTLON
	real, intent(in)               :: DELTALAT
	real, intent(in)               :: DELTALON
	real, intent(in)               :: DX
	real, intent(in)               :: DY
	real, intent(in)               :: XLONC
	real, intent(in)               :: TRUELAT1
	real, intent(in)               :: TRUELAT2
	real, intent(in)               :: NLATS

	
	integer, parameter             :: IUNIT = 10
	integer, parameter             :: OUNIT = 11
	integer                        :: ierr
	integer                        :: IFV=5
	real                           :: EARTH_RADIUS = 6367470. * .001
	logical                        :: IS_WIND_EARTH_REL = .FALSE.

	!====================================================================================!
	! 
	! You need to call the function from python and place each 2D slab here before       !
	! you can write it out to into the intermadiate file format                          !
	!                                                                                    !
	! Other information you need to know about your data:                                !
	!    Time at which data is valid                                                     !
	!    Forecast time of the data                                                       !
	!    Source of data - you can make something up, it is never used                    !
	!    Field name - NOTE THEY NEED TO MATCH THOSE EXPECTED BY METGRID                  !
	!    Units of field                                                                  !
	!    Description of data                                                             !
	!    Level of data - Pa, 200100 Pa is used for surface, and 201300 Pa is used        !
	!          for sea-level pressure                                                    !
	!    X dimension                                                                     !
	!    Y dimension                                                                     !
	!    Data projection - only recognize                                                !
	!         0:  Cylindrical Equidistant (Lat/lon) projection.                          !
	!         1:  Mercator projection.                                                   !
	!         3:  Lambert-conformal projection.                                          !
	!         4:  Gaussian projection.                                                   !
	!         5:  Polar-stereographic projection.                                        !
	!    Start location of data - "CENTER", "SWCORNER". "SWCORNER" is typical            !
	!    Start lat & long of data                                                        !
	!    Lat/Lon increment                                                               !
	!    Number of latitudes north of equator (for Gaussian grids)                       !
	!    Grid-spacing in x/y                                                             !
	!    Center long                                                                     !
	!    truelat1/2                                                                      !
	!    Has the winds been rotated                                                      !
	!====================================================================================!
    print  *, 'Writing file ', OFILE
    open(unit = IUNIT, file = OFILE, FORM = 'UNFORMATTED', position='REWIND')

	write (IUNIT, IOSTAT=IERR) 5
	
	! WRITE the second record, common to all projections:
	
	write (IUNIT) HDATE, XFCST, MAP_SOURCE, FIELD, UNITS, DESC, XLVL, NX, NY, IPROJ
	print*, HDATE//"  ", XLVL, FIELD
	
	! WRITE the third record, which depends on the projection:
	
	if (IPROJ == 0) then 
	
	   !  This is the Cylindrical Equidistant (lat/lon) projection:
	   WRITE (IUNIT) STARTLOC, STARTLAT, STARTLON, DELTALAT, DELTALON, EARTH_RADIUS
	
	elseif (IPROJ == 1) then 
	
	   ! This is the Mercator projection:
	   WRITE (IUNIT) STARTLOC, STARTLAT, STARTLON, DX, DY, TRUELAT1, EARTH_RADIUS
	
	elseif (IPROJ == 3) then
	
	   ! This is the Lambert Conformal projection:
	   WRITE (IUNIT) STARTLOC, STARTLAT, STARTLON, DX, DY, XLONC, TRUELAT1, TRUELAT2, EARTH_RADIUS
	   
	elseif (IPROJ == 4) then
	
	   ! Gaussian projection                         
	   WRITE (IUNIT) STARTLOC, STARTLAT, STARTLON, NLATS, DELTALON, EARTH_RADIUS
	    
	elseif (IPROJ == 5) then
	
	   ! This is the Polar Stereographic projection:
	   WRITE (IUNIT) STARTLOC, STARTLAT, STARTLON, DX, DY, XLONC, TRUELAT1, EARTH_RADIUS

	endif

 
	WRITE (IUNIT) IS_WIND_EARTH_REL


	WRITE (IUNIT) slab

	close(IUNIT)
	!  write(*,'(/,"End of read loop.  Program finished.")')

END SUBROUTINE writeWPS

end MODULE