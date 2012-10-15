 PROGRAM no_interp
!  Program to read wrfout data and no interpolate and providing Netalevs eta
!    levels
!  The program reads namelist.nointerp
!  Adaptation from 
!  November 2007 - Cindy Bruyere
!    by Lluis Fita. Universidad de Cantabria
!    March 2011
!
!=================================Make Executable============================
!  Make executable:
!
!  GMS Intel 11.1 with static netCDF4 and HDF5 libraries
!
!  use intel
!  export NETCDF_DIR=
!  export HDF5_DIR=
!  ifort -O3  p_interp.F90 -o p_interp -I$NETCDF_DIR/include -I$HDF5_DIR/include \
!  -Bstatic -L$NETCDF_DIR/lib -lnetcdf -lnetcdff -L$HDF5_DIR/lib -lnetcdf -lhdf5_hl -lhdf5 -lz
!
!   GMS.trueno:
!      gfortran no_interp.F90 -L//home/lluis/bin/netcdf-4.0.1/lib -lnetcdf \
!      -I/home/lluis/bin/netcdf-4.0.1/include -o no_interp
!   GMS.UC.new:
!      gfortran no_interp.F90 -L/software/CentOS/5.2/netcdf/4.1.1/gcc-gfortran4.1.2/lib -lnetcdf \
!      -I/software/CentOS/5.2/netcdf/4.1.1/gcc-gfortran4.1.2/include -o no_interp
!
!   GMS.UC.old:
!      pgf90 no_interp.F90 -L/software/ScientificLinux/4.6/netcdf/3.6.3/pgf716_gcc/lib -lnetcdf -lm \
!      -I/software/ScientificLinux/4.6/netcdf/3.6.3/pgf716_gcc/include -Mfree -o no_interp
!
!   DEC Alpha
!      f90 no_interp.F90 -L/usr/local/netcdf/lib -lnetcdf -lm  \
!      -I/usr/local/netcdf/include  -free  -o no_interp
!
!   Linux flags
!      pgf90 no_interp.F90 -L/usr/local/netcdf/lib -lnetcdf -lm  \
!      -I/usr/local/netcdf/include  -Mfree  -o no_interp
!
!   Sun flags
!      f90 no_interp.F90 -L/usr/local/netcdf/lib -lnetcdf -lm  \
!      -I/usr/local/netcdf/include  -free  -o no_interp
!
!   SGI flags
!      f90 no_interp.F90 -L/usr/local/netcdf/lib -lnetcdf -lm  \
!      -I/usr/local/netcdf/include  -freeform  -o no_interp
!
!   IBM flags 
!      xlf no_interp.F90 -L/usr/local/netcdf/lib -lnetcdf -lm  \
!      -I/usr/local/netcdf/include  -qfree=f90  -o no_interp
!
!   Mac flags (with xlf compiler)
!      xlf no_interp.F90 -L/usr/local/netcdf-xlf/lib -lnetcdf -lm  \
!      -I/usr/local/netcdf-xlf/include  -qfree=f90  -o no_interp
!
!   Mac flags (with g95 compiler)
!	g95 no_interp.F90 -L/usr/local/netcdf/lib -lnetcdf -lm  \
!	-I/usr/local/netcdf/include -ffree-form -o no_interp
!
!   Mac flags (with pgf90 compiler)
!      pgf90 no_interp.F90 -L/usr/local/netcdf/lib -lnetcdf -lm  \
!      -I/usr/local/netcdf/include  -Mfree  -o no_interp
!
!============================================================================

      IMPLICIT NONE
      
      INCLUDE 'netcdf.inc'

      REAL, PARAMETER                                    :: Rd = 287.04
      REAL, PARAMETER                                    :: Cp = 7.*Rd/2.
      REAL, PARAMETER                                    :: RCP = Rd/Cp
      REAL, PARAMETER                                    :: p0 = 100000.

      CHARACTER,         ALLOCATABLE, DIMENSION(:,:,:,:) :: text
      CHARACTER (LEN=31),ALLOCATABLE, DIMENSION(:)       :: dnamei, dnamej
      CHARACTER(LEN=250),ALLOCATABLE, DIMENSION(:)       :: input_file_names
      CHARACTER(LEN=250),ALLOCATABLE, DIMENSION(:)       :: output_file_names
      DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:,:,:,:) :: ddata1, ddata2
      REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: data1, data2, data3
      REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: pres_field, pres_out
      REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: pres_stagU, pres_stagV
      REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: height_field
      REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: height_stagU, height_stagV
      REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: ght, phb, qv, tk, rh
      REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: psfc
      REAL,              ALLOCATABLE, DIMENSION(:,:)     :: ter, znu, znw
      INTEGER,           ALLOCATABLE, DIMENSION(:)       :: dvali, dvalj
      INTEGER,           ALLOCATABLE, DIMENSION(:,:,:,:) :: idata1, idata2
      INTEGER,                        DIMENSION(4)       :: start_dims = 1
      INTEGER,                        DIMENSION(4)       :: dims_in, dims_out
      INTEGER,                        DIMENSION(8)       :: ishape, jshape

      CHARACTER (LEN=80)                                 :: cval
      CHARACTER (LEN=31)                                 :: cname, test_dim_name
      CHARACTER (LEN=300)                                :: input_file, output_file, att_text
      CHARACTER (LEN=250)                                :: path_to_input
      CHARACTER (LEN=250)                                :: path_to_output
      CHARACTER (LEN=80)                                 :: input_name
      CHARACTER (LEN=250)                                :: output_name, tmp_name
      CHARACTER (LEN=10)                                 :: option
      CHARACTER (LEN=400)                                :: command
      CHARACTER (LEN=20)                                 :: process, dummy
      CHARACTER (LEN=2000)                               :: fields, process_these_fields
      INTEGER, DIMENSION(99)                             :: eta_levels
      REAL                                               :: rval
      REAL                                               :: MISSING=1.e36
      INTEGER                                            :: LINLOG = 1
      INTEGER                                            :: interp_method=1
      INTEGER                                            :: extrapolate=0
      INTEGER                                            :: Netalevs
      INTEGER                                            :: ncid, mcid, rcode
      INTEGER                                            :: idm, ndims, nvars, natt, ngatts
      INTEGER                                            :: nunlimdimid
      INTEGER                                            :: i, ii, j, jj, ix, iy, iweg, isng,   &
        ibtg, l, jt
      INTEGER                                            :: ivar, jvar, ijk
      INTEGER                                            :: times_in_file
      INTEGER                                            :: ilen, itype, ival, na, funit, ios
      INTEGER                                            :: num_eta_levels
      INTEGER                                            :: number_of_input_files
      INTEGER                                            :: ierr, loop, loslen, strlen, lent
      INTEGER                                            :: is_there
      LOGICAL                                            :: is_used
      LOGICAL                                            :: debug=.FALSE.
      LOGICAL                                            :: interpolate=.FALSE.
      LOGICAL                                            :: unstagger_grid=.FALSE.
      LOGICAL                                            :: fix_meta_stag=.FALSE.
      LOGICAL                                            :: bit64=.FALSE.
      LOGICAL                                            :: first=.TRUE.
      LOGICAL                                            :: firsteta=.TRUE.
! GMS.UC: Lluis Dec.09
      REAL                                               :: p_top
      REAL, ALLOCATABLE, DIMENSION(:,:,:)                :: field2d_out
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:)              :: datageogrid
      CHARACTER(LEN=11), ALLOCATABLE, DIMENSION(:)       :: varnames
      INTEGER, DIMENSION(3)                              :: jshape2d, jshape1d
      INTEGER, DIMENSION(4)                              :: dims2d_out, dims1d_out
      INTEGER                                            :: grid_filt, ntimes_filt
      INTEGER                                            :: ncgid
      INTEGER                                            :: ivar1, ivar2, k 
      CHARACTER(LEN=11)                                  :: unitsIN, varnameIN
      CHARACTER(LEN=50)                                  :: section
      CHARACTER(LEN=60)                                  :: longdescIN
      CHARACTER(LEN=250)                                 :: path_to_geofile, geofile, geofilename
      CHARACTER(LEN=300)                                 :: inputlong_name
      LOGICAL                                            :: unstagvar
      INTEGER                                            :: idim
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)            :: data4
! GMS.UC: Markel March 2012 
      LOGICAL                                            :: output_netCDF4=.FALSE. 
      NAMELIST /io/ path_to_input, input_name, path_to_output, output_name,                     &
        process, fields, debug, bit64, output_netCDF4, grid_filt, ntimes_filt, path_to_geofile, geofile
      NAMELIST /interp_in/ Netalevs, firsteta, eta_levels, unstagger_grid, p_top

      path_to_input   = './'
      path_to_output  = './'
      output_name     = ' '
      eta_levels   = -99999.
      process         = 'all'

      ! Read parameters from Fortran namelist
        DO funit=10,100
           INQUIRE(unit=funit, opened=is_used)
           IF (.not. is_used) EXIT
        END DO
        OPEN(funit,file='namelist.nointerp',status='old',form='formatted',iostat=ios)
        IF ( ios /= 0 ) STOP "ERROR opening namelist.nointerp"
        READ(funit,io)
        READ(funit,interp_in)
        CLOSE(funit)

      PRINT *,'Filtering of fields will be done with grid_filt:',grid_filt,' ntimes_filt:',     &
        ntimes_filt
      ! Get all the input file names

        lent = len_trim(path_to_input)
        IF ( path_to_input(lent:lent) /= "/" ) THEN
           path_to_input = TRIM(path_to_input)//"/"
        ENDIF
        lent = len_trim(path_to_output)
        IF ( path_to_output(lent:lent) /= "/" ) THEN
           path_to_output = TRIM(path_to_output)//"/"
        ENDIF
        inputlong_name = ' '
        inputlong_name = TRIM(path_to_input)//TRIM(input_name)
! GMS.UC: Lluis Dec.09
        ! Geofile name file
        lent = len_trim(path_to_geofile)
        IF ( path_to_geofile(lent:lent) /= "/" ) THEN
           path_to_geofile = TRIM(path_to_geofile)//"/"
        ENDIF

        geofilename = TRIM(path_to_geofile)//TRIM(geofile)

        !  Build a UNIX command, and "ls", of all of the input files 
        loslen = LEN ( command )
        CALL all_spaces ( command , loslen ) 
!        WRITE ( command , FMT='("ls -1 ",A300," > .foo")' ) TRIM ( inputlong_name )
        WRITE ( command , FMT='("ls -1 ",A," > .foo")' ) TRIM ( inputlong_name )
        !  We stuck all of the matching files in the ".foo" file.  Now we place the 
        !  number of the those file (i.e. how many there are) in ".foo1". 
  
        CALL SYSTEM ( TRIM ( command ) ) 
        CALL SYSTEM ( '( cat .foo | wc -l > .foo1 )' )
  
        !  Read the number of files.
        OPEN (FILE   = '.foo1'       , &
              UNIT   = 112           , &
              STATUS = 'OLD'         , &
              ACCESS = 'SEQUENTIAL'  , &
              FORM   = 'FORMATTED'     )
  
        READ ( 112 , * ) number_of_input_files
        CLOSE ( 112 )
  
        !  If there are zero files, we are toast.
        IF ( number_of_input_files .LE. 0 ) THEN
           print*, ' Oops, we need at least ONE input file for the program to read.'
           print*, '       Make sure you have the path, and name file(s) correct,'
           print*, '       including wild characters if needed.'
           STOP
        END IF
  
        !  Allocate space for this many files.
        ALLOCATE (  input_file_names(number_of_input_files) , STAT=ierr )
        ALLOCATE ( output_file_names(number_of_input_files) , STAT=ierr )

        !  Did the allocate work OK?
        IF ( ierr .NE. 0 ) THEN
           print*, ' tried to allocate ', number_of_input_files, ' input files, (look at ./foo)'
           STOP
        END IF

        !  Initialize all of the file names to blank.
         input_file_names = '                                                  ' // &
                            '                                                  ' // &
                            '                                '
        output_file_names = '                                                  ' // &
                            '                                                  ' // &
                            '                                '
  
        !  Open the file that has the list of filenames.
        OPEN (FILE   = '.foo'        , &
              UNIT   = 111           , &
              STATUS = 'OLD'         , &
              ACCESS = 'SEQUENTIAL'  , &
              FORM   = 'FORMATTED'     )
  
        !  Read all of the file names and store them.
        DO loop = 1 , number_of_input_files
           READ ( 111 , FMT='(A)' ) input_file_names(loop)
           IF ( output_name == ' ' ) THEN
              ilen = INDEX(TRIM(input_file_names(loop)),'/',.TRUE.) 
              output_file_names(loop) = TRIM(path_to_output)//input_file_names(loop)(ilen+1:)
              output_file_names(loop) = TRIM(output_file_names(loop))//"_ETALEV"
           ELSE
              IF ( number_of_input_files == 1 ) THEN
                 output_file_names(loop) = TRIM(path_to_output)//TRIM(output_name)
              ELSE
                 write(tmp_name,'(A,A,"_",I4.4)') TRIM(path_to_output), TRIM(output_name), loop
                 output_file_names(loop) = tmp_name
              ENDIF
           ENDIF
        END DO
        CLOSE ( 112 )
        print*, " " 
  
        !   We clean up our own messes.
        CALL SYSTEM ( '/bin/rm -f .foo'  )
        CALL SYSTEM ( '/bin/rm -f .foo1' )
     
        ! Do we have a list of field that we want on output?
        process_these_fields = ','
        IF ( INDEX(process,'list') /= 0) THEN
          DO i = 1 , len(fields)
            IF (fields(i:i) /= ' ' ) THEN
              process_these_fields = trim(process_these_fields)//fields(i:i)
            ENDIF
          END DO
          process_these_fields = trim(process_these_fields)//","
        END IF

!! L. Fita March 2011. GMS-UC
!
      IF (firsteta) THEN
        eta_levels=0
        DO i=1, Netalevs     
	  eta_levels(i)=i
	END DO
      END IF

      write(6,*) 
      write(6,*) "##############################################################"
      write(6,'(A,i4,A)') " RUNNING no_interp V1.0 on ", number_of_input_files, " file(s)."
      IF (INDEX(process,'list') /= 0) PRINT *,'Fields to process:',TRIM(process_these_fields)
      LINLOG = interp_method
      write(6,*) 
      IF (.not. unstagger_grid) write(6,*)"Data will be output on C-grid" 
      IF (unstagger_grid) write(6,*)"Data will be output on unstaggered grid"

      ! Find the pressure levels to interpolate to
      write(6,*)
      write(6,*) "REDUCING to first'",Netalevs,"'ETA LEVELS:"
      num_eta_levels = Netalevs
      DO i=1, Netalevs
        PRINT *,'   ',eta_levels(i)
      END DO
      write(6,*)

      DO loop = 1, number_of_input_files

        IF (debug) write(6,*) "##############################################################"
        input_file  = input_file_names(loop)
        output_file = output_file_names(loop)

        IF ( .not. debug ) write(6,*) " Output will be written to: ",trim(output_file)

        IF (debug) THEN
          write(6,*) " INPUT FILE:         ",trim(input_file)
          write(6,*) " OUTPUT FILE:        ",trim(output_file)
          write(6,*) "  "
        ENDIF


! OPEN INPUT AND OUTPUT FILE
! output_file is input_file_new
        rcode = nf_open(input_file, 0, ncid)
        if (rcode .ne. nf_noerr) call handle_err(rcode)
        if (bit64) then
          rcode = nf_create(output_file, NF_64BIT_OFFSET, mcid)
        elseif (output_netCDF4) then
          rcode = nf_create(output_file, IOR(NF_NETCDF4, NF_CLASSIC_MODEL), mcid)
        else
          rcode = nf_create(output_file, 0, mcid)
        endif
        if (rcode .ne. nf_noerr) call handle_err(rcode)


! GET BASIC INFORMTION ABOUT THE FILE
! most important 
!   ndims:  number of dimensions
!   nvars:  number of variables
!   ngatts: number of global attributes
        rcode = nf_inq(ncid, ndims, nvars, ngatts, nunlimdimid)
        if (rcode .ne. nf_noerr) call handle_err(rcode)
        IF (debug) THEN
          write(6,*) ' INPUT file has = ',ndims, ' dimensions, '
          write(6,*) '                  ',nvars, ' variables, and '      
          write(6,*) '                  ',ngatts,' global attributes '
          write(6,*) "  "
        ENDIF
        rcode = nf_get_att_int (ncid, nf_global, 'WEST-EAST_GRID_DIMENSION', iweg)
        rcode = nf_get_att_int (ncid, nf_global, 'SOUTH-NORTH_GRID_DIMENSION', isng)
        rcode = nf_get_att_int (ncid, nf_global, 'BOTTOM-TOP_GRID_DIMENSION', ibtg)

! ALLOCATE SOME VARIABLES
        IF (ALLOCATED(dnamei)) deallocate(dnamei)
            ALLOCATE (dnamei(20))
        IF (ALLOCATED(dnamej)) deallocate(dnamej)
            ALLOCATE (dnamej(20))
        IF (ALLOCATED(dvali)) deallocate(dvali)
            ALLOCATE (dvali(20))
        IF (ALLOCATED(dvalj)) deallocate(dvalj)
            ALLOCATE (dvalj(20))


! READ ALL DIMS FROM INPUT FILE AND CREATE SOME DIMS FOR OUTPUT FILE
        j = 0
        DO i = 1, ndims
          rcode = nf_inq_dim(ncid, i, dnamei(i), dvali(i))
  
          IF ( dnamei(i) == "Time" ) THEN
            j = j + 1
            dnamej(j) = dnamei(i)
            dvalj(j) = dvali(i)
            rcode = nf_def_dim(mcid, dnamej(j), NF_UNLIMITED, j)
            times_in_file  = dvali(i)
!          ELSE IF ( dnamei(i) == "bottom_top") THEN
!	    dvalj(j) = num_eta_levels
!          ELSE IF ( dnamei(i) == "bottom_top_stag") THEN
!	    dvalj(j) = num_eta_levels + 1
          ENDIF 
        ENDDO

!! L. Fita GMS-UC March 2011
!
!        !!! Create a eta dims
!        j = j + 1
!        dnamej(j) = 'bottom_top'
!        dvalj(j) = num_eta_levels
!        rcode = nf_def_dim(mcid, dnamej(j), dvalj(j), j)
!        j = j + 1
!        dnamej(j) = 'num_eta_levels_stag'
!        dvalj(j) = num_eta_levels
!        rcode = nf_def_dim(mcid, dnamej(j), dvalj(j), j)

! DEALING WITH THE GLOBAL ATTRIBUTES
        IF (debug) THEN
          write(6,*) 
          write(6,*) " OUTPUT FILE attributes:"
        ENDIF
        do i = 1, ngatts
          rcode = nf_inq_attname(ncid, nf_global, i,    cname)
          rcode = nf_inq_atttype(ncid, nf_global, cname, itype)
          rcode = nf_inq_attlen (ncid, nf_global, cname, ilen)
  
	  if ( itype .eq. 2 ) then        ! characters
	    rcode = nf_get_att_text (ncid, nf_global, cname, cval)
	    if(cname(1:5) .eq. 'TITLE') then
               cval = cval(1:ilen)//" - ON PRES LEVELS"
               ilen = len_trim(cval)
            endif
            IF (debug) &
	      write(6,'("     i = ",i2," : ",A," = ",A)') &
                    i,cname,cval(1:ilen)
	    rcode = nf_put_att_text(mcid, nf_global, cname, ilen,&
                      cval(1:ilen))
  
	  elseif ( itype .eq. 4 ) then     ! integers
	    rcode = nf_get_att_int (ncid, nf_global, cname, ival)
            IF ( INDEX(cname,'BOTTOM-TOP_PATCH') == 0 ) THEN
	       IF (cname .eq. 'BOTTOM-TOP_GRID_DIMENSION') ival = num_eta_levels
               IF (debug) &
	         write(6,'("     i = ",i2," : ",A," = ",i7)') &
                       i,cname,ival        
               rcode = nf_put_att_int(mcid, nf_global, cname, itype,&
                         ilen, ival)
             ENDIF
  
	  elseif ( itype .eq. 5 ) then    ! real
	    rcode = nf_get_att_real (ncid, nf_global, cname, rval)
            IF (debug) &
	      write(6,'("     i = ",i2," : ",A," = ",G18.10E2)') &
                    i,cname,rval
	    rcode = nf_put_att_real(mcid, nf_global, cname, itype,&
                      ilen, rval)
	  end if
        enddo

rcode = nf_enddef(mcid)
!GMS.UC: Lluis Dec.09
! Repeating variables list formation
        ! Do we have a list of field that we want on output?
        process_these_fields = ','
        IF ( INDEX(process,'list') /= 0) THEN
          DO i = 1 , len(fields)
            IF (fields(i:i) /= ' ' ) THEN
              process_these_fields = trim(process_these_fields)//fields(i:i)
            ENDIF
          END DO
          process_these_fields = trim(process_these_fields)//","
        END IF

! WE NEED SOME BASIC FIELDS
! GMS.UC: Lluis Dec.09
         PRINT *,'Extracting basic fields...'
         PRINT *,'  P_TOP'
       	 ivar1=-1
         rcode = nf_inq_varid(ncid, "P_TOP", ivar1)
         IF (ALLOCATED(data1)) deallocate(data1)
         allocate (data1(times_in_file,1,1,1))
         IF (ivar1 == 1 ) THEN
           PRINT *,"Extracting 'P_TOP' from 'namelist.nointerp'"
           data1=p_top
         ELSE
           rcode = nf_inq_varid    ( ncid, "P_TOP", i )
           rcode = nf_get_var_real ( ncid, i, data1 )
           IF ( first ) THEN
            first = .FALSE.
           ENDIF
         ENDIF
         deallocate (data1)
         PRINT *,'  P, PB...'
         IF (ALLOCATED(pres_field)) deallocate(pres_field)
         allocate (pres_field(iweg-1, isng-1, ibtg-1, times_in_file ))
         IF (ALLOCATED(data1)) deallocate(data1)
         allocate (data1(iweg-1, isng-1, ibtg-1, times_in_file ))
         rcode = nf_inq_varid    ( ncid, "P", i )
         rcode = nf_get_var_real ( ncid, i, pres_field )
         rcode = nf_inq_varid    ( ncid, "PB", i )
         rcode = nf_get_var_real ( ncid, i, data1 )
         pres_field = pres_field + data1
         deallocate (data1)

         IF (ALLOCATED(pres_stagU)) deallocate(pres_stagU)
         IF (ALLOCATED(pres_stagV)) deallocate(pres_stagV)
         allocate (pres_stagU(iweg, isng-1, ibtg-1, times_in_file ))
         allocate (pres_stagV(iweg-1, isng, ibtg-1, times_in_file ))
         pres_stagU(1,:,:,:)        =  pres_field(1,:,:,:)
         pres_stagU(iweg,:,:,:)     =  pres_field(iweg-1,:,:,:)
         pres_stagU(2:iweg-1,:,:,:) = (pres_field(1:iweg-2,:,:,:) +    &
                                       pres_field(2:iweg-1,:,:,:))*.5
         pres_stagV(:,1,:,:)        =  pres_field(:,1,:,:)
         pres_stagV(:,isng,:,:)     =  pres_field(:,isng-1,:,:)
         pres_stagV(:,2:isng-1,:,:) = (pres_field(:,1:isng-2,:,:) +    &
                                       pres_field(:,2:isng-1,:,:))*.5
         PRINT *,'  PSFC...'
         IF (ALLOCATED(psfc)) deallocate(psfc)
         allocate (psfc(iweg-1, isng-1, times_in_file ))
         IF (ALLOCATED(data1)) deallocate(data1)
         allocate (data1(iweg-1, isng-1, 1, times_in_file ))
         rcode = nf_inq_varid    ( ncid, "PSFC", i )
         rcode = nf_get_var_real ( ncid, i, data1 )
         psfc(:,:,:) = data1(:,:,1,:)
         deallocate (data1)

         IF (ALLOCATED(pres_out)) deallocate(pres_out)
         allocate (pres_out(iweg-1, isng-1, num_eta_levels, times_in_file))
         do i = 1, num_eta_levels
           pres_out (:,:,i,:) = eta_levels(i)
         enddo
!! Lluis Feb-11. UC
        PRINT *,'  PH, PHB...'
        IF (ALLOCATED(height_field)) deallocate(height_field)
        allocate (height_field(iweg-1, isng-1, ibtg-1, times_in_file ))
        IF (ALLOCATED(phb)) deallocate(phb)
        allocate (phb(iweg-1, isng-1, ibtg, times_in_file ))
        IF (ALLOCATED(data1)) deallocate(data1)
        allocate (data1(iweg-1, isng-1, ibtg, times_in_file ))
        rcode = nf_inq_varid    ( ncid, "PH", i )
        rcode = nf_get_var_real ( ncid, i, data1 )
        rcode = nf_inq_varid    ( ncid, "PHB", i )
        rcode = nf_get_var_real ( ncid, i, phb )
        data1 = (data1 + phb) 
        height_field(:,:,1:ibtg-1,:) = (( data1(:,:,1:ibtg-1,:) + data1(:,:,2:ibtg,:) )*.5)/9.81
        deallocate (data1)

! GMS.UC: Lluis Dec.09
         PRINT *,'  PH, PHB...'
       	 ivar1=-1
         rcode = nf_inq_varid(ncid, "PH", ivar1)
         ivar2=-1
         rcode = nf_inq_varid(ncid, "PHB", ivar2)
         IF (ivar1*ivar2 > 1 ) THEN
           IF (ALLOCATED(ght)) deallocate(ght)
           allocate (ght(iweg-1, isng-1, ibtg-1, times_in_file ))
           IF (ALLOCATED(phb)) deallocate(phb)
           allocate (phb(iweg-1, isng-1, ibtg, times_in_file ))
           IF (ALLOCATED(data1)) deallocate(data1)
           allocate (data1(iweg-1, isng-1, ibtg, times_in_file ))
           rcode = nf_inq_varid    ( ncid, "PH", i )
           rcode = nf_get_var_real ( ncid, i, data1 )
           rcode = nf_inq_varid    ( ncid, "PHB", i )
           rcode = nf_get_var_real ( ncid, i, phb )
           data1 = (data1 + phb) 
           ght(:,:,1:ibtg-1,:) = ( data1(:,:,1:ibtg-1,:) + data1(:,:,2:ibtg,:) )*.5

           deallocate (data1)
           deallocate (phb)
         ENDIF

! GMS.UC: Lluis Dec.09
         PRINT *,'  HGT...'
       	 ivar1=-1
         rcode = nf_inq_varid(ncid, "HGT", ivar1)
         IF (ALLOCATED(ter)) deallocate(ter)
         allocate (ter(iweg-1, isng-1))
         IF (ivar1 == 1) THEN
           PRINT *,"Obtaining terrain height from domain file:'"//TRIM(geofilename)//"'" 
!           rcode = nf_open(geofilename, 0, ncgid)
!           rcode = nf_inq_varid    ( ncgid, "HGT_M", ivar1 )
!           rcode = nf_get_var_real ( ncgid, ivar1, ter )
!           rcode = nf_close ( ncgid )
           IF (ALLOCATED(varnames)) DEALLOCATE(varnames)
           ALLOCATE (varnames(1))
           varnames=RESHAPE((/'HGT_M'/),(/1/))
           PRINT *,'varnames allocated...'
           IF (ALLOCATED(datageogrid)) DEALLOCATE(datageogrid)
           ALLOCATE (datageogrid(iweg-1, isng-1, ibtg-1, 1))
           PRINT *,'varnames and datageogrid matrixs allocated...'
           CALL extract_from_geogrid(geofilename, varnames, 1, iweg-1, isng-1, ibtg-1, datageogrid) 
           PRINT *,'datageogrid:',datageogrid(iweg/2,isng/2,1,1)
           ter=datageogrid(1:(iweg-1),1:(isng-1),1,1)
           DEALLOCATE(varnames, datageogrid)
         ELSE
           IF (ALLOCATED(data1)) deallocate(data1)
           allocate (data1(iweg-1, isng-1, 1, times_in_file ))
           rcode = nf_inq_varid    ( ncid, "HGT", i )
           rcode = nf_get_var_real ( ncid, i, data1 )
           ter(:,:) = data1(:,:,1,1)
           deallocate (data1)
         ENDIF 

! Height is relative above specific grid topography
!
        DO i=1,ibtg-1
	  DO jt=1,times_in_file
            height_field(:,:,i,jt) = height_field(:,:,i,jt) - ter
            DO k=1,iweg-1
	      DO l=1,isng-1
                IF (height_field(k,l,i,jt) < 0.) PRINT *,'Negative height ! (',&
		  k, char(44), l, char(44), i, char(44), jt, char(44),') :',   &
		  height_field(k,l,i,jt)
	       END DO
             END DO
	  END DO
	END DO

        IF (ALLOCATED(height_stagU)) deallocate(height_stagU)
        IF (ALLOCATED(height_stagV)) deallocate(height_stagV)
        allocate (height_stagU(iweg, isng-1, ibtg-1, times_in_file ))
        allocate (height_stagV(iweg-1, isng, ibtg-1, times_in_file ))
        height_stagU(1,:,:,:)	=  height_field(1,:,:,:)
        height_stagU(iweg,:,:,:)	=  height_field(iweg-1,:,:,:)
        height_stagU(2:iweg-1,:,:,:) = (height_field(1:iweg-2,:,:,:) +    &
                                      height_field(2:iweg-1,:,:,:))*.5
        height_stagV(:,1,:,:)	=  height_field(:,1,:,:)
        height_stagV(:,isng,:,:)	=  height_field(:,isng-1,:,:)
        height_stagV(:,2:isng-1,:,:) = (height_field(:,1:isng-2,:,:) +    &
	                              height_field(:,2:isng-1,:,:))*.5
				      
!! L. Fita. UC-GMS March 2011
        PRINT *,'  ZNU...'
        IF (ALLOCATED(znu)) deallocate(znu)
        allocate (znu(ibtg-1, times_in_file ))
        rcode = nf_inq_varid    ( ncid, "ZNU", i )
        rcode = nf_get_var_real ( ncid, i, znu )

        PRINT *,'  ZNW...'
        IF (ALLOCATED(znw)) deallocate(znw)
        allocate (znw(ibtg, times_in_file ))
        rcode = nf_inq_varid    ( ncid, "ZNW", i )
        rcode = nf_get_var_real ( ncid, i, znw )
 
        PRINT *,'Will be retrieved the following eta on half/full levels'
	DO i=1, Netalevs
	  PRINT *,'lev #:',eta_levels(i),' half-eta:',znu(eta_levels(i), 1), &
	    ' full-eta:',znw(eta_levels(i), 1)
	END DO

        PRINT *,'  QVAPOR...'
        IF (ALLOCATED(qv)) deallocate(qv)
        allocate (qv(iweg-1, isng-1, ibtg-1, times_in_file ))
        rcode = nf_inq_varid    ( ncid, "QVAPOR", i )
        rcode = nf_get_var_real ( ncid, i, qv )

        PRINT *,'  TK...'
        IF (ALLOCATED(tk)) deallocate(tk)
        allocate (tk(iweg-1, isng-1, ibtg-1, times_in_file ))
        IF (ALLOCATED(data1)) deallocate(data1)
        allocate (data1(iweg-1, isng-1, ibtg-1, times_in_file ))
        rcode = nf_inq_varid    ( ncid, "T", i )
        rcode = nf_get_var_real ( ncid, i, data1 )
        tk = (data1+300.) * ( pres_field / p0 )**RCP
        deallocate (data1)

        PRINT *,'  RH...'
        IF (ALLOCATED(rh)) deallocate(rh)
        allocate (rh(iweg-1, isng-1, ibtg-1, times_in_file ))
        IF (ALLOCATED(data1)) deallocate(data1)
        IF (ALLOCATED(data2)) deallocate(data2)
        allocate (data1(iweg-1, isng-1, ibtg-1, times_in_file ))
        allocate (data2(iweg-1, isng-1, ibtg-1, times_in_file ))
        data1 = 10.*0.6112*exp(17.67*(tk-273.16)/(TK-29.65))
        data2 = 0.622*data1/(0.01 * pres_field -  (1.-0.622)*data1)
        rh = 100.*qv/data2
        rh    = 100.*AMAX1(AMIN1(qv/data2,1.0),0.0)
        deallocate (data1)
        deallocate (data2)

! TRAIN FILE 
        IF (debug) THEN
          write(6,*) 
          write(6,*) 
          write(6,*) "FILE variables:"
        ENDIF
        jvar = 0
        loop_variables : DO ivar = 1, nvars
          rcode = nf_inq_var(ncid, ivar, cval, itype, idm, ishape, natt)

          !!! Do we want this variable
          IF ( trim(cval) == 'P'  .OR. trim(cval) == 'PB' )  CYCLE loop_variables
          IF ( trim(cval) == 'PH' .OR. trim(cval) == 'PHB' ) CYCLE loop_variables
          IF ( trim(cval) == 'T' ) CYCLE loop_variables
          IF ( unstagger_grid .AND. (INDEX(cval,'_U') /= 0) ) CYCLE loop_variables !!! no sense in keeping these
          IF ( unstagger_grid .AND. (INDEX(cval,'_V') /= 0) ) CYCLE loop_variables !!! no sense in keeping these
          IF ( INDEX(process,'all') == 0 ) THEN
             !!! Only want some variables - see which
             dummy = ","//trim(cval)//","
             is_there = INDEX(process_these_fields,trim(dummy))
             IF ( is_there == 0 ) THEN
                IF ( debug ) print*,"NOTE: ", trim(cval), " - Not requested"
                CYCLE loop_variables !!! don't want this one
             ENDIF
          ENDIF

          IF ( idm >= 4 .AND. itype == 4 ) THEN
            print*,"NOTE: We cannot deal with 3D integers - maybe later"
            CYCLE loop_variables
          ENDIF
          IF ( itype == 6 ) THEN
            print*,"NOTE: We cannot deal with double precision data - maybe later"
            CYCLE loop_variables
          ENDIF
          IF ( itype == 2 .OR. itype == 4 .OR. itype == 5 ) THEN
            !!! OK I know what to do this this
          ELSE
            print*,"NOTE: Do not understand this data type ", itype, " skip field."
            CYCLE loop_variables
          ENDIF
 
          IF ( trim(cval) == 'U' )   cval = 'UU'
          IF ( trim(cval) == 'V' )   cval = 'VV'
          IF ( trim(cval) == 'TSK' ) cval = 'SKINTEMP'

          !!! OK - we want this - lets continue
          jvar = jvar + 1
          jshape = 0
          interpolate = .FALSE.
          fix_meta_stag = .FALSE.
          rcode = nf_redef(mcid)

          DO ii = 1, idm
            unstagvar = .FALSE. 
            test_dim_name = dnamei(ishape(ii))
            IF ( test_dim_name == 'bottom_top' .OR. test_dim_name == 'bottom_top_stag' ) THEN
                 IF ( test_dim_name == 'bottom_top_stag' ) fix_meta_stag = .TRUE.
!                 test_dim_name = 'num_eta_levels'
                 test_dim_name = 'bottom_top'
                 interpolate = .TRUE.
                 unstagvar = .TRUE. 
            ENDIF
            IF ( unstagger_grid .AND. test_dim_name == 'west_east_stag' )   THEN
               test_dim_name = 'west_east'
               fix_meta_stag = .TRUE.
               unstagvar = .TRUE. 
            ENDIF
            IF ( unstagger_grid .AND. test_dim_name == 'south_north_stag' ) THEN
               test_dim_name = 'south_north'
               fix_meta_stag = .TRUE.
               unstagvar = .TRUE. 
            ENDIF
            DO jj = 1,j
              IF ( test_dim_name == dnamej(jj) ) THEN
                jshape(ii) = jj
              ENDIF
            ENDDO
          IF ( jshape(ii) == 0 ) THEN
!! GMS.UC: Lluis. Jan 10 Correction in assignation of new dimensions (new dims) id
              j = j + 1
              jshape(ii) = j
              IF (unstagvar) THEN
                dvalj(j) = dvali(ishape(ii)) - 1
                dnamej(j) = test_dim_name
              ELSE
                dvalj(j) = dvali(ishape(ii))
                dnamej(j) = dnamei(ishape(ii))
              END IF
              IF (dnamej(j) == 'bottom_top' ) dvalj(j) = num_eta_levels
              IF (dnamej(j) == 'bottom_top_stag' ) dvalj(j) = num_eta_levels + 1
              rcode = nf_def_dim(mcid, dnamej(j), dvalj(j), j)
            ENDIF
          ENDDO

          rcode = nf_def_var(mcid, cval, itype, idm, jshape, jvar)
          IF (output_netCDF4) THEN
            rcode = nf_def_var_deflate(mcid, jvar, 1, 1, 3)
          ENDIF

          DO na = 1, natt
             rcode = nf_inq_attname(ncid, ivar, na, cname)
             IF ( fix_meta_stag .AND. trim(cname) == 'stagger' ) THEN
                att_text = "-"
                ilen = len_trim(att_text)
                rcode = nf_put_att_text(mcid, jvar, cname, ilen, att_text(1:ilen) )
             ELSEIF ( fix_meta_stag .AND. trim(cname) == 'coordinates' ) THEN
                att_text = "XLONG XLAT"
                ilen = len_trim(att_text)
                rcode = nf_put_att_text(mcid, jvar, cname, ilen, att_text(1:ilen) )
             ELSE
                rcode = nf_copy_att(ncid, ivar, cname, mcid, jvar)
             ENDIF
          ENDDO
          rcode = nf_put_att_real(mcid, jvar, 'missing_value', NF_FLOAT, 1, MISSING)
          rcode = nf_enddef(mcid)


! GET THE DIMS FOR INPUT AND OUTPUT FROM THE SHAPE
          dims_in  = 1
          dims_out = 1
          DO ii = 1,idm
            dims_in(ii)  = dvali(ishape(ii))
            dims_out(ii) = dvalj(jshape(ii))
          ENDDO
          IF (debug) THEN
            write(6,*) 'VAR: ',trim(cval)
            write(6,*) '     DIMS  IN: ',dims_in
            write(6,*) '     DIMS OUT: ',dims_out
          ENDIF

! ALLOCATE THE INPUT AND OUTPUT ARRAYS
! READ THE DATA FROM INPUT FILE

	  IF     (itype == 2) THEN          ! character
            allocate (text(dims_in(1), dims_in(2), dims_in(3), dims_in(4)))
	    rcode = nf_get_var_text(ncid, ivar, text)
	    rcode = nf_put_vara_text (mcid, jvar, start_dims, dims_in, text)
            IF (debug) write(6,*) '     SAMPLE VALUE = ',text(:,1,1,1)
            deallocate (text)
  
	  ELSEIF (itype == 4) THEN          ! integer
            allocate (idata1(dims_in(1), dims_in(2), dims_in(3), dims_in(4)))
	    rcode = nf_get_var_int(ncid, ivar, idata1)
	    rcode = nf_put_vara_int (mcid, jvar, start_dims, dims_in, idata1)
            IF (debug) write(6,*) '     SAMPLE VALUE = ',idata1(dims_in(1)/2,dims_in(2)/2,1,1)
            deallocate (idata1)  
  
	  ELSEIF (itype == 5) THEN          ! real
            allocate (data1(dims_in(1), dims_in(2), dims_in(3), dims_in(4)))
            allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
	    rcode = nf_get_var_real(ncid, ivar, data1)
  
            IF (idm >= 4 .AND. interpolate) THEN  
               IF (debug) write(6,*) '     THIS IS A FIELD WE NEED TO INTERPOLATE'       

               IF ( dims_in(1) == iweg .AND. .not. unstagger_grid ) THEN
                  CALL nointerp (data2, data1, eta_levels,                                 &
		                 iweg, isng-1, ibtg-1, dims_in(4),                           &
                                 num_eta_levels)
               ELSEIF ( dims_in(2) == isng .AND. .not. unstagger_grid ) THEN
                  CALL nointerp (data2, data1, eta_levels,    &
                                 iweg-1, isng, ibtg-1, dims_in(4), &
                                 num_eta_levels)
               ELSEIF ( dims_in(1) == iweg .AND. unstagger_grid ) THEN
                  allocate (data3(iweg-1, isng-1, ibtg-1, dims_in(4)))
                  data3(1:iweg-1,:,:,:) = (data1(1:iweg-1,:,:,:) + data1(2:iweg,:,:,:)) * .5
                  CALL nointerp (data2, data3, eta_levels,   &
                                 iweg-1, isng-1, ibtg-1, dims_in(4), &
                                 num_eta_levels)
                  deallocate(data3)
               ELSEIF ( dims_in(2) == isng .AND. unstagger_grid ) THEN
                  allocate (data3(iweg-1, isng-1, ibtg-1, dims_in(4)))
                  data3(:,1:isng-1,:,:) = (data1(:,1:isng-1,:,:) + data1(:,2:isng,:,:)) * .5
                  CALL nointerp (data2, data3, eta_levels,   &
                                 iweg-1, isng-1, ibtg-1, dims_in(4), &
                                 num_eta_levels)
                  deallocate(data3)
               ELSEIF ( dims_in(3) == ibtg ) THEN
                  allocate (data3(iweg-1, isng-1, ibtg-1, dims_in(4)))
                  data3(:,:,1:ibtg-1,:) = (data1(:,:,1:ibtg-1,:) + data1(:,:,2:ibtg,:)) * .5
                  CALL nointerp (data2, data3, eta_levels,   &
                                 dims_in(1), dims_in(2), ibtg-1, dims_in(4), &
                                 num_eta_levels)
                  deallocate(data3)
               ELSE
                  CALL nointerp (data2, data1, eta_levels,   &
                                 dims_in(1), dims_in(2), dims_in(3), dims_in(4), &
                                 num_eta_levels)
               END IF

               IF (debug) write(6,*) '     SAMPLE VALUE IN  = ',data1(dims_in(1)/2,dims_in(2)/2,1,1)
               IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(dims_out(1)/2,dims_out(2)/2,1,1)

            ELSEIF (idm == 3 .AND. unstagger_grid ) THEN  
               IF ( dims_in(1) == iweg ) THEN
                  data2(1:iweg-1,:,:,:) = (data1(1:iweg-1,:,:,:) + data1(2:iweg,:,:,:)) * .5
               ELSEIF ( dims_in(2) == isng ) THEN
                  data2(:,1:isng-1,:,:) = (data1(:,1:isng-1,:,:) + data1(:,2:isng,:,:)) * .5
               ELSEIF ( dims_in(3) == ibtg ) THEN
                  data2(:,:,1:ibtg-1,:) = (data1(:,:,1:ibtg-1,:) + data1(:,:,2:ibtg,:)) * .5
               ELSE
	         data2 = data1
               ENDIF
               IF (debug) write(6,*) '     SAMPLE VALUE  = ',data1(dims_in(1)/2,dims_in(2)/2,1,1)

	    ELSE
	       data2 = data1
               IF (debug) write(6,*) '     SAMPLE VALUE  = ',data1(dims_in(1)/2,dims_in(2)/2,1,1)

	    ENDIF
            
	    rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)

            deallocate (data1)
            deallocate (data2)
  
  
          ENDIF
  
        ENDDO loop_variables

!!! We have some special variables we are interested in: PRES, TT, GHT, MSLP
        IF ( debug ) print*," "
        IF ( debug ) print*,"Calculating some diagnostics..."
           
        jshape = 0
        DO ii = 1, 4
          IF ( ii == 1 ) test_dim_name = 'west_east'
          IF ( ii == 2 ) test_dim_name = 'south_north'
          IF ( ii == 3 ) test_dim_name = 'bottom_top'
          IF ( ii == 4 ) test_dim_name = 'Time'
          DO jj = 1,j
            IF ( test_dim_name == dnamej(jj) ) THEN
              jshape(ii) = jj
            ENDIF
          ENDDO

          IF ( jshape(ii) == 0 ) THEN
            j = j + 1
            jshape(ii) = j
            dnamej(j) = dnamei(ishape(ii))
            dvalj(j) = dvali(ishape(ii))
            rcode = nf_def_dim(mcid, dnamej(j), dvalj(j), j)
          ENDIF
        ENDDO
        dims_in  = 1
        dims_out = 1
        DO ii = 1,4
          dims_out(ii) = dvalj(jshape(ii))
        ENDDO

        interpolate = .TRUE.

! GMS.UC:Lluis March.11
! Allways write vector with vertical eta levels
!        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'ZNU') /= 0 ) THEN
          !!! ZNU: eta levels vector 
           jvar = jvar + 1
           jshape1d=RESHAPE((/jshape(3),jshape(2),jshape(1)/),(/3/))
           varnameIN='ZNU'
           longdescIN='eta values on half (mass) levels'
           unitsIN=''
           CALL def_var (mcid, jvar, varnameIN, 5, 1, jshape1d, "Z", longdescIN, unitsIN, "-", &
             "XLONG XLAT", output_netCDF4)
           IF (debug) THEN
             write(6,*) 'VAR: ZNU idvar:',jvar
             write(6,*) '     DIMS OUT: ',num_eta_levels
           ENDIF
           ALLOCATE (data2(1,1,num_eta_levels,1))
           data2(1,1,:,1)=znu(eta_levels(1:Netalevs),1)
           rcode = nf_put_vara_real (mcid, jvar, 1, num_eta_levels, data2(1,1,:,1))
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(1,1,num_eta_levels/2,1)
           DEALLOCATE(data2)
!         END IF

!        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'ZNW') /= 0 ) THEN
          !!! ZNU: eta levels vector 
           jvar = jvar + 1
           jshape1d=RESHAPE((/jshape(3),jshape(2),jshape(1)/),(/3/))
           varnameIN='ZNW'
           longdescIN='eta values on full (w) levels'
           unitsIN=''
           CALL def_var (mcid, jvar, varnameIN, 5, 1, num_eta_levels, "Z", longdescIN, unitsIN, "-", &
             "XLONG XLAT", output_netCDF4)
           IF (debug) THEN
             write(6,*) 'VAR: ZNW idvar:',jvar
             write(6,*) '     DIMS OUT: ',num_eta_levels
           ENDIF
           ALLOCATE (data2(1,1,dims_out(3),1))
           data2(1,1,:,1)=znw(eta_levels(1:Netalevs),1)
	   PRINT *,'data2: ',data2
           rcode = nf_put_vara_real (mcid, jvar, 1, num_eta_levels, data2(1,1,:,1))
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(1,1,num_eta_levels/2,1)
           DEALLOCATE(data2)
!         END IF

! Allways write vector with equivalent height levels
!!        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'HEIGHT') /= 0 ) THEN
           !!! GHT
           jvar = jvar + 1
           varnameIN='HEIGHT'
           longdescIN='height above surface'
           unitsIN='m'
           CALL def_var (mcid, jvar, varnameIN, 5, 4, jshape, "XZY", longdescIN, unitsIN, "-",  &
             "XLONG XLAT", output_netCDF4)
           IF (debug) THEN
             write(6,*) 'VAR: HEIGHT idvar:',jvar
             write(6,*) '     DIMS OUT: ',dims_out
           ENDIF
           allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
! GMS.UC:Lluis Dec.09
           CALL nointerpdiag (data2, height_field, eta_levels,  &
             iweg-1, isng-1, ibtg-1, dims_in(4), dims_out(4), num_eta_levels)
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',height_field(dims_out(1)/2,dims_out(2)/2,1,&
             1)
           deallocate(data2)
!!        ENDIF

        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'MSLP') /= 0 ) THEN
          PRINT *,'Computing and Writting Mean Sea Level pressure'
          !!! MSLP 
           jshape2d=RESHAPE((/jshape(1),jshape(2),jshape(4)/),(/3/))
           dims2d_out=RESHAPE((/dims_out(1),dims_out(2),dims_out(4),1/),(/4/))
           jvar = jvar + 1
           varnameIN='MSLP'
           longdescIN='Mean Sea level pressure'
           unitsIN='Pa'

           CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-", &
             "XLONG XLAT", output_netCDF4)
           IF (debug) THEN
             write(6,*) 'VAR: MSLP idvar:',jvar
             write(6,*) '     DIMS OUT: ',dims2d_out
           ENDIF
           IF ( ALLOCATED(data2) ) THEN
             DEALLOCATE(data2)
           ENDIF
           ALLOCATE (data2(dims_out(1), dims_out(2), dims_out(4), 1))
           data2=0.
           CALL mean_sealevelpress(data2, pres_field, eta_levels, psfc, ter, tk, qv,       &
             iweg-1, isng-1, ibtg-1, dims_in(4), dims_out(4),                                 &          
                          num_eta_levels)
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data2)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(dims_out(1)/2,dims_out(2)/2,1,1)
           DEALLOCATE(data2) 
         END IF
         IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'MSLPF') /= 0 ) THEN
         !!!! Filtered MSLP
           jvar = jvar + 1
           varnameIN='MSLPF'
           longdescIN='MSLP filtered 10 times 3x3'
           unitsIN='Pa'

           CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-", "XLONG XLAT", output_netCDF4)
           IF (ALLOCATED(data3)) DEALLOCATE(data3)
           ALLOCATE (data3(dims_out(1), dims_out(2), dims_out(4), 1))
           CALL spatialfiltering(data2, grid_filt, ntimes_filt, dims_out(1), dims_out(2),       &
             dims_out(4), 1, data3)
           IF (debug) THEN
             write(6,*) 'VAR: MSLPF idvar:',jvar
             write(6,*) '     DIMS OUT: ',dims2d_out
           ENDIF
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data3)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data3(dims_out(1)/2,dims_out(2)/2,1,1) 
           DEALLOCATE(data2, data3) 
         ENDIF
        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'RAINTOT') /= 0 ) THEN
        !!!! Accumulated total precipitation 
          jshape2d=RESHAPE((/jshape(1),jshape(2),jshape(4)/),(/3/))
          dims2d_out=RESHAPE((/dims_out(1),dims_out(2),dims_out(4),1/),(/4/))
          jvar = jvar + 1
          varnameIN='RAINTOT'
          longdescIN='Accumulated total precipitation'
          unitsIN='mm'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          IF (ALLOCATED(data3)) DEALLOCATE(data3)
          ALLOCATE (data3(dims_out(1), dims_out(2), dims_out(4), 1))
          IF (ALLOCATED(varnames)) DEALLOCATE(varnames)
          ALLOCATE (varnames(2))
          varnames=RESHAPE((/'RAINNC','RAINC'/),(/2/))
          CALL variablessum(ncid, varnames, 2, dims_out(1), dims_out(2), dims_out(4), 1, data3) 
          IF (debug) THEN
            write(6,*) 'VAR: RAINTOT idvar:',jvar
            write(6,*) '     DIMS OUT: ',dims2d_out
          ENDIF
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data3)
          IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data3(dims_out(1)/2,dims_out(2)/2,1,1)
          DEALLOCATE(data3)
        ENDIF

       IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'PRES') /= 0 ) THEN
           !!! PRES
           jvar = jvar + 1
           varnameIN='PRES'
           longdescIN='Pressure'
           unitsIN='Pa'

           CALL def_var (mcid, jvar, varnameIN, 5, 4, jshape, "XZY", longdescIN, unitsIN,       &
             "-", "XLONG XLAT", output_netCDF4)
           IF (debug) THEN
             write(6,*) 'VAR: PRES idvar:',jvar
             write(6,*) '     DIMS OUT: ',dims_out
           ENDIF
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, pres_out)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',pres_out(dims_out(1)/2,dims_out(2)/2&
             ,1,1)
        ENDIF

        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'TT') /= 0 ) THEN
           !!! TT
           jvar = jvar + 1
           varnameIN='TT'
           longdescIN='Temperature'
           unitsIN='K'

           CALL def_var (mcid, jvar, varnameIN, 5, 4, jshape, "XZY", longdescIN, unitsIN,       &
             "-", "XLONG XLAT", output_netCDF4)
           IF (debug) THEN
             write(6,*) 'VAR: TT idvar:',jvar
             write(6,*) '     DIMS OUT: ',dims_out
           ENDIF
           IF (ALLOCATED(data2)) DEALLOCATE(data2)
           allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
! GMS.UC:Lluis Dec.09
           CALL nointerpdiag (data2, tk, eta_levels,                                          &
             iweg-1, isng-1, ibtg-1, dims_in(4), dims_out(4), num_eta_levels)
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(dims_out(1)/2,dims_out(2)/2,1,&
             1)
           deallocate(data2)
        ENDIF

        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'GHT') /= 0 ) THEN
           !!! GHT
           jvar = jvar + 1
           varnameIN='GHT'
           longdescIN='Geopotential Height'
           unitsIN='m'

           CALL def_var (mcid, jvar, varnameIN, 5, 4, jshape, "XZY", longdescIN, unitsIN, "-",  &
             "XLONG XLAT", output_netCDF4)
           IF (debug) THEN
             write(6,*) 'VAR: GHT idvar:',jvar
             write(6,*) '     DIMS OUT: ',dims_out
           ENDIF
           allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
! GMS.UC:Lluis Dec.09
           CALL nointerpdiag (data2, ght, eta_levels,  &
             iweg-1, isng-1, ibtg-1, dims_in(4), dims_out(4), num_eta_levels)
           data2 = data2/9.81
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(dims_out(1)/2,dims_out(2)/2,1,&
             1)
           deallocate(data2)
        ENDIF

        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'RH') /= 0 ) THEN
           !!! RH
           jvar = jvar + 1
           varnameIN='RH'
           longdescIN='Relative Humidity'
           unitsIN='%'
           CALL def_var (mcid, jvar, varnameIN, 5, 4, jshape, "XZY", longdescIN, unitsIN, "-",  &
             "XLONG XLAT", output_netCDF4)
           IF (debug) THEN
             write(6,*) 'VAR: RH idvar:',jvar
             write(6,*) '     DIMS OUT: ',dims_out
           ENDIF
           allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
! GMS.UC:Lluis Dec.09
           CALL nointerpdiag (data2, rh, eta_levels,   &
             iweg-1, isng-1, ibtg-1, dims_in(4), dims_out(4), num_eta_levels)
           WHERE ( rh < 0.0 )
              rh = 0.0
           ENDWHERE
           WHERE ( rh > 100.0 )
              rh = 100.0
           ENDWHERE
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(dims_out(1)/2,dims_out(2)/2,1,&
             1)
           deallocate(data2)
        ENDIF

! Vertically integrated variables
!!
        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'VIM') /= 0 ) THEN
          !
          !   vertically integrated moisture
          ! 
          section = "'VIM'"
          IF (debug) PRINT *,'Section '//TRIM(section)//'... .. .'
          rcode = nf_inq_varndims(ncid, i, ndims)
          rcode = nf_inq_vardimid(ncid, i, ishape)
          DO idim=1, ndims
            rcode = nf_inq_dimlen(ncid, ishape(idim), dims_in(idim))
          END DO
          PRINT *,'ishape:', ishape
          PRINT *,'dims_in:', dims_in
          jshape2d=RESHAPE((/jshape(1),jshape(2),jshape(4)/),(/3/))
          dims2d_out=RESHAPE((/dims_in(1),dims_in(2),dims_in(4),1/),(/4/))
          jvar = jvar + 1
          varnameIN='VIM'
          longdescIN='Vertically integrated moisture'
          unitsIN='Kg m-2'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          IF (ALLOCATED(data3)) DEALLOCATE(data3)
          ALLOCATE (data3(dims2d_out(1), dims2d_out(2), dims2d_out(3), 1))
          CALL massvertint(ncid, (/'QVAPOR     '/), 1, debug, dims_in(1), dims_in(2), dims_in(3),&
            dims_in(4), data3) 
          IF (debug) THEN
            PRINT *,'VAR: VIM idvar:',jvar
            PRINT *,'     DIMS OUT: ',dims2d_out
          ENDIF
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data3)
          call error_nc(section, rcode)
          IF (debug) PRINT *,'     SAMPLE VALUE OUT = ',data3(dims_out(1)/2,dims_out(2)/2,1,1)
          DEALLOCATE(data3)
        ENDIF
!        STOP
        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'VIQC') /= 0 ) THEN
          !
          !   vertically integrated cloud water 
          !
          section="'VIQC'"
          IF (debug) PRINT *,'Section '//TRIM(section)//'... .. .'
          rcode = nf_inq_varndims(ncid, i, ndims)
          rcode = nf_inq_vardimid(ncid, i, ishape)
          DO idim=1, ndims
            rcode = nf_inq_dimlen(ncid, ishape(idim), dims_in(idim))
          END DO
          jshape2d=RESHAPE((/jshape(1),jshape(2),jshape(4)/),(/3/))
          dims2d_out=RESHAPE((/dims_in(1),dims_in(2),dims_in(4),1/),(/4/))
          jvar = jvar + 1
          varnameIN='VIQC'
          longdescIN='Vertically integrated cloud water'
          unitsIN='Kg m-2'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          IF (ALLOCATED(data3)) DEALLOCATE(data3)
          ALLOCATE (data3(dims_in(1), dims_in(2), dims_in(4), 1))
          CALL massvertint(ncid, (/'QCLOUD     '/), 1, debug, dims_in(1), dims_in(2), dims_in(3),&
            dims_in(4), data3)
          IF (debug) THEN
            PRINT *,'VAR: VIQC idvar:',jvar
            PRINT *,'     DIMS OUT: ',dims2d_out
          ENDIF
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data3)
          call error_nc(section, rcode)
          IF (debug) PRINT *, '     SAMPLE VALUE OUT = ',data3(dims_out(1)/2,dims_out(2)/2,1,1)
          DEALLOCATE(data3)
        ENDIF

        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'VIWIND') /= 0 ) THEN
          !
          !   vertically integrated wind 
          !
          section="'VIWIND'"
          IF (debug) PRINT *,'Section '//TRIM(section)//'... .. .'
          rcode = nf_inq_varndims(ncid, i, ndims)
          rcode = nf_inq_vardimid(ncid, i, ishape)
          DO idim=1, ndims
            rcode = nf_inq_dimlen(ncid, ishape(idim), dims_in(idim))
          END DO
          jshape2d=RESHAPE((/jshape(1),jshape(2),jshape(4)/),(/3/))
          dims2d_out=RESHAPE((/dims_in(1),dims_in(2),dims_in(4),1/),(/4/))
          IF (ALLOCATED(data3)) DEALLOCATE(data3)
          ALLOCATE (data3(dims_in(1), dims_in(2), dims_in(4), 3))
          CALL multi_massvertint(ncid, (/'U          ','V          ','W          '/), 3, debug,& 
            dims_in(1), dims_in(2), dims_in(3), dims_in(4), data3)
          IF (debug) THEN
            PRINT *, 'Vertical transport integration'
            PRINT *, 'VAR: VIU, VIV idvar:',jvar+1, char(44), jvar+2
            PRINT *, '     DIMS OUT: ',dims2d_out
          ENDIF

          jvar = jvar + 1
          varnameIN='VIU'
          longdescIN='Vertically integrated wind WE'
          unitsIN='s-1 m-1'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data3(:,:,:,1))
          CALL error_nc(section, rcode)
          IF (debug) PRINT *, '  VIU. SAMPLE VALUE OUT = ',data3(dims_out(1)/2,dims_out(2)/2,1,1)
          jvar = jvar + 1
          varnameIN='VIV'
          longdescIN='Vertically integrated wind SN'
          unitsIN='s-1 m-1'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data3(:,:,:,2))
          call error_nc(section, rcode)
          IF (debug) PRINT *, '  VIV. SAMPLE VALUE OUT = ',data3(dims_out(1)/2,dims_out(2)/2,1,2)
          jvar = jvar + 1
          varnameIN='VIW'
          longdescIN='Vertically integrated wind bottom-top'
          unitsIN='s-1 m-1'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data3(:,:,:,3))
          call error_nc(section, rcode)
          IF (debug) PRINT *, '  VIW. SAMPLE VALUE OUT = ',data3(dims_out(1)/2,dims_out(2)/2,1,3)

          DEALLOCATE(data3)
        ENDIF

        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'VIMWIND') /= 0 ) THEN
          !
          !   vertically integrated transport of moist
          !
          section="'VIMWIND'"
          IF (debug) PRINT *,'Section '//TRIM(section)//'... .. .'
          rcode = nf_inq_varndims(ncid, i, ndims)
          rcode = nf_inq_vardimid(ncid, i, ishape)
          DO idim=1, ndims
            rcode = nf_inq_dimlen(ncid, ishape(idim), dims_in(idim))
          END DO
          jshape2d=RESHAPE((/jshape(1),jshape(2),jshape(4)/),(/3/))
          dims2d_out=RESHAPE((/dims_in(1),dims_in(2),dims_in(4),1/),(/4/))
          IF (ALLOCATED(data4)) DEALLOCATE(data4)
          ALLOCATE (data4(dims_in(1), dims_in(2), dims_in(4), 1, 3))
          CALL multi_transportvertint(ncid, (/'QVAPOR     '/),1,(/'U          ','V          ',  &
            'W          '/),debug, dims_in(1), dims_in(2), dims_in(3), dims_in(4), data4) 
          IF (debug) THEN
            PRINT *,'Vertical transport integration of moist'
            PRINT *,'     VAR: VIMU, VIMV idvar:',jvar+1, char(44), jvar+2
            PRINT *,'     DIMS OUT: ',dims2d_out
          ENDIF

          jvar = jvar + 1
          varnameIN='VIMU'
          longdescIN='Vertically integrated transport WE of moist'
          unitsIN='kg s-1 m-1'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data4(:,:,:,1,1))
          IF (debug) PRINT *, '  VIMU. SAMPLE VALUE OUT = ',data4(dims_out(1)/2,dims_out(2)/2,1,1,1)
          jvar = jvar + 1
          varnameIN='VIMV'
          longdescIN='Vertically integrated transport SN of moist'
          unitsIN='kg s-1 m-1'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data4(:,:,:,1,2))
          IF (debug) PRINT *, '  VIMV. SAMPLE VALUE OUT = ',data4(dims_out(1)/2,dims_out(2)/2,1,1,2)
          jvar = jvar + 1
          varnameIN='VIMW'
          longdescIN='Vertically integrated transport bottom-top of moist'
          unitsIN='kg s-1 m-1'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data4(:,:,:,1,3))
          IF (debug) PRINT *, '  VIMV. SAMPLE VALUE OUT = ',data4(dims_out(1)/2,dims_out(2)/2,1,1,3)
          DEALLOCATE(data4)
        ENDIF

        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'VIQI') /= 0 ) THEN
          !
          !   vertically integrated cloud ice 
          !
          section="'VIQI'"
          IF (debug) PRINT *,'Section '//TRIM(section)//'... .. .'
          rcode = nf_inq_varndims(ncid, i, ndims)
          rcode = nf_inq_vardimid(ncid, i, ishape)
          DO idim=1, ndims
            rcode = nf_inq_dimlen(ncid, ishape(idim), dims_in(idim))
          END DO
          jshape2d=RESHAPE((/jshape(1),jshape(2),jshape(4)/),(/3/))
          dims2d_out=RESHAPE((/dims_in(1),dims_in(2),dims_in(4),1/),(/4/))
          jvar = jvar + 1
          varnameIN='VIQI'
          longdescIN='Vertically integrated cloud ice'
          unitsIN='Kg m-2'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          IF (ALLOCATED(data3)) DEALLOCATE(data3)
          ALLOCATE (data3(dims_in(1), dims_in(2), dims_in(4), 1))
          CALL massvertint(ncid, (/'QICE       '/), 1, debug, dims_in(1), dims_in(2), dims_in(3),&
            dims_in(4), data3)
          IF (debug) THEN
            PRINT *, 'VAR: VIQI idvar:',jvar
            PRINT *, '     DIMS OUT: ',dims2d_out
          ENDIF
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data3)
          IF (debug) PRINT *, '     SAMPLE VALUE OUT = ',data3(dims_out(1)/2,dims_out(2)/2,1,1)
          DEALLOCATE(data3)
        ENDIF

        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'CLT') /= 0 ) THEN
          !
          !   total cloud fraction
          !
          section="'CLT'"
          IF (debug) PRINT *,'Section '//TRIM(section)//'... .. .'
          rcode = nf_inq_varndims(ncid, i, ndims)
          rcode = nf_inq_vardimid(ncid, i, ishape)
          DO idim=1, ndims
            rcode = nf_inq_dimlen(ncid, ishape(idim), dims_in(idim))
          END DO
          jshape2d=RESHAPE((/jshape(1),jshape(2),jshape(4)/),(/3/))
          dims2d_out=RESHAPE((/dims_in(1),dims_in(2),dims_in(4),1/),(/4/))
          IF (ALLOCATED(data4)) DEALLOCATE(data4)
          ALLOCATE (data4(dims_in(1), dims_in(2), dims_in(4), 1, 1), STAT=ierr)
          IF (ierr /= 0) PRINT *,"Error in allocating 'data4'"
          varnameIN='CLDFRA'
          CALL clt(debug, ncid, dims_in(1), dims_in(2), dims_in(3), dims_in(4), varnameIN,     &
            data4) 
          IF (debug) THEN
            PRINT *,'total cloud fraction'
            PRINT *,'     VAR: CLT, idvar:',jvar+1
            PRINT *,'     DIMS OUT: ',dims2d_out
          ENDIF

          jvar = jvar + 1
          varnameIN='CLT'
          longdescIN='total cloud fraction'
          unitsIN='1'
          CALL def_var (mcid, jvar, varnameIN, 5, 3, jshape2d, "XY", longdescIN, unitsIN, "-"   &
            , "XLONG XLAT", output_netCDF4)
          rcode = nf_put_vara_real (mcid, jvar, start_dims, dims2d_out, data4(:,:,:,1,1))
          IF (debug) PRINT *, '  CLT. SAMPLE VALUE OUT = ',data4(dims_out(1)/2,dims_out(2)/2,1,1,1)
          DEALLOCATE(data4)
        ENDIF


!!! Close files on exit
        rcode = nf_close(ncid)
        rcode = nf_close(mcid)
        write(6,*) 
      ENDDO 

  
      write(6,*) 
      write(6,*) "##############################################################"
      write(6,*) " SUCCESS: We are out of here"      
      write(6,*) "##############################################################"

 END PROGRAM no_interp
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! GMS.UC: Lluis Dec.09
SUBROUTINE extract_from_geogrid(geoname, variables, Nvar, dx, dy, dz, geogrid_data)
! Subroutine for extracting data from geo_em.d01.nc
IMPLICIT NONE

INCLUDE 'netcdf.inc'

INTEGER                                                   :: i, ivar, rcode, ncgid
INTEGER, INTENT(IN)                                       :: Nvar, dx, dy, dz
CHARACTER(LEN=250), INTENT(IN)                            :: geoname
CHARACTER(LEN=11), DIMENSION(Nvar), INTENT(IN)            :: variables
REAL, DIMENSION(dx, dy, dz, Nvar), INTENT(OUT)            :: geogrid_data
REAL, DIMENSION(dx, dy, dz)                               :: datageo

!!!!!!!!!!! Variables
!    ncgid: unit number of 'geo_em.d01.nc' file
!    variables: vector with names of solicited variables
!    Nvar: number of solicited variables
!    dx, dy, dz: geogrid data dimension range
!    geogrid_data: 5D matrix with all solicited data

rcode = nf_open(geoname, 0, ncgid)

DO ivar=1,Nvar
  rcode = nf_inq_varid(ncgid, TRIM(variables(ivar)), i)
  rcode = nf_get_var_real ( ncgid, i, datageo )
  geogrid_data(1:dx,1:dy,1:dz,ivar) = datageo
ENDDO

rcode = nf_close ( ncgid )

END SUBROUTINE extract_from_geogrid

!---------------------------------------------------------------------

  SUBROUTINE clt(debg, ncid, dx, dy, dz, dt, cldfraname, totcloudfr)
!  Subroutine to compute total cloud cover in base 1.

  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'
!    567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567
  INTEGER, INTENT(IN)                                     :: dx, dy, dz, dt, ncid
  LOGICAL, INTENT(IN)                                     :: debg
  CHARACTER(LEN=11), INTENT(IN)                           :: cldfraname
  REAL, DIMENSION(dx,dy,dt,1,1), INTENT(OUT)              :: totcloudfr

! Local
  INTEGER                                                 :: i,j,k,it, ijk, nozero
  INTEGER                                                 :: rcode, idcldfra
  CHARACTER(LEN=50)                                       :: section
  CHARACTER(LEN=250)                                      :: message
  REAL, DIMENSION(dx,dy,dz,dt)                            :: cldfra, cldfram1, cldmax
  REAL, DIMENSION(dx,dy,dt)                               :: totcfr
  INTEGER                                                 :: halfdim
  
!!!!!!!!!!!!!! Variables
! dx, dy, dz, dt: dimensions of fields
! cldfra: cloud fraction at each level
! totcfr: total cloud fraction

  section="'subroutine_clt'"
  IF (debg) PRINT *,'Section '//TRIM(section)//'... .. .'
  
  IF (debg) PRINT *,'Dimensions: ',dx,CHAR(44), dy,CHAR(44), dz,CHAR(44), dt

  rcode = nf_inq_varid(ncid, cldfraname, idcldfra)
  rcode = nf_get_var_real(ncid, idcldfra, cldfra)

  totcfr = 1.
  
  ijk=0

  IF (debg)  PRINT *,'Computing total cloud fraction....'
  cldfram1(:,:,2:dz,:)=cldfra(:,:,1:dz-1,:)
  cldfram1(:,:,1,:)=0.
  
  cldmax=MAX(cldfram1,cldfra)
  
  WHERE (cldfram1 == 1.) cldfram1 = 0.

  vertical_levels: DO k=1, dz 
    totcfr=totcfr*((1.-cldmax(:,:,k,:))/(1.-cldfram1(:,:,k,:)))
  END DO vertical_levels
 
  totcfr=1.-totcfr
!  WHERE (totcfr > 1.) totcfr=1.

  IF (debg) THEN
    DO k=1,dz
      PRINT *,'dim/2 cloud fraction values at ',k, 'level: ', cldfra(halfdim(dx), halfdim(dy),  &
        k, halfdim(dt)),' cldframax: ',cldmax(halfdim(dx), halfdim(dy),k, halfdim(dt))
    END DO
  END IF

  IF (debg) PRINT *,'Total cloud fraction at the center dimx/2:', halfdim(dx),' dimy/2:',       &
    halfdim(dy),' dt/2:', halfdim(dt), ' =', totcfr(halfdim(dx),halfdim(dy),halfdim(dt))

  totcloudfr(:,:,:,1,1)=totcfr

  END SUBROUTINE clt


SUBROUTINE variablessum(ncid, variables, Nvar, dx, dy, dz, dt, datasum)
! Subroutine to add 'Nvar' variables
IMPLICIT NONE

INCLUDE 'netcdf.inc'

INTEGER                                                   :: i, ivar, rcode
INTEGER, INTENT(IN)                                       :: ncid, Nvar
CHARACTER(LEN=11), DIMENSION(Nvar), INTENT(IN)            :: variables
INTEGER, INTENT(IN)                                       :: dx, dy, dz, dt
REAL, DIMENSION(dx,dy,dz,dt), INTENT(OUT)                 :: datasum
REAL, DIMENSION(dx,dy,dz,dt)                              :: dataval

datasum=0.

DO ivar=1,Nvar
  rcode = nf_inq_varid(ncid, TRIM(variables(ivar)), i)
  rcode = nf_get_var_real ( ncid, i, dataval )
  datasum = datasum + dataval
ENDDO
END SUBROUTINE variablessum

SUBROUTINE multi_transportvertint(ncid, variables, Nvar, winds, dbg, nx, ny, nz, nt, integrals)
  ! Subroutine to vertically integrate a product of 'Nvar' variables
  ! in eta vertical coordinates
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER                                        :: i, iz, ivar, rcode
  INTEGER, INTENT(IN)                            :: ncid, Nvar
  CHARACTER(LEN=11), DIMENSION(Nvar), INTENT(IN) :: variables
  CHARACTER(LEN=11), DIMENSION(3), INTENT(IN)    :: winds
  INTEGER, INTENT(IN)                            :: nx, ny, nz, nt
  LOGICAL, INTENT(IN)                            :: dbg
  REAL, DIMENSION(nx,ny,nt,Nvar,3), INTENT(OUT)  :: integrals

! Local
  REAL, DIMENSION(nx,ny,nz,nt)                   :: dataval, windu, windv,   & 
    windw
  REAL, DIMENSION(nx,ny,nz,nt,3)                 :: integrand
  REAL, DIMENSION(nx,ny,nt)                      :: mu, mub
  REAL, DIMENSION(nz,nt)                         :: dz
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE          :: stagdata
  INTEGER, DIMENSION(6)                          :: vardims, lengthvardims
  INTEGER                                        :: idim, nvardims
  INTEGER                                        :: iu, iv, iw, idz, imu, imub
  REAL, PARAMETER                                :: grav=9.81
  REAL                                           :: dxm, dym
  CHARACTER(LEN=50)                              :: section

  section="'multi_transportvertint'"
  IF (dbg) PRINT *,'Vertical sigma integration of transport of ',Nvar,' variables...'
  IF (dbg) PRINT *,'Section '//TRIM(section)//'... .. .'

  integrals=0.
  rcode = nf_inq_varid(ncid, 'DNW', idz)
  rcode = nf_inq_varid(ncid, 'MU', imu)
  rcode = nf_inq_varid(ncid, 'MUB', imub)

  rcode = nf_get_var_real(ncid, idz, dz)
  rcode = nf_get_var_real(ncid, imu, mu)
  rcode = nf_get_var_real(ncid, imub, mub)
  rcode = nf_get_att_real(ncid, nf_global, 'DX', dxm)
  rcode = nf_get_att_real(ncid, nf_global, 'DY', dym)

  dz=ABS(dz)

! Adquiring winds
!!
  rcode = nf_inq_varid(ncid, TRIM(winds(1)), iu)
  IF (ALLOCATED(stagdata)) DEALLOCATE(stagdata)
  ALLOCATE (stagdata(nx+1,ny,nz,nt))
  rcode = nf_get_var_real(ncid, iu, stagdata)
  CALL error_nc(section, rcode)
  windu = 0.5 * ( stagdata(1:nx,:,:,:)+stagdata(2:nx+1,:,:,:) )

  rcode = nf_inq_varid(ncid, TRIM(winds(2)), iv)
  IF (ALLOCATED(stagdata)) DEALLOCATE(stagdata)
  ALLOCATE (stagdata(nx,ny+1,nz,nt))
  rcode = nf_get_var_real(ncid, iv, stagdata)
  CALL error_nc(section, rcode)
  windv = 0.5 * ( stagdata(:,1:ny,:,:)+stagdata(:,2:ny+1,:,:) )

  rcode = nf_inq_varid(ncid, TRIM(winds(3)), iw)
  CALL error_nc(section, rcode)
  IF (ALLOCATED(stagdata)) DEALLOCATE(stagdata)
  ALLOCATE (stagdata(nx,ny,nz+1,nt))
  rcode = nf_get_var_real(ncid, iw, stagdata)
  CALL error_nc(section, rcode)
  windw = 0.5 * ( stagdata(:,:,1:nz,:)+stagdata(:,:,2:nz+1,:) )

  DO ivar=1,Nvar
    integrand=1.
    IF (dbg) PRINT *,'  '//TRIM(variables(ivar))
    rcode = nf_inq_varid(ncid, TRIM(variables(ivar)), i)
    IF (rcode .ne. nf_noerr) call handle_err(rcode)
    rcode = nf_inq_varndims(ncid, i, nvardims)
    rcode = nf_inq_vardimid(ncid, i, vardims)
    DO idim=1, nvardims
      rcode = nf_inq_dimlen(ncid, vardims(idim), lengthvardims(idim))
    END DO
    rcode = nf_get_var_real(ncid, i, dataval)
    CALL error_nc(section, rcode)
    integrand(:,:,:,:,1) = integrand(:,:,:,:,1) * dataval * windu
    integrand(:,:,:,:,2) = integrand(:,:,:,:,2) * dataval * windv
    integrand(:,:,:,:,3) = integrand(:,:,:,:,3) * dataval * windw
    DO iz=1,nz
      integrals(:,:,:,ivar,:) = integrals(:,:,:,ivar,:) + integrand(:,:,iz,:,:)* dz(iz,1)
    ENDDO
!    integrals(:,:,:,ivar,1)=integrals(:,:,:,ivar,1)*(mu+mub)*dxm*dym/grav
!    integrals(:,:,:,ivar,2)=integrals(:,:,:,ivar,2)*(mu+mub)*dxm*dym/grav
!    integrals(:,:,:,ivar,3)=integrals(:,:,:,ivar,3)*(mu+mub)*dxm*dym/grav
    integrals(:,:,:,ivar,1)=integrals(:,:,:,ivar,1)*(mu+mub)/grav
    integrals(:,:,:,ivar,2)=integrals(:,:,:,ivar,2)*(mu+mub)/grav
    integrals(:,:,:,ivar,3)=integrals(:,:,:,ivar,3)*(mu+mub)/grav
  ENDDO

END SUBROUTINE multi_transportvertint

SUBROUTINE multi_massvertint(ncid, variables, Nvar, dbg, nx, ny, nz, nt, integrals)
  ! Subroutine to vertically integrate a product of 'Nvar' variables
  ! in eta vertical coordinates
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER                                        :: i, iz, ivar, rcode
  INTEGER, INTENT(IN)                            :: ncid, Nvar
  CHARACTER(LEN=11), DIMENSION(Nvar), INTENT(IN) :: variables
  INTEGER, INTENT(IN)                            :: nx, ny, nz, nt
  LOGICAL, INTENT(IN)                            :: dbg
  REAL, DIMENSION(nx,ny,nt,Nvar), INTENT(OUT)    :: integrals
  REAL, DIMENSION(nx,ny,nz,nt)                   :: dataval, integrand
  REAL, DIMENSION(nx,ny,nt)                      :: mu, mub
  REAL, DIMENSION(nz,nt)                         :: dz
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE          :: stagdata
  INTEGER, DIMENSION(6)                          :: vardims, lengthvardims
  INTEGER                                        :: idim, nvardims
  INTEGER                                        :: idz, imu, imub
  REAL, PARAMETER                                :: grav=9.81
  REAL                                           :: dxm, dym
  CHARACTER(LEN=50)                              :: section

  section="'multi_transportvertint'"
  IF (dbg) PRINT *,'Section '//TRIM(section)//'... .. .'

  IF (dbg) PRINT *,'Vertical sigma integration of ',Nvar,' variables...'
  integrals=0.
  rcode = nf_inq_varid(ncid, 'DNW', idz)
  CALL error_nc(section, rcode)
  rcode = nf_inq_varid(ncid, 'MU', imu)
  CALL error_nc(section, rcode)
  rcode = nf_inq_varid(ncid, 'MUB', imub)
  CALL error_nc(section, rcode)

  rcode = nf_get_var_real(ncid, idz, dz)
  CALL error_nc(section, rcode)
  rcode = nf_get_var_real(ncid, imu, mu)
  CALL error_nc(section, rcode)
  rcode = nf_get_var_real(ncid, imub, mub)
  CALL error_nc(section, rcode)
  rcode = nf_get_att_real(ncid, nf_global, 'DX', dxm)
  CALL error_nc(section, rcode)
  rcode = nf_get_att_real(ncid, nf_global, 'DY', dym)
  CALL error_nc(section, rcode)

  dz=ABS(dz)
      
  DO ivar=1,Nvar
    integrand=1.
    rcode = nf_inq_varid(ncid, TRIM(variables(ivar)), i)
    IF (dbg) PRINT *,"  '"//TRIM(variables(ivar))//"' id:", i
    CALL error_nc(section, rcode)
    if (rcode .ne. nf_noerr) call handle_err(rcode)
    IF (TRIM(variables(ivar))=='U') THEN
      IF (ALLOCATED(stagdata)) DEALLOCATE(stagdata)
      allocate (stagdata(nx+1,ny,nz,nt))
      rcode = nf_get_var_real(ncid, i, stagdata)
      CALL error_nc(section, rcode)
      dataval = 0.5 * ( stagdata(1:nx,:,:,:)+stagdata(2:nx+1,:,:,:) )
    ELSEIF (TRIM(variables(ivar))=='V') THEN
      IF (ALLOCATED(stagdata)) DEALLOCATE(stagdata)
      allocate (stagdata(nx,ny+1,nz,nt))
      rcode = nf_get_var_real(ncid, i, stagdata)
      CALL error_nc(section, rcode)
      dataval = 0.5 * ( stagdata(:,1:ny,:,:)+stagdata(:,2:ny+1,:,:) )
    ELSEIF (TRIM(variables(ivar))=='W') THEN
      IF (ALLOCATED(stagdata)) DEALLOCATE(stagdata)
      allocate (stagdata(nx,ny,nz+1,nt))
      rcode = nf_get_var_real(ncid, i, stagdata)
      CALL error_nc(section, rcode)
      dataval = 0.5 * ( stagdata(:,:,1:nz,:)+stagdata(:,:,2:nz+1,:) )
    ELSE
      rcode = nf_inq_varndims(ncid, i, nvardims)
      CALL error_nc(section, rcode)
      rcode = nf_inq_vardimid(ncid, i, vardims)
      CALL error_nc(section, rcode)
      DO idim=1, nvardims
        rcode = nf_inq_dimlen(ncid, vardims(idim), lengthvardims(idim))
        CALL error_nc(section, rcode)
      END DO
      rcode = nf_get_var_real(ncid, i, dataval)
      CALL error_nc(section, rcode)
    ENDIF
    integrand = integrand * dataval
    DO iz=1,nz
      integrals(:,:,:,ivar) = integrals(:,:,:,ivar) + integrand(:,:,iz,:)* dz(iz,1)
    ENDDO
!    integrals(:,:,:,ivar)=integrals(:,:,:,ivar)*(mu+mub)*dxm*dym/grav
    integrals(:,:,:,ivar)=integrals(:,:,:,ivar)*(mu+mub)/grav
  ENDDO

END SUBROUTINE multi_massvertint

SUBROUTINE massvertint(ncid, variables, Nvar, dbg, nx, ny, nz, nt, integral)
  ! Subroutine to vertically integrate a product of 'Nvar' variables
  ! in eta vertical coordinates
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER                                        :: i, iz, ivar, rcode
  INTEGER, INTENT(IN)                            :: ncid, Nvar
  CHARACTER(LEN=11), DIMENSION(Nvar), INTENT(IN) :: variables
  INTEGER, INTENT(IN)                            :: nx, ny, nz, nt
  LOGICAL, INTENT(IN)                            :: dbg
  REAL, DIMENSION(nx,ny,nt), INTENT(OUT)         :: integral
  REAL, DIMENSION(nx,ny,nz,nt)                   :: dataval, integrand
  REAL, DIMENSION(nx,ny,nt)                      :: mu, mub
  REAL, DIMENSION(nz,nt)                         :: dz
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE          :: stagdata
  INTEGER, DIMENSION(6)                          :: vardims, lengthvardims
  INTEGER                                        :: idim, nvardims
  INTEGER                                        :: idz, imu, imub
  REAL, PARAMETER                                :: grav=9.81
  REAL                                           :: sum=0.
  REAL                                           :: dxm, dym
  integrand=1.

  IF (dbg) PRINT *,'Vertical sigma integration of ',Nvar,' variables...'

  DO ivar=1,Nvar
    IF (dbg) PRINT *,ivar,' var: '//TRIM(variables(ivar))
    rcode = nf_inq_varid(ncid, TRIM(variables(ivar)), i)
    IF (rcode /= nf_noerr) PRINT *,nf_strerror(rcode)
    IF (TRIM(variables(ivar))=='U') THEN
      IF (ALLOCATED(stagdata)) DEALLOCATE(stagdata)
      allocate (stagdata(nx+1,ny,nz,nt))
      rcode = nf_get_var_real(ncid, i, stagdata)
      dataval = 0.5 * ( stagdata(1:nx,:,:,:)+stagdata(2:nx+1,:,:,:) )
    ELSEIF (TRIM(variables(ivar))=='V') THEN
      IF (ALLOCATED(stagdata)) DEALLOCATE(stagdata)
      allocate (stagdata(nx,ny+1,nz,nt))
      rcode = nf_get_var_real(ncid, i, stagdata)
      dataval = 0.5 * ( stagdata(:,1:ny,:,:)+stagdata(:,2:ny+1,:,:) )
    ELSE
      rcode = nf_inq_varndims(ncid, i, nvardims)
      rcode = nf_inq_vardimid(ncid, i, vardims)
      DO idim=1, nvardims
        rcode = nf_inq_dimlen(ncid, vardims(idim), lengthvardims(idim))
      END DO
      rcode = nf_get_var_real(ncid, i, dataval)
    ENDIF
    integrand = integrand * dataval
  ENDDO
  rcode = nf_inq_varid(ncid, 'DNW', idz)
  rcode = nf_inq_varid(ncid, 'MU', imu)
  rcode = nf_inq_varid(ncid, 'MUB', imub)

  rcode = nf_get_var_real(ncid, idz, dz)
  rcode = nf_get_var_real(ncid, imu, mu)
  rcode = nf_get_var_real(ncid, imub, mub)
  rcode = nf_get_att_real(ncid, nf_global, 'DX', dxm)
  rcode = nf_get_att_real(ncid, nf_global, 'DY', dym)
  dz=ABS(dz)

  integral=0.
  DO iz=1,nz
    integral = integral + integrand(:,:,iz,:)*dz(iz,1)
  ENDDO

!  integral=integral*(mu+mub)*dxm*dym/grav
  integral=integral*(mu+mub)/grav

END SUBROUTINE massvertint

 SUBROUTINE handle_err(rcode)
    INTEGER rcode
    write(6,*) 'Error number ',rcode
    stop
 END SUBROUTINE
!---------------------------------------------------------------------
 SUBROUTINE nointerp (data_out, data_in, eta_levels, ix, iy, iz, it, num_eta_levels)
! Subroutie to output a given eta-levels of a 4D field

     INTEGER                                          :: ix, iy, iz, it
     INTEGER                                          :: num_eta_levels
     REAL, DIMENSION(ix, iy, num_eta_levels, it)      :: data_out
     REAL, DIMENSION(ix, iy, iz, it)                  :: data_in
     INTEGER, DIMENSION(num_eta_levels)               :: eta_levels

     INTEGER                                          :: i, j, itt, k
!     data_out=data_in(:,:,eta_levels,:)
     do itt = 1, it
        do j = 1, iy
        do i = 1, ix
	  do k = 1, num_eta_levels
           data_out(i,j,k,itt)=data_in(i,j,eta_levels(k),itt)
	  end do
        end do
        end do
     end do

 END SUBROUTINE nointerp 
!------------------------------------------------------------------------------
! GMS.UC: Lluis
!------------------------------------------------------------------------------
SUBROUTINE spatialfiltering(datain,grid,ntimes,dx,dy,dz,dt,dataout)
! Subroutine to spatially filter a field via averaging values within boxes of 'grid' (odd) 
! number of grid points of side 'ntimes' times

IMPLICIT NONE

INTEGER                                                   :: i,j,k,it,igrid, jgrid, itime
INTEGER, INTENT(IN)                                       :: dx, dy, dz, dt, grid, ntimes
REAL, DIMENSION(dx,dy,dz,dt), INTENT(IN)                  :: datain
REAL, DIMENSION(dx,dy,dz,dt), INTENT(OUT)                 :: dataout
REAL, DIMENSION(dx,dy,dz,dt)                              :: datafilt

! Filling border values (up grid-1)
!!
DO i=0,grid-2
  dataout(1+i,:,:,:)=datain(1+i,:,:,:)
  dataout(dx-i,:,:,:)=datain(dx-i,:,:,:)
  dataout(:,1+i,:,:)=datain(:,1+i,:,:)
  dataout(:,dy-i,:,:)=datain(:,dy-i,:,:)
END DO

datafilt=datain

! Filtering
!!
DO itime=1,ntimes
  DO i=grid,dx-grid+1
    DO j=grid,dy-grid+1
      DO k=1,dz
        DO it=1,dt
          dataout(i,j,k,it)=0.
          DO igrid=-grid/2,grid/2
            DO jgrid=-grid/2, grid/2
              dataout(i,j,k,it)=dataout(i,j,k,it)+datafilt(i+igrid,j+jgrid,k,it)
            ENDDO
          ENDDO
          dataout(i,j,k,it)=dataout(i,j,k,it)/(grid*grid)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  datafilt=dataout
ENDDO

END SUBROUTINE spatialfiltering

!------------------------------------------------------------------------------
 SUBROUTINE mean_sealevelpress (data_out, pres_field, eta_levels, psfc, ter, tk, qv, ix, iy, iz, it, ito, &
        num_eta_levels)
! New subroutine to compute mean_sealevel pressure values 

     INTEGER, INTENT(IN)                                           :: ix, iy, iz, it, ito
     INTEGER, INTENT(IN)                                           :: num_eta_levels
     REAL, DIMENSION(ix, iy, ito, 1), INTENT(OUT)                  :: data_out
     REAL, DIMENSION(ix, iy, iz, ito), INTENT(IN)                  :: pres_field, tk, qv
     REAL, DIMENSION(ix, iy, ito), INTENT(IN)                      :: psfc
     REAL, DIMENSION(ix, iy), INTENT(IN)                           :: ter
     REAL, DIMENSION(num_eta_levels), INTENT(IN)                   :: eta_levels

     INTEGER                                                       :: i, j, itt, k, kk, kin

     N = 1.0

     expon=287.04*.0065/9.81

     ! Fill in missing values
     IF ( extrapolate == 0 ) RETURN       !! no extrapolation - we are out of here

     ! First find where about 400 hPa is located
     kk = 0
     find_kk : do k = 1, num_eta_levels
        kk = k
        if ( eta_levels(k) <= 40000. ) exit find_kk
     end do find_kk

     data_out=0.
     do itt = 1, ito
       do j = 1, iy
       do i = 1, ix

!                We are below both the ground and the lowest data level.

!                First, find the model level that is closest to a "target" pressure
!                level, where the "target" pressure is delta-p less that the local
!                value of a horizontally smoothed surface pressure field.  We use
!                delta-p = 150 hPa here. A standard lapse rate temperature profile
!                passing through the temperature at this model level will be used
!                to define the temperature profile below ground.  This is similar
!                to the Benjamin and Miller (1990) method, using  
!                700 hPa everywhere for the "target" pressure.

!         ptarget = (psfc(i,j,itt)*.01) - 150.
         ptarget = 700.
         dpmin=1.e4
         kupper = 0
         loop_kIN : do kin=iz,1,-1
           kupper = kin
           dp=abs( (pres_field(i,j,kin,itt)*.01) - ptarget )
           if (dp.gt.dpmin) exit loop_kIN
             dpmin=min(dpmin,dp)
           enddo loop_kIN
         ptarget=ptarget*100.
!         pbot=max(pres_field(i,j,1,itt),psfc(i,j,itt))
!         zbot=0.

!         tbotextrap=tk(i,j,kupper,itt)*(pbot/pres_field(i,j,kupper,itt))**expon
!         tvbotextrap=virtual(tbotextrap,qv(i,j,1,itt))

!         data_out(i,j,itt,1) = (zbot+tvbotextrap/.0065*(1.-(eta_levels(1)/pbot)**expon))
         tbotextrap=tk(i,j,kupper,itt)*(psfc(i,j,itt)/ptarget)**expon
         tvbotextrap=virtual(tbotextrap,qv(i,j,kupper,itt))
         data_out(i,j,itt,1) = psfc(i,j,itt)*((tvbotextrap+0.0065*ter(i,j))/tvbotextrap)**(1/expon)
!         IF (i==ix/2 .AND. j==iy/2 ) THEN
!         IF (ter(i,j) > 2500.) THEN
!           PRINT *,itt,' ptarget',ptarget,'kupper:',kupper
!           PRINT *,'tk:',tk(i,j,kupper,itt),'psfc:',psfc(i,j,itt)
!           PRINT *,'tbot:',tbotextrap,'tvbot:',tvbotextrap,'ter:',ter(i,j)
!           PRINT *,'qv:',qv(i,j,kupper,itt),'mslp:',data_out(i,j,itt,1)
!         ENDIF

       enddo
       enddo
     enddo

 END SUBROUTINE mean_sealevelpress 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! New subroutine with ito (from dims_out(4)) as a supplementary value, since it was not correctly obtained as 
! dims_in(4)
 SUBROUTINE nointerpdiag (data_out, data_in, eta_levels, ix, iy, iz, it, ito, &
                     num_eta_levels)

     INTEGER, INTENT(IN)                                           :: ix, iy, iz, it, ito
     INTEGER, INTENT(IN)                                           :: num_eta_levels
     REAL, DIMENSION(ix, iy, num_eta_levels, ito), INTENT(OUT)     :: data_out
     REAL, DIMENSION(ix, iy, iz, ito), INTENT(IN)                  :: data_in
     INTEGER, DIMENSION(num_eta_levels), INTENT(IN)                :: eta_levels

     INTEGER                                                       :: i, j, itt, k

!     data_out=data_in(:,:,eta_levels,:)
     do itt = 1, ito
        do j = 1, iy
        do i = 1, ix
          do k=1, num_eta_levels
           data_out(i,j,k,itt)=data_in(i,j,eta_levels(k),itt)
	  end do
        end do
        end do
     end do

 END SUBROUTINE nointerpdiag

!------------------------------------------------------------------------------
 FUNCTION virtual (tmp,rmix)
!      This function returns virtual temperature in K, given temperature
!      in K and mixing ratio in kg/kg.

     real                              :: tmp, rmix, virtual

     virtual=tmp*(0.622+rmix)/(0.622*(1.+rmix))

 END FUNCTION virtual

!------------------------------------------------------------------------------
 SUBROUTINE all_spaces ( command , length_of_char ) 

      IMPLICIT NONE

      INTEGER                                       :: length_of_char
      CHARACTER (LEN=length_of_char)                :: command
      INTEGER                                       :: loop

      DO loop = 1 , length_of_char
         command(loop:loop) = ' '
      END DO

 END SUBROUTINE all_spaces

!------------------------------------------------------------------------------
 SUBROUTINE def_var (mcid, jvar, cval, itype, idm, jshape, order, desc, units, stag, coord, output_netCDF4 )
      IMPLICIT NONE

      INCLUDE 'netcdf.inc'

      INTEGER                                         :: mcid, jvar
! GMS.UC: Lluis Dec.09      
      CHARACTER (LEN =  11)                           :: cval
      INTEGER                                         :: itype, idm
      REAL, DIMENSION(6)                              :: jshape
      CHARACTER (*)                                   :: order
! GMS.UC: Lluis Dec.09      
      CHARACTER (LEN = 60)                            :: desc
      CHARACTER (LEN = 10)                            :: units
      CHARACTER (LEN =  1)                            :: stag
      CHARACTER (LEN = 10)                            :: coord

      INTEGER                                         :: rcode, ilen
! GMS.UC: Lluis Dec.09
      CHARACTER (LEN=60)                              :: att_text
! GMS.UC Markel March 2012
      LOGICAL                                         :: output_netCDF4

      IF ( itype == 5 ) THEN
         rcode = nf_redef(mcid)
         rcode = nf_def_var(mcid, trim(cval), NF_REAL, idm, jshape, jvar)
         IF (output_netCDF4) THEN
           rcode = nf_def_var_deflate(mcid, jvar, 1, 1, 3)
           IF (rcode  /= nf_noerr) THEN
             PRINT *, "ERROR SETTING THE DEFLATE LEVEL"
             PRINT *, rcode
           ENDIF
         ENDIF
         rcode = nf_put_att_int(mcid, jvar, "FieldType", NF_INT, 1, 104)
      ENDIF

      att_text = order
      ilen = len_trim(att_text)
      rcode = nf_put_att_text(mcid, jvar, "MemoryOrder", ilen, att_text(1:ilen) )
      att_text = ' '

      att_text = desc
      ilen = len_trim(att_text)
      rcode = nf_put_att_text(mcid, jvar, "description", ilen, att_text(1:ilen) )
      att_text = ' '

      att_text = units
      ilen = len_trim(att_text)
      rcode = nf_put_att_text(mcid, jvar, "units", ilen, att_text(1:ilen) )
      att_text = ''

      att_text = stag
      ilen = len_trim(att_text)
      rcode = nf_put_att_text(mcid, jvar, "stagger", ilen, att_text(1:ilen) )
      att_text = ''

      att_text = coord
      ilen = len_trim(att_text)
      rcode = nf_put_att_text(mcid, jvar, "coordinates", ilen, att_text(1:ilen) )
      att_text = ''

      rcode = nf_enddef(mcid)

 END SUBROUTINE def_var

INTEGER FUNCTION halfdim(dim)
! Function to give the half value of a dimension

  IMPLICIT NONE
  
  INTEGER, INTENT(IN)                                     :: dim
  
  halfdim = dim/2
  IF (dim < 2) halfdim = 1     ! Assuming non-zero dimension range
  
END FUNCTION halfdim

SUBROUTINE error_nc(sec, rc)
! Subroutine to print error nc messages

  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'
  
  INTEGER, INTENT(IN)                                     :: rc
  CHARACTER(LEN=50), INTENT(IN)                           :: sec

! Local
  CHARACTER(LEN=50)                                       :: errmsg='ERROR - error - ERROR'//   &
    ' - error - ERROR'

  
  IF (rc /= 0) PRINT *,TRIM(errmsg)//" in "//TRIM(sec)//" "//nf_strerror(rc)

END SUBROUTINE error_nc


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
