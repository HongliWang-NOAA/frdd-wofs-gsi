module gridio_inf

  !========================================================================

  !$$$ Module documentation block
  ! 
  ! This module contains routines to read and write values in the format  
  ! from Weather Research and Forecasting (WRF) model Advanced Research
  ! WRF (ARW) 
  !
  ! Currently used for ensemble mean, prior inflation mean and sd 
  !
  ! prgmmr: jjh             org: CIMMS/OU (NSSL)       date: 2017-03-16
  !
  ! program history log:
  !   
  !   2017-03-16  Initial version.
  !
  ! attributes:
  !   language:  f90
  !
  !$$$

  !=========================================================================

  ! Define associated modules

  use gridinfo, only: dimensions, gridvarstring, npts, cross2dot, dot2cross
  use kinds,    only: r_double, r_kind, r_single
  use mpisetup, only: nproc
  use netcdf_io
  use params,   only: nlevs, nvars, nlons, nlats, cliptracers, datapath,prior_inf_file,prior_inf_sd_file,     &
       &              arw, datestring, nbackgrounds, &
                      prior_inf_lower_bound,prior_inf_upper_bound,prior_inf_sd_lower_bound !lower and upper bounds for adaptive inflation
  use constants, only: zero,one,cp,fv,rd,grav,zero

  implicit none

  !-------------------------------------------------------------------------

  ! Define all public subroutines within this module

  private
  public :: readgriddata_f
  public :: writegriddata_f

  !-------------------------------------------------------------------------

contains

  !========================================================================

  ! readgriddata_f.f90: This subroutine will receive an netcdf (inf, inf_sd)
  ! file prefix (WRF_ARW structure) and will subsequently return
  ! the respective variable interpolated to an unstaggered grid; all
  ! checks for grid staggering are contained within this subroutine
  ! only reads those varariables that will be updated in DA
  !-------------------------------------------------------------------------

  subroutine readgriddata_f(filein,vargrid)

    use constants

    !======================================================================

    ! Define array dimension variables

    integer                                                                      :: xdim, ydim, zdim

    ! Define variables passed to subroutine

    character(len=200),                                               intent(in)  :: filein

    ! Define variables returned by subroutine

    real(r_single),       dimension(npts,nvars*nlevs+1,nbackgrounds),       intent(out) :: vargrid

    ! Define variables computed within subroutine

    character(len=500)                                                   :: filename
    character(len=3)                                                     :: charnbkg
    real,       dimension(:,:,:),               allocatable              :: workgrid
    real,       dimension(:,:,:),               allocatable              :: vargrid_native
    integer                                                                      :: xdim_native
    integer                                                                      :: ydim_native
    integer                                                                      :: zdim_native
    integer                                                                      :: xdim_local
    integer                                                                      :: ydim_local
    integer                                                                      :: zdim_local

    ! Define variables requiredfor netcdf variable I/O

    character(len=12)                                                            :: varstrname
    character(len=50)                                                            :: attstr
    character(len=12)                                                            :: varstagger
    character(len=12)                                                            :: varmemoryorder

    ! Define counting variables

    integer                                                                      :: i, j, k, l, nb
    integer                                                                      :: counth, countv
    integer                                                                      :: count
    
    ! Define logic variable
    logical                                                                      :: fexist 
    !======================================================================

    ! Initialize all constants required by routine

    call init_constants(.true.)

    ! Define all local variables

    xdim = dimensions%xdim
    ydim = dimensions%ydim
    zdim = dimensions%zdim
       
    !======================================================================

    ! Begin: Loop through each (prognostic) variable (defined in
    ! gridio.F90), determine and define the spatial array
    ! dimensions, and allocate memory for ARW dynamical core

    !----------------------------------------------------------------------
    if (nbackgrounds > 1) then
       write(6,*)'gridio_inf/readgriddata_f: reading multiple prior_inf not yet supported'
       call stop2(23)
    endif
    backgroundloop: do nb=1,nbackgrounds

    ! Initialize counting variable

    countv = 1

    ! the whole path of the filename
    if (nb .lt. 10) write(charnbkg,'(i1.0)') nb
    if (nb .ge. 10 .and. nb .lt. 100) write(charnbkg,'(i2.0)') nb
    if (nb .ge. 100) write(charnbkg,'(i3.0)') nb
    filename = trim(adjustl(datapath))//trim(adjustl(filein))//'.'//charnbkg

    ! check if the file exist or not
    
    inquire( file= filename, exist=fexist )
    if ( .not. fexist ) then
        write(6,*)  "readgriddata_f: file "//filename//" does not exist!"
        return
    end if
    !----------------------------------------------------------------------

    ! Loop through all variables to be update via the EnKF

    do l = 1, nvars + 1

    !----------------------------------------------------------------------

       ! Define staggering attributes for variable grid

       attstr = 'stagger'
       call variableattribute_char(filename,gridvarstring(l),attstr,        &
            & varstagger)

       ! If variable grid is staggered in X-direction, assign array
       ! dimensions appropriately

       if(varstagger(1:1) .eq. 'X') then

          ! Assign array dimensions appropriately

          xdim_native = xdim + 1
          ydim_native = ydim
          zdim_native = zdim

       ! If variable grid is staggered in Y-direction, assign array
       ! dimensions appropriately

       else if(varstagger(1:1) .eq. 'Y') then

          ! Assign array dimensions appropriately

          xdim_native = xdim
          ydim_native = ydim + 1
          zdim_native = zdim

       ! If variable grid is staggered in Z-direction, assign array
       ! dimensions appropriately

       else if(varstagger(1:1) .eq. 'Z') then 

          ! Assign array dimensions appropriately

          xdim_native = xdim
          ydim_native = ydim
          zdim_native = zdim + 1
             
       ! If variable grid is not staggered, assign array dimensions
       ! appropriately

       else

          ! Assign array dimensions appropriately

          xdim_native = xdim
          ydim_native = ydim
          zdim_native = zdim

       end if 

    !----------------------------------------------------------------------

       ! Define memory attributes for variable grid

       attstr = 'MemoryOrder'
       call variableattribute_char(filename,gridvarstring(l),attstr,       &
            & varmemoryorder)

          ! If variable is a 2-dimensional field, rescale variables
          ! appropriately

          if(varmemoryorder(1:3) .eq. 'XY ') then

             ! Rescale grid dimension variables appropriately

             zdim_local = 1
             zdim_native = 1
       
          else

             ! Define local array dimension

             zdim_local = zdim

          end if 

          ! Define local variable dimensions

          xdim_local = xdim
          ydim_local = ydim

          ! Allocate memory for local variable arrays

          if(.not. allocated(workgrid))                                     &
               & allocate(workgrid(xdim_local,ydim_local,zdim_local))
          if(.not. allocated(vargrid_native))                               &
               & allocate(vargrid_native(xdim_native,ydim_native,           &
               & zdim_native))

          ! Ingest variable from external netcdf formatted file
       
          call readnetcdfdata(filename,vargrid_native,gridvarstring(l),     &
               & xdim_native,ydim_native,zdim_native)
       
          ! Interpolate variable from staggered (i.e., C-) grid to
          ! unstaggered (i.e., A-) grid. If variable is staggered in
          ! vertical, intepolate from model layer interfaces
          ! (including surface and top) to model layer midpoints.

          call cross2dot(vargrid_native,xdim_native,ydim_native,            &
               & zdim_native,xdim_local,ydim_local,zdim_local,workgrid)

    !----------------------------------------------------------------------

          ! Loop through vertical coordinate

          do k = 1, zdim_local

             ! Initialize counting variable
          
             counth = 1

             ! Loop through meridional horizontal coordinate
             
             do j = 1, ydim_local

                ! Loop through zonal horizontal coordinate

                do i = 1, xdim_local

                   ! Assign values to output variable array
             
                   vargrid(counth,countv,nb) = workgrid(i,j,k)

                   ! Update counting variable

                   counth = counth + 1

                end do ! do i = 1, xdim_local

             end do ! do j = 1, ydim_local

             ! Print message to user

             if (nproc .eq. 0)                                               &
                  write(6,*) 'READGRIDDATA_F: ', trim(gridvarstring(l)),   &
                  & countv, minval(vargrid(:,countv,nb)),                       &
                  & maxval(vargrid(:,countv,nb))

             ! Update counting variable

             countv = countv + 1

          end do ! do k = 1, zdim_local

    !----------------------------------------------------------------------

          ! Deallocate memory for local variables

          if(allocated(vargrid_native)) deallocate(vargrid_native)
          if(allocated(workgrid))       deallocate(workgrid)

    !----------------------------------------------------------------------

       end do ! do l = 1, nvars + 1

    end do backgroundloop ! loop over backgrounds to read in

    return

    !======================================================================

  end subroutine readgriddata_f

  ! writegriddata_f.f90: This subroutine will receive a netcdf 
  ! file prefix and will subsequently return the respective
  ! variable interpolated to the native variable grid; all checks for
  ! grid staggering are contained within this subroutine
  ! read that file first and only write the updated values
  !-------------------------------------------------------------------------

  subroutine writegriddata_f(filein,vargrid)

    use constants
    include 'netcdf.inc'      

    !----------------------------------------------------------------------

    ! Define variables passed to subroutine

    real(r_single),    dimension(npts,nvars*nlevs+1,nbackgrounds),                            intent(in)    :: vargrid
    character(len=200),                                                           intent(in)    :: filein

    !----------------------------------------------------------------------

    ! Define variables computed within subroutine

    character(len=500)                                                               :: filename
    character(len=3)                                                                 :: charnbkg
    real,    dimension(:,:,:),                allocatable                            :: vargrid_native
    real,    dimension(:,:,:),                allocatable                            :: vargridin_native
    real,    dimension(:,:,:),                allocatable                            :: workgrid
    real                                                                             :: clip
    integer iyear,imonth,iday,ihour,iminute,dh1,ierr,iw3jdn
    integer                                                                                    :: xdim_native
    integer                                                                                    :: ydim_native
    integer                                                                                    :: zdim_native
    integer                                                                                    :: xdim_local
    integer                                                                                    :: ydim_local
    integer                                                                                    :: zdim_local

    !----------------------------------------------------------------------

    ! Define array dimension variables

    integer                                                                                    :: xdim
    integer                                                                                    :: ydim
    integer                                                                                    :: zdim

    !----------------------------------------------------------------------

    ! Define variables required by for extracting netcdf variable
    ! fields

    character(len=50)                                                                          :: attstr
    character(len=12)                                                                          :: varstagger,varstrname
    character(len=12)                                                                          :: varmemoryorder
    character(len=19)  :: DateStr

    !----------------------------------------------------------------------

    ! Define counting variables

    integer                                                                                    :: i, j, k, l, nb
    integer                                                                                    :: counth, countv
    
    ! Define a logical variable

    logical                                                                                    :: fexist
    !----------------------------------------------------------------------

    ! Initialize constants required by routine

    call init_constants(.true.)

    !----------------------------------------------------------------------

    ! Define all array dimensions

    xdim = dimensions%xdim
    ydim = dimensions%ydim
    zdim = dimensions%zdim

    if (nbackgrounds > 1) then
       write(6,*)'gridio/writegriddata: writing multiple backgrounds not yet supported'
       call stop2(23)
    endif

    backgroundloop: do nb=1,nbackgrounds

    ! Allocate memory for local variable

    allocate(workgrid(xdim,ydim,zdim))

    !----------------------------------------------------------------------

    ! End: Define all local variables required by routine

    !======================================================================

    ! Begin: Loop through each prognostic variable and determine the
    ! spatial array dimensions for each variable contained within
    ! file, define appropriate array dimensions, and allocate memory;
    ! update respective analysis (e.g., prognostic model) variables

    !----------------------------------------------------------------------

    ! Initialize counting variable

    countv = 1

    !----------------------------------------------------------------------

    ! First guess file should be copied to analysis file at scripting
    ! level; only variables updated by EnKF are changed
    if (nb .lt. 10) write(charnbkg,'(i1.0)') nb
    if (nb .ge. 10 .and. nb .lt. 100) write(charnbkg,'(i2.0)') nb
    if (nb .ge. 100) write(charnbkg,'(i3.0)') nb
    filename = trim(adjustl(datapath))//trim(adjustl(filein))//'.'//charnbkg

    !----------------------------------------------------------------------

    ! Loop through all analysis variables to be updated

    do l = 1, nvars + 1

    !----------------------------------------------------------------------

       ! For WRF-ARW; analysis variables are defined on C-grid; the
       ! check for interpolation between mass and velocity points is
       ! done here

       if(arw) then

    !----------------------------------------------------------------------

          ! Define staggering attributes for variable grid
       
          attstr = 'stagger'
          call variableattribute_char(filename,gridvarstring(l),attstr,     &
               & varstagger)

    !----------------------------------------------------------------------

          ! If variable grid is staggered in X-direction, assign array
          ! dimensions appropriately

          if(varstagger(1:1) .eq. 'X') then

             ! Assign array dimensions appropriately

             xdim_native = xdim + 1
             ydim_native = ydim
             zdim_native = zdim

    !----------------------------------------------------------------------

             ! If variable grid is staggered in Y-direction, assign
             ! array dimensions appropriately

          else if(varstagger(1:1) .eq. 'Y') then

             ! Assign array dimensions appropriately

             xdim_native = xdim
             ydim_native = ydim + 1
             zdim_native = zdim

    !----------------------------------------------------------------------

             ! If variable grid is staggered in Z-direction, assign
             ! array dimensions appropriately

          else if(varstagger(1:1) .eq. 'Z') then

             ! Assign array dimensions appropriately

             xdim_native = xdim
             ydim_native = ydim
             zdim_native = zdim + 1

    !----------------------------------------------------------------------

             ! If variable grid is not staggered, assign array
             ! dimensions appropriately

          else 

             ! Assign array dimensions appropriately

             xdim_native = xdim
             ydim_native = ydim
             zdim_native = zdim

    !----------------------------------------------------------------------

          end if ! if(varstagger(1:1) .eq. 'X')

    !----------------------------------------------------------------------

       endif ! if(arw)

    !----------------------------------------------------------------------

       ! Define memory attributes for variable grid; this is done for
       ! ARW only

       if(arw) then

          attstr = 'MemoryOrder'
          call variableattribute_char(filename,gridvarstring(l),attstr,     &
               & varmemoryorder)

       end if ! if(arw)

    !----------------------------------------------------------------------

       ! If variable is a 2-dimensional field, rescale variables
       ! appropriately

       if(gridvarstring(l) .eq. 'MU' .or. gridvarstring(l) .eq. 'PD') then

          ! Rescale grid dimension variables appropriately

          zdim_local = 1
          zdim_native = 1
       
       else

          ! Define local array dimension

          zdim_local = zdim

       end if ! if(gridvarstring(l) .eq. 'MU' .or. gridvarstring(l) .eq. 
              ! 'PD')

    !----------------------------------------------------------------------

       ! Define local variable dimensions

       xdim_local = xdim
       ydim_local = ydim

    !----------------------------------------------------------------------

       ! Allocate memory local arrays (first check whether they are
       ! already allocated)

       if (allocated(vargrid_native)) deallocate(vargrid_native)
       allocate(vargrid_native(xdim_native,ydim_native,zdim_native))
       if (allocated(vargridin_native)) deallocate(vargridin_native)
       allocate(vargridin_native(xdim_native,ydim_native,zdim_native))

    !----------------------------------------------------------------------
       
       ! Read in prior values (i.e., analysis without current
       ! increments) and store in local array
       inquire( file= filename, exist=fexist )
       if ( .not. fexist ) then
          write(6,*)  "writegriddata_f: file "//filename//" does not exist!"
          return
       end if
         call readnetcdfdata(filename,vargridin_native,gridvarstring(l),      &
              & xdim_native,ydim_native,zdim_native) ! reads in old data

    !----------------------------------------------------------------------

       ! Loop through vertical coordinate

       do k = 1, zdim_local

    !----------------------------------------------------------------------

          ! Initialize counting variable

          counth = 1

    !----------------------------------------------------------------------

          ! Loop through meridional horizontal coordinate
          
          do j = 1, ydim
             
             ! Loop through zonal horizontal coordinate

             do i = 1, xdim

    !----------------------------------------------------------------------

                ! Assign values to local array

                workgrid(i,j,k) = vargrid(counth,countv,nb)

                ! Update counting variable

                counth = counth + 1

    !----------------------------------------------------------------------

             end do ! do i = 1, xdim

          end do ! do j = 1, ydim

    !----------------------------------------------------------------------

          ! Update counting variable

          countv = countv + 1

    !----------------------------------------------------------------------

       end do ! k = 1, zdim_local

    !----------------------------------------------------------------------

       ! Interpolate update value to native grid (i.e., from A-grid to
       ! C-grid; if necessary); on input, workgrid is updates on
       ! unstaggered grid; on output vargrid_native is updates on
       ! model-native (i.e., staggered grid)

       call dot2cross(xdim_local,ydim_local,zdim_local,xdim_native,          &
            ydim_native,zdim_native,workgrid,vargrid_native)

       if(arw) then

!          if (varstagger(1:1) .eq. 'Z') then ! if 'W' or 'PH' don't update surface

!             vargridin_native(:,:,2:zdim_native) = vargrid_native(:,:,2:zdim_native) 

!          else

             vargridin_native =  vargrid_native

!          endif ! if (varstagger(1:1) .eq. 'Z')

       endif ! if(arw)

       ! Clip all tracers (assume names start with 'Q')

!       if (cliptracers .and. gridvarstring(l)(1:1) .eq. 'Q') then

!          clip = tiny(vargridin_native(1,1,1))
!          where (vargridin_native < clip) vargridin_native = clip

!       end if ! if (cliptracers .and. gridvarstring(l)(1:1) .eq. 'Q')

    !----------------------------------------------------------------------

       ! Write updated inflation values.

      ! set the lower and upper bounds for prior inflation mean and sd
       if(filein==prior_inf_file) then
           where (vargridin_native < prior_inf_lower_bound ) vargridin_native = prior_inf_lower_bound
           where (vargridin_native > prior_inf_upper_bound ) vargridin_native = prior_inf_upper_bound
       end if

       if(filein==prior_inf_sd_file) then
           where (vargridin_native < prior_inf_sd_lower_bound ) vargridin_native = prior_inf_sd_lower_bound
       end if

       call writenetcdfdata(filename,vargridin_native,gridvarstring(l),       &
             xdim_native,ydim_native,zdim_native)

    end do ! do l = 1, nviars+1

    !----------------------------------------------------------------------

    ! Deallocate memory for local variable

    deallocate(workgrid)

    ! update NSTART_HOUR in NMM (HWRF) restart file.
    read(datestring(1:4),'(i4)') iyear
    read(datestring(5:6),'(i2)') imonth
    read(datestring(7:8),'(i2)') iday
    read(datestring(9:10),'(i2)') ihour
    read(datestring(11:12),'(i2)') iminute

    !  update START_DATE, SIMULATION_START_DATE, GMT, JULYR, JULDAY 
    !  global attributes.
    !
    write(DateStr,'(i4,"-",i2.2,"-",i2.2,"-",i2.2,"_",i2.2,":",i2.2)') iyear,imonth,iday,ihour,iminute,0
    ierr = NF_OPEN(trim(filename), NF_WRITE, dh1)
    IF (ierr .NE. NF_NOERR) print *, 'OPEN ',NF_STRERROR(ierr)
    ierr = NF_PUT_ATT_TEXT(dh1,NF_GLOBAL,'START_DATE',len(trim(DateStr)),DateStr)
    IF (ierr .NE. NF_NOERR) print *,'PUT START_DATE', NF_STRERROR(ierr)
    ierr = NF_PUT_ATT_TEXT(dh1,NF_GLOBAL,'SIMULATION_START_DATE',len(trim(DateStr)),DateStr)
    IF (ierr .NE. NF_NOERR) print *,'PUT SIMULATION_START_DATE', NF_STRERROR(ierr)
    ierr = NF_PUT_ATT_REAL(dh1,NF_GLOBAL,'GMT',NF_FLOAT,1,float(ihour))
    IF (ierr .NE. NF_NOERR) print *,'PUT GMT', NF_STRERROR(ierr)
    ierr = NF_PUT_ATT_INT(dh1,NF_GLOBAL,'JULYR',NF_INT,1,iyear)
    IF (ierr .NE. NF_NOERR) print *,'PUT JULYR', NF_STRERROR(ierr)
    ierr=NF_PUT_ATT_INT(dh1,NF_GLOBAL,'JULDAY',NF_INT,1,iw3jdn(iyear,imonth,iday)-iw3jdn(iyear,1,1)+1)
    IF (ierr .NE. NF_NOERR) print *,'PUT JULDAY', NF_STRERROR(ierr)
    ierr = NF_CLOSE(dh1)
    IF (ierr .NE. NF_NOERR) print *, 'CLOSE ',NF_STRERROR(ierr)

    !----------------------------------------------------------------------

    ! End: Loop through each prognostic variable and determine the
    ! spatial array dimensions for each variable contained within
    ! file, define appropriate array dimensions, and allocate memory;
    ! update respective analysis (e.g., prognostic model) variables

    !======================================================================
    end do backgroundloop ! loop over backgrounds to read in

    ! Return calculated values

    return

    !======================================================================

  end subroutine writegriddata_f

  !========================================================================


end module gridio_inf
