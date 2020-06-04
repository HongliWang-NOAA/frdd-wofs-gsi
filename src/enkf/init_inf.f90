module init_inf
!$$$  main program documentation block
!
! program: jjh              org: CIMMS/OU (NSSL)               date: 2017-03-27
!
! create an inflation file (mean and sd) according to the namelist settings
!
!

  use gridinfo, only: dimensions, gridvarstring
  use kinds,    only: r_double, r_kind, r_single
  use netcdf_io
  use params,   only: nvars,datapath, arw, nbackgrounds,fgfileprefixes
  use constants, only: zero,one,cp,fv,rd,grav,zero

  implicit none
  
  ! Define all public subroutines within this module

  private
  public :: create_netcdf_wrf


 contains

 ! create a netcdf file in wrf format and initial each variable at every grid and every lev as initv
 subroutine create_netcdf_wrf(filein,initv)
    use constants
    include 'netcdf.inc'

    character(len=200),           intent(in)    :: filein
    real(r_single),                 intent(in):: initv

    ! define some variables 
    character(len=500)                                                               :: filename
    character(len=3)                                                                 :: charnbkg
    integer               ::ncid, fstatus, dimidt, dimidt_len, dimidx, dimidy, dimidz,dimidx_stag, dimidy_stag, dimidz_stag, varidx, varidy, varidz
    character(len=20)                       ::dimstr
    character(len=50)                       :: attstr    
    character(len=12)                       :: varmemoryorder,varstagger ! from gridio, the last variable is always MU
    integer     :: fid
    integer     :: l,slen
    integer, dimension(nvars+1)     :: ndims,varid,xlen,ylen,zlen
    integer  :: dims(4) 
    real, allocatable   ::vardata4(:,:,:,:)
    real, allocatable   ::vardata3(:,:,:)
    integer  :: start(4),count(4)
    integer     :: i,j,nb

    ! Define array dimension variables
    integer                                                                                    :: xdim
    integer                                                                                    :: ydim
    integer                                                                                    :: zdim

    ! Initialize constants required by routine
    print *,'entering into create_netcdf_wrf!'
    call init_constants(.true.)

  do nb=1,nbackgrounds ! loop over background times
 
    if (nb .lt. 10) write(charnbkg,'(i1.0)') nb
    if (nb .ge. 10 .and. nb .lt. 100) write(charnbkg,'(i2.0)') nb
    if (nb .ge. 100) write(charnbkg,'(i3.0)') nb
    filename = trim(adjustl(datapath))//trim(adjustl(filein))//'.'//charnbkg   
    ! create a file 
  
    fstatus = nf_create(filename,NF_CLOBBER,ncid)
    if(fstatus .ne. nf_noerr) then
       print *,'open file:',filename,'error!'
       call stop2(1001)
    end if

    ! define the dimentions

    xdim = dimensions%xdim
    ydim = dimensions%ydim
    zdim = dimensions%zdim
    
    dimstr = 'west_east'
    fstatus = nf_def_dim(ncid,dimstr,xdim,dimidx)
    if(fstatus .ne. 0) then
      print *,'error defining dimension ', dimstr
      call stop2(1002)
    endif
    dimstr = 'south_north'
    fstatus = nf_def_dim(ncid,dimstr,ydim,dimidy)
    if(fstatus .ne. 0) then
      print *,'error defining dimension ', dimstr
      call stop2(1002)
    endif
    dimstr = 'bottom_top'
    fstatus = nf_def_dim(ncid,dimstr,zdim,dimidz)
    if(fstatus .ne. 0) then
      print *,'error defining dimension ', dimstr
      call stop2(1002)
    endif
    dimstr = 'west_east_stag'
    fstatus = nf_def_dim(ncid,dimstr,xdim+1,dimidx_stag)
    if(fstatus .ne. 0) then
      print *,'error defining dimension ', dimstr
      call stop2(1002)
    endif
    dimstr = 'south_north_stag'
    fstatus = nf_def_dim(ncid,dimstr,ydim+1,dimidy_stag)
    if(fstatus .ne. 0) then
      print *,'error defining dimension ', dimstr
      call stop2(1002)
    endif
    dimstr = 'bottom_top_stag'
    fstatus = nf_def_dim(ncid,dimstr,zdim+1,dimidz_stag)
    if(fstatus .ne. 0) then
      print *,'error defining dimension ', dimstr
      call stop2(1002)
    endif
    dimstr = 'Time'
    fstatus = nf_def_dim(ncid,dimstr,1,dimidt)
    if(fstatus .ne. 0) then
      print *,'error defining dimension ', dimstr
      call stop2(1002)
    endif
    dimstr = 'DateStrLen'
    fstatus = nf_def_dim(ncid,dimstr,19,dimidt_len)
    if(fstatus .ne. 0) then
      print *,'error defining dimension ', dimstr
      call stop2(1002)
    endif
    
    ! define the variables
    ! according to memble 001 to define the variables and their attributes
    filename = trim(adjustl(datapath))//trim(adjustl(fgfileprefixes(nb)))//'1'
    fstatus = nf_open(filename,nf_nowrite,fid)
    if(fstatus .ne. 0) then
      print *,'error opening file ', filename
      call stop2(1001)
    endif
    ! Loop through all analysis variables to be defined

    do l = 1, nvars + 1
          ! inquiry the info for each variable
          ! gridvarstring(l) variable name

          ! Memeory Order (3D or 2D)
          attstr = 'MemoryOrder'
          call variableattribute_char(filename,gridvarstring(l),attstr,     &
               & varmemoryorder)

          ! If variable is a 2-dimensional field, rescale variables
          ! appropriately
          xlen(l) = xdim
          ylen(l) = ydim
          zlen(l) = zdim
          if(varmemoryorder(1:3) .eq. 'XY ') then
             ! Rescale grid dimension variables appropriately
             ndims(l) = 3
             zlen(l) = 1  ! not necessarily set this value, but for constistency, set it to 1 for 2D variable
             dims(1) = dimidx
             dims(2) = dimidy
             dims(3) = dimidt
          else
             ndims(l) = 4
             dims(1) = dimidx
             dims(2) = dimidy
             dims(3) = dimidz
             dims(4) = dimidt
          end if          ! Memeory Order (3D or 2D)
          
          ! stagger or not
          attstr = 'stagger'
          call variableattribute_char(filename,gridvarstring(l),attstr,     &
               & varstagger)
          
          if(varstagger(1:1) .eq. 'X') then
             ! If variable grid is staggered in X-direction, assign array
             ! Assign array dimensions appropriately
             dims(1) = dimidx_stag
             xlen(l) = xlen(l)+1
             slen = 1
          else if(varstagger(1:1) .eq. 'Y') then
             ! If variable grid is staggered in Y-direction, assign array
             ! dimensions appropriately
             dims(2) = dimidy_stag
             ylen(l) = ylen(l)+1
             slen = 1
          else if(varstagger(1:1) .eq. 'Z') then
             ! If variable grid is staggered in Z-direction, assign array
             ! dimensions appropriately
             ! of course, 3D only
             dims(3) = dimidz_stag
             zlen(l) = zlen(l)+1
             slen = 1
          else
             slen = 0
          end if
         
          fstatus = nf_def_var(ncid,gridvarstring(l),NF_REAL,ndims(l),dims(1:ndims(l)),varid(l))
          fstatus = nf_put_att_text(ncid,varid(l),'MemoryOrder',3,varmemoryorder(1:3))
          if(slen .eq. 1) then
              fstatus = nf_put_att_text(ncid,varid(l),'stagger',1,varstagger(1:1))
          else
              fstatus = nf_put_att_text(ncid,varid(l),'stagger',0,'')
          end if
   end do
   
   ! close the ref file
   fstatus = nf_close(fid)

   !Change mode of netCDF operation
   fstatus = nf_enddef(ncid)
   
   ! Output the values of the variables
   do l = 1, nvars + 1
         if(ndims(l) == 3) then
             start(1) = 1
             start(2) = 1
             start(3) = 1
             count(1) = xlen(l)
             count(2) = ylen(l)
             count(3) = 1
             allocate(vardata3(count(1),count(2),count(3)))
             vardata3 = initv
             ! 2D variable
             fstatus = nf_put_vara_real(ncid,varid(l),start(1:3),count(1:3),vardata3)
             if(allocated(vardata3))   deallocate(vardata3)
         else
             start(1) = 1
             start(2) = 1
             start(3) = 1
             start(4) = 1
             count(1) = xlen(l)
             count(2) = ylen(l) 
             count(3) = zlen(l) 
             count(4) = 1
             allocate(vardata4(count(1),count(2),count(3),count(4)))
             vardata4 = initv
             ! 3D variable
             fstatus = nf_put_vara_real(ncid,varid(l),start(1:4),count(1:4),vardata4)
             if(allocated(vardata4))   deallocate(vardata4)
         end if
       
   end do 
   ! close the new file
   fstatus = nf_close(ncid)
 enddo ! end loop over background times
   print *,'complete create_netcdf_wrf'
   end subroutine create_netcdf_wrf 

end module init_inf
