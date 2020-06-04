subroutine read_goesimg_netcdf(val_img,ithin,rmesh,jsatid,gstime,&
     infile,lunout,obstype,nread,ndata,nodata,twind,sis,nobs,&
     mype_root,mype_sub,npe_sub,mpi_comm_sub, dval_use)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_goesimg                    read goes imager data
!   prgmmr: su, xiujuan      org: np23                date: 2002-02-26
!
! abstract:  This routine reads GOES imager radiance (brightness
!            temperature) files.  Optionally, the data are thinned to
!            a specified resolution using simple quality control checks.
!
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! program history log:
!   2002-02-26 su, x.  
!   2004-05-28 kleist  - update subroutine call
!   2004-06-16 treadon - update documentation
!   2004-07-23 derber - make changes to eliminate obs. earlier in thinning
!   2004-07-29  treadon - abonl  to module use, add intent in/out
!   2005-01-26  derber - land/sea determination and weighting for data selection
!   2005-09-08  derber - modify to use input group time window
!   2005-09-28  derber - modify to produce consistent surface info
!   2005-10-17  treadon - add grid and earth relative obs location to output file
!   2005-10-18  treadon - remove array obs_load and call to sumload
!   2005-11-29  parrish - modify getsfc to work for different regional options
!   2006-02-01  parrish - remove getsfc (different version called now in read_obs)
!   2006-02-03  derber  - add new obs control
!   2006-04-27  derber - clean up code
!   2006-05-19  eliu    - add logic to reset relative weight when all channels not used
!   2006-06-19  kleist - correct bug in global grid relative dlat,dlon
!   2006-07-28  derber  - add solar and satellite azimuth angles remove isflg from output
!   2007-03-01  tremolet - measure time from beginning of assimilation window
!   2008-10-14  derber - allow mpi_io
!   2009-04-21  derber  - add ithin to call to makegrids
!   2011-04-08  li      - (1) use nst_gsi, nstinfo, fac_dtl, fac_tsl and add NSST vars
!                         (2) get zob, tz_tr (call skindepth and cal_tztr)
!                         (3) interpolate NSST Variables to Obs. location (call deter_nst)
!                         (4) add more elements (nstinfo) in data array
!   2011-08-01  lueken  - added module use deter_sfc_mod, fix indentation
!   2012-03-05  akella  - nst now controlled via coupler
!   2013-01-26  parrish - change from grdcrd to grdcrd1 (to allow successful debug compile on WCOSS)
!   2016-09-20  Thomas Jones - Replace BUFR satellite data with custom netcdf data for GOES IMAGER
!   2017-12-12  Thomas Jones - Updated for GSI 3.6
!
!   input argument list:
!     mype     - mpi task id
!     val_img  - weighting factor applied to super obs
!     ithin    - flag to thin data
!     rmesh    - thinning mesh size (km)
!     jsatid   - satellite to read
!     gstime   - analysis time in minutes from reference date
!     infile   - unit from which to read BUFR data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!
!   output argument list:
!     nread    - number of BUFR GOES imager observations read
!     ndata    - number of BUFR GOES imager profiles retained for further processing
!     nodata   - number of BUFR GOES imager observations retained for further processing
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,r_double,i_kind
  use satthin, only: super_val,itxmax,makegrids,map2tgrid,destroygrids, &
      checkob,finalcheck,score_crit
  use gridmod, only: diagnostic_reg,regional,nlat,nlon,txy2ll,tll2xy,rlats,rlons
  use constants, only: deg2rad,zero,one,rad2deg,r60inv,r60
  use radinfo, only: iuse_rad,jpch_rad,nusis
  use gsi_4dvar, only: l4dvar,iwinbgn,winlen
  use deter_sfc_mod, only: deter_sfc
  use gsi_nstcouplermod, only: gsi_nstcoupler_skindepth, gsi_nstcoupler_deter, nst_gsi,nstinfo
  use mpimod, only: npe
  use netcdf

  implicit none

  include 'netcdf.inc'
  
! Declare passed variables
  character(len=*),intent(in   ) :: infile,obstype,jsatid
  character(len=20),intent(in  ) :: sis
  integer(i_kind) ,intent(in   ) :: lunout,ithin !, mype
  integer(i_kind) ,intent(inout) :: ndata,nodata
  integer(i_kind) ,intent(inout) :: nread
  real(r_kind)    ,intent(in   ) :: rmesh,gstime,twind
  real(r_kind)    ,intent(inout) :: val_img
  integer(i_kind) ,dimension(npe),intent(inout) :: nobs
  !integer(i_kind) ,dimension(1),intent(inout) :: nobs
  integer(i_kind) ,intent(in   ) :: mype_root
  integer(i_kind) ,intent(in   ) :: mype_sub
  integer(i_kind) ,intent(in   ) :: npe_sub
  integer(i_kind) ,intent(in   ) :: mpi_comm_sub
  logical         ,intent(in   ) :: dval_use

! Declare local parameters
  integer(i_kind),parameter:: nimghdr=13
  !integer(i_kind),parameter:: maxinfo=37
  real(r_kind),parameter:: r360=360.0_r_kind
  real(r_kind),parameter:: tbmin=50.0_r_kind
  real(r_kind),parameter:: tbmax=550.0_r_kind
  character(80),parameter:: hdrgoes  = &            ! goes imager header
        'SAID YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SAZA SOZA BEARAZ SOLAZI'

! Declare local variables
  logical outside,iuse,assim

  character(8) subset

  integer(i_kind) nchanl,ilath,ilonh,ilzah,iszah,irec,next
  integer(i_kind) nmind,lnbufr,idate,ilat,ilon, sthin, maxinfo
  integer(i_kind) ireadmg,ireadsb,iret,nreal,nele,itt
  integer(i_kind) itx,i,k,isflg,kidsat,n,iscan,idomsfc, ii
  integer(i_kind) idate5(5)
  integer(i_kind),allocatable,dimension(:)::nrec

  real(r_kind) dg2ew,sstime,tdiff,t4dv,sfcr
  real(r_kind) dlon,dlat,timedif,crit1,dist1
  real(r_kind) dlon_earth,dlat_earth, thislon, thislat
  real(r_kind) pred, boxcwp3,  boxcwp5
  real(r_kind),dimension(0:4):: rlndsea
  real(r_kind),dimension(0:3):: sfcpct
  real(r_kind),dimension(0:3):: ts
  real(r_kind) :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10
  real(r_kind) :: zob,tref,dtw,dtc,tz_tr
  real(r_kind),allocatable,dimension(:,:):: data_all

  real(r_double),dimension(nimghdr) :: hdrgoesarr       !  goes imager header
  real(r_double),dimension(3,6) :: dataimg              !  goes imager data

  real(r_kind) cdist,disterr,disterrmax,dlon00,dlat00
  integer(i_kind) ntest

!------------------
!  NETCDF-RELATED
!------------------
   INTEGER(i_kind)   :: ncdfID
   INTEGER(i_kind)   :: status
   INTEGER(i_kind)   :: datestlen, varID, DimID
   INTEGER(i_kind)   :: vardim, natts
   INTEGER(i_kind)   :: vartype, id_time
   INTEGER(i_kind)   :: numdim, numvars, numatt
   INTEGER(i_kind)   :: xx,yy,mx,my

!------------------
!  SATELLITE DATA
!------------------
   REAL(r_kind), ALLOCATABLE, DIMENSION( :, : )       :: nir16
   REAL(r_kind), ALLOCATABLE, DIMENSION( :, : )       :: sir39
   REAL(r_kind), ALLOCATABLE, DIMENSION( :, : )       :: wv68
   REAL(r_kind), ALLOCATABLE, DIMENSION( :, : )       :: ir108
   REAL(r_kind), ALLOCATABLE, DIMENSION( :, : )       :: ir133
   REAL(r_kind), ALLOCATABLE, DIMENSION( :, : )       :: cwp
   INTEGER(i_kind), ALLOCATABLE, DIMENSION( :, : )    :: phase 
   REAL(r_kind), ALLOCATABLE, DIMENSION( :, : )       :: vza
   REAL(r_kind), ALLOCATABLE, DIMENSION( :, : )       :: sza
   REAL(r_kind), ALLOCATABLE, DIMENSION( :, : )       :: lon
   REAL(r_kind), ALLOCATABLE, DIMENSION( :, : )       :: lat

!**************************************************************************
! Initialize variables
  sthin = 2 !THIN OBSERVATIONS:  1 = no thin, 2 = every other grid point, etc

  maxinfo=35
  lnbufr = 10
  disterrmax=zero
  ntest=0
  dg2ew = r360*deg2rad

  ilon=3
  ilat=4

  if (nst_gsi > 0 ) then
     call gsi_nstcoupler_skindepth(obstype, zob)         ! get penetration depth (zob) for the obstype
  endif

  rlndsea(0) = zero
  rlndsea(1) = 15._r_kind
  rlndsea(2) = 10._r_kind
  rlndsea(3) = 15._r_kind
  rlndsea(4) = 30._r_kind

  ndata=0
  nodata=0
  nchanl=4       ! the channel number
  ilath=8        ! the position of latitude in the header
  ilonh=9        ! the position of longitude in the header
  ilzah=10       ! satellite zenith angle
  iszah=11       ! solar zenith angle

! If all channels of a given sensor are set to monitor or not
! assimilate mode (iuse_rad<1), reset relative weight to zero.
! We do not want such observations affecting the relative
! weighting between observations within a given thinning group.

  assim=.false.
  search: do i=1,jpch_rad
     if ((nusis(i)==sis) .and. (iuse_rad(i)>0)) then
        assim=.true.
        exit search
     endif
  end do search
  if (.not.assim) val_img=zero


! Make thinning grids
  call makegrids(rmesh,ithin)


! OPEN NETCDF FILE
status = nf90_open(TRIM(infile), NF90_NOWRITE, ncdfID)
print*, '*** OPENING GOES RADIANCE NETCDF FILE', status
!------------------------
! Get global attributes
!-------------------------
status = nf90_get_att( ncdfID, nf90_global, 'sat_year', idate5(1) )
status = nf90_get_att( ncdfID, nf90_global, 'sat_month', idate5(2) )
status = nf90_get_att( ncdfID, nf90_global, 'sat_day', idate5(3) )
status = nf90_get_att( ncdfID, nf90_global, 'sat_hour', idate5(4) )
status = nf90_get_att( ncdfID, nf90_global, 'sat_minute', idate5(5) )

print*, idate5

!------------------------
! Get Dimension Info (2-D)
!-------------------------
status = NF_INQ_DIMID(ncdfID, 'Lon', varID)
status = NF_INQ_DIMLEN(ncdfID, varID, mx)

status = NF_INQ_DIMID(ncdfID, 'Lat', varID)
status = NF_INQ_DIMLEN(ncdfID, varID, my)

!------------------------
! Allocate data arrays
!-------------------------
ALLOCATE( lat( mx,my ) )
ALLOCATE( lon( mx,my ) )
ALLOCATE( sza( mx,my ) )
ALLOCATE( vza( mx,my ) )
ALLOCATE( cwp( mx,my ) )
ALLOCATE( phase( mx,my ) )
!ALLOCATE( nir16( mx,my ) )
ALLOCATE( sir39( mx,my ) )
ALLOCATE( wv68( mx,my ) )
ALLOCATE( ir108( mx,my ) )
ALLOCATE( ir133( mx,my ) )

!------------------------
! Get useful data arrays
!-------------------------
! LAT
status = nf90_inq_varid( ncdfID, 'Lat2D', varID )
status = nf90_get_var( ncdfID, varID, lat )

! LON
status = nf90_inq_varid( ncdfID, 'Lon2D', varID )
status = nf90_get_var( ncdfID, varID, lon )

! VZA
status = nf90_inq_varid( ncdfID, 'vza', varID )
status = nf90_get_var( ncdfID, varID, vza )

! SZA
status = nf90_inq_varid( ncdfID, 'sza', varID )
status = nf90_get_var( ncdfID, varID, sza )

! CWP 
status = nf90_inq_varid( ncdfID, 'cwp', varID )
status = nf90_get_var( ncdfID, varID, cwp )

! Cloud Phase (use to determine clear-sky pixels)
status = nf90_inq_varid( ncdfID, 'phase', varID )
status = nf90_get_var( ncdfID, varID, phase )

! NEAR IR 1.6 um
!status = nf90_inq_varid( ncdfID, 'mrms_nir16', varID )
!status = nf90_get_var( ncdfID, varID, nir16 )

! SOLAR IR 3.9 um (2)
status = nf90_inq_varid( ncdfID, 'mrms_sir39', varID )
status = nf90_get_var( ncdfID, varID, sir39 )

! WATER VAPOR 6.8 um (3)
status = nf90_inq_varid( ncdfID, 'mrms_wv68', varID )
status = nf90_get_var( ncdfID, varID, wv68 )

! IR 10.8 um (4)
status = nf90_inq_varid( ncdfID, 'mrms_ir108', varID )
status = nf90_get_var( ncdfID, varID, ir108 )

! IR 13.3 um (5 or 6)
status = nf90_inq_varid( ncdfID, 'mrms_ir133', varID )
status = nf90_get_var( ncdfID, varID, ir133 )

! CLOSE NETCDF FILE
status = nf90_close( ncdfID )


! SAT ID
  if(jsatid == 'g08') kidsat = 252
  if(jsatid == 'g09') kidsat = 253
  if(jsatid == 'g10') kidsat = 254
  if(jsatid == 'g11') kidsat = 255
  if(jsatid == 'g12') kidsat = 256
  if(jsatid == 'g13') kidsat = 257
  if(jsatid == 'g14') kidsat = 258
  if(jsatid == 'g15') kidsat = 259
!  if(jsatid == 'g16') kidsat = 2XX

! Allocate arrays to hold all data for given satellite
  if(dval_use) maxinfo = maxinfo + 2
  nreal = maxinfo + nstinfo
  nele  = nreal   + nchanl

  print*, nreal, nele, itxmax

  allocate(data_all(nele,itxmax),nrec(itxmax))

  next=0
  nrec=999999
  irec=0
  itx=1
  ii=1
  iuse = .true.

! Big loop over data array
  do xx=3, mx-3, sthin
    do yy=3, my-3, sthin

     irec=irec+1
     next=next+1
     
        nread=nread+nchanl
        
        !print*, xx,yy,irec, wv68(xx,yy), cwp(xx,yy), phase(xx,yy)
        
!       Convert obs location from degrees to radians
        thislon=lon(xx,yy)
        thislat=lat(xx,yy)

        if (thislon >= r360) thislon=thislon-r360
        if (thislon < zero)  thislon=thislon+r360

        dlon_earth=thislon*deg2rad
        dlat_earth=thislat*deg2rad

!       If regional, map obs lat,lon to rotated grid.
        if(regional)then

!          Convert to rotated coordinate.  dlon centered on 180 (pi), 
!          so always positive for limited area
           call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)

           if(diagnostic_reg) then
              call txy2ll(dlon,dlat,dlon00,dlat00)
              ntest=ntest+1
              cdist=sin(dlat_earth)*sin(dlat00)+cos(dlat_earth)*cos(dlat00)* &
                   (sin(dlon_earth)*sin(dlon00)+cos(dlon_earth)*cos(dlon00))
              cdist=max(-one,min(cdist,one))
              disterr=acos(cdist)*rad2deg
              disterrmax=max(disterrmax,disterr)
           end if

!          Check to see if in domain.  outside=.true. if dlon_earth,
!          dlat_earth outside domain, =.false. if inside
!           if(outside) cycle read_loop

!       Global case
        else
           dlon=dlon_earth
           dlat=dlat_earth
           call grdcrd1(dlat,rlats,nlat,1)
           call grdcrd1(dlon,rlons,nlon,1)
        endif

        !ONLY CONTINUE IF OBS IN DOMAIN
        if (outside .eqv. .false.) then

        if (l4dvar) then
           crit1=0.01_r_kind
        else
           timedif = 6.0_r_kind*abs(tdiff)        ! range:  0 to 18
           crit1=0.01_r_kind+timedif
        endif
        
!        call map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,sis)
!       if(.not. iuse)cycle read_loop


!       Locate the observation on the analysis grid.  Get sst and land/sea/ice
!       mask.  
!     isflg    - surface flag
!                0 sea      LARC: 1 or 6
!                1 land     LARC: 4
!                2 sea ice  LARC: 2 or 7
!                3 snow     LARC: 0
!                4 mixed                 
       ! SET DEFAULT TO LAND
       isflg = 1        
       if ( phase(xx,yy).eq.1 .or. phase(xx,yy).eq.6 ) isflg = 0
       if ( phase(xx,yy).eq.4 ) isflg = 1 
       if ( phase(xx,yy).eq.2 .or. phase(xx,yy).eq.7 ) isflg = 2
       if ( phase(xx,yy).eq.0 ) isflg = 3


      ! TEMPORARY...set isflg to 1 (all land)
      !isflg = 1
        call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg,idomsfc,sfcpct, &
            ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)

!       Set common predictor parameters
        crit1=crit1+rlndsea(isflg)
!        call checkob(dist1,crit1,itx,iuse)
!        if(.not. iuse)cycle read_loop

!       Set data quality predictor 
!        pred =(10.0_r_kind-dataimg(2,1)/10.0_r_kind)+dataimg(3,3)*10.0_r_kind  ! clear sky and
                                                                 ! bt std as quality indicater
                                                                                                                    
!       Compute "score" for observation.  All scores>=0.0.  Lowest score is "best"
!        crit1 = crit1+pred 
!        call finalcheck(dist1,crit1,itx,iuse)
!        if(.not. iuse) cycle read_loop

!       Map obs to grids
!        iscan = nint(hdrgoesarr(ilzah))+1.001_r_kind ! integer scan position
        
!
!       interpolate NSST variables to Obs. location and get dtw, dtc, tz_tr
        if ( nst_gsi > 0 ) then
           tref  = ts(0)
           dtw   = zero
           dtc   = zero
           tz_tr = one
           if ( sfcpct(0) > zero ) then
              call gsi_nstcoupler_deter(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr)
           endif
        endif

!       Calculate CWP for 3x3 or 5x5 box surrounding observation. Define obs as cloudy if sum(cwp) > 0
        boxcwp3 = zero
        boxcwp5 = zero

        boxcwp3=cwp(xx-1,yy-1)+cwp(xx,yy-1)+cwp(xx+1,yy-1)+cwp(xx-1,yy)+cwp(xx,yy)+cwp(xx+1,yy)+cwp(xx-1,yy+1)+cwp(xx,yy+1)+cwp(xx+1,yy+1)
        boxcwp5=cwp(xx-2,yy-2)+cwp(xx-1,yy-2)+cwp(xx,yy-2)+cwp(xx+1,yy-2)+cwp(xx+2,yy-2)+ &
                    cwp(xx-2,yy-1)+cwp(xx-1,yy-1)+cwp(xx,yy-1)+cwp(xx+1,yy-1)+cwp(xx+2,yy-1)+ &
                    cwp(xx-2,yy)+cwp(xx-1,yy)+cwp(xx,yy)+cwp(xx+1,yy)+cwp(xx+2,yy)+ &
                    cwp(xx-2,yy+1)+cwp(xx-1,yy+1)+cwp(xx,yy+1)+cwp(xx+1,yy+1)+cwp(xx+2,yy+1)+ &
                    cwp(xx-2,yy+2)+cwp(xx-1,yy+2)+cwp(xx,yy+2)+cwp(xx+1,yy+2)+cwp(xx+2,yy+2)

!       Only Assimilate clear sky pixels 
        if ( cwp(xx,yy) >= zero .and. boxcwp3 >= 0 .and. boxcwp3 < 0.0001 .and. .not.outside .and. iuse) then
  
           !print*, xx,yy, dlat, dlon, irec, wv68(xx,yy), cwp(xx,yy), sza(xx,yy), boxcwp5
 
           !if (wv68(xx,yy).gt.0.1) print*, xx,yy, dlat, dlon, irec, wv68(xx,yy), cwp(xx,yy), phase(xx,yy), sza(xx,yy), outside, iuse, sfcpct(:)
           !print*, xx, yy, sfcpct(:)        
 
!       Transfer information to work array
        data_all( 1,itx) = kidsat                    ! satellite id
        data_all( 2,itx) = t4dv                       ! analysis relative time
        data_all( 3,itx) = dlon                       ! grid relative longitude
        data_all( 4,itx) = dlat                       ! grid relative latitude
        data_all( 5,itx) = vza(xx,yy)*deg2rad         ! satellite zenith angle (radians)
        data_all( 6,itx) = zero             ! satellite azimuth angle (radians)
        data_all( 7,itx) = phase(xx,yy)               ! clear sky amount
        data_all( 8,itx) = iscan                      ! integer scan position
        data_all( 9,itx) = sza(xx,yy)                 ! solar zenith angle
        data_all(10,itx) = zero             ! solar azimuth angle
        data_all(11,itx) = sfcpct(0)                  ! sea percentage of
        data_all(12,itx) = sfcpct(1)                  ! land percentage
        data_all(13,itx) = sfcpct(2)                  ! sea ice percentage
        data_all(14,itx) = sfcpct(3)                  ! snow percentage
        data_all(15,itx)= ts(0)                       ! ocean skin temperature
        data_all(16,itx)= ts(1)                       ! land skin temperature
        data_all(17,itx)= ts(2)                       ! ice skin temperature
        data_all(18,itx)= ts(3)                       ! snow skin temperature
        data_all(19,itx)= tsavg                       ! average skin temperature
        data_all(20,itx)= vty                         ! vegetation type
        data_all(21,itx)= vfr                         ! vegetation fraction
        data_all(22,itx)= sty                         ! soil type
        data_all(23,itx)= stp                         ! soil temperature
        data_all(24,itx)= sm                          ! soil moisture
        data_all(25,itx)= sn                          ! snow depth
        data_all(26,itx)= zz                          ! surface height
        data_all(27,itx)= idomsfc + 0.001_r_kind      ! dominate surface type
        data_all(28,itx)= sfcr                        ! surface roughness
        data_all(29,itx)= ff10                        ! ten meter wind factor
        data_all(30,itx)= dlon_earth*rad2deg          ! earth relative longitude (degrees)
        data_all(31,itx)= dlat_earth*rad2deg          ! earth relative latitude (degrees)

        if(dval_use)then
           data_all(36,itx) = val_img
           data_all(37,itx) = itt
        end if

        if ( nst_gsi > 0 ) then
           data_all(maxinfo+1,itx) = tref         ! foundation temperature
           data_all(maxinfo+2,itx) = dtw          ! dt_warm at zob
           data_all(maxinfo+3,itx) = dtc          ! dt_cool at zob
           data_all(maxinfo+4,itx) = tz_tr        ! d(Tz)/d(Tr)
        endif

!       Transfer observation location and other data to local arrays
!        do k=1,nchanl
!           data_all(k+31,itx) = 2.0
           !data_all(k+31,ii)=dataimg(3,k+1)		!STD DEVS
           !data_all(k+nreal,itx)=dataimg(1,k+1)     !TBs
!        end do
        
        !Transfer Observations to local array
        data_all(1+nreal,itx) = sir39(xx,yy)
        data_all(2+nreal,itx) = wv68(xx,yy)
        data_all(3+nreal,itx) = ir108(xx,yy)
        data_all(4+nreal,itx) = ir133(xx,yy)

        score_crit(itx) = 1.0
        
        nrec(itx)=irec
        ndata = ndata+1
        itx = itx+1 
        ii = ii+1

      endif
     endif !domain check
    enddo  !yy loop
  enddo  !xx loop 

!print*, 'READ NUM CHECK', npe, ii, mype_sub, mype_root, nreal, ndata, nst_gsi
!  call combine_radobs(mype_sub,mype_root,npe_sub,mpi_comm_sub,&
!     nele,itxmax,nread,ndata,data_all,score_crit,nrec)

! Write final set of "best" observations to output file
! Allow single task to check for bad obs, update superobs sum,
! and write out data to scratch file for further processing.
  if (mype_sub==mype_root.and.ndata>0) then

   call count_obs(ndata,nele,ilat,ilon,data_all,nobs)
   print*, ndata,nreal,ilat,ilon, nobs, k
   write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
   write(lunout) ((data_all(k,n),k=1,nele),n=1,ndata)

    !do ii=1, ndata
    !  print*, npe, ii, data_all(2+nreal,ii), data_all(11,ii), data_all(12,ii), data_all(13,ii), data_all(14,ii)
    !enddo
  endif

! Deallocate local arrays
  deallocate(data_all,nrec)

  DEALLOCATE(lat)
  DEALLOCATE(lon)
  DEALLOCATE(vza)
  DEALLOCATE(sza)
  DEALLOCATE(cwp)
  DEALLOCATE(phase)
  !DEALLOCATE(nir16)
  DEALLOCATE(sir39)
  DEALLOCATE(wv68)
  DEALLOCATE(ir108)
  DEALLOCATE(ir133)
  
   print*, '*** FINISHED READING RADIANCE NETCDF FILE'
   
! Deallocate satthin arrays
900 continue
  call destroygrids

!  if(diagnostic_reg.and.ntest>0) write(6,*)'READ_GOESIMG_NETCDF:  ',&
!     'mype,ntest,disterrmax=',mype,ntest,disterrmax

  return
end subroutine read_goesimg_netcdf
