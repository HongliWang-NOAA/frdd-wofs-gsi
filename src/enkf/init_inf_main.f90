program init_inf_main
!$$$  main program documentation block
!
! program:  init_inf_main             high level driver program for 
!                                      prior adaptive inflation
!
! prgmmr: jjh          org: CIMMS/OU (NSSL)               date: 2017-03-16
!
! abstract: This is the main program for prior adaptive inflation. It does the following:
!           a) initialize MPI, read namelist from enkf.nml on each task.
!           b) read horizontal grid information (lat/lon of each grid point) and
!              pressure at each grid point/vertical level. 
!              to each task (one task one member)
!           c) initialize inflation file 
!
! program history log:
!   2017-03-16  Initial version.
!
! attributes:
!   language: f90
!
!$$$
 use kinds, only: r_kind,r_double,i_kind,r_single
 ! reads namelist parameters.
 use params, only : prior_inf_initial,prior_inf_sd_initial,read_namelist,nvars,prior_inf,prior_inf_lower_bound,prior_inf_upper_bound,prior_inf_sd_lower_bound,prior_inf_file,prior_inf_sd_file
 ! grid information
 use gridinfo, only: getgridinfo, gridinfo_cleanup, gridvarstring
 ! mpi functions and variables.
 use mpisetup, only:  mpi_initialize, mpi_cleanup, nproc
 
 use init_inf, only: create_netcdf_wrf

 implicit none
 real(r_double) t1,t2
 real(r_single) inf_init,inf_sd_init

 ! initialize MPI.
 call mpi_initialize()
 
 ! read namelist.
 call read_namelist()

 ! read horizontal grid information and pressure fields from
 ! 6-h forecast ensemble mean file.
 call getgridinfo()
 
 if(nproc .eq. 0) then
     inf_init = prior_inf_initial
     inf_sd_init = prior_inf_sd_initial
     if(prior_inf == 1) then
        if(prior_inf_initial < prior_inf_lower_bound) inf_init=prior_inf_lower_bound
        if(prior_inf_initial > prior_inf_upper_bound) inf_init=prior_inf_upper_bound
        if(prior_inf_sd_initial < prior_inf_sd_lower_bound) inf_sd_init=prior_inf_sd_lower_bound
        
        call create_netcdf_wrf(prior_inf_file, inf_init)
        call create_netcdf_wrf(prior_inf_sd_file, inf_sd_init)
     end if
 end if

 call gridinfo_cleanup()
 call mpi_cleanup()

end program init_inf_main
