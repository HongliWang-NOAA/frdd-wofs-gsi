program prior_inf_main
!$$$  main program documentation block
!
! program:  prior_inf_main             high level driver program for 
!                                      prior adaptive inflation
!
! prgmmr: jjh          org: CIMMS/OU (NSSL)               date: 2017-03-16
!
! abstract: This is the main program for prior adaptive inflation. It does the following:
!           a) initialize MPI, read namelist from enkf.nml on each task.
!           b) read horizontal grid information (lat/lon of each grid point) and
!              pressure at each grid point/vertical level. 
!           c) read state variables for each ensemble member, distribute 
!              to each task (one task one member)
!           d) read inflation value, ensemble mean on each task
!           e) inflate the state variables
!           f) write out inflated analysis ensembls
!           i) deallocate all allocatable arrays, finalize MPI.
!
! program history log:
!   2017-03-16  Initial version.
!
! usage:
!   input files:
!     prior_inf.nc  - prior inflation file
!     firstguess.ensmean        - ensemble mean before running GSI forward operator
!     firstguess.mem*  - ensemble members before running GSI forward operator
!
!   output files: 
!     analysis.mem*  - updated (inflated) ensemble members, can call GSI forward operator afterwards
!     comments: the ensemble mean file should be recomputed before calling GSI forward operator afterwards                    

! This program is run before the forward operator code. The forward operator is run on each ensemble
! member to create the diag*mem* input files.
!
! attributes:
!   language: f90
!
!$$$
 use kinds, only: r_kind,r_double,i_kind
 ! reads namelist parameters.
 use params, only : read_namelist,prior_inf,nanals
 ! mpi functions and variables.
 use mpisetup, only:  mpi_initialize,mpi_initialize_io, mpi_cleanup, nproc, numproc, mpi_wtime
 ! grid information
 use gridinfo, only: getgridinfo, gridinfo_cleanup
 ! prior inflation
 use priorinf, only:prior_inf_ens, prior_inf_cleanup

 implicit none
 real(r_double) t1,t2
 
 ! initialize MPI.
 call mpi_initialize()

 ! read namelist.
 call read_namelist()

 ! initialize MPI communicator for IO tasks.
 call mpi_initialize_io(nanals)
 if (nproc==0) call w3tagb('ENKF_ANL',2011,0319,0055,'NP25')

 ! read horizontal grid information and pressure fields from
 ! 6-h forecast ensemble mean file.
 call getgridinfo()

 ! do prior inflation
 if (prior_inf == 1) then
     ! read inflation, ensemble mean and do prior inflation
     call prior_inf_ens()
     call prior_inf_cleanup()
 end if

 call gridinfo_cleanup()

 ! write log file (which script can check to verify completion).
 if (nproc .eq. 0) then
    call write_logfile()
 endif

 if (nproc==0) call w3tage('PRIOR_INF')
 call mpi_cleanup()

end program prior_inf_main
