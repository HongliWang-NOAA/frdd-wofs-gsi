module priorinf
!$$$  main program documentation block
!
! program:  priorinf              
!
! prgmmr: jjh          org: CIMMS/OU (NSSL)               date: 2017-03-16
!
! abstract: This module does the prior adaptive inflation.
!
! program history log:
!   2017-03-16  Initial version.
!
! attributes:
!   language: f90

use mpisetup
use params, only: nlevs,nvars,ndim,nanals,nbackgrounds,prior_inf_file,prior_inf_sd_file, &
                  prior_inf_damping, prior_inf_initial, prior_inf_sd_initial,prior_inf_lower_bound,prior_inf_upper_bound,prior_inf_sd_lower_bound !
use kinds, only: r_kind, i_kind, r_double, r_single
use gridio, only: readgriddata, writegriddata
use gridio_inf, only: readgriddata_f,writegriddata_f
use gridinfo, only: npts
use loadbal, only: npts_max,indxproc,numptsperproc

implicit none
private
public prior_inf_ens,prior_inf_cleanup
real(r_single),public, allocatable, dimension(:,:,:) :: gridin
real(r_single),public, allocatable, dimension(:,:,:) :: gridin_prior_inf,gridin_ensmean
real(r_double),public, allocatable, dimension(:,:,:) :: gsatin

contains

subroutine prior_inf_ens()
implicit none
real(r_double) t1,t2
integer(i_kind) nanal,nn,i,nb,n
integer(i_kind) ierr, np
real(r_single), allocatable, dimension(:,:,:) :: newgridin
! New variables needed for the deterministic approach
real(r_single), allocatable, dimension(:) :: recv_buffer
integer(i_kind) :: root_proc
! --- END DECLARATION SECTION ---

! --- BEGIN EXECUTABLE SECTION ---
root_proc = 0

! must at least nanals tasks allocated.
if (numproc < nanals) then
  print *,'need at least nanals =',nanals,'MPI tasks, exiting ...'
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_finalize(ierr)
end if

! each process is reading one member 
if (nproc <= nanals-1) then
   ! read an ensemble member
   allocate(gridin(npts,ndim,nbackgrounds))
   allocate(gsatin(npts,nlevs,nbackgrounds))
   nanal = nproc + 1
   t1 = mpi_wtime()
   call readgriddata(nanal,gridin,gsatin)
   if (nproc == 0) then
     t2 = mpi_wtime()
     print *,'time in readgridata on root',t2-t1,'secs'
   end if

   ! ensemble mean
   allocate(gridin_ensmean(npts,ndim,nbackgrounds))
   gridin_ensmean = gridin
   t1 = mpi_wtime()

   ! Allocate a buffer on the root process to receive data from all other processes.
   ! The size is npts (for one variable) * nanals (number of processes).
   if (nproc == root_proc) then
      allocate(recv_buffer(npts * nanals))
   end if

   ! Loop over backgrounds and dimensions as before
   do nb=1, nbackgrounds
       do nn=1, ndim
           ! 1. GATHER data from all processes to the root process.
           ! Each process sends its slice of the gridin array for the current variable.
           call mpi_gather(gridin(1,nn,nb), npts, mpi_real4, &
                           recv_buffer, npts, mpi_real4, &
                           root_proc, mpi_comm_io, ierr)

           ! 2. SUM SERIALLY on the root process
           if (nproc == root_proc) then
               ! Clear the ensemble mean array for summation
               gridin_ensmean(1:npts,nn,nb) = 0.0_r_single
               ! Loop through the gathered data from each process and add it up
               do n=1, nanals
                   gridin_ensmean(1:npts,nn,nb) = gridin_ensmean(1:npts,nn,nb) + &
                                                  recv_buffer((n-1)*npts + 1 : n*npts)
               end do

               ! Calculate the final mean
               gridin_ensmean(1:npts,nn,nb) = gridin_ensmean(1:npts,nn,nb) / real(nanals)
           end if

           ! 3. BROADCAST the deterministic result from the root to all other processes
           call mpi_bcast(gridin_ensmean(1,nn,nb), npts, mpi_real4, root_proc, mpi_comm_io, ierr)

       end do ! end nn loop
   end do ! end nb loop

   ! Deallocate the receive buffer on the root process
   if (nproc == root_proc) then
       deallocate(recv_buffer)
   end if

! Begin old broken code
!   ! convert gridin_ensmean from gridin to sum(gridin).
!   do nb=1,nbackgrounds ! loop over time levels in background
!      do nn=1,ndim
!         call mpi_allreduce(mpi_in_place,gridin_ensmean(1,nn,nb),npts, mpi_real4, mpi_sum,mpi_comm_io,ierr)
!      enddo
!   enddo
!
!   ! calculate ensemble mean
!   gridin_ensmean = gridin_ensmean/real(nanals) 
! End old broken code

   !print *,'nproc=',nproc,',gridin(45000,:,1)=',gridin(45000,1:5,1),',gridin_ensmean(45000,:,1)=',gridin_ensmean(45000,1:5,1)

   if (nproc == 0) then
     t2 = mpi_wtime()
     print *,'time in calculating ensmean on root',t2-t1,'secs'
   end if

   ! read the inflation value (every task reads one copy)
   allocate(gridin_prior_inf(npts,ndim,nbackgrounds))
   gridin_prior_inf = prior_inf_initial
   t1 = mpi_wtime()
   call readgriddata_f(prior_inf_file,gridin_prior_inf)
   if (nproc == 0) then
       t2 = mpi_wtime()
       print *,'time in readgridata_f for prior inflation on root',t2-t1,'secs'
   end if

   ! do the inflation for each state variable
   allocate(newgridin(npts,ndim,nbackgrounds))

   do nb=1,nbackgrounds
      if(prior_inf_damping>0.0_r_single .and. prior_inf_damping<1.0_r_single) then
         ! damping first
         !print *,'nproc=',nproc,',init_inf=',gridin_prior_inf(45000,1:5,nb)
         gridin_prior_inf(:,:,nb) = 1.0_r_single + prior_inf_damping*(gridin_prior_inf(:,:,nb)-1.0_r_single) 
         ! set lower and upper bounds
         where (gridin_prior_inf(:,:,nb) < prior_inf_lower_bound ) gridin_prior_inf(:,:,nb) = prior_inf_lower_bound 
         where (gridin_prior_inf(:,:,nb) > prior_inf_upper_bound ) gridin_prior_inf(:,:,nb) = prior_inf_upper_bound
      end if 
      newgridin(:,:,nb) = (gridin(:,:,nb)-gridin_ensmean(:,:,nb))*sqrt(gridin_prior_inf(:,:,nb)) + gridin_ensmean(:,:,nb)
      !print *,'nproc=',nproc,',damp_inf=',gridin_prior_inf(45000,1:5,nb),',gridin(:,:,nb)=',gridin(45000,1:5,nb),'newgridin(:,:,nb)=',newgridin(45000,1:5,nb)
      gridin(:,:,nb) = newgridin(:,:,nb) - gridin(:,:,nb) ! the writegriddata in gridio.F90 writes the increment (not the updated values)  
      !print *,'nproc=',nproc,',diff=',gridin(45000,1:5,nb)
    enddo ! loop over nbackgrounds

    ! deallocate allocatable arrays.
    deallocate(newgridin)

! write the inflated ensemble members back
! each process is writing one member 
   t1 = mpi_wtime()
   call writegriddata(nanal,gridin)
   if (nproc == 0) then
     t2 = mpi_wtime()
     print *,'time in writegridata on root',t2-t1,'secs'
   end if
endif

! write the inflation values back only on root processor ( if damping performed)
if (nproc == 0) then
   t1 = mpi_wtime()
   call writegriddata_f(prior_inf_file,gridin_prior_inf) ! write the updated inflation value
   t2 = mpi_wtime()
   print *,'time in writegridata_f on root',t2-t1,'secs'
end if


end subroutine prior_inf_ens

subroutine prior_inf_cleanup()
! deallocate module-level allocatable arrays.
 if (nproc <= nanals-1 .and. allocated(gridin))  deallocate(gridin)
 if (nproc <= nanals-1 .and. allocated(gsatin))  deallocate(gsatin)
 if (nproc <= nanals-1 .and. allocated(gridin_prior_inf))  deallocate(gridin_prior_inf)
 if (nproc <= nanals-1 .and. allocated(gridin_ensmean))  deallocate(gridin_ensmean)
end subroutine prior_inf_cleanup
end module priorinf

