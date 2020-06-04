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
integer(i_kind), allocatable, dimension(:) :: scounts, sdispls, rcounts, rdispls

contains

subroutine prior_inf_ens()
implicit none
real(r_double) t1,t2
integer(i_kind) nanal,nn,i,nb,n
integer(i_kind) ierr, np
real(r_single), allocatable, dimension(:,:,:) :: newgridin
real(r_single), allocatable, dimension(:) :: sendbuf,recvbuf

! must at least nanals tasks allocated.
if (numproc < nanals) then
  print *,'need at least nanals =',nanals,'MPI tasks, exiting ...'
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_finalize(ierr)
end if

! scounts is number of data elements to send to processor np.
! rcounts is number of data elements to recv from processor np.
! displs is displacement into send array for data to go to proc np
allocate(scounts(0:numproc-1))
allocate(sdispls(0:numproc-1))
allocate(rcounts(0:numproc-1))
allocate(rdispls(0:numproc-1))
sdispls = 0 

if (nproc <= nanals-1) then
   scounts = npts*ndim ! every task sends corresponding member to root task
else
   scounts = 0
endif
if (nproc == 0) then
   do np=0,numproc-1
      rcounts(np) = npts*ndim ! only root task receives member info and calculate the ensemble mean
      rdispls(np) = np*npts*ndim
   enddo
else
   rcounts = 0
   rdispls = 0 
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
endif

call mpi_barrier(mpi_comm_world, ierr)
! ensemble mean
allocate(gridin_ensmean(npts,ndim,nbackgrounds))
gridin_ensmean = 0
t1 = mpi_wtime()
allocate(sendbuf(npts*ndim)) ! every task send info for one member to root task
if(nproc == 0) then
   allocate(recvbuf(nanals*npts*ndim)) ! root task receives info of every member
else
   allocate(recvbuf(1))
endif
write(6,*) 'nproc,nbackgrounds,numproc,nanals,ndim,npts=',nproc,nbackgrounds,numproc,nanals,ndim,npts
! send and receive buffers.
do nb=1,nbackgrounds ! loop over time levels in background

if (nproc <= nanals-1) then
   write(6,*) 'nproc,npts,ndim=',nproc,npts,ndim
   ! fill up send buffer.
   do nn=1,ndim
     do i=1,npts
       n = ((nn-1)*ndim + i
       sendbuf(n) = gridin(i,nn,nb)
     enddo
   enddo
end if
print *,'----1111111111111111111'
call mpi_alltoallv(sendbuf, scounts,sdispls, mpi_real4, recvbuf, rcounts, rdispls,&
                   mpi_real4, mpi_comm_world, ierr)
print *,'########1111111111111111111'

!==> compute ensemble mean of first guesses on root task.
if(nproc == 0) then
  do nn=1,ndim
    do i=1,npts
      do nanal=1,nanals
         n = (nanal-1)*npts*ndim+((nn-1)*ndim + i
         gridin_ensmean(i,nn,nb) = gridin_ensmean(i,nn,nb)+recvbuf(n)
      enddo
      gridin_ensmean(i,nn,nb) = gridin_ensmean(i,nn,nb)/float(nanals)
    enddo
  end do
  print *,'********1111111111111111111'
end do

!broadcast 

enddo ! loop over nbackgrounds
print *,'1111111111111122222'
deallocate(sendbuf, recvbuf)
if (nproc == 0) then
     t2 = mpi_wtime()
     print *,'time in calculating ensmean on root',t2-t1,'secs'
end if

! read the inflation value
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
       gridin_prior_inf(:,:,nb) = 1.0_r_single + prior_inf_damping*(gridin_prior_inf(:,:,nb)-1.0_r_single) 
      ! set lower and upper bounds
       where (gridin_prior_inf(:,:,nb) < prior_inf_lower_bound ) gridin_prior_inf(:,:,nb) = prior_inf_lower_bound 
       where (gridin_prior_inf(:,:,nb) > prior_inf_upper_bound ) gridin_prior_inf(:,:,nb) = prior_inf_upper_bound
   end if 
   newgridin(:,:,nb) = (gridin(:,:,nb)-gridin_ensmean(:,:,nb))*sqrt(gridin_prior_inf(:,:,nb)) + gridin_ensmean(:,:,nb)
   gridin(:,:,nb) = newgridin(:,:,nb) - gridin(:,:,nb) ! the writegriddata in gridio.F90 writes the increment (not the updated values)  
enddo ! loop over nbackgrounds

! write the inflated ensemble members back
t1 = mpi_wtime()
call writegriddata(nanal,gridin)
if (nproc == 0) then
   t2 = mpi_wtime()
   print *,'time in writegridata on root',t2-t1,'secs'
end if

! write the inflation values back only on root processor ( if damping performed)
t1 = mpi_wtime()
if (nproc == 0) then
   call writegriddata_f(prior_inf_file,gridin_prior_inf) ! write the updated inflation value
   t2 = mpi_wtime()
   print *,'time in writegridata on root',t2-t1,'secs'
end if

   ! deallocate allocatable arrays.
deallocate(newgridin)

end subroutine prior_inf_ens

subroutine prior_inf_cleanup()
! deallocate module-level allocatable arrays.
 if (nproc <= nanals-1 .and. allocated(gridin))  deallocate(gridin)
 if (nproc <= nanals-1 .and. allocated(gsatin))  deallocate(gsatin)
 if (nproc <= nanals-1 .and. allocated(gridin_prior_inf))  deallocate(gridin_prior_inf)
 if (nproc <= nanals-1 .and. allocated(gridin_ensmean))  deallocate(gridin_ensmean)
 if (nproc <= nanals-1 .and. allocated(scounts))  deallocate(scounts)
 if (nproc <= nanals-1 .and. allocated(displs))  deallocate(sdispls)
 if (nproc <= nanals-1 .and. allocated(rcounts))  deallocate(rcounts)
 if (nproc <= nanals-1 .and. allocated(displs))  deallocate(rdispls)
end subroutine prior_inf_cleanup
end module priorinf

