module intdbzmod
!$$$ module documentation block
!           .      .    .                                       .
! module:   intdbzmod    module for intdbz and its tangent linear intdbz_tl
!   prgmmr: carley         2011-04-06
!
! abstract: module for intdbz and its tangent linear intdbz_tl
!
! program history log:
!   2005-05-13  Yanqiu zhu - wrap intrw and its tangent linear intrw_tl into one module
!   2005-11-16  Derber - remove interfaces
!   2008-11-26  Todling - remove intrw_tl; add interface back
!   2009-08-13  lueken - update documentation
!
! subroutines included:
!   sub intdbz_
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
use m_obsNode, only: obsNode
use m_dbzNode, only: dbzNode
use m_dbzNode, only: dbzNode_typecast
use m_dbzNode, only: dbzNode_nextcast

implicit none

PRIVATE
PUBLIC intdbz
interface intdbz; module procedure &
          intdbz_
end interface
contains

subroutine intdbz_(dbzhead,rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intdbz       apply nonlin qc operator for radar reflectivity
!   prgmmr: derber           org: np23                date: 1991-02-26
!
! abstract: apply observation operator for radar winds
!             with nonlinear qc operator
!
! program history log:
!   1991-02-26  derber
!   1999-11-22  yang
!   2004-08-02  treadon - add only to module use, add intent in/out
!   2004-10-08  parrish - add nonlinear qc option
!   2005-03-01  parrish - nonlinear qc change to account for inflated obs error
!   2005-04-11  treadon - merge intrw and intrw_qc into single routine
!   2005-08-02  derber  - modify for variational qc parameters for each ob
!   2005-09-28  derber  - consolidate location and weight arrays
!   2006-07-28  derber  - modify to use new inner loop obs data structure
!                       - unify NL qc
!   2007-02-15  rancic  - add foto
!   2007-03-19  tremolet - binning of observations
!   2007-06-05  tremolet - use observation diagnostics structure
!   2007-07-09  tremolet - observation sensitivity
!   2008-01-04  tremolet - Don't apply H^T if l_do_adjoint is false
!   2008-11-28  todling  - turn FOTO optional; changed ptr%time handle
!   2010-05-13  todlng   - update to use gsi_bundle; update interface
!   2011-04-05  carley   - based upon nonlinear branch of intspd.f90, but
!                          modified for radar reflectivity 
!   2011-05-15  carley   - Updated to include tangent linear forward model
!   2011-05-19  carley  - Complete rewrite based upon intrad - routine is entirely linear
!                         
!   input argument list:
!     dbzhead   - obs type pointer to obs structure     
!     sqr       - current qr solution increment
!     sqli      - current qli solution increment
!     rqr
!     rqli
!
!   output argument list:
!     rqr        - qr results from observation operator
!     rqli       - qli results from observation operator
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: zero,half,one,two,tiny_r_kind,cg_term,r3600
  use obsmod, only: luse_obsdiag,lsaveobsens,l_do_adjoint
  use qcmod, only: nlnqc_iter,varqc_iter
  use gridmod, only: latlon1n
  use mpimod, only: mype
  use jfunc, only: jiter
  use gsi_bundlemod, only: gsi_bundle,gsi_bundlegetpointer
  use gridmod, only: wrf_mass_regional,wrf_nmm_regional,&
                     nems_nmmb_regional,cmaq_regional

  class(obsNode), pointer, intent(in   ) :: dbzhead
  type(gsi_bundle),        intent(in   ) :: sval
  type(gsi_bundle),        intent(inout) :: rval
! Declare local parameters
  real(r_kind),parameter:: ten    = 10.0_r_kind
! Declare local varibles
  integer(i_kind) j1,j2,j3,j4,j5,j6,j7,j8,ier,istatus
  real(r_kind) dbz,qrtl,qlitl, qstl, qgtl
  real(r_kind) dbztl,term
  real(r_kind) val,valqr,valqs,valqi,valqg,w1,w2,w3,w4,w5,w6,w7,w8 ,valqli !,valnr,valni modified
  real(r_kind) cg_dbz,p0,grad,wnotgross,wgross,time_dbz,pg_dbz
  real(r_kind),pointer,dimension(:) :: sqr,sqli, sqs, sqg !modified
  real(r_kind),pointer,dimension(:) :: rqr,rqli, rqs, rqg !,rnr,rni !,rqli
  type(dbzNode), pointer :: dbzptr
! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'qr',sqr,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'qr',rqr,istatus);ier=istatus+ier
 
  if( nems_nmmb_regional )then
    call gsi_bundlegetpointer(sval,'qli',sqli,istatus);ier=istatus+ier !modified
    call gsi_bundlegetpointer(rval,'qli',rqli,istatus);ier=istatus+ier
  else if ( wrf_mass_regional )then
    call gsi_bundlegetpointer(sval,'qs',sqs,istatus);ier=istatus+ier !modified
    call gsi_bundlegetpointer(rval,'qs',rqs,istatus);ier=istatus+ier
    call gsi_bundlegetpointer(sval,'qg',sqg,istatus);ier=istatus+ier !modified
    call gsi_bundlegetpointer(rval,'qg',rqg,istatus);ier=istatus+ier    
  end if
 
!  if(l_foto) then
!     write(6,*) 'Radar reflectivity assimilation is incompatible with l_foto=true'
!     return
!  endif
  if(ier/=0)return
  !dbzptr => dbzhead
  dbzptr => dbzNode_typecast(dbzhead)
  do while (associated(dbzptr))
     j1=dbzptr%ij(1)
     j2=dbzptr%ij(2)
     j3=dbzptr%ij(3)
     j4=dbzptr%ij(4)
     j5=dbzptr%ij(5)
     j6=dbzptr%ij(6)
     j7=dbzptr%ij(7)
     j8=dbzptr%ij(8)
     w1=dbzptr%wij(1)
     w2=dbzptr%wij(2)
     w3=dbzptr%wij(3)
     w4=dbzptr%wij(4)
     w5=dbzptr%wij(5)
     w6=dbzptr%wij(6)
     w7=dbzptr%wij(7)
     w8=dbzptr%wij(8)


!DEBUGGING:
!w1=1/8.0
!w2=1/8.0
!w3=1/8.0
!w4=1/8.0
!w5=1/8.0
!w6=1/8.0
!w7=1/8.0
!w8=1/8.0                   
 
     qrtl  = w1* sqr(j1)+w2* sqr(j2)+w3* sqr(j3)+w4* sqr(j4)+      &
             w5* sqr(j5)+w6* sqr(j6)+w7* sqr(j7)+w8* sqr(j8) !modified:

     if( nems_nmmb_regional )then
        qlitl = w1* sqli(j1)+w2* sqli(j2)+w3* sqli(j3)+w4* sqli(j4)+  &
                w5* sqli(j5)+w6* sqli(j6)+w7* sqli(j7)+w8* sqli(j8)                                     

        dbztl = (dbzptr%jqr)*qrtl + (dbzptr%jqli)*qlitl !modified

     else if ( wrf_mass_regional ) then
        qstl  = w1* sqs(j1)+w2* sqs(j2)+w3* sqs(j3)+w4* sqs(j4)+  &
                w5* sqs(j5)+w6* sqs(j6)+w7* sqs(j7)+w8* sqs(j8)

        qgtl  = w1* sqg(j1)+w2* sqg(j2)+w3* sqg(j3)+w4* sqg(j4)+  &
                w5* sqg(j5)+w6* sqg(j6)+w7* sqg(j7)+w8* sqg(j8)

        dbztl = (dbzptr%jqr)*qrtl + (dbzptr%jqs)*qstl + (dbzptr%jqg)*qgtl

     end if


!for comments here, D be delta, x is statespace variable(s) [Qr Qli]
!dbztl:  HDx [scalar]


     if (lsaveobsens) then
        dbzptr%diags%obssen(jiter)=dbzptr%raterr2*dbzptr%err2*dbztl
     else
        if (dbzptr%luse) dbzptr%diags%tldepart(jiter)=dbztl
     endif
     if (l_do_adjoint) then
        if (lsaveobsens) then
           grad=dbzptr%diags%obssen(jiter)
        else                         
           dbz=dbztl-dbzptr%ddiff              !from setupdbz (ob-NLges)
!let Yo be (Hxb - Y)
!dbz:  HDx - Yo
           grad=dbzptr%raterr2*dbzptr%err2*dbz
!grad = R^-1  * (HDx - Yo)
        endif

!!!apply adjoint here
![DQr DQs DQg] = H^T * grad

!!variable transform part of adjoint
        valqr=dbzptr%jqr*grad !modified
        if( nems_nmmb_regional ) then
           valqli=dbzptr%jqli*grad       
        else if ( wrf_mass_regional ) then
           valqs=dbzptr%jqs*grad
           valqg=dbzptr%jqg*grad
        end if

!print *, "dbztl",dbztl
!print *, "intdbzline197:",j1,j2,j3,j4,j5,j6,j7,j8
!print *, "intdbzline198:",w1,w2,w3,w4,w4,w6,w7,w8
!!print *, "rqr:",rqr(j1),rqr(j2),rqr(j3),rqr(j4),rqr(j5),rqr(j6),rqr(j7),rqr(j8)
!print *, "intdbzline199:",valqs,dbzptr%err2
!!print *, "intdbzline200:",sqr(j1),sqr(j2),sqr(j3),sqr(j4),sqr(j5),sqr(j6),sqr(j7),sqr(j8)
!print *, "intdbzline203:",dbzptr%ddiff

!spatial interpolation part of adjoint
        rqr(j1)=rqr(j1)+w1*valqr
        rqr(j2)=rqr(j2)+w2*valqr
        rqr(j3)=rqr(j3)+w3*valqr
        rqr(j4)=rqr(j4)+w4*valqr
        rqr(j5)=rqr(j5)+w5*valqr
        rqr(j6)=rqr(j6)+w6*valqr
        rqr(j7)=rqr(j7)+w7*valqr
        rqr(j8)=rqr(j8)+w8*valqr  !modified  
 !print *, 'valqr: ',valqr
        
     if( nems_nmmb_regional ) then

        rqli(j1)=rqli(j1)+w1*valqli
        rqli(j2)=rqli(j2)+w2*valqli
        rqli(j3)=rqli(j3)+w3*valqli
        rqli(j4)=rqli(j4)+w4*valqli
        rqli(j5)=rqli(j5)+w5*valqli
        rqli(j6)=rqli(j6)+w6*valqli
        rqli(j7)=rqli(j7)+w7*valqli
        rqli(j8)=rqli(j8)+w8*valqli

     else if( wrf_mass_regional ) then

        rqs(j1)=rqs(j1)+w1*valqs
        rqs(j2)=rqs(j2)+w2*valqs
        rqs(j3)=rqs(j3)+w3*valqs
        rqs(j4)=rqs(j4)+w4*valqs
        rqs(j5)=rqs(j5)+w5*valqs
        rqs(j6)=rqs(j6)+w6*valqs
        rqs(j7)=rqs(j7)+w7*valqs
        rqs(j8)=rqs(j8)+w8*valqs

        rqg(j1)=rqg(j1)+w1*valqg
        rqg(j2)=rqg(j2)+w2*valqg
        rqg(j3)=rqg(j3)+w3*valqg
        rqg(j4)=rqg(j4)+w4*valqg
        rqg(j5)=rqg(j5)+w5*valqg
        rqg(j6)=rqg(j6)+w6*valqg
        rqg(j7)=rqg(j7)+w7*valqg
        rqg(j8)=rqg(j8)+w8*valqg

     end if ! nems_nmmb_regional or wrf_mass_regional 
     !print*,'qs-intdbz', minval(rqs),maxval(rqs)
       
     end if ! l_adjoint
     dbzptr => dbzptr%llpoint
  end do
  return
end subroutine intdbz_
end module intdbzmod
