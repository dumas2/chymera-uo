Module RelaxMod
!=============================================================================
! Name        : poissonRelax.F90
! Version     : 1.0
! Description : Solves the Poisson equation in two dimensions.
!
! Method      :
!
!   Iterative (over)relaxation. 
!==============================================================================
Contains

Subroutine poissonRelax(Nj,Nk,Nl,rho,dr,dz,V1h)

Use MultiGrid, only : Residual,RelaxB
Use io       , only : readBoundary, readDensity, writeData

Implicit None

integer, Intent(in) :: Nj,Nk,Nl

real   , parameter  :: tol = 1e-8
integer   :: i,m,k,mn,ir,iz,p,je,ke
integer   :: nsteps = 1
integer   :: msteps = 5000
integer   :: diag = 500  

Real, Intent(inout) :: rho(1:Nj,1:Nk), V1h(0:Nj+1,0:Nk+1)
Real, allocatable   :: Tmp(:,:),Resid(:,:)

real      :: dr, dz, w, errmax

character(len=6) :: iter
character(len=6) :: numrlx

!write nsteps to string for textual output
write(numrlx,"(i6.6)") nsteps
print *, "size(rho) =",size(rho)
! Define dr and dz, these should be passed in eventually
dr = 0.09316442463639991
dz = 0.09316442463639991
!stop
! We are considering the m=1 Fourier mode for testing.
m = Nl

! Allocate arrays
Allocate(Tmp(-1:Nj+1,-1:Nk+1))
Allocate(Resid(-1:Nj+1,-1:Nk+1))

!! Initialize
!V1h =  -0.14
!rho =  0.0
resid= 0.0

!! Read in source term and boundary data from file
!call readDensity(rho,Nj,Nk)
! Assigns boundary data for potential to rho(:,Nk) and rho(Nj,:).
!call readBoundary(rho,Nj,Nk)

!! CHANGE THIS so that we start on the coursest mesh, and go up!!!
call writeData(0,Nj+1,Nk+1,V1h, numrlx // ".000")
call writeData(1,Nj,Nk,rho, numrlx // ".rho.000")

errmax = 1e6
mn = 0
!do while(mn.lt.msteps)
do while(errmax.gt.tol)
!... Relax solution on 1h mesh

!print *, "mn= ",mn
!    -------------------------
 do i = 1, nsteps
mn = mn + 1
!!! Update boundary before relaxing.
!! Top boundary, calculated directly before call to potential solve
  do ir = 1,Nj-1
    V1h(ir,Nk) = rho(ir,Nk)
 !   V1h(ir,1) = rho(ir,0)
  end do

!! Left and right boundary.
! Left boundary uses symmetry, right boundary is calculated directly.
  do iz = 1,Nk-1
    V1h(Nj,iz) = rho(Nj,iz)
 !   V1h(0 ,iz) = -V1h(1,iz)  ! Left boundary actually not required for relaxation.
  end do
! Dirichlet boundary conditions in the midplane.
do ir = 1,Nj-1
   V1h(ir,0) = V1h(ir,1) 
end do
 

! Relax on V1h using rho as source
  call RelaxB(Nj, Nk, m, V1h, Tmp, rho, dr, dz)

 end do

!! Call residual
! do i = 1, nsteps
!mn=mn+1
!!! Update boundary before relaxing.
!! Top boundary, calculated directly before call to potential solve
!  do ir = 1,Nj-1
!    V1h(ir,Nk) = rho(ir,Nk)
 !   V1h(ir,1) = rho(ir,0)
!  end do

!! Left and right boundary.
! Left boundary uses symmetry, right boundary is calculated directly.
!  do iz = 1,Nk-1
!    V1h(Nj,iz) = rho(Nj,iz)
 !   V1h(0 ,iz) = -V1h(1,iz)  ! Left boundary actually not required for relaxation.
!  end do
! Dirichlet boundary conditions in the midplane.
!do ir = 1,Nj-1
!   V1h(ir,0) = V1h(ir,1) 
!end do
 

! Relax on V1h using rho as source
! call RelaxB(Nj, Nk, m, V1h, Tmp, rho, dr, dz)

!end do


 call Residual(Nj, Nk, m, V1h, Resid, rho, dr, dz)

 if(mod(mn,diag)==0)then
   write(iter,"(i6.6)") mn
   call writeData(0,Nj+1,Nk+1,V1h, "sol." // iter )
   print *, "Wrote to file at iteration " // iter // "."

  call writeData(1,Nj,Nk, Resid, "resid." // iter ) 
 end if


  errmax = maxval(abs(Resid))
  write(75,"(i6.6,1x,f30.24)") mn,errmax



end do ! While loop

! Write results
call writeData(1,Nj-1,Nk-1,V1h, "final" )

deallocate(Tmp,Resid)

end subroutine PoissonRelax

end module relaxmod
