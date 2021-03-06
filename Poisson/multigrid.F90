Module MultiGrid
!
! old grid with HALO(-1:*:1), DIMENSION(0:N)
!
!  x |          o          | x    N/4-1 interior  N/4+1  total  [N=8]
!  x |    o     o     o    | x    N/2-1 interior  N/2+1  total  [N=8]
!  x | o  o  o  o  o  o  o | x    N-1   interior  N  +1  total  [N=8]
!
!    |       x| ^ |x      ???       N interior N+2 total for one parallel node (boundary shared)
!         proc boundary   ???       proc boundary is shared point
!

! new grid with HALO(-2:*:2), DIMENSION(-1:N+1)
!
! N-1, N/2-1,N/4-1,... interior cells, N+3,N/2+3,N/4+3,... total cells
!  ____________________________________________________________________________
!  x           x |          o          | o |          o          | x           x
!        x     x |    o     o     o    | o |    o     o     o    | x     x
!           x  x | o  o  o  o  o  o  o | o | o  o  o  o  o  o  o | x  x
!                | 2  3  4                                    256| 257  258
!          -1  0 | 1  2  3                                    255| 256 257        
!
!  N=8, NP=2
!
!  N-1 interior, 4 halo cells with N+3 total points
!    - one cell on each boundary shared with other processors
!    - one halo cell at each end
!    - could treat first border cell as interior for   
!
  
Contains

Subroutine ExchangeHalo(Nj, Nk, A)
!
! Exchange halo information between neighboring processes
!
use MPI_F08, only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
use MPI_F08, only : MPI_Send, MPI_Recv, MPI_DOUBLE_PRECISION, MPI_Status
   Implicit None
   Integer, intent(in   ) :: Nj,Nk
   Real,    intent(inout) :: A(-1:Nj+1,-1:Nk+1)

   Integer :: above, below, tag,bufSize
   Real    :: topBC, bottomBC   
   integer :: numRanks, rank
   type(MPI_Status) :: status
   real  , allocatable :: bufFromAbove(:),bufFromBelow(:)

Call MPI_Comm_size( MPI_COMM_WORLD, numRanks)
Call MPI_Comm_rank( MPI_COMM_WORLD, rank)

   allocate(bufFromAbove(Nj+3),bufFromBelow(Nj+3))
   bufSize=NJ+3
   tag    = 156
   above  = rank + 1
   below  = rank - 1
    
   !! MPI halo exchange for parallel version
!Rank=0, bottom of grid, no below
if (rank == 0) then
  Call MPI_Send(A(:,Nk-1), bufSize , MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD)
  Call MPI_Recv(bufFromAbove, bufSize, MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD, status)
  A(:,Nk+1) = bufFromAbove

!Rank=numRanks-1, top of grid, no above
else if (rank == numRanks-1) then

  Call MPI_Send(A(:,1), bufSize , MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD)
  Call MPI_Recv(bufFromBelow, bufSize, MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD, status)
  A(:,-1) = bufFromBelow
!Middle of grid
else
  Call MPI_Send(A(:,Nk-1), bufSize , MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD)

  Call MPI_Send(A(:,1), bufSize , MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD)
  Call MPI_Recv(bufFromBelow, bufSize, MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD, status)
  Call MPI_Recv(bufFromAbove, bufSize, MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD, status)
  A(:,Nk+1) = bufFromAbove
  A(:, -1) = bufFromBelow

end if   

deallocate(bufFromAbove,bufFromBelow)
End Subroutine ExchangeHalo


Subroutine Prolongate(J, K, V1h, V2h)
!
!  Interpolation (prolongation) operator R^(J/2+1,K/2+1) => R^(J+1,K+1)
!
!  J-1 is the number of interior fine grid cells in radial direction
use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
!
  implicit none
  integer, intent(in) :: J, K
  real :: V1h(-1:J+1,-1:K+1), V2h(-1:J/2+1,-1:K/2+1)
  integer :: ir, iz, jj, kk, m, n, rank, numRanks

call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)

  m = J/2 - 1     ! # interior coarse cells in radial direction
  n = K/2 - 1     ! # interior coarse cells in vertical direction

  do iz = 0, n+1
    kk = 2*iz
    do ir = 0, m
      jj = 2*ir

      V1h(jj,kk) = V2h(ir,iz)
      V1h(jj+1,kk  ) =  .5*(V2h(ir,iz) + V2h(ir+1,iz  ))
      V1h(jj  ,kk+1) =  .5*(V2h(ir,iz) + V2h(ir  ,iz+1))
      V1h(jj+1,kk+1) = .25*(V2h(ir,iz) + V2h(ir,iz+1) + V2h(ir+1,iz) + V2h(ir+1,iz+1))
    end do
  end do

 if(rank.eq.numRanks-1) then
  do ir = 0, m
     V1h(2*ir,K) = V2h(ir,n+1)      ! top   halo cells
     V1h(2*ir+1,K) = .5*(V2h(ir,n+1)+V2h(ir+1,n+1))      ! top   halo cells
  end do
 end if


  do iz = 0, n+1
     V1h(J,2*iz) = V2h(m+1,iz)      ! right halo cells
     V1h(J,2*iz+1) = .5*(V2h(m+1,iz)+V2h(m+1,iz+1))      ! right halo cells
  end do
End Subroutine Prolongate

Subroutine Restrict(J, K, V1h, V2h)
!
!  Restriction operator R^(J+1,K+1) => R^(J/2+1,K/2+1)
!
!  J-1 is the number of interior fine grid cells in radial direction
!  K-1 is the number of interior fine grid cells in vertical direction
!
  implicit none
  integer, intent(in) :: J, K
  real :: V1h(-1:J+1,-1:K+1), V2h(-1:J/2+1,-1:K/2+1)
  integer :: ir, iz, jj, kk, m, n

  m = J/2 - 1     ! # interior coarse cells in radial direction
  n = K/2 - 1     ! # interior coarse cells in vertical direction

  do iz = 0, n+1
    kk = 2*iz
    do ir = 0, m+1
      jj = 2*ir

      V2h(ir,iz) = .25*(.25*V1h(jj-1,kk+1) + .5*V1h(jj,kk+1) + .25*V1h(jj+1,kk+1)   &
                      +  .5*V1h(jj-1,kk  ) +    V1h(jj,kk  ) +  .5*V1h(jj+1,kk  )   &
                      + .25*V1h(jj-1,kk-1) + .5*V1h(jj,kk-1) + .25*V1h(jj+1,kk-1))
    end do
  end do

End Subroutine Restrict

Subroutine RestrictRho(J, K, V1h, V2h)
!
!  Restriction operator R^(J+1,K+1) => R^(J/2+1,K/2+1)
!
!  J-1 is the number of interior fine grid cells in radial direction
!  K-1 is the number of interior fine grid cells in vertical direction
!
use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
!
  implicit none
  integer, intent(in) :: J, K
  real :: V1h(-1:J+1,-1:K+1), V2h(-1:J/2+1,-1:K/2+1)
  integer :: ir, iz, jj, kk, m, n,rank,numRanks

call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)
  m = J/2 - 1     ! # interior coarse cells in radial direction
  n = K/2 - 1     ! # interior coarse cells in vertical direction

  do iz = 0, n
    kk = 2*iz
    do ir = 0, m
      jj = 2*ir

      V2h(ir,iz) = .25*(.25*V1h(jj-1,kk+1) + .5*V1h(jj,kk+1) + .25*V1h(jj+1,kk+1)   &
                      +  .5*V1h(jj-1,kk  ) +    V1h(jj,kk  ) +  .5*V1h(jj+1,kk  )   &
                      + .25*V1h(jj-1,kk-1) + .5*V1h(jj,kk-1) + .25*V1h(jj+1,kk-1))
    end do
  end do

  do iz=0,n+1
     V2h(m+1,iz) = 0.25*(V1h(J,2*iz-1)+2.0*V1h(J,2*iz)+V1h(J,2*iz+1))
  end do
if(rank.eq.numRanks-1)then
  do ir=0,m
     V2h(ir,n+1) = 0.25*(V1h(2*ir-1,K)+2.*V1h(2*ir,K)+V1h(2*ir+1,K))
  end do
     V2h(m+1,n+1) = V1h(J,K)
else
  do ir=0,m
    jj = 2*ir
     kk = 2*(n+1)
      V2h(ir,n+1) = .25*(.25*V1h(jj-1,kk+1) + .5*V1h(jj,kk+1) + .25*V1h(jj+1,kk+1)   &
                      +  .5*V1h(jj-1,kk  ) +    V1h(jj,kk  ) +  .5*V1h(jj+1,kk  )   &
                      + .25*V1h(jj-1,kk-1) + .5*V1h(jj,kk-1) + .25*V1h(jj+1,kk-1))
  end do
end if

End Subroutine RestrictRho

Subroutine Relax(Nj, Nk, m, A, Tmp, rho, dr, dz)
!
! Relax_2D on the interior and the two halo cells shared with the left and right neighbors
!   - shared halo cells are computed twice and are not exchanged
!   - the outside halo cells are from neighbors and cannot be not computed
!
use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
   Implicit None
   Integer, intent(in   ) :: Nj, Nk, m
   Real,    intent(inout) :: A(-1:Nj+1,-1:Nk+1)
   Real,    intent(inout) :: Tmp(-1:Nj+1,-1:Nk+1)
   Real,    intent(in)    :: rho(-1:Nj+1,-1:Nk+1),dr,dz
   Real                   :: r_var(-1:Nj+1)
   Real                   :: pi,dtheta,m1,w
   Integer                :: ir,jk,i,numRanks,rank,iz
call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)

   pi=acos(-1.0)
   dtheta=2.*pi/dble(512)
   
   w = 1.3d0
!   m1 = (cos((m-1)*dtheta)-1.)/dtheta/dtheta
   m1 = float(m)
   do i = 0,Nj+1
     r_var(i) = (float(i)-0.5)*dr
   end do

 if (rank==0) then
  if (numRanks>1) then
   do jk = Nk,1,-1
    do ir = Nj-1,1,-1
    A(ir,jk) =  (1.0-w)*A(ir,jk) + w*1.0/(2.0/dr/dr+2.0/dz/dz+m1*m1/r_var(ir)/r_var(ir))*(                     &
                 (1.0/dr/dr-1.0/2.0/r_var(ir)/dr)*A(ir-1,jk) &
              +  (1.0/dr/dr+1.0/2.0/r_var(ir)/dr)*A(ir+1,jk) &
              +  1.0/dz/dz                     *A(ir,jk+1) &
              +  1.0/dz/dz                     *A(ir,jk-1) &
              -  rho(ir,jk)    )
       end do
   end do
  else
   do jk = Nk-1,1,-1
    do ir = Nj-1,1,-1
    A(ir,jk) =  (1.0-w)*A(ir,jk) + w*1.0/(2.0/dr/dr+2.0/dz/dz+m1*m1/r_var(ir)/r_var(ir))*(                     &
                 (1.0/dr/dr-1.0/2.0/r_var(ir)/dr)*A(ir-1,jk) &
              +  (1.0/dr/dr+1.0/2.0/r_var(ir)/dr)*A(ir+1,jk) &
              +  1.0/dz/dz                     *A(ir,jk+1) &
              +  1.0/dz/dz                     *A(ir,jk-1) &
              -  rho(ir,jk)    )
       end do
   end do
  end if 
! compute over extended region including boundary cells
else if (rank==numRanks-1) then 
   do jk = Nk-1,0,-1
    do ir = Nj-1,1,-1
    A(ir,jk) =  (1.0-w)*A(ir,jk) + w*1.0/(2.0/dr/dr+2.0/dz/dz+m1*m1/r_var(ir)/r_var(ir))*(                     &
                 (1.0/dr/dr-1.0/2.0/r_var(ir)/dr)*A(ir-1,jk) &
              +  (1.0/dr/dr+1.0/2.0/r_var(ir)/dr)*A(ir+1,jk) &
              +  1.0/dz/dz                     *A(ir,jk+1) &
              +  1.0/dz/dz                     *A(ir,jk-1) &
              -  rho(ir,jk)    )
    end do
   end do
else
   do jk = Nk,0,-1
    do ir = Nj-1,1,-1
    A(ir,jk) =  (1.0-w)*A(ir,jk) + w*1.0/(2.0/dr/dr+2.0/dz/dz+m1*m1/r_var(ir)/r_var(ir))*(                     &
                 (1.0/dr/dr-1.0/2.0/r_var(ir)/dr)*A(ir-1,jk) &
              +  (1.0/dr/dr+1.0/2.0/r_var(ir)/dr)*A(ir+1,jk) &
              +  1.0/dz/dz                     *A(ir,jk+1) &
              +  1.0/dz/dz                     *A(ir,jk-1) &
              -  rho(ir,jk)    )
       end do
   end do
end if
End Subroutine Relax


Subroutine Residual(Nj, Nk, m, A, Resid, rho, dr, dz)
!
! Relax_2D on the interior and the two halo cells shared with the left and right neighbors
!   - shared halo cells are computed twice and are not exchanged
!   - the outside halo cells are from neighbors and cannot be not computed
!
use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
   Implicit None
   Integer, intent(in   ) :: Nj, Nk, m
   Real,    intent(inout) :: A(-1:Nj+1,-1:Nk+1)
   Real,    intent(inout) :: Resid(-1:Nj+1,-1:Nk+1)
   Real,    intent(in)    :: rho(-1:Nj+1,-1:Nk+1),dr,dz
   Real                   :: r_var(-1:Nj+1)
   Real                   :: pi,w,m1,dtheta
   Integer                :: ir,jk,i,numRanks,rank
call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)


   pi=acos(-1.0)
   dtheta = 2.0*pi/dble(512)
   m1 = float(m)
!  m1 = (cos((m-1)*dtheta)-1.)/dtheta/dtheta
   do i = 0,Nj
     r_var(i) = (float(i)-0.5)*dr
   end do

if (numRanks>1) then
 if (rank==numRanks-1) then
!! Calculate the residual 
   do jk = 0, Nk-1
     do ir = 1, Nj-1
      Resid(ir,jk) =    (                                                 &
                 (1.0/dr/dr-1.0/2.0/r_var(ir)/dr)*A(ir-1,jk)              &
              +  (1.0/dr/dr+1.0/2.0/r_var(ir)/dr)*A(ir+1,jk)              &
              +  1.0/dz/dz                     *A(ir,jk+1)                &
              +  1.0/dz/dz                     *A(ir,jk-1)                &
              - (2.0/dr/dr+2.0/dz/dz+m1*m1/r_var(ir)/r_var(ir))*A(ir,jk)     &  
              -  rho(ir,jk)    )
     end do
   end do
 else
   do jk = 0, Nk
     do ir = 1, Nj-1
      Resid(ir,jk) =    (                                                 &
                 (1.0/dr/dr-1.0/2.0/r_var(ir)/dr)*A(ir-1,jk)              &
              +  (1.0/dr/dr+1.0/2.0/r_var(ir)/dr)*A(ir+1,jk)              &
              +  1.0/dz/dz                     *A(ir,jk+1)                &
              +  1.0/dz/dz                     *A(ir,jk-1)                &
              - (2.0/dr/dr+2.0/dz/dz+m1*m1/r_var(ir)/r_var(ir))*A(ir,jk)     &  
              -  rho(ir,jk)    )
     end do
   end do
 end if
else
   do jk = 1, Nk-1
     do ir = 1, Nj-1
      Resid(ir,jk) =    (                                                 &
                 (1.0/dr/dr-1.0/2.0/r_var(ir)/dr)*A(ir-1,jk)              &
              +  (1.0/dr/dr+1.0/2.0/r_var(ir)/dr)*A(ir+1,jk)              &
              +  1.0/dz/dz                     *A(ir,jk+1)                &
              +  1.0/dz/dz                     *A(ir,jk-1)                &
              - (2.0/dr/dr+2.0/dz/dz+m1*m1/r_var(ir)/r_var(ir))*A(ir,jk)     &  
              -  rho(ir,jk)    )
     end do
   end do
end if
End Subroutine Residual

End Module MultiGrid
