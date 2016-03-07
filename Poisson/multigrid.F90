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
!
!  x           x |          o          | o |          o          | x           x
!        x     x |    o     o     o    | o |    o     o     o    | x     x
!           x  x | o  o  o  o  o  o  o | o | o  o  o  o  o  o  o | x  x
!  
!  N=8, NP=2
!
!  N-1 interior, 4 halo cells with N+3 total points
!    - one cell on each boundary shared with other processors
!    - one halo cell at each end
!    - could treat first border cell as interior for   
!
  
Contains

Subroutine Prolongate(J, K, V1h, V2h)
!
!  Interpolation (prolongation) operator R^(J/2+1,K/2+1) => R^(J+1,K+1)
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

  do iz = 0, n
    kk = 2*iz
    do ir = 0, m
      jj = 2*ir

      V1h(jj,kk) = V2h(ir,iz)

      V1h(jj+1,kk  ) =  .5*(V2h(ir,iz) + V2h(ir+1,iz  ))
      V1h(jj  ,kk+1) =  .5*(V2h(ir,iz) + V2h(ir  ,iz+1))

      V1h(jj+1,kk+1) = .25*(V2h(ir,iz) + V2h(ir,iz+1) + V2h(ir+1,iz) + V2h(ir+1,iz+1))
    end do
  end do

  do iz = 0, n
     V1h(J,iz) = V2h(n+1,iz)      ! right halo cells
  end do
  do ir = 0, m
     V1h(ir,K) = V2h(ir,m+1)      ! top   halo cells
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

End Module MultiGrid