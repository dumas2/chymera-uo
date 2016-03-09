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

Pure Subroutine Relax_2D(N, A, Tmp, rho, r_var)
!
! Relax_1D on the interior and the two halo cells shared with the left and right neighbors
!   - shared halo cells are computed twice and are not exchanged
!   - the outside halo cells are from neighbors and cannot be not computed
!
   Implicit None
   Integer, intent(in   ) :: N, m
   Real,    intent(inout) :: A  (-1:N+1,-1:N+1)
   Real,    intent(inout) :: Tmp(-1:N+1,-1:N+1)
   Real,    intent(in)    :: rho(-1:N+1,-1:N+1)
   Real,    intent(in)    :: r_var(-1:N+1)
   Real,                  :: dz,dr,pi
   Integer                :: i,j

   ! compute over extended region including boundary cells
   do i = 0, N
     do j = 0, N
      Tmp(i) = (1.0 - w)*A(i) + 0.5*w*( 0.5*(2.0*r_var(j)-dr)*dz/(r_var(j)*(dr*dr+dz*dz)+(m*dz*dr)**2)*A(i-1,j) &
                                     +  0.5*(2.0*r_var(j)+dr)*dz/(r_var(j)*(dr*dr+dz*dz)+(m*dz*dr)**2)*A(i+1,j) &
                                     +  r_var(j)*dr**2          /(r_var(j)*(dr*dr+dz*dz)+(m*dz*dr)**2)*A(i,j+1) &
                                     +  r_var(j)*dr**2          /(r_var(j)*(dr*dr+dz*dz)+(m*dz*dr)**2)*A(i,j-1) &
                                   +  r_var(j)*(dr*dz)**2*4.0*pi/(r_var(j)*(dr*dr+dz*dz)+(m*dz*dr)**2)*rho(i,j) )
     end do
   end do

   ! Do this in exchange halo...
   !    - probably should have rank information so that physical boundaries aren't changed
   !

   ! compute over just the interior
   do i = 1, N-1
      A(i) = (1.0 - w)*A(i) + 0.5*w*( 0.5*(2.0*r_var(j)-dr)*dz/(r_var(j)*(dr*dr+dz*dz)+(m*dz*dr)**2)*Tmp(i-1,j) &
                                  +  0.5*(2.0*r_var(j)+dr)*dz/(r_var(j)*(dr*dr+dz*dz)+(m*dz*dr)**2)*Tmp(i+1,j) &
                                  +  r_var(j)*dr**2          /(r_var(j)*(dr*dr+dz*dz)+(m*dz*dr)**2)*Tmp(i,j+1) &
                                  +  r_var(j)*dr**2          /(r_var(j)*(dr*dr+dz*dz)+(m*dz*dr)**2)*Tmp(i,j-1) &
                                  +  r_var(j)*(dr*dz)**2*4.0*pi/(r_var(j)*(dr*dr+dz*dz)+(m*dz*dr)**2)*rho(i,j) )
   end do


End Subroutine Relax_2D


End Module MultiGrid
