      subroutine Fluid_Setup
      use Fluid_Mod
      implicit  none
      logical,save :: first=.true.
      integer :: I,J,K,L,N,SEED
      real(DOUBLE) :: dev_r,dev_z,dev_phi,ran4
      real(DOUBLE) :: annulus,dummy,limiter
      print *,first
      if (first) then
        limiter=den*phylim*10.
        dr=rof3n
        dz=zof3n
        dphi=dtheta
        TSTEP=0
        print *,rof3n
        if (restart_fluid) then
          call Fluid_Restart
        else
          if (NPART>1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(L,K,J,dev_phi,dev_z,dev_r,N,I,   &
!$OMP& annulus)                                                         &
!$OMP&  SHARED(limiter,dphi,dz,dr,pi)
!$OMP MASTER
            print *,                                                    &
     &      '->FLUID_SETUP OUT: placing fluid elements randomly.'
!$OMP END MASTER
            annulus=r(R_END)-r(R_START)
!!$OMP DO SCHEDULE(STATIC)
!$OMP MASTER
              SEED=0
              do I=0,NPART-1,1
              N=0
                DO WHILE ((N==0))
                SEED=SEED+1
                dev_r= sqrt(ran4(SEED)); SEED=SEED+1
                dev_z= ran4(SEED); SEED=SEED+1
                dev_phi= ran4(SEED)
                r_e(0,I)=dev_r*annulus+r(R_START)
                z_e(0,I)=dev_z*z(Z_TOP)
                angle_e(0,I)=dev_phi*2*pi
                J=r_e(0,I)/dr+2
                K=z_e(0,I)/dz+2
                L=angle_e(0,I)/dphi+1
                if (rho(J,K,L)>limiter) then
                  N=1
                endif
                partskip(I)=.false.
                enddo
              enddo
!$OMP END MASTER
!!$OMP END DO nowait
!$OMP END PARALLEL
          else
          endif
        endif
        print *,'->FLUID_SETUP OUT: exiting setup. '
      else
! else statement for leading if statement
      endif
      return
      endsubroutine

