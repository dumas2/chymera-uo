      subroutine fluidtrace(ITSTEP,ITSTOP,half)
      use fluid_mod
      implicit  none
      integer I,J,K,L,partcheck,ITSTEP,ITSTOP
      real*8 limiter,sum
      logical half
      limiter=den*phylim
      delta_t=delt
      partcheck=0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)                               &
!$OMP&  SHARED(dtheta,limiter,dz,dr,TSTEP) REDUCTION(+:partcheck)
!     allocate workspace for interp in gsl library ! not yet implemented
!      call interp_alloc(5)
!$OMP DO
        do I=0,NPART-1,1
        if (.not.partskip(I)) then
          partcheck=partcheck+1
          call fluid_interp(I,vr_e(TSTEP,I),omega_e(TSTEP,I),vz_e(TSTEP,&
     &    I),.true.)
          call fluid_advance(I,half)
          if (((r_e(TSTEP,I)<r(JMIN+2)).or.(r_e(TSTEP,I)>r(JMAX-2))).or.&
     &    (z_e(TSTEP,I)>z(KMAX))) then
            print *,' Particle',I,' went out of the ballpark.'
            partskip(I)=.true.
          elseif (z_e(TSTEP,I)<0.d0) then
            print *,                                                    &
     &      ' Oh, shit.  Something really bad happened with particle ',I
            print *,'             It went below the midplane.'
            partskip(I)=.true.
          elseif (rho_e(TSTEP,I)<=limiter) then
            print *,' For particle',I,                                  &
     &      ', code cannot handle density. rho = ',rho_e(TSTEP,I)
            partskip(I)=.true.
            print *,int(r_e(TSTEP-2,I)/dr+2),int(angle_e(TSTEP-2,I)/    &
     &      dtheta+1),int(z_e(TSTEP-2,I)/dz+2),vz_e(TSTEP-2,I)
            print *,int(r_e(TSTEP-1,I)/dr+2),int(angle_e(TSTEP-1,I)/    &
     &      dtheta+1),int(z_e(TSTEP-1,I)/dz+2),vz_e(TSTEP-1,I)
            print *,int(r_e(TSTEP,I)/dr+2),int(angle_e(TSTEP,I)/dtheta+1&
     &      ),int(z_e(TSTEP,I)/dz+2),vz_e(TSTEP,I)
          else
          endif
        endif
        enddo
!$OMP ENDDO nowait
!     free workspace used in gsl library calls ! not yet implemented
!      call interp_free
!$OMP END PARALLEL
      if(.not.half)TSTEP=TSTEP+1
      if ((TSTEP>=TSTOP-1).or.(ITSTEP==ITSTOP).and..not.half) then
        TSTEP=TSTEP-1
        call fluid_writeout(time-delt,ITSTOP,ITSTEP)
        TSTEP=TSTEP+1
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(J,I)                          &
!$OMP&  SHARED(TSTEP)
          do I=0,NPART-1,1
          r_e(0,I)=r_e(TSTEP,I)
          z_e(0,I)=z_e(TSTEP,I)
          angle_e(0,I)=angle_e(TSTEP,I)
          vz_e(0,I)=vz_e(TSTEP,I)
          vr_e(0,I)=vr_e(TSTEP,I)
          omega_e(0,I)=omega_e(TSTEP,I)
          rho_e(0,I)=rho_e(TSTEP,I)
          eps_e(0,I)=eps_e(TSTEP,I)
          temp_e(0,I)=temp_e(TSTEP,I)
          gamma1_e(0,I)=gamma1_e(TSTEP,I)
            do J=1,TSTOP-1,1
            r_e(J,I)=0.
            z_e(J,I)=0.
            angle_e(J,I)=0.
            vr_e(J,I)=0.
            vz_e(J,I)=0.
            omega_e(J,I)=0.
            rho_e(J,I)=0.
            eps_e(J,I)=0.
            temp_e(J,I)=0.
            gamma1_e(J,I)=0.
            enddo
          enddo
!$OMP END PARALLEL DO
        TSTEP=0
      endif
      if(.not.half)print *,'->FLUID_START OUT: partcheck = ',partcheck
      if (partcheck==0) then
        stop
      endif
      return
      endsubroutine

