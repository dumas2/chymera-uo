      subroutine Fluid_Restart
      use Fluid_Mod
      implicit  none
      integer :: I,N
      open(unit=945,file="fluid_restart.dat",err=911,form="UNFORMATTED")
      N=TSTEP
        do I=0,NPART-1,1
        read(945)r_e(N,I),angle_e(N,I),z_e(N,I),rho_e(N,I),temp_e(N,I), &
     &  gamma1_e(N,I),vr_e(N,I),omega_e(N,I),vz_e(N,I)
        enddo
      close(945)
      print *,'->FLUID_RESTART OUT: exiting setup. '
      return
  911 print *,                                                          &
     &'->ERROR REPORT: failed to get file in FLUID_RESTART. Stopping...'
      stop
      endsubroutine

