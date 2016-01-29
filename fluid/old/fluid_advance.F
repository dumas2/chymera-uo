      subroutine fluid_advance(i)
      use fluid_mod
      implicit  none
      integer :: i
      real(DOUBLE) :: dvr,domega,dvz,ro,zo,phio
      ro=r_e(TSTEP,i)
      phio=angle_e(TSTEP,i)
      zo=z_e(TSTEP,i)
      r_e(TSTEP+1,i)=ro+vr_e(TSTEP,I)*delta_t 
      angle_e(TSTEP+1,i)=phio+omega_e(TSTEP,I)*delta_t
      z_e(TSTEP+1,i)=zo+vz_e(TSTEP,I)*delta_t
      if (angle_e(TSTEP+1,i)>2.d0*pi) then
        angle_e(TSTEP+1,i)=angle_e(TSTEP+1,i)-2.d0*pi
      ENDIF
      return
      endsubroutine
