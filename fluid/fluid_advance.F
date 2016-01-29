      subroutine fluid_advance(i,half)
      use fluid_mod
      implicit  none
      integer :: i
      logical half
      if(half)then
       ro(I)=r_e(TSTEP,i)
       phio(I)=angle_e(TSTEP,i)
       zo(I)=z_e(TSTEP,i)
       r_e(TSTEP,i)=ro(I)+vr_e(TSTEP,I)*delta_t 
       angle_e(TSTEP,i)=phio(I)+omega_e(TSTEP,I)*delta_t
       z_e(TSTEP,i)=zo(I)+vz_e(TSTEP,I)*delta_t
       if (angle_e(TSTEP,i)>2.d0*pi) then
         angle_e(TSTEP,i)=angle_e(TSTEP,i)-2.d0*pi
       ENDIF
      else
       r_e(TSTEP+1,i)=ro(I)+vr_e(TSTEP,I)*delta_t  
       angle_e(TSTEP+1,i)=phio(I)+omega_e(TSTEP,I)*delta_t
       z_e(TSTEP+1,i)=zo(I)+vz_e(TSTEP,I)*delta_t
       if (angle_e(TSTEP+1,i)>2.d0*pi) then
         angle_e(TSTEP+1,i)=angle_e(TSTEP+1,i)-2.d0*pi
       ENDIF
       r_e(TSTEP,I)=ro(I)
       angle_e(TSTEP,I)=phio(I)
       z_e(TSTEP,I)=zo(I)
      endif
      return
      endsubroutine
