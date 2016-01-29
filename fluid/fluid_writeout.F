      subroutine fluid_writeout(ttime,ITSTOP,ITSTEP)
      use fluid_mod
      implicit  none
      integer :: I,N,ITSTOP,ITSTEP
      real*8 :: ttime
      character*25 :: file
      logical,save :: first=.true.
      integer,save :: ibuf=0,j,k
      real*8, save :: rbuffer(12,100,0:NPART-1)
!90  format('## Columns: r, phi, z, rho, temp, gamma_1, eps, tau, vr, omega, vz, time, fluid element')
!100 format(2X,e22.14,1X,e22.14,1X,e22.14,1X,e22.14,1X,e22.14,1X, &
!      e22.14,1X,e22.14,1X,e22.14,1X,e22.14,1X,e24.14,1X,e24.14,1X,e24.14,1X,I6.6) 
!
!     store in buffer as file I/O is very expensive
!
      ibuf=ibuf+1
      do I=0,NPART-1
        rbuffer(1,ibuf,I)=r_e(TSTEP,I)
        rbuffer(2,ibuf,I)=angle_e(TSTEP,I)
        rbuffer(3,ibuf,I)=z_e(TSTEP,I)
        rbuffer(4,ibuf,I)=rho_e(TSTEP,I)
        rbuffer(5,ibuf,I)=temp_e(TSTEP,I)
        rbuffer(6,ibuf,I)=gamma1_e(TSTEP,I)
        rbuffer(7,ibuf,I)=eps_e(TSTEP,I)
        rbuffer(8,ibuf,I)=tau_e(TSTEP,I)
        rbuffer(9,ibuf,I)=vr_e(TSTEP,I)
        rbuffer(10,ibuf,I)=omega_e(TSTEP,I)
        rbuffer(11,ibuf,I)=vz_e(TSTEP,I)
        rbuffer(12,ibuf,I)=ttime
      enddo

      if (ibuf.eq.100.or.ITSTEP.EQ.ITSTOP) then
!       buffer is full or this is the final time step so output 
        N=TSTEP
        if (first) then
          first=.false.
            do I=0,NPART-1,1
            write(file,'(A,I6.6,A,I6.6)')'trace/trace.',I,'.',ITSTOP
            open(unit=956,file=file,form='UNFORMATTED')
            write(956)
              do j=1,ibuf
              write(956)(rbuffer(k,j,I),k=1,12),I
              enddo
            close(956)
            enddo
        else
            do I=0,NPART-1,1
            write(file,'(A,I6.6,A,I6.6)')'trace/trace.',I,'.',ITSTOP
            open(unit=956,file=file,form='UNFORMATTED',position='APPEND'&
     &      )
              do j=1,ibuf
              write(956)(rbuffer(k,j,I),k=1,12),I
              enddo
            close(956)
            enddo
        endif
!     go back to start of buffer
      ibuf=0
      endif
      return
      endsubroutine


      subroutine Fluid_WriteRestart(ITSTOP)
      use fluid_mod
      implicit  none
      integer :: I,N,ITSTOP
      character*20 :: file
      write(file,'(A,I6.6)')'fluid_restart.',ITSTOP
      open(unit=945,file=file,err=911,form="UNFORMATTED")
      N=TSTEP
        do I=0,NPART-1,1
        write(945)r_e(N,I),angle_e(N,I),z_e(N,I),rho_e(N,I),temp_e(N,I),&
     &  gamma1_e(N,I),vr_e(N,I),omega_e(N,I),vz_e(N,I)
        enddo
      close(945)
      print *,'->FLUID_WRITERESTART OUT: exiting setup. '
      return
  911 print *,                                                          &
     &'->ERROR REPORT: failed to get file in FLUID_WRITERESTART. Stop'//&
     &'ping...'
      stop
      endsubroutine


