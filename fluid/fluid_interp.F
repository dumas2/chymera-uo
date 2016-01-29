      subroutine fluid_interp(I,out1,out2,out3,do_rho)
      use Fluid_Mod
      implicit  none
      integer :: I,IYP,IYPP,IYM,IYMM,IZP,IZPP,IZM,IZMM
      integer,dimension(0:2) :: INT_A
      integer,dimension(0:4) :: ptix,ptiy,ptiz
      real(KIND=8) :: xx,yy,zz,ret,out1,out2,out3,zsignnm,offset
      real(KIND=8),dimension(0:2) :: pos
      real(KIND=8),dimension(0:4) :: zsign,nosign
      real(KIND=8) :: interp_project
      logical :: do_rho
! find cells
      xx=r_e(TSTEP,I)/dr
      yy=angle_e(TSTEP,I)/dphi
      zz=z_e(TSTEP,I)/dz
      INT_A(0)=xx
      INT_A(1)=yy
      INT_A(2)=zz
      pos(0)=xx
      pos(1)=yy
      pos(2)=zz
      IYPP=INT_A(1)+2
! set array indices for L direction
      IYP=IYPP-1
      IYM=IYP-2
      IYMM=IYM-1
      if (IYPP>LMAX-1) then
        IYPP=IYPP-LMAX
      endif
! take into account periodic bc
      if (IYP>LMAX-1) then
        IYP=IYP-LMAX
      endif
      if (IYMM<0) then
        IYMM=IYMM+LMAX
      endif
      if (IYM<0) then
        IYM=IYM+LMAX
      endif
      if (INT_A(2)>=KMAX) then
        INT_A(2)=INT_A(2)-2
      endif
! if at edge, back off
      IZPP=INT_A(2)+2
! set array indices for K direction
      IZP=IZPP-1
      IZM=IZP-2
      IZMM=IZM-1
      nosign=1.d0
      if (INT_A(2)==1) then
        IZMM=0
! take care of reflection symmetry
      endif
      if (INT_A(2)==0) then
        IZM=0
        IZMM=1
      endif
      ptix(0)=INT_A(0)-2
      ptix(1)=INT_A(0)-1
      ptix(2)=INT_A(0)
      ptix(3)=INT_A(0)+1
      ptix(4)=INT_A(0)+2
      ptiy(0)=IYMM
      ptiy(1)=IYM
      ptiy(2)=INT_A(1)
      ptiy(3)=IYP
      ptiy(4)=IYPP
      ptiz(0)=IZMM
      ptiz(1)=IZM
      ptiz(2)=INT_A(2)
      ptiz(3)=IZP
      ptiz(4)=IZPP
      if (do_rho) then
        ret=interp_project(ptix,ptiy,ptiz,pos,JMAX,KMAX,LMAX,rho,INT_A, &
     &  nosign,0.5,0.5,0.5)
        rho_e(TSTEP,I)=ret
        ret=interp_project(ptix,ptiy,ptiz,pos,JMAX,KMAX,LMAX,eps,INT_A, &
     &  nosign,0.5,0.5,0.5)
        eps_e(TSTEP,I)=ret
        ret=interp_project(ptix,ptiy,ptiz,pos,JMAX,KMAX,LMAX,TempK,INT_A&
     &  ,nosign,0.5,0.5,0.5)
        temp_e(TSTEP,I)=ret
        ret=interp_project(ptix,ptiy,ptiz,pos,JMAX,KMAX,LMAX,gamma1,    &
     &  INT_A,nosign,0.5,0.5,0.5)
        gamma1_e(TSTEP,I)=ret
        ret=interp_project(ptix,ptiy,ptiz,pos,JMAX,KMAX,LMAX,tau(:,:,:,1&
     &  ),INT_A,nosign,0.5,0.5,0.5)
        tau_e(TSTEP,I)=ret
      endif
      ret=interp_project(ptix,ptiy,ptiz,pos,JMAX,KMAX,LMAX,u,INT_A,     &
     &nosign,0.0,0.5,0.5)
      out1=ret
      ret=interp_project(ptix,ptiy,ptiz,pos,JMAX,KMAX,LMAX,omega,INT_A, &
     &nosign,0.5,0.5,0.5)
      out2=ret
      zsign=1.d0
      if (INT_A(2)==1) then
        IZMM=1
! take care of reflection symmetry
        zsign(0)=-1.d0
      endif
      if (INT_A(2)==0) then
        IZM=1
        IZMM=2
        zsign(1)=-1.d0
        zsign(0)=-1.d0
      endif
      ptiz(0)=IZMM
      ptiz(1)=IZM
      ptiz(2)=INT_A(2)
      ptiz(3)=IZP
      ptiz(4)=IZPP
      ret=interp_project(ptix,ptiy,ptiz,pos,JMAX,KMAX,LMAX,w,INT_A,zsign&
     &,0.5,0.5,0.0)
      out3=ret
      return
      endsubroutine
