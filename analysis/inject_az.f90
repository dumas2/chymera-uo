program chym2part
  implicit none
  
  integer, parameter :: jmax=256,kmax=64,lmax=512
  integer, parameter :: jmax2=256,kmax2=32,lmax2=1024
  integer :: j,k,l,I
  
  real*8, parameter :: msun=1.989d33,au=1.496d13,rgas=8.254d7,grav=6.67e-8
  real*8, dimension(-1:jmax,-1:kmax,0:lmax-1) :: s,t,a,rho,eps
  real*8, dimension(-1:jmax2,-1:kmax2,0:lmax2-1) :: s2,t2,a2,rho2,eps2

  real*8 :: den,dummy,dx,rholmt,tmass,rhoconv,vconv,tconv,dtheta
  real*8 :: rof3n,zof3n,delt,time,elost,ommax,sound
  real*8 :: tmassini,tmassadd,tmassout,tmassacc,totcool,totdflux,totheat,totirr,etotfl,eflufftot

  integer::jreq

  open(unit=11,file="../run/saved.150000.save",form='unformatted')

  read(11)s
  read(11)t
  read(11)a
  read(11)rho
  read(11)eps
  read(11)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
  read(11) tmassini,tmass,tmassadd,tmassout,tmassacc,totcool,totdflux,totheat,totirr,etotfl, &
              eflufftot  !ACB
   close(11)

  if (lmax2/=2*lmax)then
    print *, "At the moment, the user must specify an azimuth for lmax2 that is twice lmax."
    stop
  endif
  do L=0,LMAX-1
   do K=-1,kmax2
     do J=-1,jmax2
       s2(J,K,L*2)=s(J,K,L)
       s2(J,K,L*2+1)=s(J,K,L)
       t2(J,K,L*2)=t(J,K,L)
       t2(J,K,L*2+1)=t(J,K,L)
       a2(J,K,L*2)=a(J,K,L)
       a2(J,K,L*2+1)=a(J,K,L)
       rho2(J,K,L*2)=rho(J,K,L)
       rho2(J,K,L*2+1)=rho(J,K,L)
       eps2(J,K,L*2)=eps(J,K,L)
       eps2(J,K,L*2+1)=eps(J,K,L)
     enddo
    enddo
  enddo

   open(unit=11,file="../run/saved.interpazimuth.150000",form='unformatted')

  write(11)s2
  write(11)t2
  write(11)a2
  write(11)rho2
  write(11)eps2
  write(11)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
  write(11) tmassini,tmass,tmassadd,tmassout,tmassacc,totcool,totdflux,totheat,totirr,etotfl, &
              eflufftot  !ACB

   close(11)


  stop
end

