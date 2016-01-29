! make Toomre Q profile data
! take saved file and creates output for gnuplot.
! note: good only for fixed gamma
program qprofiler_fixed
 implicit none
!
! SWITCH: the following parameters should be self-explanatory for
! any CHYMERA user.  YMAX and XMAX represent the x-y grid for
! creating the map.  rconv converts from cells to AU.
! dx is zof3n
!
 integer :: JMAX, KMAX, LMAX, YMAX, XMAX
 integer :: ISTART, IEND, ISKIP
 integer :: II, JJ, IRR, IAN, IAN2, I, L, J, K

 real(KIND=8), parameter :: rgas = 8.254d7

 real(KIND=8), allocatable, dimension(:,:,:) :: rho,p,readdum,a,temp,gamma1
 real(KIND=8), allocatable, dimension(:) :: q,sig,kappa,omega,sound,gam,tk,qalt

 real(KIND=8) :: dr,dz,rconv,torp,dummy,average_q,innerav,outerav,count
 real(KIND=8) :: y1, y2, y3, y4, t, u, angle, dphi, pi, ir, jr, rr,time,engconv,vconv,muc

 logical :: avq = .true.

 character :: filein(3)*72, fileout*72,filenum*6,dum(16)*72

  do I = 1,16
   dum(I) = "                                                                       "
   call GetArg(I,dum(I))
  enddo 

  if ( len(trim(dum(1))) == 0 .or. &
       len(trim(dum(2))) == 0 .or. &
       len(trim(dum(3))) == 0 .or. &
       len(trim(dum(4))) == 0 .or. &
       len(trim(dum(5))) == 0 .or. &
       len(trim(dum(6))) == 0 .or. &
       len(trim(dum(7))) == 0 .or. &
       len(trim(dum(8))) == 0 .or. &
       len(trim(dum(9))) == 0 .or. &
       len(trim(dum(10))) == 0 .or. &
       len(trim(dum(11))) == 0 .or. &
       len(trim(dum(12))) == 0 .or. &
       len(trim(dum(13))) == 0 .or. &
       len(trim(dum(14)))== 0        ) then
    print *, "QPROF_F USAGE: QPROF_F {stepbegin} {stepend} {stepskip   }"
    print *, "                       {jmax} {kmax} {lmax} {torp} {mmw  }"
    print *, "                       {cell-to-size unit conversion     }"
    print *, "                       {savedfile prefix} {temperature file prefix}"
    print *, "                       {gamma1 file prefix} {fileout prefix}" 
    print *, "                       { specific energy conversion }      }"
    print *, "                       [ [ lower av limit] [high limit]  ]" 
    stop
  endif

  dum(1) = trim(dum(1)) 
  dum(2) = trim(dum(2)) 
  dum(3) = trim(dum(3)) 
  dum(4) = trim(dum(4)) 
  dum(5) = trim(dum(5)) 
  dum(6) = trim(dum(6)) 
  dum(7) = trim(dum(7)) 
  dum(8) = trim(dum(8)) 
  dum(9) = trim(dum(9)) 
  dum(14) = trim(dum(14))
  dum(15) = trim(dum(15))
  dum(16) = trim(dum(16)) 

  read(dum(1),"(i8)"   )ISTART
  read(dum(2),"(i8)"   )IEND
  read(dum(3),"(i8)"   )ISKIP
  read(dum(4),"(i8)"   )JMAX
  read(dum(5),"(i8)"   )KMAX
  read(dum(6),"(i8)"   )LMAX
  read(dum(7),"(e20.8)")torp
  read(dum(8),"(e20.8)")muc
  read(dum(9),"(e20.8)")rconv
  read(dum(14),"(e20.8)")engconv

  if (len(trim(dum(15))) == 0 ) avq = .false.
  if (avq) then
    read(dum(15),"(f15.8)")innerav
    read(dum(16),"(f15.8)")outerav
  endif

  print *, " QPLOT OUT: ISTART -> ", ISTART
  print *, " QPLOT OUT: IEND -> ", IEND
  print *, " QPLOT OUT: ISKIP -> ", ISKIP
  print *, " QPLOT OUT: JMAX -> ", JMAX
  print *, " QPLOT OUT: KMAX -> ", KMAX
  print *, " QPLOT OUT: LMAX -> ", LMAX
  print *, " QPLOT OUT: torp -> ", torp
  print *, " QPLOT OUT: rconversion -> ", rconv
  print *, " QPLOT OUT: econversion -> ", engconv
  print *, " QPLOT OUT: mmw -> ", muc
  print *, " QPLOT OUT: saved prefix -> ", trim(dum(10))
  print *, " QPLOT OUT: temp  prefix -> ", trim(dum(11))
  print *, " QPLOT OUT: gamma prefix -> ", trim(dum(12))
  print *, " QPLOT OUT: fileout prefix-> ", trim(dum(13))
  if (avq) then
    print *, " QPLOT OUT: Inner bound for fit -> ",innerav
    print *, " QPLOT OUT: Inner bound for fit -> ",outerav
  endif

  allocate(readdum(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(a(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(rho(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(temp(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(gamma1(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(p(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(q(-1:JMAX)      )
  allocate(qalt(-1:JMAX)      )
  allocate(gam(-1:JMAX)      )
  allocate(sig(-1:JMAX)      )
  allocate(kappa(-1:JMAX)      )
  allocate(omega(-1:JMAX)      )
  allocate(sound(-1:JMAX)      )
  allocate(tk(-1:JMAX)      )

 vconv = sqrt(engconv)

 I = ISTART
 do while ( I <= IEND )
    write (filenum,'(I6.6)')I
    filein(1)=trim(dum(10))//filenum//' ' 
    filein(2)=trim(dum(11))//filenum//' ' 
    filein(3)=trim(dum(12))//filenum//' ' 
    fileout=trim(dum(13))//filenum//' ' 
    print "(a,1x,a)", trim(filein(1)),trim(filein(2)),trim(filein(3)), trim(fileout) 
    OPEN(UNIT=12, FILE=trim(filein(1)),FORM='UNFORMATTED', STATUS="OLD")
    OPEN(UNIT=13, FILE=trim(filein(2)),FORM='UNFORMATTED', STATUS="OLD")
    OPEN(UNIT=14, FILE=trim(filein(3)),FORM='UNFORMATTED', STATUS="OLD")
    open(unit=15,file=trim(fileout))

 read(12)readdum
 read(12)readdum
 read(12)a
 read(12)rho
 read(12)p ! actually internal energy density
 read(12)dr,dz,dummy,time
 
 read(13)temp

 read(14)gamma1

 close(12)
 close(13)
 close(14)

 pi = acos(-1.d0)
 q = 0.
 sig = 0.
 omega = 0.
 kappa = 0.
 sound = 0.
 gam=0.
 tk = 0.

 print "(a,I9,a,I6.6,a)", " THE TIME IS ", int(time/torp), " FOR FILE ",I,"."

 do L = 0, LMAX-1
  do K = 0, KMAX-1
    do J = 0, JMAX-1
      sig(J) = sig(J)+ rho(J,K,L)*dz*2.d0/dble(LMAX) 
             ! take into account both sides of disk
      omega(J) = omega(J)+a(J,K,L)/( (dble(J)+.5d0)*dr )**2 *dz*2./dble(LMAX)
      sound(J) = sound(J)+rho(J,K,L)*sqrt(gamma1(J,K,L)*temp(J,K,L)/(muc)*rgas) &
               * dz*2.d0/(dble(LMAX))
      gam(J) = gam(J) + gamma1(J,K,L)*rho(J,K,L)*dz*2.d0/dble(LMAX)
    enddo
  enddo
 enddo

 omega(-1)=omega(0)

 do J = 0, JMAX-1
   omega(J) = omega(J)/sig(J)
   sound(J) = sound(J)/sig(J)
   gam(J) = gam(J)/sig(J)
   tk(J)= sound(J)**2/gam(J)/rgas*muc
 enddo
 do J = 0, JMAX-1
   if (omega(J) == 0.) then
    kappa(J) = 0.
   else
    kappa(J) = (dble(J)+.5d0)*omega(J)*(omega(J+1)-omega(J-1))+4.d0*omega(J)**2
    if (kappa(J) < 0.) kappa(J) = 0.
    kappa(J) = sqrt(kappa(J))
   endif
 
   q(J) = kappa(J)*sound(J)/(pi*sig(J)*vconv) ! vconv required because c is in cgs
   qalt(J) = omega(J)*sound(J)/(pi*sig(J)*vconv) 
 enddo

 if (avq) then
  count = 0.
  average_q = 0.
  do J = 0, JMAX-1
   if ( (J*1.+.5)*rconv >= innerav .and. (J*1.+.5)*rconv <= outerav) then
     average_q= average_q + q(J)
     count = count+1
   endif
  enddo
  average_q = average_q/count
  print *, " QPLOT OUT: Average Q between [innerav,outerav] -> ",average_q
 endif
 
 do J = 0, JMAX-1
   write(15,'(10(1X,1pe15.8),1X,I6)'),(J*1.+.5)*rconv,q(J),omega(J),sig(J)*8.9e6,sound(J),gam(J),temp(J,0,0),tk(J),&
        qalt(J),kappa(J),J
 enddo
 close(15)

 I = I+ISKIP
 enddo ! end while loop

 deallocate(readdum,a,rho,p,sig,omega,sound,q,temp,gamma1)

 stop
end
  

    
