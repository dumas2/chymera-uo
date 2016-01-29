C****************************************************************************/
c     computes the hydrostatic equilibrium structures of rigidly 
c     rotating polytropes and white dwarfs - 2d.
c     to build tori, kout should be a negative integer.
c     to build spheroids, kout should be a positive integer.
c     when building polytropes, the output appears in polyout.
c     the units in polyout are the dimesionless units given in
c     hachisu, i. 1986, apjsuppl., 61, 479. the output for white
c     dwarf models is in wdout, where the units are cgs and are
c     akin to the normalizations in table 4 of hachisu, 1986.
c*****************************************************************************
c
c     compile double precision on workstations!!!
c
c     creating output file for hydro code called fort.2;
c     as written, minor changes in hydro code are required in order to
c     read the input file.  they are:
c             1.) change the dimensions of the denny and anggy arrays in hydro
c                 to denny(jmax2,kmax2) and anngy(jmax1,kmax1).  both
c                 parameters are currently included in the prmtr.h file.
c             2.) the input data are written from indices starting from
c                 j,k=1, not  2 as in the hydro code.
c             3.) note that anngy and denny are both already cell centered.
c

      program scf

      implicit real*8 (a-h,o-z)

c     prmtr.h contains the grid dimensions
 
      include "prmtr.h"

      parameter (delta=1.0e-10, rho0=1.0d0,
     1           maxit=10000, jmax1=jmax+1, kmax1=kmax+1, 
     2           jmax2=jmax+2, kmax2=kmax+2, g=6.674e-8,
     3           a=6.0e22, b=1.964e6)

      integer    case,rinit

      dimension phi(jmax2,kmax2),phi_star(jmax2,kmax2), 
     &           psi(jmax1),
     1           rhf(jmax1),rho(jmax2,kmax2), h(jmax2,kmax2), 
     2           c(maxit), h0sq(maxit), hmax(maxit), f(jmax2,kmax2),
     3           fmax(maxit),x(jmax2,kmax2),
     4           xmas(jmax1),omega(jmax1),r(jmax1)

      real*8 jr(jmax1),lgrhomax,np

      dimension rho1(jmax2,kmax2),rho2(jmax2,kmax2)

      dimension omegac(jmax1)
      dimension omegap(jmax1)
      real*8 jtot,jtotp,jtotc,jfac,jxa,jxb,jxc
      real*8 ixa,ixb,ixc,ifac
      dimension anngy(jmax1,kmax1),cen_pot(jmax1)
      dimension sumi(jmax1),sums(jmax1),
     1                 rhfpoly(jmax1)

      common /vec/ phi, rho, rhf, pi, dr
      common /jofm/ hp(2,1000)

      dimension denny(jmax2,jmax2)

c     open output files

        open(unit=19,file='rho.dat',form='unformatted')
        open(unit=94,file='anal.dat',form='formatted')
        open(unit=9,file='fact')
        open(unit=13,file='wdout',form='formatted')
        open(unit=12,file='polyout',form='formatted')
        open(unit=17,file='outfile',form='formatted')
c

      write(6,*) 'enter 1 to construct a polytrope or'
      write(6,*) 'enter 2 to construct a white dwarf:'
      read(5,*) case
c
c     to start building models start with a fairly spherical model 
c     and use a density array with all 1's. it turns out it really doesn't 
c     matter. we tried initial arrays with 1/r^2 and others and the code 
c     still convereges pretty quickly. after you have a rho.dat file you can 
c     use successive files to build the flatter models. if you want to save 
c     the density arrays you need to rename rho.dat each time.
c
      write(*,*)'enter 1 to read old rho.dat'
      write(*,*)'enter 0 to use all 1.0'
      read(*,*)rinit
      if (case.eq.1) then
        write(6,*) 'enter polytropic index, np,jold,kmax,jout,kout'
        write(6,*) 'and log maximum density (in cgs units),del:  '
        read(5,*) pin,np,jold,kold,jout,kout,lgrhomax,del,starnot
        starm    = starnot
        rhomax   = 10.0**lgrhomax 
c
c     np is nprime (for  infinity, enter a value .gt. 3.0);
c     jold is the previous jout;
c     jout and kout set the axis ratio. note that the ratio is 
c       (kout-2)/(jout-2); 
c     the log max dens we use 0;
c     del is the percentage of the previous iteration to keep. The flatter
c       flatter you get the closer to 1.o del needs to be.  Typically, 0.8
c       or 0.9 works well. 
c
      else
        write(6,*) 'enter jout, kout and '
        write(6,*) 'log maximum density (in cgs units):'
        read(5,*) jout,kout,lgrhomax
        rhomax = 10.0**lgrhomax 
        xmax = (rhomax/b)**(1.0/3.0)
      end if     

      do j_model =  1 , 1

      rb = (kout/abs(kout))*float(abs(kout)-2)/float(jout-2)
      pi = acos(-1.0d0)

      koutest = kout
      kin = 2
      if (koutest.gt.0) then
        jin = 2
      else
        jin = -kout 
        kout = (jout-jin)
      end if

c     dr is the grid spacing
 
      dr = 1.0/(float(jout)-2.0)

c     rhf is the radial coordinate at the grid center

         rhf(1) = -0.5*dr
         do j=2,jmax1
            rhf(j) = rhf(j-1) + dr 
         end do

         r(2) = 0.d0
         do j = 2,jout-1
            r(j+1) = r(j) + dr
         end do

      if(j_model.eq.1) then 

c     zeros out the density on the entire grid

      do 20 k=1,kmax2
        do 20 j=1,jmax2
          rho(j,k) = 0.0
 20     continue
c
c  modification to make mass of the * + disk = 1
c
      cff = 4.0*pi*dr*dr 
c
c
c
c       sets the density of the model to the initial guess
c       and calculates the total mass, tmss

        if(rinit .eq. 0) then
          tmss = 0.d0
          do 40 k=kin,kmax
            rhoz = rho0
            do 40 j=jin,jout
              rho(j,k) = rhoz*(rhf(j)/rhf(jout))**2
              tmss = tmss + rhf(j)*rho(j,k)
 40       continue
          tmss = cff*tmss
        else
          tmss = 0.d0
            do j = 2,jold
              do k = 2,kold
                read(19) rho(j,k)
                tmss = tmss + rhf(j)*rho(j,k)
              end do
            end do
            tmss = cff*tmss
 1111       format(1x,e16.8)
 1112       continue 
        end if

      end if
c
c     Big iteration loop begins here.
c
c     Construct model based on specific angular momentum j(m)
c
c     psi is the effective centrifugal potential

      do 30 n=1,maxit

         do j = 2,jmax2
            do k = 2,kmax2
               rho1(j,k) = rho(j,k)
            end do
         end do

c       the subroutines called here determine the gravitational
c       potential phi by solving poisson's equation

         jout1=jout-1
         jin1=jin-1
         kout1=kout-1

         call jhelp(jout,kout,jin,rho,rhf,dr,psi,cen_pot,xmas,
     &             jr,omega,psia,psib,r,pi,tmss,
     &             jmax2,jmax1,kmax2,np)

         dpsiold = dpsi
         dpsi = psia - psib

         call bndry
         call pot2d

         starm = starnot*tmss

         do j = 2,jmax1
           do k = 2,kmax1
             radius = sqrt (rhf(j)**2 + rhf(k)**2)
             phi_star(j,k) = - starm / radius
             phi(j,k)      =  phi(j,k) + phi_star(j,k)
           end do
         end do

c       phia and phib are the boundary values of phi

        phia = 0.5*(phi(jout,2) + phi(jout1,2))
        if (koutest.gt.0) then
          phib = 0.5*(phi(2,kout) + phi(2,kout1))
        else
          phib = 0.5*(phi(jin,2) + phi(jin1,2))
        end if

        dphiold = dphi
        dphi = phia - phib

c       for rigid rotation, h0sq is the square of the angular velocity

        h0sq(n) = -dphi/dpsi
        if(n.gt.10000) then 
          h0sqold = -dphiold/dpsiold
          h0sq(n) = h0sqold + 0.1*rand(0)*(h0sq(n)-h0sqold)
        end if

c       c is the integration constant

        if (case.eq.1) then
          c(n) = phia + h0sq(n)*psia
        else
          do 43 k=kin,kmax1
            do 43 j=jin,jout
              f(j,k) = - phi(j,k) - h0sq(n)*psi(j)
   43     continue
          fmax(n) = f(jin,kin)
          do 45 k=kin,kmax1
            do 45 j=jin,jout
              if (fmax(n).lt.f(j,k)) fmax(n) = f(j,k)
   45     continue
          c(n) = (fmax(n) - ((1.0 + xmax**2.0)**0.5)
     1           *f(jout,2))/((1.0 + xmax**2.0)**0.5
     2           -1.0)
        end if

c       calculates the enthalpy h and updates the density and the
c       total mass

        do 50 k=kin,kmax2
          do 50 j=jin,jout
            h(j,k) = c(n) - phi(j,k) - h0sq(n)*psi(j)           
   50   continue

c       finds the maximum value of the enthalpy 

        hmax(n) = h(jin,kin)
        do 60 k=kin,kmax1
          do 60 j=jin,jout
            if (hmax(n).lt.h(j,k)) hmax(n) = h(j,k)
   60   continue

        tmss = 0.0
        do 65 k=kin,kmax1
          do 65 j=jin,jout
            if (h(j,k).gt.0.0) then
              if (case.eq.1) then
                rho(j,k) = rhomax*(h(j,k)/hmax(n))**pin
              else
                x(j,k) = (1.0 + xmax**2.0)*
     1            ((h(j,k)/hmax(n))**2.0) - 1.0
                if (x(j,k).lt.0.0) then
                  x(j,k) = 0.0
                else
                  x(j,k) = x(j,k)**0.5
                endif
                rho(j,k) = (x(j,k)/xmax)**3.0
              endif         
              tmss = tmss + rhf(j)*rho(j,k)
            else
              rho(j,k) = 0.0
            end if
   65   continue    
        tmss = cff*tmss

         do j = 2,jmax2
            do k = 2,kmax2
              rho2(j,k) = rho(j,k)
            end do
         end do

         do j = 2,jmax2
           do k =2,kmax2
             rho(j,k) = rho1(j,k) + (1.0-del)*(rho2(j,k) - rho1(j,k))
           end do
         end do

        if (n.eq.1)  goto 30

c       checks for convergence of hmax, h0sq, and c

        hmaxd = abs((hmax(n) - hmax(n-1))/hmax(n))
        h0sqd = abs((h0sq(n) - h0sq(n-1))/h0sq(n))
        cd    = abs((c(n)    - c(n-1)   )/c(n))

        write(*,*)'hmax(n) =',hmax(n)
        write(*,*)'h0sq(n) =',h0sq(n)
        write(*,*)'c(n)    =',c(n)
        write(*,*)' '
        write(*,*)'dh      =',hmaxd
        write(*,*)'dh02    =',h0sqd
        write(*,*)'cd      =',cd
        write(*,*)' '

        if (hmaxd.le.delta) then
          if (h0sqd.le.delta) then
            if (cd.le.delta) then
              write(6,*) 'converged after ',n,' iterations'
              nfin = n
              goto 35
            end if
          end if
        else
          if (n.eq.maxit) write(6,*) 'iteration did'
     1    ,' not converge'
          nfin =  maxit
        end if
   30 continue

c     ss is the thermal energy, tt is the rotational kinetic
c     energy, and ww is the gravitational energy

   35 ss = 0.0
      tt = 0.0
      ww = 0.0
      ws = 0.0 
      if ((case.eq.1).and.(pin.ne.0.0)) then
        re = ((1.0+pin)*(1.0/g)
     1      *(rhomax**(-1.0+1.0/pin))/hmax(nfin))**0.5
      else
        re = (8.0*a*(1.0 + xmax**2.0)**0.5/(b*g*rhomax*
     1       hmax(nfin)))**0.5
      end if
        do 70 k=kin,kmax2
        do 70 j=jin,jout
          dm = rho(j,k)
          if (dm.gt.0.0) then
            dm = dm*rhf(j)
            if (case.eq.1) then 
              ss = ss + dm*h(j,k)
            else
              ss = ss + (x(j,k)*((2.0*x(j,k)**2.0)-3.0)
     1         *((x(j,k)**2.0+1.0)**0.5)+3.0*log(x(j,k)+
     2         sqrt(1.0+x(j,k)**2.0)))*rhf(j)
            end if 
            tt = tt + dm*(omega(j)*rhf(j))**2.d0
            ww = ww + dm*(phi(j,k)-phi_star(j,k))
            ws = ws + dm*phi_star(j,k)
          end if
   70 continue
      oc = 1.0/(1.0+pin)
      enrm = g*(re/1.0e10)**5.0*rhomax**2.0
      if (case.eq.1) then
        ss = cff*ss*oc
      else
        ss = cff*ss*a/(g*re**2.0*rhomax**2.0)
      end if
      tt = 0.5*cff*h0sq(nfin)*tt
      ww = 0.5*cff*ww
      ws =     cff*ws
      aww = -1.0/(ww+ws)
      sow = 1.5*ss*aww
      tow = tt*aww

c     vc is the virial correction, a measure of the numerical
c     accuracy of the solution

      vc = -1.0/(ww+ws)*abs(2.0*tt + (ww+ws) + 3.0*ss)

c     write converged model to fort.71

        rewind(71)
        do j=2,jold
          do k=2,kold
            write(71) rho(j,k)
          end do
        end do

c     converts to cgs units

      hmax(nfin) = hmax(nfin)*g*(re**2.0)*rhomax
      do 80 k=kin,kmax2
        do 80 j=jin,jout
          rho(j,k) = rhomax*rho(j,k)
          phi(j,k) = phi(j,k)*g*(re**2.0)*rhomax
          phi_star(j,k) = phi_star(j,k)*g*(re**2.0)*rhomax
          h(j,k) = h(j,k)*g*(re**2.0)*rhomax
   80 continue

      if (case.eq.2) then
        tmss = tmss*((re/1.258e11)**3.0)*rhomax
        ww = ww*enrm
        ws = ws*enrm
        tt = tt*enrm
        ss = ss*enrm
      end if

      if (case.eq.1) then

      write(12,101) pin
      write(12,1001) np
      write(12,*) 'converged after ', nfin,' iterations'
      write(12,100) hmaxd,h0sqd,cd,
     1              sow,tow,vc
      write(12,110) rb,h0sq(nfin),tmss,starm,
     &              tt,-ww,-ws,3.0*ss,starm/tmss

      pi = acos(-1.0)
      do j=2,jmax1
         if(rho(j,2).gt.0.0) then
           column=(xmas(j)-xmas(j-1))/(r(j)**2-r(j-1)**2)/pi
         else
           column=0.0
         end if
         write(17,99)j,rhf(j),jr(j),omega(j),xmas(j),rho(j,2),column
      end do

      else
        write(13,*)
        write(13,102) lgrhomax
        write(13,*) ' converged after ', nfin,' iterations'
        write(13,*)
        write(13,100) hmaxd,h0sqd,cd,sow,tow,vc
        write(13,111) rb,tmss,starm,
     &                tt,-ww,-ws,3.*ss,re/1.0e8
      end if

   99 format(1x,i3,1x,1p6e10.3)
  100 format(1x,'relative diff. in hmax  :  ',e12.4,/,
     1  1x, 'relative diff. in h0sq  :  ',e12.4,/,
     2  1x, 'relative diff. in c     :  ',e12.4,/,
     3  1x, '1.5*s/|w|               :  ',e12.4,/,
     4  1x, 't/|w|                   :  ',e12.4,/,
     5  1x, 'virial error            :  ',e12.4,/)
  101 format(1x,'polytropic index  :  ',f5.1,/)
 1001 format(1x,'n prime :  ',f8.1,/)
  102 format(1x,'rhomax:  10^',f5.1,/)
  110 format(1x,'axis ratio:  ',f6.3,/,
     1       1x,'omega0^2:  ',1pe12.4,/,
     2       1x,'mass:  ',1p2e12.4,/,
     3       1x,'t:  ',1pe12.4,/,
     4       1x,'-w:  ',1pe12.4,/,
     4       1x,'-ws  ',1pe12.4,/,
     5       1x,'3.0*s:  ',1pe12.4,/,
     6       1x,'star/disk:   ',1pe12.4,/) 
  111 format(1x,'axis ratio  :  ',f6.3,/,
     1       1x,'mass:  ',1p2e12.4,/,
     2       1x,'t   :  ',1pe12.4,/,
     3       1x,'-w  :  ',1pe12.4,/,
     4       1x,'3*s :  ',1pe12.4,/,
     5       1x,'re  :  ',1pe12.4)      
 1113 format(1x,e16.8)


c*************************************************************************
c     output in polytrope units (k=1,m=1,g=1)                            *
c*************************************************************************

c     xk is polytropic constant

      gamma = 1.d0 + 1.d0/pin

      rhomax = (hmax(nfin)/(1+pin))**(1.d0/(gamma-1.d0))
      xk = (1.d0/(1.d0+pin))*h(jin,kin)*(rho(jin,kin)**(-1.d0/pin))
      xnum = (3.d0 * gamma) - 4.d0

c     rotational inertia and total angular momentum

      do j =2,jmax1
         dum = 0.d0
         do k =2,kmax1
            dum = dum + rho(j,k)*(rhf(j)**3.d0)
         end do
         sumi(j) = dum
         sums(j) = dum*omega(j)*sqrt(h0sq(nfin))
      end do

      dummy1 = 0.d0
      dummy2 = 0.d0
      do  j= 2,jmax1
            dummy1 = dummy1 + sumi(j)
            dummy2 = dummy2 + sums(j)
      end do
      roti = cff*dummy1
      jtot = cff*dummy2

c*************************************************************************
c     convert quantites to cgs units. h, rho, phi already done           
c*************************************************************************

      tmcgs = tmss*(re**3.d0)*rhomax

      dummy = h0sq(nfin)
      do j=2,jmax1
         omegac(j) = sqrt(dummy)*omega(j)*(g**0.5d0)*(rhomax**0.5d0)
      enddo

      do j =2,jmax1
         r(j) = r(j)*re
      enddo

      ttcgs = tt * g * (re**5.d0) * (rhomax**2.d0)
      wwcgs = ww * g * (re**5.d0) * (rhomax**2.d0)
      sscgs = ss * g * (re**5.d0) * (rhomax**2.d0)
      recgs = r(jout)*re
      drcgs = dr*re
      rotic = roti *(re**5.d0)*(rhomax)
      h0sqc = h0sq(nfin) * g * (re**5.d0) * (rhomax**2.d0)
      jtotc = jtot * (g**0.5d0) * (re**5.d0) * (rhomax**(3.d0/2.d0))

c     factor for rho, h

      rfac = (tmcgs**2.d0*g**3.d0*xk**(-3.d0))**(1.d0/xnum)
      hfac = g*(re**2.d0)*rhomax

      do j=2,jmax2
         do k=2,kmax2
            rho(j,k) = rho(j,k)/rfac
            h(j,k) = h(j,k)/hfac
         end do
      end do

c     interpolate to get rhocen

      rslope = (rho(3,2) - rho(2,2))/(rhf(3) - rhf(2))
      br = rho(2,2) - rslope*rhf(2)

c     zero out rho's not considered in hydro code

      rhocut = 0.0
      do j = 2,jmax2
        if(rho(j,2).ge.rhocut) rhocut = rho(j,2)
      end do

      do j = 2,jmax2
         do k = 2,kmax2
           if(rho(j,k) .lt. (rhocut*1.e-10)) rho(j,k) = 0.d0
         end do
      end do

c     factor for omega

      opxa = 1.d0/xnum
      opxb = (3.d0*gamma-1.d0)/(2.d0*xnum)
      opxc = (-3.d0/2.d0)/xnum
      ofac = (tmcgs**opxa)*(g**opxb)*(xk**opxc)

      omegap(1) = 0.d0
      do j=2,jmax1
         omegap(j) = omegac(j)/ofac  
         if(omegap(j) .gt. omegap(j-1)) then
            omax = omegap(j)
            omaxj = j
         endif   
      enddo


c     interpolate to get omegacen
c     should really set to zero for protostars

      oslope = (omegap(3) - omegap(2))/(rhf(3) - rhf(2))
      ob = omegap(2) - oslope*rhf(2)

c     factor for phi

      pxa = (2.d0*(gamma-1.d0))/xnum
      pxb = (3.d0*(gamma-1.d0))/xnum
      pxc = -1.d0/xnum
      pfac = (tmcgs**pxa)*(g**pxb)*(xk**pxc)

      do j =2,jmax1
         phi(j,2) = phi(j,2)/pfac
         phi_star(j,2) = phi_star(j,2)/pfac
      enddo

c     factor for ww,tt,ss,h0sq

      wxa = (5.d0*gamma-6.d0)/xnum
      wxb = (3.d0*gamma-3.d0)/xnum
      wxc = -1.d0/xnum
      wfac = (tmcgs**wxa)*(g**wxb)*(xk**wxc)
      ttpoly = ttcgs/wfac
      wwpoly = wwcgs/wfac
      sspoly = sscgs/wfac
      upoly = 1.5d0*sspoly
      etot = ttpoly+wwpoly+upoly
      h0sqp = h0sqc/wfac

c     factor for radius

      rxa = (gamma-2.d0)/xnum
      rxb = -1.d0/xnum
      rxc = 1.d0/xnum
      radfac = (tmcgs**rxa)*(g**rxb)*(xk**rxc)
      repoly = recgs/radfac
      drpoly = drcgs/radfac

      do j = 2,jmax
         rhfpoly(j) = drpoly*(j-1.5d0)
         r(j) = r(j)/radfac
      enddo

c     factor for roti

      ixa = (5.d0*gamma-8.d0)/xnum
      ixb = -2.d0/xnum
      ixc = 2.d0/xnum
      ifac = (tmcgs**ixa)*(g**ixb)*(xk**ixc)
      rotip = rotic/ifac

c     factor for jtot

      jxa = (5.d0*gamma-7.d0)/xnum
      jxb = (3.d0*gamma-5.d0)/(2.d0*xnum)
      jxc = 1.d0/(2.d0*xnum)
      jfac = (tmcgs**jxa)*(g**jxb)*(xk**jxc)

      jtotp = jtotc/jfac

      rat = upoly/abs(wwpoly)
      
      write(9,12)'its       =',nfin        !# of iterations
      write(9,11)'pin       =',pin         !polytropic index
      write(9,11)'nprime    =',np          !n'
      write(9,12)'jmax      =',jmax        !max extent of hach grid in r
      write(9,12)'kmax      =',kmax        !"  "  "  "  "  "  "  "  "  z
      write(9,12)'jout      =',jout        !jreq for hydro
      write(9,12)'kout      =',kout        !kzpol for hydro
      write(9,11)'omega max =',omax  
      write(9,11)'omaxj     =',omaxj       !where the max occurs in j
      write(9,11)'omegacen  =',ob
      write(9,11)'rhocen    =',br
      write(9,11)'rhomax    =',rhomax
      write(9,11)'xk        =',xk          !polytropic constant k
      write(9,11)'wwpoly    =',wwpoly      !w in pu's
      write(9,11)'ttpoly    =',ttpoly      !t "  "
      write(9,11)'sspoly    =',sspoly      !s "  "
      write(9,11)'u/|w|     =',rat
      write(9,11)'t/|w|     =',tow
      write(9,11)'etot      =',etot        !total energy 
      write(9,11)'ve        =',vc          !virial test
      write(9,11)'repoly    =',r(jout)     !r equatorial in pu's
      write(9,11)'re/rp     =',1.d0/rb
      write(9,11)'jtotp     =',jtotp       !total spec ang mom in pu's
      write(9,11)'jtotc     =',jtotc
      write(9,11)'jfac      =',jfac
      write(9,11)'jtot      =',jtot
      write(9,11)'jtoth     =',jtot        !"  "  "  "  "  "   in hach units
      write(9,11)'rotip     =',rotip       !moment of inertia in pu's
      write(9,11)'rotih     =',roti        !" " " " " " " "   in hach units
      write(9,11)'rhf(jout) =',rhf(jout)
      write(9,11)'drpoly    =',drpoly      !rof3n

 11   format(1x,a15,1x,e14.7)
 12   format(1x,a15,1x,i8)

c     calculate the anngy array

      rhomax=0.0
      do j = 2,jmax1
         if(rho(j,2).gt.rhomax) then
           rhomax=rho(j,2)
           jrhomax=j
         end if
         do k = 2,kmax1
            anngy(j,k) = omegap(j)*rhfpoly(j)**2
         end do
      end do
c
c
      gamma=(5.0/3.0)
      deltaz=drpoly
      do j=2,jmax1
        epi=sqrt(2.0*(2.0+np))*omegap(j)
        radnew=0.5*(r(j)+r(j+1))
        frequency=anngy(j,2)/radnew**2
        epi1=2.0*(frequency/radnew)*
     &   ((anngy(j+1,2)-anngy(j,2))/drpoly)
        surface=0.0
        do k=2,kmax1
          surface=surface+2.0*rho(j,k)*deltaz
        end do
        csq=(gamma*surface**(gamma-1.0))
        if(surface.ne.0.0) then
          toomre=sqrt(csq)*epi/(3.14*surface)
          toomre1=sqrt(csq)*sqrt(epi1)/(3.14*surface)
          vortensity=(2.0+np)*omegap(j)/surface
          write(6,1061) j,rhfpoly(j),toomre,toomre1,epi,sqrt(epi1),
     &                   epi/sqrt(epi1)
 1061     format(i6,1p6e12.4)
        end if
      end do
c
c
      rewind(2)
      write(2,*) jout,kout
      write(2,*) r(jout),starm/tmss
      write(2,96) ((rho(j,k),j=2,jmax1),k=2,kmax1)
      write(2,96) (anngy(j,2),j=2,jmax1)
 94   format(2i5)
 95   format(1pe16.7)
 96   format(1p6e13.5)
      write(*,*) ' ' 
      write(13,*) 'dr,drpoly,drcgs,re,repoly,recgs,jrhomax,r(jrhomax)'

      write(13,*) dr,drpoly,drcgs,re,repoly,recgs,jrhomax,r(jrhomax)

      rewind(93)
      ainewz = starm/tmss                          ! used as input for chymera
      write(93,195) pin,jtotp,r(jout),omegap(jrhomax),rho(jrhomax,2)
     &            ,tow,drpoly,drpoly,ainewz,jout,kout
      write(93,196) rho
      write(93,196) anngy
      
      rewind(47)
      do k = 2,kmax
        do j = 2,jmax 
          write(47,499) j,k,rho(j,k)
        end do
        write(47,*)
      end do

 499  format(2i5,1p1e12.5)

 195  FORMAT(3X,1PE22.15,2X,8(1PE22.15,2X),2I4)
 196  FORMAT(8(1PE22.15,2X))

c
c  Self-gravity parameters (Andalib, Tohline, & Christodoulou 1997)
c
        pi     =  acos(-1.0d0)
        r_not  =  r(jrhomax)
c    Kepler frequency
        womega = sqrt(starm/tmss/r_not**3) 
c    Frequency at density maximum
        wpeak  =  omegap(jrhomax)
c    "Kepler" disk parameter
        eta    = (womega/wpeak)**2
c    Disk thickness parameter, epsilon
        epsilon= (r(jout)-r(jin))
        g2     = (2.0*pi)*epsilon**2
        p2     = (2.0*pi)*r_not**2
c    Jeans frequency
        wjeans =  2.0*sqrt(pi*rho(jrhomax,2))
c    Disk fragmentation parameter, p_J
        p2a    = (wjeans/wpeak)**2

        write(41,411) starm/tmss,pin,np
        write(41,*) 'Self-gravity parameters'
        write(41,*) 'eta,p,g,epsilon,r_-/r_+,t/|W|: '
        write(41,412) eta,sqrt(p2),sqrt(g2),epsilon/r(jout),
     &                r(jin)/r(jout),tow
        write(42,413) starm/tmss,eta,sqrt(p2a),sqrt(g2),epsilon/r(jout),
     &                r(jin)/r(jout),tow
        write(41,412) eta,sqrt(p2a),sqrt(g2),epsilon/r(jout),
     &                r(jin)/r(jout),tow
        write(41,*) ' '
 411    format(' M,n,q:',1p6e11.3)
 412    format(' eta,p,g:',1p6e11.3)
 413    format(1p7e11.3)
c
c
c     output for equilibrium analysis
c     we've used this data  to locate linblad resonances, calculate
c     toomre's q, keplerian rotation, etc. and to make plots of
c     the equilibrium models.  glen is working on a verison that
c     incororates sm directly for automatic plotting.

      do j = 2,jout
        write(94,1114)j,r(j),rho(j,2),xmas(j),phi(j,2),phi_star(j,2),
     &                  omegap(j)
      enddo

 1114 format(1x,i4,1x,e10.3,1x,e12.5,1x,e10.3,1x,e12.5,1x,e12.5,1x,
     7                e10.3)

        jold = jout
c_torus        kold = kout
c_torus        kout = koutest

        if(j_model.le.20) then
          kout = kout - 8
        else
c_torus         kout = kout - 2
        end if

 7100    continue

      end do

      end


      subroutine jhelp(jout,kout,jin,rho,rhf,dr,psi,cen_pot,xmas,
     &                 jr,omega,psia,psib,r,pi,tmss,
     &                 jmax2,jmax1,kmax2,np)

      implicit real*8 (a-h,o-z)

c     calculate omegas for the chosen j(m)

      dimension xmas(jmax1),omega(jmax1),
     &       psi(jmax1),cen_pot(jmax1)
      real*8 jr(jmax1)
      dimension r(jmax1)
      dimension rho(jmax2,kmax2),rhf(jmax1)
      real*8 np

c     zero out mass fraction array

      do j = 2,jout
         xmas(j) = 0.0
      end do

      sum1 = 0.d0
      cf = 2.0*pi*dr

c     calculate the mass fraction for 1st radius, leave off 2*pi*dr

      sum = 0.0
      vol1 = rhf(2)**2.0
            do k =2,kmax2
               sum = sum + rho(2,k)*vol1
            end do
         xmas(2) = sum

c     calc the rest of the mr by adding up shell masses. no 2*pi*dr

         sumz = 0.0
         do j = 3,jout-1
               do k = 2,kmax2
            sumz = sumz + rho(j-1,k) * (r(j)**2.d0 - rhf(j-1)**2.d0)+
     &            rho(j,k)* (rhf(j)**2.d0 - r(j)**2.d0)
               end do
            xmas(j) = sumz + xmas(j-1) 
            sumz = 0.d0
         end do

c     last 1/2 cell to get total mass. no 2*pi*dr

         tms = 0.d0
         emas = 0.d0
         do k = 2,kmax2
            emas = emas + rho(jout-1,k)*(r(jout)**2.d0 - 
     &           rhf(jout-1)**2.d0)
         enddo

         tms = emas + xmas(jout-1)
         xmas(jout)=tms

c     normalize the mass

      xnorm = 1.0d0/tms
      do j =2,jout
         xmas(j) = xmas(j)*xnorm
      enddo

      if (np .gt. 5.0d0) then
        call inomega(np,omega,xmas,rhf,jout)
      elseif (np .le. 5.0d0) then
        call jfunc(np,xmas,jmax2,jmax1,rhf,jout,omega,jr)
      endif

c     integrate to get psi(j)

      sum1 = -0.5*dr*(omega(2)**2.0*rhf(2))
      power = 2.0
      czero = 0.5*(omega(2)*rhf(2))**2/(power-1.0)

      do j =3,jout
         sum1 = sum1 - 0.5*dr*((omega(j-1)**2.0)*rhf(j-1) +
     1                          (omega(j)**2.0)*rhf(j))
         psi(j) = sum1
         cen_pot(j) = 0.5*(omega(j)*rhf(j))**2/(power-1.0)
     &                       - czero
      enddo

c     boundary values for psia and psib

      psib = 0.d0
      psib = -(omega(jin)**2.0*rhf(jin))*dr*0.5d0
     &           + psi(jin-1)

      psia = -(2.5d0**2.d0)*dr*0.5d0 
      psia = -(omega(jout)**2.d0*rhf(jout))*dr*0.5d0
     &           + psi(jout-1)

      psi(jout) = psia

      return
      end

      subroutine jfunc(np,xmas,jmax2,jmax1,rhf,jout,omega,jr)

      implicit real*8 (a-h,o-z)

      dimension xmas(jmax1),omega(jmax1),rhf(jmax1)
      real*8 np,jr(jmax1)

      do j = 2,jout
         q = xmas(j)
         rad = rhf(j)
         jr(j) = f(q,np,rad)
c
c   Specific angular momentum distribution described in Pickett et al.
c
cpickett         diskfraction = 0.10
cpickett         jr(j)=((q+(1.0-diskfraction))/(1.0+(1.0-diskfraction)))**2
c
         omega(j) = jr(j)/rad**2
         if(np.lt.0.0) then 
           omega(j) = rad**np
           jr(j)    = rad**(2+np)
         end if
      end do
c
      return
c
      xcut   = 0.4
      xwidth = 0.02
      xnorm  = 2.5e-01

      do j =2,jout
        q = xmas(j)
        if (q.le.xcut) jcut=j
      end do

      do j =2,jout
         q = xmas(j)
         rad = rhf(j)
         fermi = exp(-(q-xcut)/xwidth)
         omega(j) = xnorm
         omega(j) = omega(j)+rhf(jcut)**2/(fermi+rad**2)
      end do


      return
      end

      function f(x,np,pomega)

      implicit real*8 (a-h,o-z)
      real*8 np
      common/fermi/xcut,width,transition
      data xcut,width/0.75,0.025/

      if (np .eq. 0.0) then
         f = 2.5d0 -2.5d0*(1-x)**(2.d0/3.d0)

      elseif (np .eq. 0.5) then
         f = 3.068133d0+0.203667d0*(1-x)**0.801297d0
     &           -3.271800d0*(1-x)**0.5d0

      elseif (np .eq. 1.0) then
         f = 3.825819d0+0.857311d0*(1-x)**0.650981d0
     &           -4.68313d0*(1-x)**0.4d0

      elseif (np .eq. 1.5) then
         f = 4.887588d0+2.345310d0*(1-x)**0.525816d0
     &           -7.232898d0*(1-x)**(1.d0/3.d0)

      elseif (np .eq. 2.0) then
          f = 6.457897d0+6.018111d0*(1-x)**0.417472d0
     &           -12.476007d0*(1-x)**0.285714d0

      elseif (np .eq. 2.5) then
          f = 8.944150d0+18.234305d0*(1-x)**0.321459d0
     &           -27.178455d0*(1-x)**0.25d0

      elseif (np .eq. 3.0) then
          xcut = 0.99
          x    = xcut*x
          f = 13.270061d0+163.26149d0*(1-x)**0.235287d0
     &           -176.53154d0*(1-x)**0.222222d0
          f = f - 0.000011

      elseif (np .eq. 5.0) then
          starm = 50.0
          f = ((x+starm)/(1.0+starm))

      end if      

      return
      end


      subroutine inomega(prime,omega,xmas,rhf,jout)
      implicit real*8 (a-h,o-z)
      include "prmtr.h"
      
      parameter (jmax2=jmax+2,jmax1=jmax+1,itmax=500, zzz=1.0e-6)
      dimension val(jmax2)
      dimension omega(jmax1),xmas(jmax1),rhf(jmax1)
c
c     parametrized function for the n'=infinity case.  see shelby
c     yang's dissertaion or the '89 paper for a discussion.
c      
      fin(x,y)= x + (1.d0 - sqrt(1.d0 - x**2.d0) - y)/acos(x)
c
c...zero out arrays
c
      do 3 j=2,jmax1
      	omega(j)=0.d0
        val(j)=0.d0
 3    continue

c
c...main iteration loop.  converge j(m(ksi)) and spit out
c...the omega's. 
c
      do 10 j=2,jout-1

         q = xmas(j)
         guess = rhf(j)/rhf(jout)
         ii = 1
         val(1) = fin(guess,q)

         do 12 ii = 2, itmax
           val(ii) = val(ii-1) - fin(val(ii-1),q)
           check = -val(ii)*(val(ii)-1.d0)
           if(check.lt.0.d0) then 
              write(6,*) 'ksi out of range=',val(ii)
              write(6,*) 'step # =',ii
              write(6,*) 'guess =',guess
c_stop        stop
           endif

           diff = abs(val(ii) - val(ii-1))
           if(diff.le.zzz) go to 100
 12      continue

         write(6,*) 'ksi failed @',ii,j,diff
         stop

 100     omega(j)=4.5d0*(val(ii)/rhf(j))**2.d0

 10   continue

 999  return
      end




      subroutine bndry

      implicit real*8 (a-h,o-z)

      include "prmtr.h"
      parameter     (c0 = 0.5, c1 = 1.5, c2 = 0.375, c3 = 3.75,
     1		     c4 = 4.375, c5 = 0.3125, c6 = 6.5625,
     2		     c7 = 19.6875, c8 = 14.4375, c9 = 0.2734375,
     3		     ca = 9.84375, cb = 54.140625, cc = 93.84375,
     4		     cd = 50.2734375, ce = 0.2460938, cf = 13.5351563,
     5               cg = 117.3046875, ch = 351.9140625,
     6               ci = 427.3242188, cj = 180.4257813,
     7		     jmax1 = jmax + 1, kmax1 = kmax + 1,
     8		     jmax2 = jmax + 2, kmax2 = kmax + 2)
      dimension phi(jmax2,kmax2), rho(jmax2,kmax2), rhf(jmax1)
      common /vec /  phi, rho, rhf, pi, dr

      pg2 = -(pi + pi)
      dd = 2.0 * dr * dr
      q0 = 0.0
      q2 = 0.0
      q4 = 0.0
      q6 = 0.0
      q8 = 0.0
      qa = 0.0
      do   100   k = 2, kmax1
      zz = rhf(k)
      z2 = zz * zz
      do   100   j = 2, jmax1
      if(rho(j,k) .gt. 0.0) then
      rr = rhf(j)
      s2 = rr * rr + z2
      x2 = z2 / s2
      p2 = x2 * c1 - c0
      p4 = x2 * (x2 * c4 - c3) + c2
      p6 = x2 * (x2 * (x2 * c8 - c7) + c6) - c5
      p8 = x2 * (x2 * (x2 * (x2 * cd - cc) + cb) - ca) + c9
      pa = x2 * (x2 * (x2 * (x2 * (x2 * cj - ci) + ch) - cg) + cf) - ce
      dm = rr * rho(j,k)
      sm = s2 * dm
      s4 = s2 * sm
      s6 = s2 * s4
      s8 = s2 * s6
      q0 = q0 + dm
      q2 = q2 + p2 * sm
      q4 = q4 + p4 * s4
      q6 = q6 + p6 * s6
      q8 = q8 + p8 * s8
      qa = qa + pa * s8 * s2
      endif
 100  continue
      pd = dd * pg2
      q0 = pd * q0
      q2 = pd * q2
      q4 = pd * q4
      q6 = pd * q6
      q8 = pd * q8
      qa = pd * qa
      j  = jmax1
      rr = rhf(j)
      r2 = rr * rr
      do   200   k = 2, kmax1
      zz = rhf(k)
      z2 = zz * zz
      s2 = 1.0 / (r2 + z2)
      s1 = sqrt(s2)
      x2 = z2 * s2
      p2 = x2 * c1 - c0
      p4 = x2 * (x2 * c4 - c3) + c2
      p6 = x2 * (x2 * (x2 * c8 - c7) + c6) - c5
      p8 = x2 * (x2 * (x2 * (x2 * cd - cc) + cb) - ca) + c9
      pa = x2 * (x2 * (x2 * (x2 * (x2 * cj - ci) + ch) - cg) + cf) - ce
      phi(j,k) = s1 * (q0 + s2 * (p2 * q2 + s2 * (p4 * q4 + s2
     1              * (p6 * q6 + s2 * (p8 * q8 + s2 * pa * qa)))))
 200  continue
      k  = kmax1
      zz = rhf(k)
      z2 = zz * zz
      do   300   j = 2, jmax1
      rr = rhf(j)
      r2 = rr * rr
      s2 = 1.0 / (r2 + z2)
      s1 = sqrt(s2)
      x2 = z2 * s2
      p2 = x2 * c1 - c0
      p4 = x2 * (x2 * c4 - c3) + c2
      p6 = x2 * (x2 * (x2 * c8 - c7) + c6) - c5
      p8 = x2 * (x2 * (x2 * (x2 * cd - cc) + cb) - ca) + c9
      pa = x2 * (x2 * (x2 * (x2 * (x2 * cj - ci) + ch) - cg) + cf) - ce
      phi(j,k) = s1 * (q0 + s2 * (p2 * q2 + s2 * (p4 * q4 + s2
     1              * (p6 * q6 + s2 * (p8 * q8 + s2 * pa * qa)))))
 300  continue
      phi(jmax1,1) = phi(jmax1,2)
      phi(1,kmax1) = phi(2,kmax1)
      return
      end
      subroutine pot2d

      implicit real*8 (a-h,o-z)
      include "prmtr.h"

      parameter      (mw = 50000, jm2 = jmax - 2,
     1               jmax0 = jmax - 1, kmax0 = kmax - 1,
     2   	     jmax1 = jmax + 1, kmax1 = kmax + 1,
     3   	     jmax2 = jmax + 2, kmax2 = kmax + 2)
      dimension phi(jmax2, kmax2), rho(jmax2, kmax2),
     1               an(kmax0), bn(kmax0), cn(kmax0), am(jmax0),
     2		     bm(jmax0), cm(jmax0), y(jmax0, kmax0), wfw(mw),
     3		     c(jmax0), rd3(jmax), rhf(jmax1)
      common /vec /  phi, rho, rhf, pi, dr

      pp           = 0.25 / (pi * pi)
      pg4          = 4.0 * pi
      odr          = 1.0 / dr
      od2          = odr * odr
      tdr          = 2.0 * odr
      tm2          = -2.0* od2
      hdr          = 0.5 * odr
      do   210   j = 1, jmax
      rd3(j)       = 1.0 / rhf(j)
 210  continue
      do   230   j = 2, jmax
      i            = j - 1
      rr           = rd3(j)
      am(i)        = hdr * (tdr - rr)
      cm(i)        = hdr * (rr + tdr)
      c(i)         = rr  * rr * pp
 230  continue
      do   240   k = 1, kmax0
      cn(k)        = od2
      an(k)        = od2
      bn(k)        = tm2
 240  continue
      cmmax        = cm(jmax0)
      cm(jmax0)    = 0.0
      cnmax        = cn(kmax0)
      cn(kmax0)    = 0.0
      am(1)        = 0.0
      bn(1)        = bn(1) + an(1)
      an(1)        = 0.0
      do   310   k = 2, kmax
      ik           = k - 1
      do   310   j = 2, jmax
      ij           = j - 1
      y(ij,ik)     = pg4 * rho(j,k)
 310  continue
      do   340   j = 2, jm2
      bm(j)        = - (cm(j) + am(j))
 340  continue
      do   330   k = 2, kmax
      ik           = k - 1
      y(jmax0,ik)  = y(jmax0,ik) - cmmax * phi(jmax1,k)
 330  continue
      bm(jmax0)    = - (cmmax + am(jmax0))
      do   320   j = 2, jmax
      ij           = j - 1
      y(ij,kmax0)  = y(ij,kmax0) - cnmax * phi(j,kmax1)
 320  continue
      bm(1)        = - cm(1)
      j0           = jmax0
      k0           = kmax0
      call blktri(0,1,k0,an,bn,cn,1,j0,am,bm,cm,j0,y,ier0,wfw)
      call blktri(1,1,k0,an,bn,cn,1,j0,am,bm,cm,j0,y,ier1,wfw)
      do   350   k = 2, kmax
      ik           = k - 1
      do   350   j = 2, jmax
      ij           = j - 1
      phi(j,k)     = y(ij,ik)
 350  continue
      do   570   k = 1, kmax1
      phi(1,k)     = phi(2,k)
 570  continue
      do   575   j = 1, jmax2
 575  phi(j,1)     = phi(j,2)
      return
      end

      subroutine blktr1(n,an,bn,cn,m,am,bm,cm,idimy,y,b,w1,w2,w3,wd,
     1                ww,wu)

      implicit real*8 (a-h,o-z)

      dimension an(1),bn(1),cn(1),am(1),bm(1),cm(1),b(1),
     1                w1(1),w2(1),w3(1),wd(1),ww(1),wu(1),y(idimy,1)
      common /cblkt/  eps,cnv,npp,k,nm,ncmplx,ik

      kdo = k-1
      do 109 l=1,kdo
      ir = l-1
      i2 = 2**ir
      i1 = i2/2
      i3 = i2+i1
      i4 = i2+i2
      irm1 = ir-1
      call indxb (i2,ir,im2,nm2)
      call indxb (i1,irm1,im3,nm3)
      call indxb (i3,irm1,im1,nm1)
      call prod (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,y(1,i2),w3,
     1               m,am,bm,cm,wd,ww,wu)
      ig = 2**k
      do 108 i=i4,ig,i4
      if (i-nm) 101,101,108
 101  ipi1 = i+i1
      ipi2 = i+i2
      ipi3 = i+i3
      call indxc (i,ir,idxc,nc)
      if (i-ig) 102,108,108
 102  call indxa (i,ir,idxa,na)
      call indxb (i-i1,irm1,im1,nm1)
      call indxb (ipi2,ir,ip2,np2)
      call indxb (ipi1,irm1,ip1,np1)
      call indxb (ipi3,irm1,ip3,np3)
      call prod (nm1,b(im1),0,dum,0,dum,na,an(idxa),w3,w1,m,am,
     1                  bm,cm,wd,ww,wu)
      if (ipi2-nm) 105,105,103
 103  do 104 j=1,m
      w3(j) = 0.0
      w2(j) = 0.0
 104  continue
      go to 106
 105  call prod (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum,
     1                  y(1,ipi2),w3,m,am,bm,cm,wd,ww,wu)
      call prod (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w3,w2,m,am,
     1                  bm,cm,wd,ww,wu)
 106  do 107 j=1,m
      y(j,i) = w1(j)+w2(j)+y(j,i)
 107  continue
 108  continue
 109  continue
 132  do 144 ll=1,k
      l = k-ll+1
      ir = l-1
      irm1 = ir-1
      i2 = 2**ir
      i1 = i2/2
      i4 = i2+i2
      ifd = ig-i2
      do 143 i=i2,ifd,i4
      if (i-nm) 133,133,143
 133  imi1 = i-i1
      imi2 = i-i2
      ipi1 = i+i1
      ipi2 = i+i2
      call indxa (i,ir,idxa,na)
      call indxc (i,ir,idxc,nc)
      call indxb (i,ir,iz,nz)
      call indxb (imi1,irm1,im1,nm1)
      call indxb (ipi1,irm1,ip1,np1)
      if (i-i2) 134,134,136
 134  do 135 j=1,m
      w1(j) = 0.0
 135  continue
      go to 137
 136  call prod (nm1,b(im1),0,dum,0,dum,na,an(idxa),y(1,imi2),
     1                  w1,m,am,bm,cm,wd,ww,wu)
 137  if (ipi2-nm) 140,140,138
 138  do 139 j=1,m
      w2(j) = 0.0
 139  continue
      go to 141
 140  call prod (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),y(1,ipi2),
     1                  w2,m,am,bm,cm,wd,ww,wu)
 141  do 142 j=1,m
      w1(j) = y(j,i)+w1(j)+w2(j)
 142  continue
      call prod (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,y(1,i),
     1                  m,am,bm,cm,wd,ww,wu)
 143  continue
 144  continue
      return
      end

      subroutine blktri(iflg,np,n,an,bn,cn,mp,m,am,bm,cm,idimy,y,
     1                ierror,w)

      implicit real*8 (a-h,o-z)

      dimension an(1),bn(1),cn(1),am(1),
     1                bm(1),cm(1),y(idimy,1) ,w(1)
      common /cblkt/  eps,cnv,npp,k,nm,ncmplx,ik

      nm = n
      ierror = 0
      if (m-5) 101,102,102
 101  ierror = 1
      go to 119
 102  if (nm-3) 103,104,104
 103  ierror = 2
      go to 119
 104  if (idimy-m) 105,106,106
 105  ierror = 3
      go to 119
 106  nh = n
      npp = np
 107  nh = nh+1
 108  ik = 2
      k = 1
 109  ik = ik+ik
      k = k+1
      if (nh-ik) 110,110,109
 110  nl = ik
      ik = ik+ik
      nl = nl-1
      iwah = (k-2)*ik+k+6
 111  iw1 = iwah
      iwbh = iw1+nm
      w(1) = float(iw1-1+max0(2*nm,6*m))
      go to 113
 112  iwbh = iwah+nm+nm
      iw1 = iwbh
      w(1) = float(iw1-1+max0(2*nm,6*m))
      nm = nm-1
 113  if (ierror) 119,114,119
 114  iw2 = iw1+m
      iw3 = iw2+m
      iwd = iw3+m
      iww = iwd+m
      iwu = iww+m
      if(iflg) 116,115,116
 115  call compb(nl,ierror,an,bn,cn,w(2),w(iwah),w(iwbh))
      go to 119
 116  continue 
 117  call blktr1(nl,an,bn,cn,m,am,bm,cm,idimy,y,w(2),w(iw1),w(iw2),
     1             w(iw3),w(iwd),w(iww),w(iwu))
 119  continue
      return
      end
      subroutine compb(n,ierror,an,bn,cn,b,ah,bh)

      implicit real*8 (a-h,o-z)

      dimension an(1),bn(1),cn(1),b(1),ah(1),bh(1)
      common /cblkt/  eps,cnv,npp,k,nm,ncmplx,ik

      t08 = 256.0
      t16 = t08 * t08
      t32 = t16 * t16
      t64 = t32 * t32
      eps = 100.0 / t64
      bnorm = abs(bn(1))
      do 102 j=2,nm
      bnorm = dmax1(bnorm,abs(bn(j)))
      arg = an(j)*cn(j-1)
      if (arg) 119,101,101
 101  b(j) = sign(sqrt(arg),an(j))
 102  continue
      cnv = eps*bnorm
      if = 2**k
      kdo = k-1
      do 108 l=1,kdo
      ir = l-1
      i2 = 2**ir
      i4 = i2+i2
      ipl = i4-1
      ifd = if-i4
      do 107 i=i4,ifd,i4
      call indxb (i,l,ib,nb)
      if (nb) 108,108,103
 103  js = i-ipl
      jf = js+nb-1
      ls = 0
      do 104 j=js,jf
      ls = ls+1
      bh(ls) = bn(j)
      ah(ls) = b(j)
 104  continue
      call tevls (nb,bh,ah,ierror)
      if (ierror) 118,105,118
 105  lh = ib-1
      do 106 j=1,nb
      lh = lh+1
      b(lh) = -bh(j)
 106  continue
 107  continue
 108  continue
      do 109 j=1,nm
      b(j) = -bn(j)
 109  continue
 117  return
 118  ierror = 4
      return
 119  ierror = 5
      return
      end
      subroutine indxa(i,ir,idxa,na)

      implicit real*8 (a-h,o-z)
      common /cblkt/eps,cnv,npp,k,nm,ncmplx,ik

      na = 2**ir
      idxa = i-na+1
      if (i-nm) 102,102,101
 101  na = 0
 102  return
      end
      subroutine indxb(i,ir,idx,idp)

      implicit real*8 (a-h,o-z)

      common /cblkt/eps,cnv,npp,k,nm,ncmplx,ik

      idp = 0
      if (ir) 107,101,103
 101  if (i-nm) 102,102,107
 102  idx = i
      idp = 1
      return
 103  izh = 2**ir
      id = i-izh-izh
      idx = id+id+(ir-1)*ik+ir+(ik-i)/izh+4
      ipl = izh-1
      idp = izh+izh-1
      if (i-ipl-nm) 105,105,104
 104  idp = 0
      return
 105  if (i+ipl-nm) 107,107,106
 106  idp = nm+ipl-i+1
 107  return
      end
      subroutine indxc(i,ir,idxc,nc)

      implicit real*8 (a-h,o-z)

      common /cblkt/eps,cnv,npp,k,nm,ncmplx,ik

      nc = 2**ir
      idxc = i
      if (idxc+nc-1-nm) 102,102,101
 101  nc = 0
 102  return
      end
      subroutine prod(nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,w,u)

      implicit real*8 (a-h,o-z)

      dimension a(1),b(1),c(1),x(1),y(1),d(1),w(1),bd(1),
     1                bm1(1),bm2(1),aa(1),u(1)

      do 101 j=1,m
      w(j) = x(j)
      y(j) = w(j)
 101  continue
      mm = m-1
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
 102  if (ia) 105,105,103
 103  rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
      do 104 j=1,m
      y(j) = rt*w(j)
 104  continue
 105  if (id) 125,125,106
 106  rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
      d(m) = a(m)/(b(m)-rt)
      w(m) = y(m)/(b(m)-rt)
      do 107 j=2,mm
      k = m-j
      den = b(k+1)-rt-c(k+1)*d(k+2)
      d(k+1) = a(k+1)/den
      w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
 107  continue
      den = b(1)-rt-c(1)*d(2)
      w(1) = 1.
      if (den) 108,109,108
 108  w(1) = (y(1)-c(1)*w(2))/den
 109  do 110 j=2,m
      w(j) = w(j)-d(j)*w(j-1)
 110  continue
      if (na) 113,113,102
 111  do 112 j=1,m
      y(j) = w(j)
 112  continue
      ibr = 1
      go to 102
 113  if (m1) 114,114,115
 114  if (m2) 111,111,120
 115  if (m2) 117,117,116
 116  if (abs(bm1(m1))-abs(bm2(m2))) 120,120,117
 117  if (ibr) 118,118,119
 118  if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 111,119,119
 119  rt = rt-bm1(m1)
      m1 = m1-1
      go to 123
 120  if (ibr) 121,121,122
 121  if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 111,122,122
 122  rt = rt-bm2(m2)
      m2 = m2-1
 123  do 124 j=1,m
      y(j) = y(j)+rt*w(j)
 124  continue
      go to 102
 125  return
      end
      subroutine tevls(n,d,e2,ierr)

      implicit real*8 (a-h,o-z)
      real*8 machep

      dimension d(n),e2(n)
      integer*4       i,j,l,m,n,ii,l1,mml,ierr
      common /cblkt/  machep,cnv,npp,k,nm,ncmplx,ik

      ierr = 0
      if (n .eq. 1) go to 115
      do 101 i=2,n
      e2(i-1) = e2(i)*e2(i)
 101  continue
      f = 0.0
      b = 0.0
      e2(n) = 0.0
      do 112 l=1,n
      j = 0
      h = machep*(abs(d(l))+sqrt(e2(l)))
      if (b .gt. h) go to 102
      b = h
      c = b*b
 102  do 103 m=l,n
      if (e2(m) .le. c) go to 104
 103  continue
 104  if (m .eq. l) go to 108
 105  if (j .eq. 30) go to 114
      j = j+1
      l1 = l+1
      s = sqrt(e2(l))
      g = d(l)
      p = (d(l1)-g)/(2.0*s)
      r = sqrt(p*p+1.0)
      d(l) = s/(p+sign(r,p))
      h = g-d(l)
      do 106 i=l1,n
      d(i) = d(i)-h
 106  continue
      f = f+h
      g = d(m)
      if (g .eq. 0.0) g = b
      h = g
      s = 0.0
      mml = m-l
      do 107 ii=1,mml
      i = m-ii
      p = g*h
      r = p+e2(i)
      e2(i+1) = s*r
      s = e2(i)/r
      d(i+1) = h+s*(h+d(i))
      g = d(i)-e2(i)/g
      if (g .eq. 0.0) g = b
      h = g*p/r
 107  continue
      e2(l) = s*g
      d(l) = h
      if (h .eq. 0.0) go to 108
      if (abs(e2(l)) .le. abs(c/h)) go to 108
      e2(l) = h*e2(l)
      if (e2(l) .ne. 0.0) go to 105
 108  p = d(l)+f
      if (l .eq. 1) go to 110
      do 109 ii=2,l
      i = l+2-ii
      if (p .ge. d(i-1)) go to 111
      d(i) = d(i-1)
 109  continue
 110  i = 1
 111  d(i) = p
 112  continue
      if (abs(d(n)) .ge. abs(d(1))) go to 115
      nhalf = n/2
      do 113 i=1,nhalf
      ntop = n-i
      dhold = d(i)
      d(i) = d(ntop+1)
      d(ntop+1) = dhold
 113  continue
      go to 115
 114  ierr = l
 115  return
      end
