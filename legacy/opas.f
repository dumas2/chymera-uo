c	::::::::::::::
c	opas.f
c	::::::::::::::
	subroutine opas(t,p,rosstab,plancktab)
	implicit double precision (a-h,o-z)

c	Para ser usada en programas de DISCO 
c	Interpola en tablas de Rosseland mean opacity,
c	de planck, y de mu
c	 en funcion de presion y temperatura 

c	Archivo (binario) de entrada  
c.	limites: 10 a 1.e7 K, 1.e-21 a 1.e9 din/cm^2
c	Tablas en log T, log p, log chi
	parameter(jmax=100,kmax=100,np=jmax,nt=kmax)
	character*100 fileros,fileplan
	dimension o1(2),o2(2)
	dimension p1(2),p2(2)
	common/rosstabla/chitab(jmax,kmax),ttab(kmax),
     .	ptab(jmax)
	common/plnktabla/chitabp(jmax,kmax),ttabp(kmax),
     .  ptabp(jmax)
	common/ntablas/fileros,fileplan
	
	rosstab=0.
	plancktab=0.

	if(chitab(1,1).eq.0.) then

	pmin=-12.
	pmax=9.
	tmin=0.5
	tmax=7.

	dp=(pmax-pmin)/(np-1)
	dt=(tmax-tmin)/(nt-1)
		do it=1,nt
		ttab(it)=tmin+(it-1)*dt
		enddo
			do ip=1,np
			ptab(ip)=pmin+(ip-1)*dp
			enddo
	open(unit=50,file=fileros,status='old', form='formatted')
	do j=1,jmax
	do k=1,kmax
	read(50,*) chitab(j,k)
	enddo
	enddo
	close(50)

	open(unit=52,file=fileplan, status='old', form='formatted')
	do j=1,jmax
	do k=1,kmax
        read(52,*) chitabp(j,k)
	enddo
	enddo
        close(52)

	endif


c
	xlt=dlog10(t)
	xlp=dlog10(p)



c
c	Fuera de limites
c	!T en orden decreciente

	if(xlt.gt.ttab(nt)) then
	write(*,*) 'opas T>T(max)' 
	write(*,*)'10.**xlt=',10.**xlt
	write(*,*)'10.**ttab(nt)=',10.**ttab(nt)
	xlt=ttab(nt)
c	stop
	endif

c	if(xlp.lt.ptab(1).or.xlp.gt.ptab(np)) then
	write(*,*) 'opas P fuera,t,p=',t,p
		if (xlp.lt.ptab(1))xlp=ptab(1)
		if (xlp.gt.ptab(np))xlp=ptab(np)
c	stop
c	endif
	
c	asigna valores a T(1) para t<T(1)
	if(xlt.lt.ttab(1)) xlt=ttab(1)

c	prueba
c	if (xlt.lt.1.)xlt=1.
c

	it1=dint((xlt-ttab(1))/(ttab(nt)-ttab(1))*
     *  dfloat(nt-1)+1.)
	it2=it1+1
		
		if (it2.gt.nt) it2=it1

c.	prueba
c	if (it1.gt.30) then
c	it1=30
c	it2=30
c	endif

	ip1=dint((xlp-ptab(1))/(ptab(np)-ptab(1))*
     *  dfloat(np-1)+1.)
	ip2=ip1+1

		if (ip2.gt.np) ip2=ip1

c.	ubico las presiones y temperaturas
c
	l=1
	do it=it1,it2
		j=1
		do ip=ip1,ip2
		o1(j)=chitab(ip,it)
		p1(j)=chitabp(ip,it)
		j=j+1
c		de p
		enddo 

	if (ip2.ne.ip1) then
	o2(l)=o1(1)+(o1(2)-o1(1))*(xlp-ptab(ip1))/(ptab(ip2)-ptab(ip1))
	p2(l)=p1(1)+(p1(2)-p1(1))*(xlp-ptab(ip1))/(ptab(ip2)-ptab(ip1))
	else
	o2(l)=o1(1)
	p2(l)=p1(1)
	endif

	l=l+1
c	de t
	enddo 

	if (it2.ne.it1) then
	rosstab=o2(1)+(o2(2)-o2(1))*(xlt-ttab(it1))/
     *	(ttab(it2)-ttab(it1))
	plancktab=p2(1)+(p2(2)-p2(1))*(xlt-ttab(it1))/
     *	(ttab(it2)-ttab(it1))
	else
	rosstab=o2(1)
	plancktab=p2(1)
	endif

	rosstab=10.**rosstab
	plancktab=10.**plancktab


	if (xlt.ge.6.) then
        rosstab=0.3429
        endif
	return
	end

c	::::::::::::::
c	pmm.f
c	::::::::::::::
	function pmm(t,p)
	implicit double precision (a-h,o-z)
c	real*4 chitab,ttab,ptab

c	Para ser usada en programas de DISCO 
c	Interpola en tablas de Rosseland mean opacity
c	 en funcion de presion y temperatura 

c	Archivo (binario) de entrada  
c.	limites: 10**(0.5) a 1.e7 K, 1.e-12 a 1.e9 din/cm^2
c	Tablas en log T, log p, log chi
	parameter(jmax=100,kmax=100,np=jmax,nt=kmax)
	dimension o1(2),o2(2)
	common/pmmtabla1/chitab(jmax,kmax),ttab(kmax),
     .	ptab(jmax)
c
	rosstab=0.
	opa1=0.
        opa2=0.

	if(chitab(1,1).eq.0.) then
	open(unit=50,file=
     *'/alp7_2/dalessio/COMUN/pmm.dat',
     *status='old',
     *	form='formatted')
	do j=1,jmax
	do k=1,kmax
	read(50,*)chitab(j,k)
	enddo
	enddo
	close(50)

	pmin=-12.
        pmax=9.
        tmin=0.5
        tmax=7.
        dp=(pmax-pmin)/(np-1)
        dt=(tmax-tmin)/(nt-1)
                do it=1,nt
                ttab(it)=tmin+(it-1)*dt
                enddo
                        do ip=1,np
                        ptab(ip)=pmin+(ip-1)*dp
                        enddo
c
	endif

	xlt=dlog10(t)
	xlp=dlog10(p)

c	Fuera de limites
c	!T en orden decreciente


	if(xlt.gt.ttab(nt)) stop 'pmm T>T(max)' 
c	if(xlp.lt.ptab(1).or.xlp.gt.ptab(np)) then
c	write(*,*) 'pmm P fuera,t,p=',t,p
		if (xlp.lt.ptab(1))xlp=ptab(1)
		if (xlp.gt.ptab(np))xlp=ptab(np)
c	stop
c	endif
	
c	asigna la pmm a T(1) para t<T(1)
	if(xlt.lt.ttab(1)) xlt=ttab(1)
c

        it1=dint((xlt-ttab(1))/(ttab(nt)-ttab(1))*
     *  dfloat(nt-1)+1.)
        it2=it1+1

                if (it2.gt.nt) it2=it1
        ip1=dint((xlp-ptab(1))/(ptab(np)-ptab(1))*
     *  dfloat(np-1)+1.)
        ip2=ip1+1

                if (ip2.gt.np) ip2=ip1

c
	l=1
	
	
	do it=it1,it2
		j=1
		do ip=ip1,ip2
		o1(j)=chitab(ip,it)
		j=j+1
c		de p
		enddo 

	if (ip2.ne.ip1) then
	o2(l)=o1(1)+(o1(2)-o1(1))*(xlp-ptab(ip1))/(ptab(ip2)-ptab(ip1))
	else
	o2(l)=o1(1)
	endif

	l=l+1
c	de t
	enddo 

	if (it2.ne.it1) then
	pmm=o2(1)+(o2(2)-o2(1))*(xlt-ttab(it1))/(ttab(it2)-ttab(it1))
	else
	pmm=o2(1)
	endif

c	write(*,*)'xlt,xlp,pmm=',xlt,xlp,pmm
	return
	end

