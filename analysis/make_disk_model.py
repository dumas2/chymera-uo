#!/usr/bin/python

import math

PRINTOPAC=1

# some important constants
N_grid = 256  # Number of radial grid points
disk_size=100. # disk size in AU
mstar=2e33 # Stellar mass cgs
T0=140. # Irradiation temperature at 1 AU
Q=1.5 # Target Q
innerBound=10. # inner region to consider for mass calculation.  AU
opac_frac=1.0 #
itable = 100 # Opacity table.  Don't change
plcktablefile="planck.abpoll.p3p5.amax1.dat"
rosstablefile="rosseland.abpoll.p3p5.amax1.dat"

#declare arrays
p_tab=[]
t_tab=[]
ros_table=[]
plk_table=[]

#constants
AU=1.5e13
rgas=8.254e7
mu=2.33
GN=6.67e-8

def read_opacity_table():
	f=open(rosstablefile,"r")
	for i in xrange(itable):
		ros_table.append([])
		for j in xrange(itable):
			line=f.readline()
			ros_table[i].append(float(line))
	f.close()
        f=open(plcktablefile,"r")
        for i in xrange(itable):
                plk_table.append([])
                for j in xrange(itable):
                        line=f.readline()
                        plk_table[i].append(float(line))
        f.close()
	Pmin=-12.
	Pmax=9.0
	Tmin=0.5
	Tmax=7.0
	dP=(Pmax-Pmin)/(float(itable-1))
	dT=(Tmax-Tmin)/(float(itable-1))
	for i in xrange(itable):
		t_tab.append(Tmin+float(i)*dT)
		p_tab.append(Pmin+float(i)*dT)
	return 1

def dalessio(Tk,Pcgs):
	xlt=math.log10(Tk)
	xlp=math.log10(Pcgs)
	itable1=itable-1

	o1=[]
	o2=[]
	p1=[]
	p2=[]
	for i in xrange(2):
		o1.append(0.)
		o2.append(0.)
		p1.append(0.)
		p2.append(0.)

        if xlt>t_tab[itable1]:xlt=t_tab[itable1]
        if xlp>p_tab[itable1]:xlp=p_tab[itable1]
        if xlp<p_tab[0]:xlp=p_tab[0]
        if xlt<t_tab[0]:xlt=t_tab[0]

        it1=int((xlt-t_tab[0])/(t_tab[itable1]-t_tab[0])*float(itable1))
        it2=it1+1
        if it2>itable1:it2=it1

        ip1=int((xlp-p_tab[0])/(p_tab[itable1]-p_tab[0])*float(itable1))
        ip2=ip1+1
        if ip2>itable1:ip2=ip1
	
        l=0
	for it in xrange(it1,(it2+1)):
		ij=0
		for ip in xrange(ip1,(ip2+1)):
			o1[ij]=ros_table[ip][it]
			p1[ij]=plk_table[ip][it]
			ij+=1

		if  ip2 != ip1:
			o2[l]=o1[0]+(o1[1]-o1[0])*(xlp-p_tab[ip1])/(p_tab[ip2]-p_tab[ip1])
			p2[l]=p1[0]+(p1[1]-p1[0])*(xlp-p_tab[ip1])/(p_tab[ip2]-p_tab[ip1])
                else: 
			o2[l]=o1[0]
			p2[l]=p1[0]

		l+=1

	if it2 != it1:
		Ors=o2[0]+(o2[1]-o2[0])*(xlt-t_tab[it1])/(t_tab[it2]-t_tab[it1])
		Opl=p2[0]+(p2[1]-p2[0])*(xlt-t_tab[it1])/(t_tab[it2]-t_tab[it1])
	else:
		Ors=o2[0]
		Opl=p2[0]

	Ors=10.**Ors
	Opl=10.**Opl		
	return Ors*opac_frac,Opl*opac_frac

def main():
	OK=read_opacity_table()
	if PRINTOPAC:
		p=rgas/mu*2e-11*100.
		print dalessio(100.,p)
		return
	surfden=[]
	t_mid=[]
	t_irr=[]
	tau = []
	csound=[]
	for i in xrange(N_grid):
		r=(float(i)+.5)/float(N_grid)*disk_size
		omega=math.sqrt(GN*mstar/(r*AU)**3)
		t_irr.append(T0/r**0.5+10.)
		csound.append(math.sqrt(rgas/mu*t_irr[i]*1.6))
		surfden.append(csound[i]*omega/math.pi/GN/Q)
		t_mid.append(t_irr[i])
		tau.append(0.)

	it=0
	while True:
		frac=0.
		for i in xrange(N_grid):
			r=(float(i)+.5)/float(N_grid)*disk_size
			p=rgas/mu*t_mid[i]
			omega=math.sqrt(GN*mstar/(r*AU)**3)
			oross,oplck=dalessio(t_mid[i],p)
			tau_temp=oplck*surfden[i]
			tau_temp=oross*surfden[i]*(1.-math.exp(-tau_temp*4.))+oplck*surfden[i]*math.exp(-tau_temp*4.)
			tau[i]=tau[i]*.9+tau_temp*.1
			t_mid[i]=t_irr[i]*(3./4.*(tau[i]/2.+4./3.))**0.25
			csound[i]=math.sqrt(rgas/mu*t_mid[i])
			olds=surfden[i]
			surfden[i]=0.9*surfden[i]+0.1*(csound[i]*omega/math.pi/GN/Q)
			frac+=math.sqrt((surfden[i]-olds)**2)/olds/float(N_grid)
		if frac<1e-6:break
		it+=1
	
	mass=0.
	for i in xrange(N_grid):
		r=(float(i)+.5)/float(N_grid)*disk_size
                omega=math.sqrt(GN*mstar/(r*AU)**3)
		q=csound[i]*omega/math.pi/GN/surfden[i]
		if r>innerBound:mass+=disk_size/float(N_grid)*r*surfden[i]/2e33*AU**2
		print i,r,t_mid[i],surfden[i],tau[i],t_irr[i],q

	print "#DISK MASS ",mass*2.*math.pi
	return 1

OK=main()
if OK: print "# Program exited normally."
else: print "# Error encountered in main."
	
