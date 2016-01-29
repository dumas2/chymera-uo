#!/usr/bin/python

import math as M
import numpy as N
import pylab as P

f=open("sig.053000","r")
Size=16
od=3.

x=[]
y=[]
v=[]

n=0
mass = 0.
for line in f:
	if (line[0]=="#"): continue
	row=line.split()
	if (len(row)==0): continue
	x.append([])
	y.append([])
	v.append([])
	x[n]=float(row[0])
	y[n]=float(row[1])
	v[n]=float(row[2])
	n+=1

dx = max(x[1]-x[0],y[1]-y[0])

print n, "lines read.", M.sqrt(n)
print n, "Grid dimension: ", dx

m=int(M.sqrt(n))

v=N.reshape(v,(m,m))#,order="FORTRAN")
flag=N.zeros((m,m),"int")

mass=0.
for i in xrange(m):
      for j in xrange(m):
           mass+=10.**v[i,j]
print "Total mass:",mass*(1.5e13*dx)**2/2e33

for iii in xrange(2):
	mass=0.
	maxv=0.
	for i in xrange(m):
		for j in xrange(m):
			mass+=10.**v[i,j]
			if(maxv<v[i,j] and not flag[i,j] ):
				maxv=v[i,j]
				ix=i
				iy=j
	lm=m/Size
	lv=N.zeros((lm,lm),"float")
	rd=N.zeros((lm*lm),"float")
	rr=N.zeros((lm*lm),"float")
        flag[ix-16:ix+16,iy-16:iy+16]=True
#	print mass/2e33*9e26,maxv,ix,iy
	mass = 0.
	for i in xrange(lm):
		for j in xrange(lm):
			lv[i,j]=v[ix-lm/2+i,iy-lm/2+j]
			rd[i*j+j]=lv[i,j]
			rr[i*j+j]=M.sqrt(float( (lm/2-i)**2+(lm/2-j)**2))*dx
			mass += 10.**lv[i,j]
        avgDen=mass/float(lm)**2
        maxv=10.**maxv
        mass = 0.
        for i in xrange(lm):
                for j in xrange(lm):
			den=lv[i,j]
                        lv[i,j]=M.log10(max(10.**v[ix-lm/2+i,iy-lm/2+j]-od*avgDen,1e-10))
                        rd[i*j+j]=den
                        rr[i*j+j]=M.sqrt(float( (lm/2-i)**2+(lm/2-j)**2))*dx
			#if lv[i,j] > -10. : mass+=den/2e33*(1.5e13*dx)**2
                        #mass += 10.**lv[i,j]/2e33*(1.5e13*dx)**2
                        if 10.**lv[i,j]>od*avgDen: mass += 10.**lv[i,j]/2e33*(1.5e13*dx)**2

#	if mass*2e33/(float(lm)**2*9e26) < (10.**maxv)/10.:continue
        if mass*1e3<1.:continue
        print mass*1e3,ix,iy,maxv,avgDen*od,M.sqrt((ix-256.)**2+(iy-256.)**2)*2.

	P.figure()
	P.contourf(lv,100)
	P.figure()
	P.plot(rr,rd,"x")


P.figure()
P.imshow(v,aspect=1)
P.show()
