#!/usr/bin/python2.5

import math as M
import numpy as N
import pylab as P

f=open("image/loadm5_c_t_w2l_sig.320000","r")
g=open("image/loadm5_c_t_w2l_sig.321000","r")
Size=16
od=1.5
au=1.49e13
Msun=1.989e33

x=[]
y=[]
v=[]
x2=[]
y2=[]
v2=[]

n=0
mass = 0.
for line in f:
	if line[0]=="#":
		junk,time=line.split()
		time=float(time)
		continue
	row=line.split()
	if (len(row)==0): continue
	x.append(float(row[0]))
	y.append(float(row[1]))
	v.append(float(row[2]))
	n+=1

for line in g:
        if line[0]=="#":
                junk,time2=line.split()
                time2=float(time2)
                continue
	row=line.split()
        if (len(row)==0): continue
	x2.append(float(row[0]))
	y2.append(float(row[1]))
	v2.append(float(row[2]))


n=len(v)
n2=len(v2)
if n!=n2:
	print "Files do not have same data structure."
	sys.exit(-1)

dx = max(x[1]-x[0],y[1]-y[0])

print n, "lines read.", M.sqrt(n)
print n, "Grid dimension: ", dx
print "Time difference: ",time2-time

m=int(M.sqrt(n))

v=N.reshape(v,(m,m))#,order="FORTRAN")
v2=N.reshape(v2,(m,m))#,order="FORTRAN")
flag=N.zeros((m,m),"int")


mass=0.
mass2=0.
for i in xrange(m):
	for j in xrange(m):
		mass+=dx*au*dx*au/Msun*10.**v[i][j]
		mass2+=dx*au*dx*au/Msun*10.**v2[i][j]
print "Mass in file 1, file 2, mdot", mass,mass2,(mass2-mass)/(time2-time)


mdot=[]
mdotsum=0.
for i in xrange(m):
	mdot.append([])
        for j in xrange(m):
		mdot[i].append((dx*dx*au*au*(10.**v2[i][j]-10.**v[i][j])/(time2-time))/Msun)
		mdotsum+=mdot[i][j]

print "Mdotsum", mdotsum

P.figure()
C=P.imshow(mdot,aspect=1)
P.colorbar(C,format='%5.2e')
P.show()
