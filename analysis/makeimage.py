'''Script for making multiple images from sigplot or tkplot output. User should adjust colorscale limits
minVal and maxVal, as well as the TICKS for value indication.'''
#!/usr/bin/python
#Use your appropriate python
# ASSUMES AU UNITS

import math as M
import numpy as N
import pylab as P

# USER DEFINED STUFF
START=10
STOP =20
SKIP =1000
outputdir='./png/'     # YOU MAY NEED TO CREATE THESE DIRECTORIES.
inputdir='./images/'
prefix='sig_'

minVal=-1.
maxVal=3.0
TICKS=-.5,0.,.5,1.,1.5,2.,2.5,3.0
FORMAT='%3.1f'
x0=-100.
x1=100.
y0=-100.
y1=100.
XLABEL="AU"
YLABEL="AU"
CMAP=P.cm.jet

def main():
	for thisFile in xrange(START,STOP+1,SKIP):
		filein=inputdir+prefix+repr(thisFile).zfill(6)
		f=open(filein,"r")


		x=[]
		y=[]
		v=[]

		n=0
		mass = 0.
		for line in f:
			if (line[0]=="#"): continue
			row=line.split()
			if (len(row)==0): continue
			x.append(float(row[0]))
			y.append(float(row[1]))
			v.append(float(row[2]))
			n+=1
		f.close()
		savefile=outputdir+prefix+repr(thisFile).zfill(6)+".png"

		dx = max(x[1]-x[0],y[1]-y[0])

		m=int(M.sqrt(n))
		print n, "lines read.", m, ' by ',m
		print n, "Grid dimension: ", dx

		mass=0.
		vg=[]
		vgflip=[]
		xg=[]
		yg=[]
		iter=0
		jiter=0
		for i in xrange(m):
			vg.append([])
			vgflip.append([])
			xg.append(dx*(float(i)+.5)-dx*m/2.)
			for j in xrange(m):
				vgflip[i].append(0.)
				if(jiter==j):
					yg.append(dx*(float(j)+.5)-dx*m/2.)
					jiter+=1
				vg[i].append(v[iter])
				iter+=1
				mass+=10.**vg[i][j]
		for i in xrange(m):
			for j in xrange(m):
				vgflip[j][i]=vg[i][j]

		print "Total mass:",mass*(1.5e13*dx)**2/2e33

		dz=(maxVal-minVal)/256.
		cont=[]
		for i in xrange(256):
			cont.append(minVal+(i+1)*dz)

		P.figure()
		IM=P.contourf(xg,yg,vgflip,cont,cmap=CMAP)
		CB=P.colorbar(IM,format=FORMAT,ticks=TICKS)
		IM=P.xlim(x0,x1)
		IM=P.ylim(y0,y1)
		IM=P.xlabel(XLABEL)
		IM=P.ylabel(YLABEL)
		P.savefig(savefile)

if __name__ == "__main__":main()
