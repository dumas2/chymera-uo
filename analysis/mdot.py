#!/usr/bin/python2.5
''' This program plots an mdot profile.  The data are assumed to be float64, but that can
be changed by altering the IMGSIZE variable.  Record marker length is also taken to be 4, but that to
can be changed with the RM variable.'''

import numpy
import math
import array
import struct
import pylab
import sys

# file-dependent variables.
RM=4
IMGSIZE='d' # change to 'f' for float32
SWAP = 0 # swap bytes

i = len(sys.argv)

if i < 9:
	print "Incorrect Usage. Received ", i-1, " arguments."
	print "$: mdot.py file1 file2 time_conversion(yr) mass_conv(Msun) dx JMAX KMAX LMAX"
        sys.exit()
if i > 9:
	print "Incorrect Usage. Received ", i-1, " arguments."
	print "$: mdot.py file1 file2 time_conversion(yr) mass_conv(Msun) dx JMAX KMAX LMAX"
	sys.exit()

jmax = int(sys.argv[6])
kmax = int(sys.argv[7])
lmax = int(sys.argv[8])
jmax2 = jmax+2
kmax2 = kmax+2

dx =float( sys.argv[5])
mconv = float(sys.argv[4])
dphi = 2.*math.pi/float(lmax)

f = open(sys.argv[1],'rb')
h = open(sys.argv[2],'rb')

s = f.read(RM)
if SWAP == 1:
	par = struct.unpack(">i",s)
else:
	par = struct.unpack("i",s)

print "Data size according to record in file ", sys.argv[1],".", par[0]

binvalues = array.array(IMGSIZE)
binvalues.read(f,jmax2*kmax2*lmax)
if SWAP == 1: binvalues.byteswap()
data0=numpy.array(binvalues,'float64')
data0 = numpy.reshape(data0,(lmax,kmax2,jmax2))

s = f.read(2*RM+struct.calcsize(IMGSIZE))

if SWAP == 1:
	par,par,time0 = struct.unpack(">2id",s)
else:
	par,par,time0 = struct.unpack("2id",s)

time0 *= float(sys.argv[3])
print time0

s = h.read(RM)
if SWAP == 1:
        par = struct.unpack(">i",s)
else:
        par = struct.unpack("i",s)

print "Data size according to record in file ", sys.argv[2],".", par[0]

binvalues = array.array(IMGSIZE)
binvalues.read(h,jmax2*kmax2*lmax)
if SWAP == 1: binvalues.byteswap()
data1 = numpy.array(binvalues,'float64')
data1 = numpy.reshape(data1,(lmax,kmax2,jmax2))

s = h.read(2*RM+struct.calcsize(IMGSIZE))
if SWAP == 1:
        par,par,time1 = struct.unpack(">2id",s)
else:
        par,par,time1 = struct.unpack("2id",s)

time1 *= float(sys.argv[3])

print time1

# clean up file handles.
f.close()
h.close()

# define some new arrays

mass0    = numpy.zeros((jmax2),'float')
mass1    = numpy.zeros((jmax2),'float')
diff     = numpy.zeros((jmax2),'float')
summass0  = numpy.zeros((jmax2),'float')
summass1  = numpy.zeros((jmax2),'float')
rad  = numpy.zeros((jmax2),'float')
tmass0 = 0.
tmass1 = 0.

for j in xrange(1,jmax2-1):
	fj = float(j)
	for k in xrange(1,kmax2-1):
		for l in xrange(lmax):
			dum0 = data0[l,k,j]*dx**3*(fj-.5)*2.*mconv*dphi
			dum1 = data1[l,k,j]*dx**3*(fj-.5)*2.*mconv*dphi
			tmass0   += dum0
			tmass1   += dum1
			mass0[j] += dum0
			mass1[j] += dum1

print "Total mass in file ", sys.argv[1], " is ",tmass0, "."
print "Total mass in file ", sys.argv[1], " is ",tmass1, "."

for j in xrange(1,jmax2):
	summass0[j] = summass0[j-1]+mass0[j]
	summass1[j] = summass1[j-1]+mass1[j]
	diff[j] = (summass1[j]-summass0[j])/(time1-time0)
	rad[j]=(float(j)-.5)*dx


pylab.figure()
pylab.plot(rad,mass0,'r')
pylab.plot(rad,mass1),'b'
pylab.figure()
pylab.plot(rad,summass0,'r')
pylab.plot(rad,summass1,'b')
pylab.figure()
pylab.plot(rad,diff)
pylab.yaxis(format="%5.2e")
pylab.show()
#pylab.savefig(sys.argv[1]+".png")

#
