#!/usr/bin/env python

omega_dwk=[]
datafile=open('dwk_omega_dwk.out');
data=datafile.readlines()
num =len(data)

import numpy
returnMat = numpy.zeros((num,4))
index = 0
for line in data:
	line = line.strip()
	linelist = line.split()
	returnMat[index,:] = linelist[0:4]
	index +=1


import matplotlib.pyplot as plt
import matplotlib
from matplotlib.pyplot import rc, figure, axes, plot, xlabel, ylabel, title, \
	grid, savefig, show

omega_r =returnMat[:,0]
omega_i =returnMat[:,1]
dwk_r =returnMat[:,2]
dwk_i =returnMat[:,3]


font ={ 'size' : '18'}
rc('font',**font)
rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman']})

rc('text', usetex=True)
fig1=plt.figure(figsize=(8,8.5),)

rect = fig1.patch
rect.set_facecolor('white')
plt.subplot(2,1,1)

plt.plot(omega_r,dwk_r,'ro',markersize=6)
plt.xlabel(r'$\omega$')
plt.ylabel(r'$real(dwk)$')


ax=plt.subplot(2,1,2,axisbg='white')
plt.plot(omega_r,dwk_i,'ro',markersize=6)

plt.xlabel(r'$\omega$')
plt.ylabel(r'$imag(dwk)$')

plt.show()
