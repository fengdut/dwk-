#!/usr/bin/env python


from sys import argv
import csv
import numpy as np

if (len(argv)<2):
	print('The default data file is omega_dwk.out')
	filename='omega_dwk.out'
else:
	filename=argv[1]


returnMat	=	np.loadtxt(filename)

	
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.pyplot import rc, figure, axes, plot, xlabel, ylabel, title,\
	grid, savefig, show

omegar	=	returnMat[:,0]
omegai	=	returnMat[:,1]
dwkr	=	returnMat[:,2]
dwki	=	returnMat[:,3]

font = {'size' : '18'}
rc('font',**font)
rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman']})
rc('text', usetex=True)
#fig1=plt.figure(figsize=(6,6),)	
fig1=plt.figure()	
rect = fig1.patch
rect.set_facecolor('white')

plt.plot(omegar,dwkr)
plt.xlabel(r'$\omega$')
plt.ylabel(r'$\delta W_k$')

plt.show()

