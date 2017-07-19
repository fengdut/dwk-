#!/usr/bin/env python


from sys import argv
import numpy as np

if (len(argv)<2):
        print('The default data file is dwk.cfg')
        filename='dwk.cfg'
else:
        filename=argv[1]


F       =open(filename,'r')
FN      =open('dwk_new.cfg','w');

domegai=0.001


for line in F:
        ti=line.find('omega_i')
        if(ti==-1):
                FN.write(line)
        else:
                te=line.find('=')
                tend=line.find(';');
                stromegai=line[te+1:tend]
                omegaiold=float(stromegai)
                omegainew=omegaiold +domegai
                print stromegai
                newomegai='\tomega_i='+str(omegainew)+';\n'
                FN.write(newomegai)

