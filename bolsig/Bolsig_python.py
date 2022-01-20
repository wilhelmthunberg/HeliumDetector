
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 12:04:10 2021

@author: aronj
! IF FILENOTFOUND ERRORS ARE ENCOUNTERED, MAKE SURE YOU RUN THE MAIN METHOD FROM THE 'CODE' FOLDER!
"""

from scipy import interpolate

import subprocess as sp
import numpy as np
import time
import matplotlib.pyplot as plt
import os
import sys
os.chdir(sys.path[0]+'/bolsig')
#!IF FILENOTFOUND ERRORS ARE ENCOUNTERED, MAKE SURE YOU RUN THE MAIN METHOD FROM THE 'CODE' FOLDER!
path=os.getcwd()+'/bolsigminus'

run_file_original='input-master.dat'

query_E=('Electric field')
query_save=('testrun_0.dat')

def mobilityCalc(E_start, E_end, steps):
    mobility=[]
    E_vec=np.linspace(E_start,E_end,steps)
    for i in range(0,steps):
        run_file_copy='input-runseries_'+str(i)+'_E_min_'+str(E_start)+'_E_max'+str(E_end)+'.dat'

        if os.path.exists(run_file_copy):
            pass
        else:
            fin = open(run_file_original)
            fout = open(run_file_copy, "wt")
            for line in fin:
                if query_E in line:
                    fout.write( line.replace(line, str(E_vec[i])+"       / Electric field / N (Td) \n") )
                elif 'testrun_0.dat' in line:
                    fout.write( line.replace(line,'testrun_'+str(i)+'_E_min_'+str(E_start)+'_E_max'+str(E_end)+'.dat        / File \n'))
                else:
                    fout.write(line)
            fin.close()
            fout.close()
            #calculate mobility using bolsig in this specific section
            args=[path, run_file_copy]
            process = sp.Popen(args, stdin=sp.PIPE, stdout=sp.PIPE)
            process.stdin.close()
            time.sleep(0.61)
            os.remove(run_file_copy)# remove file

    for i in range(0,steps):
        result_file='testrun_'+str(i)+'_E_min_'+str(E_start)+'_E_max'+str(E_end)+'.dat'
        with open(result_file,'r') as fin:
            mobility_string='E/N (Td)	Mobility *N (1/m/V/s)'
            lines=fin.readlines()
            for j in range(0,len(lines)):
                line=lines[j]
                if mobility_string in line:
                    mobility.append(float(lines[j+1].split('\t')[1].strip()))
            time.sleep(0.1)
            os.remove(result_file) #remove result_file after being used
    mob_intpltd = interpolate.interp1d(E_vec,mobility,kind='slinear',  fill_value="extrapolate") #linear interpolated function for 
    
    return mob_intpltd, E_vec
