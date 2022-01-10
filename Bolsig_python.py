# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 12:04:10 2021

@author: aronj
"""

from scipy import interpolate
import subprocess as sp
import numpy as np
import time
import matplotlib.pyplot as plt
import os
import sys

# path=r'C:/Users/aronj/OneDrive/Dokument/Skola/AR5/PROJEKTITILLAMPADFYSIK/runs/bolsigminus.exe'
path=sys.path[0]+'/bolsig/bolsigminus'

# run_file_original=r'C:/Users/aronj/OneDrive/Dokument/Skola/AR5/PROJEKTITILLAMPADFYSIK/runs/input-master.dat'
run_file_original=sys.path[0]+'/bolsig/input-master.dat'

query_E=('Electric field')
query_save=('bolsig/testrun_0.dat')

def mobilityCalc(E_start, E_end, steps):
    mobility=[]
    E_vec=np.linspace(E_start,E_end,steps)
    for i in range(0,steps):
        # run_file_copy=r'C:/Users/aronj/OneDrive/Dokument/Skola/AR5/PROJEKTITILLAMPADFYSIK/runs/input-runseries_'+str(i)+'_E_min_'+str(E_start)+'_E_max'+str(E_end)+'.dat'
        run_file_copy=sys.path[0]+'/bolsig/input-runseries_'+str(i)+'_E_min_'+str(E_start)+'_E_max'+str(E_end)+'.dat'

        if os.path.exists(run_file_copy):
            pass
        else:
            fin = open(run_file_original)
            fout = open(run_file_copy, "wt")
            for line in fin:
                if query_E in line:
                    fout.write( line.replace(line, str(E_vec[i])+"       / Electric field / N (Td) \n") )
                elif query_save in line:
                    fout.write( line.replace(line,sys.path[0]+'/bolsig/testrun_'+str(i)+'_E_min_'+str(E_start)+'_E_max'+str(E_end)+'.dat        / File \n'))
                else:
                    fout.write(line)
            fin.close()
            fout.close()
            print(run_file_copy)
            args=[path, run_file_copy]
            process = sp.Popen(args, stdin=sp.PIPE, stdout=sp.PIPE)
            process.stdin.close()
            time.sleep(2)
    for i in range(0,steps):
        # result_file=r'C:/Users/aronj/OneDrive/Dokument/Skola/AR5/PROJEKTITILLAMPADFYSIK/runs/testrun_'+str(i)+'_E_min_'+str(E_start)+'_E_max'+str(E_end)+'.dat'
        result_file=sys.path[0]+'/bolsig/testrun_'+str(i)+'_E_min_'+str(E_start)+'_E_max'+str(E_end)+'.dat'
        with open(result_file,'r') as fin:
            mobility_string='E/N (Td)	Mobility *N (1/m/V/s)'
            lines=fin.readlines()
            for j in range(0,len(lines)):
                line=lines[j]
                if mobility_string in line:
                    mobility.append(float(lines[j+1].split('\t')[1].strip()))
    mob_intpltd = interpolate.interp1d(E_vec,mobility,kind='slinear',  fill_value="extrapolate")
    return mob_intpltd, E_vec
