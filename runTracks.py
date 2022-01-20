    #track in detector 
import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import expon
from chargeTrack import chargeTrack
import random
from scipy import interpolate

def runTracks(det,filename,rateFile,cutOff,t_start,detCoords, det_of_interest, rate=-1):

  with open(os.path.join(sys.path[0], rateFile), "r") as rf:
    df_rate = pd.read_csv(rf)[filename.split('_tracks')[0].split('/')[-1]+' [cps]']
    dt=[]
    if rate>0:
      for i in range(0,len(df_rate)): 
        dt.append(expon.rvs(scale=1/rate, loc=0, size=1000))
    else:
      for i in range(0,len(df_rate)): 
        dt.append(expon.rvs(scale=1/df_rate[i], loc=0, size=1000))
      print('Specific rate not given or less than 0, gathered from file instead\n')
  with open(os.path.join(sys.path[0], filename), "r") as f:
    df=pd.read_csv(f)
    #array of timesteps from exponetial dist based on rate param.
  #Generate signal from gamma track 0
  t=t_start
  x=0
  event = df['NPS'][x]
  allSig=[]
  ts= []
  detector=[]
  Is = []
  NPS =[]
  t0=[]
  for x in range(0,int(len(df['Start_E'])/cutOff)):
    if (df['NPS'][x]!=event):
      event=df['NPS'][x]
      t+= dt[int(df['Detector_no'][x])][random.randint(0,len(dt)-1)]
    print(x,' Calculating pulses from ',filename,' ... ', int(cutOff*x/len(df['Start_E'])*100), '%',end= '\r')
    if (df['Delta_ray_flag'][x]==1) or (df['Detector_no'][x] not in det_of_interest) :
      print(end= '\r')
      #do nothing
    else:
      start_x =  (df['Start_x'][x]-detCoords[int(df['Detector_no'][x])][0])*1e-2;  end_x =  (df['End_x'][x]-detCoords[int(df['Detector_no'][x])][0])*1e-2
      start_y =  (df['Start_y'][x]-detCoords[int(df['Detector_no'][x])][1])*1e-2;  end_y =  (df['End_y'][x]-detCoords[int(df['Detector_no'][x])][1])*1e-2
      I, sigs, t_sigs, num_seg=chargeTrack(det,df['Start_E'][x]*1e6,df['End_E'][x]*1e6,start_x, end_x, start_y, end_y, df['Start_z'][x]*1e-2, df['End_z'][x]*1e-2,t,t)
      allSig +=sigs
      ts+=t_sigs
      Is+=I
      for i in range(0,len(sigs)):
        detector.append(int(df['Detector_no'][x]))
        NPS.append(event)
        t0.append(t)
  print('\r\nDone')
  
  return NPS, t0 ,detector,dt,Is, ts, allSig

