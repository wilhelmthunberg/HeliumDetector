import matplotlib.pyplot as plt
import pandas as pd
import os 
import sys
import numpy as np
with open(os.path.join(sys.path[0], 'results_tau=1.0775076213251382e-05_gamma_3months.csv'), "r") as f:
    df = pd.read_csv(f)
    t = np.array(df['Time'])
    for det in range(0,18):
        plt.figure(det)
        plt.plot(t*1e6, df['Detector_'+str(det)])
        plt.title('Accumulated signal from gamma radiation in detector '+str(det))
        plt.xlabel('Time [micro seconds]')
        plt.ylabel('Voltage of pulse [V]')
        plt.savefig(os.path.join(sys.path[0],'plots/plot_detector_'+str(det)+'_gamma'))


with open(os.path.join(sys.path[0], 'results_tau=1.0775076213251382e-05_neutron_3months.csv'), "r") as f:
    df = pd.read_csv(f)
    t = np.array(df['Time'])
    for det in range(0,18):
        plt.figure(det+100)
        plt.plot(t*1e6, df['Detector_'+str(det)])
        plt.title('Accumulated signal from neutron radiation in detector '+str(det))
        plt.xlabel('Time [micro seconds]')
        plt.ylabel('Voltage of pulse [V]')
        plt.savefig(os.path.join(sys.path[0],'plots/plot_detector_'+str(det)+'_neutron'))
