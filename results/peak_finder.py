import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import special
import pandas as pd
from scipy.signal import find_peaks
import scipy.stats as stats
import os

os.chdir('/Users/w/Project_AppliedPhys/Code/results_dec09')


## Get the raw energy deposition (per event per detector)
df_raw = pd.read_csv('../data/SF_5yr_tracks.csv')
df_raw['Delta_E'] = df_raw['Start_E'] - df_raw['End_E']
df_raw_grouped = df_raw.groupby(['NPS', 'Detector_no']).sum()
Edep_raw_n = df_raw_grouped['Delta_E'].to_numpy()

df_raw = pd.read_csv('../data/Gamma_5yr_tracks.csv')
df_raw['Delta_E'] = df_raw['Start_E'] - df_raw['End_E']
df_raw_grouped = df_raw.groupby(['NPS', 'Detector_no']).sum()
Edep_raw_g = df_raw_grouped['Delta_E'].to_numpy()


# Everything done for CT=5yr at the moment

peak_width = 20             ## required to use the scipy peak finder and not trigger on rising-edge noise


# Read the neutron pulse train:
df_n = pd.read_csv('SF/tau=1.0775076213251382e-05/count_rate=15925.0')

t_n = 1e6*df_n['Time']

V_n = [None]*18
peaks_n = [None]*18
peak_V_n = [None]*18

for i in range(0, 18):
    V_n[i] = df_n['Detector_' + str(i)]
    peaks_n[i], _ = find_peaks(V_n[i], height=0, width=peak_width)        ## width=20 to avoid identifying jagged rising edge as multiple peaks

    peak_V_n[i] = V_n[i][peaks_n[i]]


# Read the gamma pulse train:
# I processed multiples of np.linspace(1, 50, 10) of the original rate
neutron_rate = 15925
#gamma_rates = ['1179947.0', '15634293.0', '30088640.0', '44542987.0', '58997333.0']
gamma_rates = ['1179947.0', '39724871.0','46149025.0', '52573179.0', '58997333.0']
N_gamma_rates = len(gamma_rates)

df_g = [None]*N_gamma_rates
t_g = [ [ [] for i in range(18) ] for i in range(N_gamma_rates) ]
V_g = [ [ [] for i in range(18) ] for i in range(N_gamma_rates) ]
peaks_g = [ [ [] for i in range(18) ] for i in range(N_gamma_rates) ]
peak_V_g = [ [ [] for i in range(18) ] for i in range(N_gamma_rates) ]
V_dist_g = [ [ [] for i in range(18) ] for i in range(N_gamma_rates) ]              ## "all voltages" (from element 200 to 1400) (including dips) in the gamma pulse train


for i in range(N_gamma_rates):
    df_g[i] = pd.read_csv('gamma/tau=1.0775076213251382e-05/count_rate=' + gamma_rates[i]) 
    t_g[i] = 1e6*df_g[i]['Time']

    for j in range(0, 18):
        V_g[i][j] = df_g[i]['Detector_' + str(j)]
        peaks_g[i][j], _ = find_peaks(V_g[i][j], height=0, width=peak_width)        ## width=20 to avoid identifying jagged rising edge as multiple peaks
        peak_V_g[i][j] = V_g[i][j][peaks_g[i][j]]
        V_dist_g[i][j] = V_g[i][j][200:1400]





# Currently, there is one sublist for each detector. I want a list of all pulse heights (equivalent to summing over all detectors). To do that, flatten lists.
peak_V_g_flat_list = [ [ ] for i in range(N_gamma_rates) ]
V_dist_g_flat_list = [ [ ] for i in range(N_gamma_rates) ]
peak_V_n_flat_list = []

for i in range(N_gamma_rates):
    for sublist in peak_V_g[i]:
        if sublist is not None:
            for item in sublist:
                peak_V_g_flat_list[i].append(item)

for i in range(N_gamma_rates):
    for sublist in V_dist_g[i]:
        if sublist is not None:
            for item in sublist:
                V_dist_g_flat_list[i].append(item)

for sublist in peak_V_n:
    if sublist is not None:
        for item in sublist:
            peak_V_n_flat_list.append(item)


for i in range(N_gamma_rates):
    np.savetxt("gamma_pulse_heights_" + gamma_rates[i] + ".csv", peak_V_g_flat_list[i], delimiter=",")

np.savetxt("neutron_pulse_heights_15925.0.csv", peak_V_n_flat_list, delimiter=",")

f, ax = plt.subplots(nrows=N_gamma_rates+1, sharex=True)
n, bins, patches = ax[0].hist(peak_V_n_flat_list, bins=200, range=[0, 700], density=True, facecolor='r', alpha=0.75, label='SF neutrons')
ax[0].set_yscale('log')
ax[0].text(0.9, 0.5,'Neutrons', horizontalalignment='center', verticalalignment='center', transform = ax[0].transAxes)



for i in range(N_gamma_rates):
    n, bins, patches = ax[1+i].hist(peak_V_g_flat_list[i], bins=200, range=[0, 3000], density=True, facecolor='g', alpha=0.75, label='Gammas, rate ' + str(i))
    n, bins, patches = ax[1+i].hist(V_dist_g_flat_list[i], bins=200, range=[0, 3000], density=True, facecolor='y', alpha=0.75, label= '"All voltages"')
    ax[1+i].set_yscale('log')
    ax[1+i].text(0.78, 0.5,'Gammas (hit rate: ' + "{:.2f}".format(float(gamma_rates[i])/1e6) + ' MHz)', horizontalalignment='center', verticalalignment='center', transform = ax[1+i].transAxes)



plt.xlabel('Pulse height [a. u.]')
f.text(0.02, 0.5, 'Counts per bin (not normalised)', va='center', rotation='vertical')
plt.suptitle('Pulse height per detector')

### Draw histo of raw energy deposits
plt.figure()
n_raw, bins_raw, patches_raw = plt.hist(Edep_raw_g, bins=200, range=[0, 1.5], density=True, facecolor='g', alpha=0.75, label='Gammas')
n_raw, bins_raw, patches_raw = plt.hist(Edep_raw_n, bins=200, range=[0, 1.5], density=True, facecolor='r', alpha=0.75, label='SF neutrons')
plt.xlabel('Energy deposition [MeV]')
plt.ylabel('Counts per bin (not normalised)')
plt.title('Energy deposition per detector')
plt.legend()


f, ax = plt.subplots(18)
for i in range(0, 18):
    ax[i].plot(t_g[0], V_g[0][i], linestyle='-', color='red', linewidth=2)
    ax[i].plot(t_g[0][peaks_g[0][i]], V_g[0][i][peaks_g[0][i]], "x")



fig = plt.figure()
plt.plot(t_g[-1], V_g[-1][1], linestyle='-', color='red', linewidth=2)
plt.plot(t_g[-1][peaks_g[-1][1]], V_g[-1][1][peaks_g[-1][1]], "x")
plt.plot(t_g[-1][200:1400], V_g[-1][1][200:1400], linestyle='-', color='black', linewidth=2)

fig = plt.figure()
plt.plot(t_n, V_n[1], linestyle='-', color='red', linewidth=2)
plt.plot(t_n[peaks_n[1]], V_n[1][peaks_n[1]], "x")






##
#normalize
# How to generate continuos prob dist from hist?

from collections import Iterable
def flatten(lis):
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, str):
             for x in flatten(item):
                 yield x
         else:
             yield item

data = np.array(list(flatten(peak_V_g_flat_list)))

# generate a list of kde estimators for each bw
loc1, scale1 = stats.norm.fit(data) 

dist = stats.norm.pdf(data)

plt.figure()
plt.plot(data,dist)

plt.show()