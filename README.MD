# Helium DEtector
A model of gamma-ray interference in a $^3$He neutron detector in high count-rate environments
README
Authors:
Aron Johanson & Wilhelm Thunberg

The following document describes how to use the code written for project 16
of the applied physics course.

The code is structured as a set of files each containing a necessary function.
"detector.py" contains a class object for the detector.
The necessary files and their home directory are as follows:


## Run simulation
To run the simulation launch main.py
The main file is a script. Input data and detector parameters etc. are defined inside the code.
In the script, the input data is gathered from the directory *data* along with a **track-fraction** for how many percent of the rows in the input data should be used.
**det_of_interest** defines which detectors to simulate. For example det_of_interest=range(1,18)
or det_of_interest=[4,7,9].
The characteristics of the detector is also provided here. All units are in *SI units*.

**main.py** calls upon **runTracks.py** that reads the input files and launches **chargeTrail.py**
which in turn calculates the ionisation along a track (a single line in the input data file) and divides it into segements.
For each segment, **pulse.py** is called upon to calculate the induced current and voltage.
**pulse.py** solves a differential equation using the Euler backwards method.

The results from main are saved under the directory **main/results/** as .csv files. The file names indicate
what data they contain. Using **plotResults.py** the data can be plotted.

## Important parameters

- *specify input data* - > **main.py**

- *detector parameters* (Geometry, Resistance, Capacitance and Voltage) - > **main.py**

- *tau* -> USer can change tau in **main.py**

- *num_seg* -> **chargeTrail.py** (Calculated by dividing total energy deposited along track with 40KeV
which gives roughly 240 ions per segment. Maximum num_seg is 10, if track is relatively long (0.1*L_det) the minium is 3.)

- *dr* -> pulse.py (The legnth of the segments that the detector radius is divided into.
dr governs the time step according to dt=dr/v and is decreases close to the anode.)

- radiation detection rate -> **runTracks.py** from **data/Count_rates.csv**

Multiplication factor -> **detector.py**

t_global (for gamma and neutron) -> **main.py** 


## File structure
/main
main.py
chargeTrail.py
driftVel.py
field.py
inducedCurrent.py
main.py
multFactor.py
numDens.py
processResults.py
runTracks.py
pulse.py 

/main/bolsig
Bolsig_python.py
bolsigminus
SigloDataBase-LXCat-04Jun2013.txt
LXCat-June2013

/main/data
Count_rates.csv
Gamma_3months_tracks.csv
Gamma_5yr_tracks
SF_3months_tracks.csv
SF_5yr_tracks.csv

/main/results/

