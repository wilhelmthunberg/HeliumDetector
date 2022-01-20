import numpy as np
import math
from field import field ; from multFactor import multFactor ; from driftVel import driftVel; from inducedCurrent import inducedCurrent; from numDens import numDens
from pulse import pulse
def chargeTrack(det,E_in, E_out, x_in, x_out, y_in, y_out, z_in, z_out, t_in, t_out):
    """
    Function that models the ionisation of the gas along a photoelectron's path inside the volume.
    The path is divided into a number of segments and the electrons created in each segment is summed into
    a single charge. The modelled charge is placed in the middle of the segment and its pulse response 
    is calculated using the function pulse
    
    Param: 
    det= (Class object) Detector
    E_in/E_out= Energy at entry/exit of volume [eV]
    x,y,z, coordinates at entry/exit [m]
    t_in/t_out timestamp at entry/exit [s]
    num_seg = number of segments that the path should be divided into [-]
    
    """
    q=1.60217662*10**-19    
    ion_E=42.3 #Ionisation energy
    seg_E=40e3 #Energy deposited per segment 
    deltaE=abs(E_in-E_out) #Photoelectron energy difference at entry/exit
    if t_out == t_in:
      t_out = t_in+1e-11
    T=t_out-t_in #Time inside volume
    L=np.sqrt((x_out-x_in)**2+(y_out-y_in)**2+(z_out-z_in)**2) #Length of path inside volume
    num_seg=(int(math.ceil((deltaE/seg_E))))
    if num_seg>10:
      num_seg=10
    if L>det.L*0.05 and num_seg<3:
      num_seg=3
    num_ion=deltaE/ion_E #Number of ions created along path
    q_tot_L=num_ion*q #Total charge
    q_dL=q_tot_L/(num_seg) #Charge per segment
    #print('deltaE/numseg: ', deltaE/num_seg)
    #print('q_dL')
    #step size in x,y and z dir
    dx= -(x_in - x_out)/num_seg  
    dy= -(y_in - y_out)/num_seg
    dz= -(z_in - z_out)/num_seg
    #initial step
    x=x_in + dx/2
    y=y_in +dy/2
    z=z_in+dz/2
    #print('numseg: ',num_seg ,'\nx_in:', x_in, '\nx_out: ',x_out,'\nx: ',x, '\ndx',dx)
    #input()
    vel = L/T #Velocity in
    dt=0 #Time between the segments being initialised

    sigs_vec=[]
    t_sig_vec=[]
    Is = []
    #from 1 because initial step has already been made
    for i in range(0,num_seg):
        t=t_in+i*(dt/num_seg)
        #print('\nx_in:', x_in, '\nx_out: ',x_out,'\nx: ',x, '\ndx',dx)
        #print('\ny_in:', y_in, '\ny_out: ',y_out,'\nx: ',y, '\ndx',dy)
        #print('\nr: ',np.sqrt(x**2))
        #input()
        I,sigs,t_sig=pulse(det,mu_e=det.Mob,mu_i=det.mu_i,x=x,y=y,z=z,dr=1e-4,q=q_dL,t_start=t)
        sigs_vec.append(sigs)
        t_sig_vec.append(t_sig)
        Is.append(I)
        x+=dx;y+=dy;z+=dz

  
    return Is,sigs_vec, t_sig_vec, num_seg
