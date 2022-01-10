import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt
from detector import detector
from field import field ; from multFactor import multFactor ; from driftVel import driftVel; from inducedCurrent import inducedCurrent; from numDens import numDens
def root_eulerBack(sig_next,I,C,dt,tau,sig):
  """
  Function used for newton raphson used in the euler backward method in signal()
  sig_next is unknown and the variable of interest.

  See "Signal Formation in Various Detectors" eq.22  for derivation of equation
  Param: 
        sig_next: (int,double, float) next value of signal [V]
        I       : (int,double, float) induced current [C/s]
        C       : (int,double, float) Capactiance of detector [F]
        dt      : (int,double, float) time step [s]
        tau     : (int,double, float) R*C [ohm*F]
        sig     : (int,double, float) Current value of signal [V]
  
  """
  return sig_next-dt*((I/C) - sig_next/tau)-sig
def f(I,sig, C, tau):
  return (I/C) - (sig/tau)

def signal(det,mu_e,mu_i,x,y,z,dr,q=1.60217662 *10**-19,t_start=0):
  """
  Function that returns array containing signal (pulse) values.
  
  Param: 
        det    : (Class object) detector
        x      : (int,double, float) x-coord of particle in det
        y      : (int,double, float) y-coord of particle in det
        z      : (int,double, float) z-coord of particle in det
        dr     : (int,double, float) spatial step used in numerical calculation of current and output signal.
        q      : (int,double, float) charge of particle that produces signal
        t_start: (int,double, float) time when charge is introduced/created in detector
  Out: 
      (list) [I,sig_intpltd, ts]: 
  
  """
  r= np.sqrt(x**2 + y**2)
  if r > det.r_c:
    raise ValueError('The particle is outside detector. \nr is: ', r, '\n det.r_c is: ',det.r_c)
  tau = det.R*det.C
  V0=det.V
  #initial values
  E0 = field(det.r_a,det.r_c,det.V,r)
  v0=0
  t0 = t_start
  I0 = 0 
  sig0 =  0
  #initiate lists to store values
  ts= [t0] #time
  rs = [r] #radii
  vs=[v0] #velocity
  Es=[E0] # E-field
  Is = [I0] # induced Current
  sigs = [sig0] # Signal

  r_prev=r  
  avalanche=False # used when shortening dr near center of anode
  last_e =False # used to determine final contribution from electrons
  #Move electron towards center
  # Calculate Signal change using Euler-backward method
  # Stop first loop when eletrons reach anode
  while last_e==False:
    if r-dr<det.r_a:
      last_e =True #final contribution always in same pos since input-dr does not changexw
      r_final =det.r_a +dr*0.1; dr_final = abs(r_prev-r_final)
      E = field(det.r_a,det.r_c,det.V,r_final)
      v = driftVel(E,(mu_e(E)/det.N))
      dt = dr_final/v
    else:
      E = field(det.r_a,det.r_c,det.V,r)
      v = driftVel(E,(mu_e(E)/det.N))
      dt=dr/v
    t=ts[-1]+dt  
    I = inducedCurrent(V0, q , E,v)
    if ((E/det.P)>det.K):#avalanche effect
      if avalanche==False:
        #r_av=r
        avalanche=True
        dr/=100
      v_i =driftVel(E,mu_i)      #M=det.M*np.exp((det.r_a-r)/r_av) # M exponetially increases to max M as closer to r_a
      #M_next = det.M*np.exp((det.r_a-(r-dr))/r_av)
      I=det.M*(I + inducedCurrent(V0, q , E,v_i)) #add current from ions 
    sig = newton(root_eulerBack, x0=sigs[-1], args= (I,det.C,dt,tau,sigs[-1]),tol=1e-8 )# Euler-backward method, See "Signal Formation in Various Detectors" eq.22 for equation  
    rs.append(r);Es.append(E);vs.append(v);Is.append(I); ts.append(t); sigs.append(sig)
    r_prev=r
    r-=dr
  r+=dr
  # Calculate induced current from ions moving towards cathode
  while r<det.r_c*0.4:
    E = field(det.r_a,det.r_c,det.V,r)
    v_i = driftVel(E,mu_i)
    dt = dr/v_i
    t=ts[-1]+dt 
    I = det.M*inducedCurrent(V0, q , E,v_i)
    sig = newton(root_eulerBack, x0=sigs[-1], args= (I,det.C,dt,tau,sigs[-1]),tol=1e-10 ) 
    Is.append(I); ts.append(t); sigs.append(sig)
    dr*=1.2
    r+=dr
  return [Is,sigs,ts]
  
det = detector(r_a=0.0055e-3,r_c=(2.54e-2)/2,L=3e-2,V=1.1e3,P=4e5,T=300, R=1.5e6,tau=6.5e-6)
Is,sigs,ts = signal(det, det.Mob, det.mu_i, det.r_c/2, det.r_c/2, det.r_c/1.1, 1e-5)
plt.plot(np.array(ts)*1e6, sigs)
plt.xlabel('Time '+r'[$\mu s$]')
plt.ylabel('Output voltage [V]')
plt.title('Pulse genrated by a single charge. Tau='+str(round(det.R*det.C*1e6,2) )+r'$\mu s$')
plt.show()
