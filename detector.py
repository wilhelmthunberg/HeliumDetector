import numpy as np
from field import field ; from multFactor import multFactor ; from driftVel import driftVel; from inducedCurrent import inducedCurrent; from numDens import numDens

from Bolsig_python import mobilityCalc
class detector:
  def __init__(self,r_a,r_c,L,V,P,T,R,tau,K=1.48):
    """
    Initialize the detector

    Param: 
          r_a: (int, double,float) radius of anode [m]
          r_c: (int, double,float) radius of cathode [m]
          V  : (int, double,float) Applied voltage [V]
          P  : (int, double,float) Pressure  [Pa]
          R  : (int, double,float) Resistance [Ohm]  
          T  : (int, double,float) Temperature [K]
          K  : (int, double,float) E/P threshold for avalanche effect [Ohm]  
    """
    self.r_a = r_a
    self.r_c = r_c
    self.L = L
    self.V =V
    self.P = P
    self.K = 1.48
    self.C = tau/R
    self.R = R
    self.N = numDens(P,T)
    self.M = multFactor(r_a,r_c,V,P)
    self.mu_i = 1e-4*10.4* numDens(P=1.1e5,T=293)/numDens(P=P,T=T)
    self.E_min=field(r_a, r_c, V, r_c)/(numDens(P=P, T=T)*1e-21)
    self.E_max=field(r_a, r_c, V, r_a)/(numDens(P=P, T=T)*1e-21)
    self.Mob=mobilityCalc(self.E_min,self.E_max, 10)[0]    
