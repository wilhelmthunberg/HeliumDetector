import numpy as np
def multFactor(a,b,V,P,dV=27.6,K=1.48):
    """
    Function that calculates the multiplication factor of a detector. 
    Presetes of dV and K taken from Knoll.
    Param:
          a   : (int, double,float) radius of anode [m] 
          b   : (int, double,float) radius of cathode [m]
          V   : (int, double,float) applied voltage [V]
          P   : (int, double,float) Pressure inside detector [Pa]
          dV  : (int, double,float) average voltage increase between ionizing events. [V]

    Out: 
          M: (float) Multiplication factor
    """
    lnM = V/(np.log(b/a)) * np.log(2)/dV * (np.log(V/(P*a*np.log(b/a))) - np.log(K))
    M= np.exp(lnM)
    return M