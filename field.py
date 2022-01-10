import numpy as np
def field(a,b,V,r):
  """
  Function that returns electric field at specific radius in a cylindrical detector

  Param:
        r_a: (int, double,float) radius of anode [m]
        r_c: (int, double,float) radius of cathode [m]
        V  : (int, double,float) Applied voltage [V]
        r  : (int, double,float) radial position of electron [m]

  Out: 
      E: (int, double,float) Electric field at r

  """
  E= V/(r*np.log(b/a))
  return E