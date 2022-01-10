def numDens(P,T):
  """
  Function that returns the number density.
  
  Param: 
        P: (int,double, float) Pressue [Pa]
        T: (int,double, float) Temperature [K]
  
  Out: Particles per volume [number/m^3]
  
  """
  R=8.31446261815324
  N_A = 6.022*10**23	
  N=N_A*P/(R*T)
  return N