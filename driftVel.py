def driftVel(E,mu=1):
  """
  Function that calculates drift velocity of charged particle

  Param:
        E: (int,double, float) Value of electric field.
        mu: (int,double, float) electron mobility in given position 

  Out: 
        E*mu:(int,double, float) electron drift velocity

  """
  return E*mu