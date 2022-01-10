
def  inducedCurrent(V, q ,E,v): 
    """
    Function that calculates induced current from free moving chrge in E-field using shockley

    Param:
          V: (int,double, float) Applied voltage over entire E-field [V]
          q: (int,double, float) charge of particle [Coulomb]
          E: (int,double, float) value of electric field where the charge is. [V/m]
          v: (int,double, float) drift velocity of charge. [m/s]  

    """
    return (1/V)* q*E*v 