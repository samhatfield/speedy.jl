const γ      = RealType(6.0)       # Reference temperature lapse rate (-dT/dz in deg/km)
const hscale = RealType(7.5)       # Reference scale height for pressure (in km)
const hshum  = RealType(2.5)       # Reference scale height for specific humidity (in km)
const refrh1 = RealType(0.7)       # Reference relative humidity of near-surface air
const thd    = RealType(2.4)       # Max damping time (in hours) for horizontal diffusion
                                      # (del^6) of temperature and vorticity
const thdd   = RealType(2.4)       # Max damping time (in hours) for horizontal diffusion
                                      # (del^6) of divergence
const thds   = RealType(12.0)      # Max damping time (in hours) for extra diffusion (del^2)
                                      # in the stratosphere
const tdrs   = RealType(24.0*30.0) # Damping time (in hours) for drag on zonal-mean wind
                                      # in the stratosphere