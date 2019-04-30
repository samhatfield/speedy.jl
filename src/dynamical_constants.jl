const γ      = Real(6.0)       # Reference temperature lapse rate (-dT/dz in deg/km)
const hscale = Real(7.5)       # Reference scale height for pressure (in km)
const hshum  = Real(2.5)       # Reference scale height for specific humidity (in km)
const refrh1 = Real(0.7)       # Reference relative humidity of near-surface air
const thd    = Real(2.4)       # Max damping time (in hours) for horizontal diffusion
                                      # (del^6) of temperature and vorticity
const thdd   = Real(2.4)       # Max damping time (in hours) for horizontal diffusion
                                      # (del^6) of divergence
const thds   = Real(12.0)      # Max damping time (in hours) for extra diffusion (del^2)
                                      # in the stratosphere
const tdrs   = Real(24.0*30.0) # Damping time (in hours) for drag on zonal-mean wind
                                      # in the stratosphere