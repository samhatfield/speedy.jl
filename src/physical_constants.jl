# Physical constants for dynamics
const rearth = RealType(6.371e+6)    # Radius of Earth
const Ω      = RealType(7.292e-05)   # Angular frequency of Earth's rotation
const g      = RealType(9.81)        # Gravitational acceleration
const akap   = RealType(2.0/7.0)     # Ratio of gas constant to specific heat of dry air at constant
                                     # pressure
const R      = RealType(akap*1004.0) # Gas constant