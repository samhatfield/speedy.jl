# Model geometry parameters
const trunc   = 30        # Spectral truncation
const nlon    = 96        # Number of longitudes
const nlat    = 48        # Number of latitudes
const nlev    = 8         # Number of vertical levels
const nx      = trunc + 2 # Number of total wavenumbers
const mx      = trunc + 1 # Number of zonal wavenumbers
const n_trace = 1         # Number of tracers

# Time stepping parameters
const n_steps_day     = 36                            # Number of time steps in one day
const Δt              = Real(86400.0/n_steps_day) # Time step in seconds
const rob             = Real(0.05)                # Damping factor in Robert time filter
const wil             = Real(0.53)                # Parameter of Williams filter
const α               = Real(0.5)                 # Coefficient for semi-implicit computations
                                                      # 0 -> forward step for gravity wave terms,
                                                      # 1 -> backward implicit
                                                      # 0.5 -> centered implicit

# Output parameters
const nstdia = 36*5  # Period (number of steps) for diagnostic print-out
const nsteps_out = 1 # Number of time steps between outputs