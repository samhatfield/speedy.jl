using Dates

mutable struct Model
    constants::Constants
    geometry::Geometry
    current_datetime::DateTime
    end_datetime::DateTime
    spectral_trans::SpectralTrans
end

function Model(;
    # Model grid point resolution
    nlon, nlat,
    # Number of model levels
    nlev,
    # Model spectral resolution
    trunc,
    # Real number type
    real_type = Float64,
    # Radius of Earth
    Rₑ = 6.371e+6,
    # Angular frequency of Earth's rotation
    Ω = 7.292e-05,
    # Gravitational acceleration
    g = 9.81,
    # Start and end datetimes
    start_datetime = DateTime(1982,1,1),
    end_datetime = DateTime(1982,1,2)
    )

    # Default gas parameters
    akap = 2.0/7.0
    R    = akap*1000.4

    # Default dynamical constant parameters
    γ      = 6.0       # Reference temperature lapse rate (-dT/dz in deg/km)
    hscale = 7.5       # Reference scale height for pressure (in km)
    hshum  = 2.5       # Reference scale height for specific humidity (in km)
    refrh1 = 0.7       # Reference relative humidity of near-surface air
    thd    = 2.4       # Max damping time (in hours) for horizontal diffusion (del^6) of temperature
                       # and vorticity
    thdd   = 2.4       # Max damping time (in hours) for horizontal diffusion (del^6) of divergence
    thds   = 12.0      # Max damping time (in hours) for extra diffusion (del^2) in the stratosphere
    tdrs   = 24.0*30.0 # Damping time (in hours) for drag on zonal-mean wind in the stratosphere

    current_datetime = start_datetime

    constants = Constants(real_type, Rₑ, Ω, g, akap, R, γ, hscale, hshum, refrh1, thd, thdd, thds,
                          tdrs)
    geometry = Geometry(real_type, constants, nlon, nlat, nlev, trunc)
    spectral_trans = SpectralTrans(real_type, geometry, constants.Rₑ)

    Model(constants, geometry, current_datetime, end_datetime, spectral_trans)
end
