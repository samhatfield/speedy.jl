mutable struct Model
    geometry::Geometry
end

function Model(;
    # Model grid point resolution
    nlon, nlat,
    # Number of model levels
    nlev,
    # Model spectral resolution
    trunc,
    # Real number type
    real_type = Float64
    )

    geometry = Geometry(real_type, nlon, nlat, nlev, trunc)

    Model(geometry)
end