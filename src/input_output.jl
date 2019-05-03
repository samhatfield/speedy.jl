# Load boundary condition from given file
function load_boundary_file(file_name, field_name)
    # Read variable from boundary file
    raw_input = ncread("data/boundaries/t30/clim/$file_name", field_name)

    # Flip latitudes
    raw_input = raw_input[:,end:-1:1]

    # Fix undefined values
    raw_input[raw_input .<= -999.0] .= 0.0
    convert(Array{RealType}, raw_input)
end

function output(timestep)
    # Construct file name
    file_name = Dates.format(current_datetime, "yyyymmddHHMM.nc")

    # Construct time units string
    time_string = "hours since $(Dates.format(current_datetime, "yyyy-mm-dd HH:MM:0.0"))"

    # Define Time
    timedim = NcDim("time", 0, unlimited=true)
    timevar = NcVar("time", timedim, t=Int32, atts=Dict("units"=>time_string))

    # Define Space
    londim  = NcDim("lon", nlon)
    latdim  = NcDim("lat", nlat)
    levdim  = NcDim("lev", nlev)
    lonvar  = NcVar("lon",  londim,  t=Float32,
        atts=Dict("long_name"=>"longitude"))
    latvar  = NcVar("lat",  latdim,  t=Float32,
        atts=Dict("long_name"=>"latitude"))
    levvar  = NcVar("lev",  levdim,  t=Float32,
        atts=Dict("long_name"=>"atmosphere_sigma_coordinate"))

    # Define prognostics fields
    uvar   = NcVar("u", [londim, latdim, levdim, timedim], t=Float32,
        atts=Dict("long_name"=>"eastward_wind", "units"=>"m/s"))
    vvar   = NcVar("v", [londim, latdim, levdim, timedim], t=Float32,
        atts=Dict("long_name"=>"northward_wind", "units"=>"m/s"))
    tvar   = NcVar("t", [londim, latdim, levdim, timedim], t=Float32,
        atts=Dict("long_name"=>"air_temperature", "units"=>"K"))
    qvar   = NcVar("q", [londim, latdim, levdim, timedim], t=Float32,
        atts=Dict("long_name"=>"specific_humidity", "units"=>"1"))
    ϕvar = NcVar("phi", [londim, latdim, levdim, timedim], t=Float32,
        atts=Dict("long_name"=>"geopotential_height", "units"=>"m"))
    pₛvar = NcVar("ps", [londim, latdim, timedim], t=Float32,
        atts=Dict("long_name"=>"surface_air_pressure", "units"=>"Pa"))

    # Create output file
    nc = NetCDF.create("rundir/$file_name",
        [timevar, lonvar, latvar, levvar, uvar, vvar, tvar, qvar, ϕvar, pₛvar])

    # Write dimensions to file
    NetCDF.putvar(nc, "time", [(timestep - 1)*24.0/Float32(n_steps_day)])
    NetCDF.putvar(nc, "lon", [3.75*k for k in 0:95])
    NetCDF.putvar(nc, "lat", [radang[k]*90.0/asin(1.0) for k in 1:nlat])
    NetCDF.putvar(nc, "lev", fsg)

    # Convert prognostic fields from spectral space to grid point space
    u_grid = zeros(RealType, nlon, nlat, nlev)
    v_grid = zeros(RealType, nlon, nlat, nlev)
    t_grid = zeros(RealType, nlon, nlat, nlev)
    q_grid = zeros(RealType, nlon, nlat, nlev)
    ϕ_grid = zeros(RealType, nlon, nlat, nlev)
    pₛ_grid = zeros(RealType, nlon, nlat)
    for k in 1:nlev
        ucos = zeros(RealType, mx, nx)
        vcos = zeros(RealType, mx, nx)
        uvspec!(vorU[:,:,k,1], divU[:,:,k,1], ucos, vcos)
        u_grid[:,:,k] = spec_to_grid(ucos, scale=true)
        v_grid[:,:,k] = spec_to_grid(vcos, scale=true)
        t_grid[:,:,k] = spec_to_grid(tem[:,:,k,1], scale=false)
        q_grid[:,:,k] = spec_to_grid(tr[:,:,k,1,1], scale=false)
        ϕ_grid[:,:,k] = spec_to_grid(ϕ[:,:,k], scale=false)
    end
    pₛ_grid = spec_to_grid(pₛ[:,:,1], scale=false)

    println("Writing output file for $(Dates.format(current_datetime, "YYYY-mm-dd HH:MM:00"))")

    # Write prognostic variables to file
    NetCDF.putvar(nc, "u",   Float32.(u_grid),  start=[1,1,1,1], count=[-1,-1,-1,1])
    NetCDF.putvar(nc, "v",   Float32.(v_grid),  start=[1,1,1,1], count=[-1,-1,-1,1])
    NetCDF.putvar(nc, "t",   Float32.(t_grid),  start=[1,1,1,1], count=[-1,-1,-1,1])
    NetCDF.putvar(nc, "q",   Float32.(q_grid),  start=[1,1,1,1], count=[-1,-1,-1,1])
    NetCDF.putvar(nc, "phi", Float32.(ϕ_grid),  start=[1,1,1,1], count=[-1,-1,-1,1])
    NetCDF.putvar(nc, "ps",  Float32.(pₛ_grid), start=[1,1,1],   count=[-1,-1,1])

    NetCDF.close(nc)
end
