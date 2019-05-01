using NetCDF

# Load boundary condition from given file
function load_boundary_file(file_name, field_name)
    # Read variable from boundary file
    raw_input = ncread("data/boundaries/t30/clim/$file_name", field_name)

    # Flip latitudes
    raw_input = raw_input[:,end:-1:1]

    # Fix undefined values
    raw_input[raw_input .<= -999.0] .= 0.0
    convert(Array{Real}, raw_input)
end