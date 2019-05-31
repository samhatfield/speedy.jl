struct Params{T<:AbstractFloat}
    # Frequency at which to output diagnostics (in number of time steps)
    n_diag::Int
    # Number of time steps in a day
    n_steps_day::Int
    # Time step in seconds
    Δt::T
    # Coefficient for semi-implicit computations
    # 0 -> forward step for gravity wave terms,
    # 1 -> backward implicit, 0.5 -> centered implicit
    α::T
end
