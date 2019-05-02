# Prints global means of eddy kinetic energy and temperature.
# Also stops the integration if the computed diagnostics are outside of
# allowable ranges.
function check_diagnostics(vor, div, tem, step)
    diagnostics = zeros(RealType, nlev, 3)

    # 1. Get global-mean temperature and compute eddy kinetic energy
    for k in 1:nlev
        diagnostics[k,3] = √(half)*real(tem[1,1,k])

        temp = ∇⁻²(vor[:,:,k])

        for m in 2:mx
            for n in 1:nx
                diagnostics[k,1] = diagnostics[k,1] - real(temp[m,n]*conj(vor[m,n,k]))
            end
        end

        temp = ∇⁻²(div[:,:,k])

        for m in 2:mx
            for n in 1:nx
                diagnostics[k,2] = diagnostics[k,2] - real(temp[m,n]*conj(div[m,n,k]))
            end
        end
    end

    # 2. Print results to screen
    if mod(step, nstdia) == 0
        println("Step = $step")
        println(diagnostics[:,1])
        println(diagnostics[:,2])
        println(diagnostics[:,3])
    end

    # 3. Stop integration if model variables are out of range
    for k in 1:nlev
        if diagnostics[k,1] > 500.0 || diagnostics[k,2] > 500.0 ||
            diagnostics[k,3] < 180.0 || diagnostics[k,3] > 320.0

            println(diagnostics[:,1])
            println(diagnostics[:,2])
            println(diagnostics[:,3])

            throw("Model variables out of accepted range")
        end
    end
end