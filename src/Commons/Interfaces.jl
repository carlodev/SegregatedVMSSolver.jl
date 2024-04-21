module Interfaces

export evaluate_convergence


function evaluate_convergence(ratio::Float64, field::String)
    if ratio<1
        @error "Simulation not converging, $field is diverging, error =$ratio"
    else
        println("$field error = $ratio")
    end
end


end