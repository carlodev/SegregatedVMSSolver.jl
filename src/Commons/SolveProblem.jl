module SolveProblem
using Gridap
using GridapDistributed
using GridapPETSc
using SparseArrays
using PartitionedArrays
using Parameters
using LinearAlgebra
using Gridap.FESpaces
using TimerOutputs

using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using SegregatedVMSSolver.CreateProblem
using SegregatedVMSSolver.MatrixCreation
using SegregatedVMSSolver.VectorsOperations
using SegregatedVMSSolver.ExportUtility
using SegregatedVMSSolver.Interfaces

export solve_case
export CBStepSetupStart, CBStepSetupEnd, CBSolveStart, CBSolveEnd

struct CBStepSetupStart end
struct CBStepSetupEnd end
struct CBSolveStart end
struct CBSolveEnd end


# ---------------------------------------------------------------------------
# Indices used to destructure the matrices/vectors tuples in a single place.
# The tuples are kept (rather than introducing a struct) to minimise the
# blast radius of the refactor. If/when a struct is introduced, only these
# helpers need to change.
# ---------------------------------------------------------------------------

@inline function _unpack_matrices(matrices::Tuple)
    (
        Tuu     = matrices[1],
        Tpu     = matrices[2],
        Auu     = matrices[3],
        Aup     = matrices[4],
        Apu     = matrices[5],
        App     = matrices[6],
        ML      = matrices[7],
        inv_ML  = matrices[8],
        S       = matrices[9],
        Vec_Auu = matrices[10],
        Vec_Aup = matrices[11],
        Vec_Apu = matrices[12],
        Vec_App = matrices[13],
        Vec_Au  = matrices[14],
        Vec_Ap  = matrices[15],
    )
end

@inline function _unpack_vectors(vectors::Tuple)
    (
        pm       = vectors[1],
        um       = vectors[2],
        am       = vectors[3],
        sum_pm   = vectors[4],
        Δa_star  = vectors[5],
        Δpm1     = vectors[6],
        Δa       = vectors[7],
        b1       = vectors[8],
        b2       = vectors[9],
        ũ_vector = vectors[10],
    )
end


function initialize_solve(simcase::SimulationCase, params::Dict{Symbol,Any})
    @unpack trials, tests = params
    U, P = trials

    @sunpack t0, dt, save_sim_dir = simcase

    Ut0   = U(t0)
    Pt0   = P(t0)
    Ut0_1 = U(t0 + dt)
    Pt0_1 = P(t0 + dt)

    merge!(params, Dict(:Utn => Ut0, :Ptn => Pt0, :Utn1 => Ut0_1, :Ptn1 => Pt0_1))

    uh0, ph0 = create_initial_conditions(simcase, params)
    @info "Initial Conditions Created"

    matrices = initialize_matrices(uh0, params, simcase)
    vectors  = initialize_vectors(matrices, uh0, ph0)

    initialize_export_nodes(params, simcase)
    mkpath(save_sim_dir)

    uh_avg = FEFunction(Ut0, vectors[2])
    set_zeros!(uh_avg.fields)
    ph_avg = FEFunction(Pt0, vectors[1])
    set_zeros!(ph_avg.fields)

    return matrices, vectors, (uh_avg, ph_avg)
end


"""
    solve_case(params::Dict{Symbol,Any}, simcase::SimulationCase)

It solves iteratively the velocity and pressure system using the LS-VMS
segregated scheme described in equations (27) and (30) of the paper.
"""
function solve_case(params::Dict{Symbol,Any}, simcase::SimulationCase; callback)
    @unpack trials, tests = params
    U, P = trials

    @sunpack t0, dt, tF = simcase.simulationp.timep
    @sunpack petsc_options, matrix_freq_update, M,
             a_err_threshold, θ, Number_Skip_Expansion = simcase.simulationp.solverp

    time_step = collect(t0 + dt : dt : tF)

    matrices, vectors, (uh_avg, ph_avg) = initialize_solve(simcase, params)
    @unpack Utn, Utn1, Ptn, Ptn1 = params

    @info "Starting up solver"

    GridapPETSc.with(args = split(petsc_options)) do
        any(kw -> occursin(kw, petsc_options), ["cuda", "aijcusparse"]) && run(`nvidia-smi`)

        M_ = _unpack_matrices(matrices)
        V_ = _unpack_vectors(vectors)

        to = TimerOutput()

        @timeit to "PETSc setup" begin
            @timeit to "vel" ns1 = create_PETSc_setup(M_.ML, vel_kspsetup)
            @timeit to "pres" ns2 = create_PETSc_setup(M_.S,  pres_kspsetup)
        end

        uh_tn_updt = FEFunction(Utn, V_.um)

        for (ntime, tn) in enumerate(time_step)
            @info "outer iteration $ntime, t = $tn"

            if mod(ntime, matrix_freq_update) == 0
                @info "updating matrices, vectors and PETSc setup"
                @timeit to "Step setup" begin
                    cbstate = callback(CBStepSetupStart())
                    @timeit to "matrices and vectors" update_all_matrices_vectors!(matrices, uh_tn_updt, params, simcase)
                    @timeit to "ns vel" numerical_setup!(ns1, M_.ML)
                    @timeit to "ns pres" numerical_setup!(ns2, M_.S)
                    callback(CBStepSetupEnd(), cbstate)
                end
            end

            @timeit to "Solution" begin
                cbstate = callback(CBSolveStart())
                # Reset accumulators for this time step. Use fill! to avoid
                # allocating a new PVector just to copy zeros into it.
                fill!(V_.am,     0.0)
                fill!(V_.sum_pm, 0.0)

                m            = 0
                norm_Δa0     = 10.0
                norm_Δp0     = 10.0
                err_norm_Δa0 = 1.0
                err_norm_Δp0 = 1.0

                while (m <= M) && (err_norm_Δa0 < a_err_threshold)
                    fill!(V_.Δpm1,    0.0)
                    fill!(V_.Δa_star, 0.0)

                    @timeit to "solve velocity" solve_velocity!(ns1, M_, V_, dt, θ)
                    @timeit to "solve pressure" solve_pressure!(ns2, M_, V_, dt)

                    @timeit to "change ghost" Δpm1 = GridapDistributed.change_ghost(V_.Δpm1, M_.Aup)

                    # Δa = Δa* - θ * inv(ML) * (Aup * Δpm1)
                    # Use V_.Δa as scratch storage for Aup * Δpm1 to avoid
                    # an extra allocation; then update in place.
                    @timeit to "Aup * Δpm1" mul!(V_.Δa, M_.Aup, Δpm1)
                    @timeit to "Δa" V_.Δa .= V_.Δa_star - θ * M_.inv_ML .* V_.Δa

                    # In-place updates of velocity and pressure free dofs.
                    @timeit to "um" axpy!(dt, V_.Δa, V_.um)   # um += dt * Δa
                    @timeit to "pm" V_.pm .+= Δpm1            # pm += Δpm1   (no temporary)

                    @timeit to "copy and norm" if m == 0
                        copy!(V_.sum_pm, Δpm1)
                        copy!(V_.am,     V_.Δa)
                        norm_Δa0 = norm(V_.Δa)
                        norm_Δp0 = norm(Δpm1)
                    else
                        V_.sum_pm .+= Δpm1
                        V_.am     .+= V_.Δa
                    end

                    err_norm_Δa0 = norm_Δa0 / norm(V_.Δa)
                    err_norm_Δp0 = norm_Δp0 / norm(Δpm1)

                    evaluate_convergence(err_norm_Δa0, "velocity")
                    evaluate_convergence(err_norm_Δp0, "pressure")

                    m += 1
                end
                callback(CBSolveEnd(), cbstate)
            end

            #GridapPETSc.GridapPETSc.gridap_petsc_gc()

            @timeit to "update u vector" update_ũ_vector!(V_.ũ_vector, V_.um)

            Utn  = Utn1
            Utn1 = U(tn + dt)
            Ptn  = Ptn1
            Ptn1 = P(tn + dt)
            params[:Ptn1] = Ptn1
            params[:Utn1] = Utn1

            uh_tn_updt = FEFunction(Utn1, V_.um)
            if ntime > Number_Skip_Expansion
                @timeit to "uh_tn update" uh_tn_updt = FEFunction(Utn1, update_ũ(V_.ũ_vector))
            end

            uh_tn = FEFunction(Utn, V_.um)
            ph_tn = FEFunction(Ptn, V_.pm)

            uh_avg = update_time_average(uh_tn, uh_avg, Utn, tn, ntime, time_step, simcase.simulationp.timep)
            ph_avg = update_time_average(ph_tn, ph_avg, Ptn, tn, ntime, time_step, simcase.simulationp.timep)

            @timeit to "write solution" writesolution(params, simcase, ntime, tn, (uh_tn, ph_tn), (uh_avg, ph_avg))
            @timeit to "export fields" export_fields(params, simcase, tn, uh_tn, ph_tn)
        end
        show(to) # show timer output
    end
end


# ---------------------------------------------------------------------------
# RHS assembly for the velocity sub-problem, paper equation (27):
#
#   (Tvu + θ Δt Avu) Δa* =  −Auu·uᵐ − Aup·pᵐ − ML·aᵐ
#                          + Δt·Auu·aᵐ + (1−θ)·Aup·Σ Δp_i   + Vec_Au
#
# where ML = Tvu + θ Δt Avu. Vec_Au carries the Dirichlet contribution that
# Gridap returns as the RHS of the corresponding AffineFEOperator.
#
# The previous implementation built this expression with five matrix-vector
# products composed via `*`, each of which allocates a fresh PVector and
# triggers a halo exchange. Using the 5-argument `mul!(C, A, B, α, β)`,
# which computes  C .= α·A·B + β·C  in place, we accumulate the result
# directly into the pre-allocated buffer `b1` with zero extra allocations.
#
# Auu and ML share their column partition (trial space Utn1), so a single
# `change_ghost` per operand suffices.
# ---------------------------------------------------------------------------

function _assemble_velocity_rhs!(b1, M_, um, pm, am, sum_pm, dt::Float64, θ::Float64)
    um     = GridapDistributed.change_ghost(um,     M_.Auu)   # cols: Utn1
    pm     = GridapDistributed.change_ghost(pm,     M_.Aup)   # cols: Ptn1
    am     = GridapDistributed.change_ghost(am,     M_.ML)    # cols: Utn1
    sum_pm = GridapDistributed.change_ghost(sum_pm, M_.Aup)   # cols: Ptn1

    copy!(b1, M_.Vec_Au)                              # b1  =  Vec_Au
    mul!(b1, M_.Auu, um,     -1.0,    1.0)            # b1 += -Auu·um
    mul!(b1, M_.Aup, pm,     -1.0,    1.0)            # b1 += -Aup·pm
    mul!(b1, M_.ML,  am,     -1.0,    1.0)            # b1 += -ML·am
    mul!(b1, M_.Auu, am,      dt,     1.0)            # b1 +=  dt·Auu·am
    mul!(b1, M_.Aup, sum_pm,  1 - θ,  1.0)            # b1 += (1-θ)·Aup·Σ Δp
    return b1
end

function solve_velocity!(ns1, M_::NamedTuple, V_::NamedTuple, dt::Float64, θ::Float64)
    _assemble_velocity_rhs!(V_.b1, M_, V_.um, V_.pm, V_.am, V_.sum_pm, dt, θ)
    @info "solving velocity"
    solve!(V_.Δa_star, ns1, V_.b1)
end

# Backwards-compatible signature that takes the raw tuples.
function solve_velocity!(ns1, matrices::Tuple, vectors::Tuple, dt::Float64, θ::Float64)
    solve_velocity!(ns1, _unpack_matrices(matrices), _unpack_vectors(vectors), dt, θ)
end


# ---------------------------------------------------------------------------
# RHS assembly for the pressure sub-problem, paper equation (30):
#
#   ((Tqu + Δt·Aqu)(M̃L)⁻¹ θ Avp − Aqp) Δp^{m+1}
#       = Tqu·Δa* + Aqu·(uᵐ + Δt·Δa*) + Aqp·pᵐ + Tqu·aᵐ
#
# The sign on Vec_Ap is flipped because the continuity equation uses the
# opposite sign convention from the momentum equation (the original code
# had a comment to that effect).
# ---------------------------------------------------------------------------

function _assemble_pressure_rhs!(b2, M_, um, pm, am, Δa_star, dt::Float64)
    um      = GridapDistributed.change_ghost(um,      M_.Apu)
    pm      = GridapDistributed.change_ghost(pm,      M_.App)
    am      = GridapDistributed.change_ghost(am,      M_.Tpu)
    Δa_star = GridapDistributed.change_ghost(Δa_star, M_.Tpu)

    # b2 starts as −Vec_Ap (continuity sign convention).
    copy!(b2, M_.Vec_Ap)
    rmul!(b2, -1.0)

    mul!(b2, M_.Tpu, Δa_star, 1.0, 1.0)               # b2 += Tpu·Δa*
    mul!(b2, M_.Apu, um,      1.0, 1.0)               # b2 += Apu·um
    mul!(b2, M_.Apu, Δa_star,  dt, 1.0)               # b2 += dt·Apu·Δa*
    mul!(b2, M_.App, pm,      1.0, 1.0)               # b2 += App·pm
    mul!(b2, M_.Tpu, am,      1.0, 1.0)               # b2 += Tpu·am
    return b2
end


function solve_pressure!(ns2, M_::NamedTuple, V_::NamedTuple, dt::Float64)
    _assemble_pressure_rhs!(V_.b2, M_, V_.um, V_.pm, V_.am, V_.Δa_star, dt)
    @info "solving pressure"
    solve!(V_.Δpm1, ns2, V_.b2)
end

# Backwards-compatible signature. θ is accepted for API stability but unused:
# the pressure RHS in eq. (30) does not depend on θ.
function solve_pressure!(ns2, matrices::Tuple, vectors::Tuple, dt::Float64, θ::Float64=1.0)
    solve_pressure!(ns2, _unpack_matrices(matrices), _unpack_vectors(vectors), dt)
end

end # module
