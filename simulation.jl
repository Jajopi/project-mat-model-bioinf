#!/usr/bin/env julia

# Functions to run the simulation and plot the results

using Parameters
using Plots
using LaTeXStrings


# Simulation

struct SimulationState
        t::Float64 # Time
    # Non-mucosal compartment
      N_η::Float64 # Naive T cells
    T_1_η::Float64 # Th1 cells
    T_2_η::Float64 # Th2 cells
    T_r_η::Float64 # T_reg cells
    # Mucosal compartment
    T_1_μ::Float64 # Th1 cells
    T_2_μ::Float64 # Th2 cells
    T_r_μ::Float64 # T_reg cells
      A_1::Float64 # Naive APCs
      A_2::Float64 # Mature APCs
        I::Float64 # IL-6
end

function Base.:+(s1::SimulationState, disturbances::SimulationState)
    d = disturbances
    return SimulationState(s1.t, s1.N_η + d.N_η,
        s1.T_1_η + d.T_1_η, s1.T_2_η + d.T_2_η, s1.T_r_η + d.T_r_η,
        s1.T_1_μ + d.T_1_μ, s1.T_2_μ + d.T_2_μ, s1.T_r_μ + d.T_r_μ,
        s1.A_1 + d.A_1, s1.A_2 + d.A_2, s1.I + d.I)
end

@kwdef struct SimulationParams
     α::Float64 = 0.02    # The production rate of naive cells [27]
     v::Float64 = 1       # The proliferation rate of stimulated T cells [14]
     c::Float64 = 10^-4   # Cytokines produced by the immune system [14]
    μ2::Float64 = 0.1     # The strength of the suppression rate of Th1 by Th2 [14]
    μ1::Float64 = 0.2     # The strength of suppression of Th2 by Th1 [14]
    μr::Float64 = 0.25    # The strength of the suppression rate by Treg [14]
     ϕ::Float64 = 0.05    # The differences in the autocrine action of the three subsets in Equation (3) [14]
     κ::Float64 = 0.1     # The differences in the autocrine action of the three subsets in Equation (4) [14]
    γ1::Float64 = 0.3     # The death rate of naive APCs [27]
    γ2::Float64 = 0.1     # The death rate of mature APCs [27]
    γ3::Float64 = 0.03    # The death rate of naive T cells [27]
    γ4::Float64 = 0.4152  # The natural decay of IL-6 [9]
    η1::Float64 = 0.60    # The inhibition rate of Th1 cells by IL-6 [28]
    η2::Float64 = 0.007   # The stimulation rate of Th2 cells by IL-6 [28]
    ηr::Float64 = 0.39    # The inhibition rate of Treg cells by IL-6 [28]
     Λ::Float64 = 1       # The amount of the injected drug dose during desensitization [13]
    τ1::Float64 = 0.25    # Delay to production of IL-6 [28]
    τ2::Float64 = 0.0794  # The second time delay [18]
    τ3::Float64 = 0.08    # The third time delay [11]
     θ::Float64 = 0.1     # The inter compartment migration rate [14]
    k1::Float64 = 0.055   # The production of IL-6 by other cells [29]
     λ::Float64 = 0.3     # The birth rate of naive APCs [27]
     β::Float64 = 0.5     # The rate of APC activation by the antigen [13]
     μ::Float64 = 10^-2   # The rate of APC inhibition by regulatory T cells [13]
end

mutable struct Simulation
    params::SimulationParams
    Δ::Float64
    d_τ1::Int
    d_τ2::Int
    d_τ3::Int
    actual_states_stored::Int
    actual_states::Array{SimulationState, 1}
    save_every_nth_step::Int
    saved_states::Array{SimulationState}

    function Simulation(params::SimulationParams, initial_state::SimulationState; Δ::Float64 = 1e-5, save_every_nth_step::Int = 1)
        d_τ1 = Int(floor(params.τ1 / Δ))
        d_τ2 = Int(floor(params.τ2 / Δ))
        d_τ3 = Int(floor(params.τ3 / Δ))
        actual_states_stored = max(d_τ1, d_τ2, d_τ3) + 1
        actual_states = [initial_state for _ in 1:actual_states_stored]
        new(params, Δ, d_τ1, d_τ2, d_τ3, actual_states_stored, actual_states, max(1, save_every_nth_step), [])
    end
end
function perform_step!(sim::Simulation, step::Int; force_push::Bool=false)
    p = sim.params
    Δ = sim.Δ
    s = sim.actual_states[step % sim.actual_states_stored + 1]
    s_τ1 = sim.actual_states[(step + sim.actual_states_stored - sim.d_τ1) % sim.actual_states_stored + 1]
    s_τ2 = sim.actual_states[(step + sim.actual_states_stored - sim.d_τ2) % sim.actual_states_stored + 1]
    s_τ3 = sim.actual_states[(step + sim.actual_states_stored - sim.d_τ3) % sim.actual_states_stored + 1]

    if step % sim.save_every_nth_step == 0 || force_push push!(sim.saved_states, s) end

    new_state = SimulationState(
        s.t + Δ,
        s.N_η   + Δ * (p.α - p.γ3 * s.N_η - s.N_η * s_τ2.A_2 * s.T_1_η / (1 + p.μ2 * s.T_2_η) - p.ϕ * s.N_η * s_τ2.A_2 * s.T_2_η - p.κ * s.N_η * s_τ2.A_2 * s.T_r_η),
        s.T_1_η + Δ * (-(1 + p.θ) * s.T_1_η + p.θ * s_τ3.T_1_μ +       p.v * s.N_η * s_τ2.A_2 / (1 + p.μr * s.T_r_η) * s.T_1_η / (1 + p.μ2 * s.T_2_η)                       - p.η1 * s.I * s.T_1_η / (1 + s.I)),
        s.T_2_η + Δ * (-(1 + p.θ) * s.T_2_η + p.θ * s_τ3.T_2_μ + p.ϕ * p.v * s.N_η * s_τ2.A_2 / (1 + p.μr * s.T_r_η) * s.T_2_η / (1 + p.μ1 * s.T_1_η / (1 + p.μ2 * s.T_2_η)) + p.η2 * s.I * s.T_2_η / (1 + s.I)),
        s.T_r_η + Δ * (-(1 + p.θ) * s.T_r_η + p.θ * s_τ3.T_r_μ + p.κ * p.v * s.N_η * s_τ2.A_2 * s.T_r_η                                                                   - p.ηr * s.I * s.T_r_η / (1 + s.I)),
        s.T_1_μ + Δ * (-(1 + p.θ) * s.T_1_μ + p.θ * s_τ3.T_1_η),
        s.T_2_μ + Δ * (-(1 + p.θ) * s.T_2_μ + p.θ * s_τ3.T_2_η),
        s.T_r_μ + Δ * (-(1 + p.θ) * s.T_r_μ + p.θ * s_τ3.T_r_η),
        s.A_1   + Δ * (p.λ - p.γ1 * s.A_1 - p.β * p.Λ * s.A_1),
        s.A_2   + Δ * (p.β * p.Λ * s.A_1 - p.γ2 * s.A_2 - p.μ * s.A_2 * s.T_r_η),
        s.I     + Δ * (-p.γ4 * s.I + p.k1 * (s_τ1.A_2 + s_τ1.N_η + s_τ1.T_1_η + s_τ1.T_2_η + s_τ1.T_r_η))
    )
    sim.actual_states[(step + 1) % sim.actual_states_stored + 1] = new_state
end
function run!(sim::Simulation, steps::Int)
    for step in 1:steps perform_step!(sim, step) end
end
function run!(sim::Simulation, time::Float64)
    steps = Int(floor(time / sim.Δ)) + 1
    run!(sim, steps)
end
function run_until_convergence!(sim::Simulation; convergence_threshold::Float64=1e-8, divergence_threshold::Float64=1e12, max_time::Float64=1e3, max_steps::Int=0)
    time_steps = Int(floor(max_time / sim.Δ)) + 1
    max_steps = (max_steps > 0) ? min(time_steps, max_steps) : time_steps
    Δ = sim.Δ * sim.actual_states_stored
    for step in 1:sim.actual_states_stored perform_step!(sim, step) end
    for step in sim.actual_states_stored:max_steps
        perform_step!(sim, step)
        prev_state = sim.actual_states[(step - sim.actual_states_stored + 2) % sim.actual_states_stored + 1]
        curr_state = sim.actual_states[(step + 1) % sim.actual_states_stored + 1]
        if all(abs(getfield(curr_state, field) - getfield(prev_state, field)) < convergence_threshold * Δ for field in fieldnames(SimulationState)[2:end])
            println("Simulation converged at step $step., time $(curr_state.t).")
            return
        end
        if any(abs(getfield(curr_state, field)) > divergence_threshold * Δ for field in fieldnames(SimulationState)[2:end])
            println("Simulation diverged at step $step, time $(curr_state.t).")
            return
        end
    end
    println("$max_steps steps (time $max_time), reached without convergence or divergence.")
end


# Plotting

function plot_variable(sim::Simulation, variable::Symbol; steps::Int=0, max_points_roughly::Int=1000)
    if steps == 0 steps = length(sim.saved_states) end
    step_size = max(1, Int(floor(steps / max_points_roughly)))
    range = sim.saved_states[1:step_size:end]
    t = [s.t for s in range]
    y = [getfield(s, variable) for s in range]
    plot!(t, y, label = L"%$variable",
        xlabel = "Time", ylabel = L"%$variable", title = L"%$variable / t",
        dpi = 300, lw=2, color = :blue, legend=:topright, grid=true)
end
function plot_multiple_variables(sim::Simulation, variables::Vector{Symbol}; steps::Int=0, max_points_roughly::Int=1000)
    if steps == 0 steps = length(sim.saved_states) end
    step_size = max(1, Int(floor(steps / max_points_roughly)))
    range = sim.saved_states[1:step_size:end]
    t = [s.t for s in range]
    variable_names = join([string(variable) for variable in variables], ", ")
    plot!(xlabel = "Time", title = L"%$variable_names / t",
        dpi = 300, legend=:topright, grid=true)
    for variable in variables
        y = [getfield(s, variable) for s in range]
        plot!(t, y, label = L"%$variable", lw=2)
    end
end
function save_plot(variable::Symbol; prefix::String = "plot_", plot_type::String = "png")
    file_name = "plots/$(prefix)$(variable).$(plot_type)"
    mkpath(dirname(file_name))
    savefig(file_name)
end
function save_plot(name::String; prefix::String = "plot_", plot_type::String = "png")
    file_name = "plots/$(prefix)$(name).$(plot_type)"
    mkpath(dirname(file_name))
    savefig(file_name)
end


# Running basic experiments

function perform_experiment(S::SimulationState, name::String; E::Union{SimulationState, Nothing}=nothing, time::Float64=0.0)
    sim = Simulation(SimulationParams(), S, save_every_nth_step=1000)
    if time == 0 run_until_convergence!(sim, max_time=2e3) else run!(sim, time) end
    for variable in [:N_η, :T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ, :A_1, :A_2, :I]
        plot()
        if E !== nothing hline!([getfield(E, variable)], label="Equilibrium", color=:red, lw=2, ls=:dash) end
        plot_variable(sim, variable)
        save_plot(variable, prefix=name)
    end
    plot()
    plot_multiple_variables(sim, [:T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ])
    save_plot("T_cells", prefix=name)
    plot()
    plot_multiple_variables(sim, [:A_1, :A_2])
    save_plot("APCs", prefix=name)
    println("Done")
end

function perform_experiment_with_repeated_addition(S::SimulationState, name::String; E::Union{SimulationState, Nothing}=nothing, time::Float64=0.0, addition::SimulationState=SimulationState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), addition_interval::Float64=1.0)
    sim = Simulation(SimulationParams(), S, save_every_nth_step=1000)
    steps_per_addition = Int(floor(addition_interval / sim.Δ))
    total_steps = Int(floor(time / sim.Δ))
    for step in 1:total_steps
        perform_step!(sim, step, force_push=(step % steps_per_addition == 0))
        if step % steps_per_addition == 0
            current_state = sim.actual_states[(step + 1) % sim.actual_states_stored + 1]
            new_state = current_state + addition
            sim.actual_states[(step + 1) % sim.actual_states_stored + 1] = new_state
        end
    end
    for variable in [:N_η, :T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ, :A_1, :A_2, :I]
        plot()
        if E !== nothing hline!([getfield(E, variable)], label="Equilibrium", color=:red, lw=2, ls=:dash) end
        plot_variable(sim, variable, max_points_roughly=100000)
        save_plot(variable, prefix=name)
    end
    plot()
    plot_multiple_variables(sim, [:T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ], max_points_roughly=100000)
    save_plot("T_cells", prefix=name)
    plot()
    plot_multiple_variables(sim, [:A_1, :A_2], max_points_roughly=100000)
    save_plot("APCs", prefix=name)
    println("Done")
end
