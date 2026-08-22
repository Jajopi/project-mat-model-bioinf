#!/usr/bin/env julia

# import Pkg;
# Pkg.add("Parameters")
# Pkg.add("Plots")
# Pkg.add("LaTeXStrings")

using Parameters
using Plots
using LaTeXStrings

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

struct Simulation
    params::SimulationParams
    Δ::Float64
    d_τ1::Int
    d_τ2::Int
    d_τ3::Int
    max_delay::Int
    states::Array{SimulationState, 1}
    function Simulation(params::SimulationParams, Δ::Float64, state::SimulationState)
        d_τ1 = Int(floor(params.τ1 / Δ))
        d_τ2 = Int(floor(params.τ2 / Δ))
        d_τ3 = Int(floor(params.τ3 / Δ))
        max_delay = max(d_τ1, d_τ2, d_τ3)
        states = [state for _ in 1:(max_delay + 1)]
        new(params, Δ, d_τ1, d_τ2, d_τ3, max_delay, states)
    end
end
function perform_step!(sim::Simulation)
    p = sim.params
    Δ = sim.Δ
    s = sim.states[end]
    s_τ1 = sim.states[end - sim.d_τ1]
    s_τ2 = sim.states[end - sim.d_τ2]
    s_τ3 = sim.states[end - sim.d_τ3]

    new_state = SimulationState(
        s.t + Δ,
        s.N_η  + Δ * (p.α - p.γ3 * s.N_η - s.N_η * s_τ2.A_2 * s.T_1_η / (1 + p.μ2 * s.T_2_η) - p.ϕ * s.N_η * s_τ2.A_2 * s.T_2_η - p.κ * s.N_η * s_τ2.A_2 * s.T_r_η),
        s.T_1_η + Δ * (-(1 + p.θ) * s.T_1_η + p.θ * s_τ3.T_1_μ +       p.v * s.N_η * s_τ2.A_2 / (1 + p.μr * s.T_r_η) * s.T_1_η / (1 + p.μ2 * s.T_2_η)                       - p.η1 * s.I * s.T_1_η / (1 + s.I)),
        s.T_2_η + Δ * (-(1 + p.θ) * s.T_2_η + p.θ * s_τ3.T_2_μ + p.ϕ * p.v * s.N_η * s_τ2.A_2 / (1 + p.μr * s.T_r_η) * s.T_2_η / (1 + p.μ1 * s.T_1_η / (1 + p.μ2 * s.T_2_η)) + p.η2 * s.I * s.T_2_η / (1 + s.I)),
        s.T_r_η + Δ * (-(1 + p.θ) * s.T_r_η + p.θ * s_τ3.T_r_μ + p.κ * p.v * s.N_η * s_τ2.A_2 * s.T_r_η                                                                   - p.ηr * s.I * s.T_r_η / (1 + s.I)),
        s.T_1_μ + Δ * (-(1 + p.θ) * s.T_1_μ + p.θ * s_τ3.T_1_η),
        s.T_2_μ + Δ * (-(1 + p.θ) * s.T_2_μ + p.θ * s_τ3.T_2_η),
        s.T_r_μ + Δ * (-(1 + p.θ) * s.T_r_μ + p.θ * s_τ3.T_r_η),
        s.A_1  + Δ * (p.λ - p.γ1 * s.A_1 - p.β * p.Λ * s.A_1),
        s.A_2  + Δ * (p.β * p.Λ * s.A_1 - p.γ2 * s.A_2 - p.μ * s.A_2 * s.T_r_η),
        s.I    + Δ * (-p.γ4 * s.I + p.k1 * (s_τ1.A_2 + s_τ1.N_η + s_τ1.T_1_η + s_τ1.T_2_η + s_τ1.T_r_η))
    )
    push!(sim.states, new_state)
end
function run!(sim::Simulation, time::Float64)
    steps = Int(time / sim.Δ)
    for _ in 1:steps perform_step!(sim) end
end
