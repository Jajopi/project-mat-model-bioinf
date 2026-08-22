using Parameters

struct SimulationState
       t::Float64
     N_η::Float64
    T_1η::Float64
    T_2η::Float64
    T_rη::Float64
    T_1μ::Float64
    T_2μ::Float64
    T_rμ::Float64
     A_1::Float64
     A_2::Float64
       I::Float64
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
    p::SimulationParams
    Δ::Float64
    d_τ1::Int
    d_τ2::Int
    d_τ3::Int
    states::Array{SimulationState, 1}
    function Simulation(params::SimulationParams, Δ::Float64, state::SimulationState)
        new(params, Δ, Int(params.τ1 / Δ), Int(params.τ2 / Δ), Int(params.τ3 / Δ), [state])
    end
    function perform_step()
        s = states[end]
        s_τ1 = states[end - d_τ1]
        s_τ2 = states[end - d_τ2]
        s_τ3 = states[end - d_τ3]

        states.push!(
            s.t + Δ,
            s.N_η  + Δ * (p.α - p.γ3 * s.N_η - s.N_η * s_τ2.A_2 * s.T_1η / (1 + p.μ2 * s.T_2η) - p.ϕ * s.N_η * s_τ2.A_2 * s.T_2η - p.κ * s.N_η * s_τ2.A_2 * s.T_rη),
            s.T_1η + Δ * (-(1 + p.θ) * s.T_1η + p.θ * s_τ3.T_1μ +       p.v * s.N_η * s_τ2.A_2 / (1 + p.μr * s.T_rη) * s.T_1η / (1 + p.μ2 * s.T_2η)                       - p.η1 * s.I * s.T_1η / (1 + s.I)),
            s.T_2η + Δ * (-(1 + p.θ) * s.T_2η + p.θ * s_τ3.T_2μ + p.ϕ * p.v * s.N_η * s_τ2.A_2 / (1 + p.μr * s.T_rη) * s.T_2η / (1 + p.μ1 * s.T_1η / (1 + p.μ2 * s.T_2η)) + p.η2 * s.I * s.T_2η / (1 + s.I)),
            s.T_rη + Δ * (-(1 + p.θ) * s.T_rη + p.θ * s_τ3.T_rμ + p.κ * p.v * s.N_η * s_τ2.A_2 * s.T_rη                                                                   - p.ηr * s.I * s.T_rη / (1 + s.I)),
            s.T_1μ + Δ * (-(1 + p.θ) * s.T_1μ + p.θ * s_τ3.T_1η),
            s.T_2μ + Δ * (-(1 + p.θ) * s.T_2μ + p.θ * s_τ3.T_2η),
            s.T_rμ + Δ * (-(1 + p.θ) * s.T_rμ + p.θ * s_τ3.T_rη),
            s.A_1  + Δ * (p.λ - p.γ1 * s.A_1 - p.β * p.Λ * s.A_1),
            s.A_2  + Δ * (p.β * p.Λ * s.A_1 - p.γ2 * s.A_2 - p.μ * s.A_2 * s.T_rη),
            s.I    + Δ * (-p.γ4 * s.I + p.k1 * (s_τ1.A_2 + s_τ1.N_η + s_τ1.T_1η + s_τ1.T_2η + s_τ1.T_rη))
        )
    end
end
