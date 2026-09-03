#!/usr/bin/env julia

include("simulation.jl")

# Replicating the plots from the paper

E1 = SimulationState(0, 0.6667, 0, 0, 0, 0, 0, 0, 0.3750, 1.8750, 0.3367)
E2 = SimulationState(0, 0.6623, 1.0549e-4, 0, 0, 9.5901e-6, 0, 0, 0.3750, 1.8750, 0.3361)
# Disturbance of + 1%
D001 = SimulationState(0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)

sim = Simulation(SimulationParams(), E1 + D001)
println("Running simulation 1 from paper for 100 steps...")
run!(sim, 100.0)
println("Plotting results...")
for variable in [:N_η, :T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ, :A_1, :A_2, :I]
    plot()
    hline!([getfield(E1, variable)], label = "E1", color = :red, lw=2, ls=:dash)
    plot_variable(sim, variable)
    save_plot(variable, prefix="paper_1/")
end

sim = Simulation(SimulationParams(), E2 + D001)
println("Running simulation 2 from paper for 100 steps...")
run!(sim, 100.0)
println("Plotting results...")
for variable in [:N_η, :T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ, :A_1, :A_2, :I]
    plot()
    hline!([getfield(E2, variable)], label = "E2", color = :red, lw=2, ls=:dash)
    plot_variable(sim, variable)
    save_plot(variable, prefix="paper_2/")
end
