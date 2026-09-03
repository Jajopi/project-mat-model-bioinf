#!/usr/bin/env julia

include("simulation.jl")

# Searching for the new equilibria

# Testing convergence for E1 and E2

E1 = SimulationState(0, 0.6667, 0, 0, 0, 0, 0, 0, 0.3750, 1.8750, 0.3367)
E2 = SimulationState(0, 0.6623, 1.0549e-4, 0, 0, 9.5901e-6, 0, 0, 0.3750, 1.8750, 0.3361)

sim = Simulation(SimulationParams(), E1, save_frequency=1000)
println("Running stability test of E1 from paper...")
run_until_convergence!(sim, max_time=1e3)
println("Plotting results...")
for variable in [:N_η, :T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ, :A_1, :A_2, :I]
    plot()
    hline!([getfield(E1, variable)], label = "E1", color = :red, lw=2, ls=:dash)
    plot_variable(sim, variable)
    save_plot(variable, prefix="paper_E1/")
end

sim = Simulation(SimulationParams(), E2, save_frequency=1000)
println("Running stability test of E2 from paper ...")
run_until_convergence!(sim, max_time=1e3)
println("Plotting results...")
for variable in [:N_η, :T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ, :A_1, :A_2, :I]
    plot()
    hline!([getfield(E2, variable)], label = "E2", color = :red, lw=2, ls=:dash)
    plot_variable(sim, variable)
    save_plot(variable, prefix="paper_E2/")
end


# More precise equilibria solved using agent

E1_precise = SimulationState(0, 0.6666666666666667, 0, 0, 0, 0, 0, 0, 0.37499999999999994, 1.8749999999999996, 0.3366851316634553)
E2_precise = SimulationState(0, 0.6623193028718382, 0.0001050215816082241, 0, 0, 9.547416509838555e-6, 0, 0, 0.37499999999999994, 1.8749999999999996, 0.33612316436642464)

sim = Simulation(SimulationParams(), E1_precise, save_frequency=1000)
println("Running stability test of E1_precise (higher-precision equilibrium)...")
run!(sim, 100.0)
println("Plotting results...")
for variable in [:N_η, :T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ, :A_1, :A_2, :I]
    plot()
    hline!([getfield(E1_precise, variable)], label = "E1_precise", color = :red, lw=2, ls=:dash)
    plot_variable(sim, variable)
    save_plot(variable, prefix="paper_E1_precise/")
end

sim = Simulation(SimulationParams(), E2_precise, save_frequency=1000)
println("Running stability test of E2_precise (higher-precision equilibrium)...")
run!(sim, 100.0)
println("Plotting results...")
for variable in [:N_η, :T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ, :A_1, :A_2, :I]
    plot()
    hline!([getfield(E2_precise, variable)], label = "E2_precise", color = :red, lw=2, ls=:dash)
    plot_variable(sim, variable)
    save_plot(variable, prefix="paper_E2_precise/")
end


# Disturbance of 1%, waiting for convergence

D001 = SimulationState(0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)

sim = Simulation(SimulationParams(), E1_precise + D001; save_frequency=1000)
println("Running simulation 1 from paper until convergence...")
run_until_convergence!(sim, max_time=1e3)
println("Plotting results...")
for variable in [:N_η, :T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ, :A_1, :A_2, :I]
    plot()
    hline!([getfield(E1, variable)], label = "E1", color = :red, lw=2, ls=:dash)
    plot_variable(sim, variable)
    save_plot(variable, prefix="paper_1_conv/")
end

sim = Simulation(SimulationParams(), E2_precise + D001, save_frequency=1000)
println("Running simulation 2 from paper until convergence...")
run_until_convergence!(sim, max_time=1e3)
println("Plotting results...")
for variable in [:N_η, :T_1_η, :T_2_η, :T_r_η, :T_1_μ, :T_2_μ, :T_r_μ, :A_1, :A_2, :I]
    plot()
    hline!([getfield(E2, variable)], label = "E2", color = :red, lw=2, ls=:dash)
    plot_variable(sim, variable)
    save_plot(variable, prefix="paper_2_conv/")
end
