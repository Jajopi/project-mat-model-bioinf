#!/usr/bin/env julia

include("simulation.jl")

# Searching for the new equilibria

# Testing convergence for E1 and E2

E1 = SimulationState(0, 0.6667, 0, 0, 0, 0, 0, 0, 0.3750, 1.8750, 0.3367)
E2 = SimulationState(0, 0.6623, 1.0549e-4, 0, 0, 9.5901e-6, 0, 0, 0.3750, 1.8750, 0.3361)

println("Running stability test of E1 from paper...")
perform_experiment(E1, "paper_E1/", E1)

println("Running stability test of E2 from paper ...")
perform_experiment(E2, "paper_E2/", E2)


# More precise equilibria solved using agent

E1_precise = SimulationState(0, 0.6666666666666667, 0, 0, 0, 0, 0, 0, 0.37499999999999994, 1.8749999999999996, 0.3366851316634553)
E2_precise = SimulationState(0, 0.6623193028718382, 0.0001050215816082241, 0, 0, 9.547416509838555e-6, 0, 0, 0.37499999999999994, 1.8749999999999996, 0.33612316436642464)

println("Running stability test of E1_precise (higher-precision equilibrium)...")
perform_experiment(E1_precise, "paper_E1_precise/"; time=100.0)

println("Running stability test of E2_precise (higher-precision equilibrium)...")
perform_experiment(E2_precise, "paper_E2_precise/"; time=100.0)


# Disturbance of 1%, waiting for convergence

D001 = SimulationState(0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)

println("Running simulation 1 from paper (with higher-precision equilibrium) until convergence...")
perform_experiment(E1_precise + D001, "paper_1_conv/", E1_precise)

println("Running simulation 2 from paper (with higher-precision equilibrium) until convergence...")
perform_experiment(E2_precise + D001, "paper_2_conv/", E2_precise)
