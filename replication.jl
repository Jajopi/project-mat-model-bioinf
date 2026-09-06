#!/usr/bin/env julia

include("simulation.jl")

# Replicating the plots from the paper

E1 = SimulationState(0, 0.6667, 0, 0, 0, 0, 0, 0, 0.3750, 1.8750, 0.3367)
E2 = SimulationState(0, 0.6623, 1.0549e-4, 0, 0, 9.5901e-6, 0, 0, 0.3750, 1.8750, 0.3361)
# Disturbance of + 1%
D001 = SimulationState(0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)

println("Running simulation 1 from paper for 100 steps...")
perform_experiment(E1 + D001, "paper_1/", E1; time=100.0)

println("Running simulation 2 from paper for 100 steps...")
perform_experiment(E2 + D001, "paper_2/", E2; time=100.0)
