---
title: "Replication and numerical analysis of allergic reaction model"
author: "Ján Plachý"
affiliation: "The Faculty of Science, Charles University"
project: "Project for Mathematical modeling in bioinformatics"
abstract: |
  Complexity of immune system requires rather complex models to be developed to make it possible to model even simple and well-known interactions, such as allergic response.
  In this project, we replicate a recent multiple-compartment model [1] using delayed differential equations to simulate interactions of four types of T cells, two types of APCs and one type of interleukins.
  We develop a practical Julia implementation for numerical testing of behaviour of the model, and test several cases of initial conditions, including equilibria proposed in the original paper, and corrected equilibria with more precise values.
  We search for combinations of inital conditions that would lead to interesting behaviour, and discuss the practical usability of the model in the context of results we obtain.
---

# Introduction

Allergic reactions are a result of immune system overreacting in a contact with certain substances, producing high amounts of immunoglobulin E antibodies (IgE).
Those antibodies then cause common allergic symptoms, which are usually annoying and harmless, but in some cases can lead to anaphylactic shock and death.
One of those cases is a re-exposure of a patient to various chemotherapy drugs, which are widely used today for cancer threatment.

The authors of paper Stability Analysis in a Mathematical Model for Allergic Reactions [1] use this as a primary motivation for proposing a new complex model of immune system interactions.
They model the concentrations of four types of T cells, namely naive T cells, Th1 and Th2 helper cells, and regulatory T cells, naive and mature antigen-presenting cells (APCs) and interleukins IL-6.
Th1 and Th2 cells compete via reciprocal inhibition, while regulatory T cells supress both Th1 and Th2 without being inhibited themselves.
Imbalance towards Th2 may lead to allergies, as the production of Th2-type cytokines corellates with IgE, which are associated with allergic reactions.
Naive T cells differentiate into one of the three other types via interaction with APCs. 
Interleukin IL-6, which promotes the growth of Th2 cells over Th1 cells, is produced by both helper T cells and APCs.
APCs start as naive cells and mature over time to be able to interact with T cells.

According to the previous works, the immune system needs to be modeled considering two separate compartments: the inner one, such as blood and well-perfused organs, and the peripheral one.
However, this is considered not enough, as peripheral tissues need to be further divided into mucosal and non-mucosal.
The authors mention that the model considers three compartments, but in the equations, only concentrations of T cells in mucosal and non-mucosal peripheral tissues is considered.
The third compartment is therefore only considered as the source of allergens, which slowly difund into the non-mucosal compartment.
The non-mucosal compartment $\eta$ is where the T cells are initially recruited, and they can transit to mucosal compartment $\mu$.

The transition of T cells between compartments requires a delay, as well as IL-6 production and allergen diffusion from the central compartment.
The model therefore utilizes delayed differential equations, which make theoretical analysis harder; however, numerically simulating the model stays pretty simple.
In total, the model consists of 10 equations: for naive T cells in $\eta$ compartment, for each of three types of differentiated T cells in both compartments, for naive and mature APCs and for IL-6 interleukins.

The authors state that the this model can be used to model drug desentization, and provide values for all 24 parameters used in the model, backed up by sources from which the corresponding values were taken.

# Implementation

We implemented the simulation in a standalone Julia file `simulation.jl` using a struct of high precision floats to hold the values of concentrations of tracked elements for every time step. The simulation begins with a given state and then computes a new one using a small time difference, according to the equations and parameters from the paper. The simulation then runs for the number of steps given by the time which should be reached, or until convergence (when the change in all variables within a given time frame is less than $10^{-8}$ per time unit) or divergence (when any of the variables reaches threshold $10^{12}$).

To be able to model delayed differential equations, we need to keep track of past states. However, working with small time steps neccesary to keep reasonable precision would lead to a high number of states being saved even for short simulations, leading to extremely high RAM usage even for short experiments. We circumvent this limitation by only storing a specific number of states which would still be needed for future computations, in a rotating queue. Once a new state is created, it overwrites a formerly stored state, which was just enough steps apart from the new state as the longest delay would require. For the given model, the longest delay is $0.25$, requiring $2500$ steps for simulation we did.

To be able to see the progress of the simulation, some of the states are periodically sampled and saved with configurable frequency. Also, if in a given step of the simulation, some values are artificially changed, to, for example, increase a concentration of specific cell or chemical type, the simulation can be forced to store the given state. As each state stores its time from beginning of the simulation, this does not introduce any errors in the plots.

Results of simulations are plotted separately for every simulated variable, with equilibria shown in the plot. Also, one plot with all T cell concentrations and one plot with all APCs is generated. The plots can be limited in the number of points they show, usually $1000$ points is enough to get precisely-looking results without extensive overhead; however, for experiments when we artificially increase concentrations, we need to plot all captured data points to keep the right proportions of the peaks.

# Experiments and results

The experiments we implemented are of two main directions: replicating the original results and searching for interesting behaviour that can be studied using the model.

## Replication of original results

After implementing the model, we checked the correctness of the implementation by running the simulation with conditions from the paper and plotting the results (in the file `replication.jl`). The results can be found in the supplementary plots, labeled `paper_1` and `paper_2`, for the first and the second experiment from the paper. In both experiment, the authors tried small disturbances ($0.01$) from equilibria they computed. The results were similar, with no visible differences.

To be sure the implementation runs correctly, we also checked the behaviour of the model when ran from equilibria without disturbances. Interestingly, the plots were showing convergence to slightly different equilibria (supplementary plots `paper_E1` and `paper_E2`, for the respective equilibria).

We had our code checked with an AI agent against the equations from the paper. It found out the equilibria in the paper were not precise enough, and recalculated the right values from the equations provided by the authors. We ran the simulation with the new equlibria to see that they do not change (supplementary plots `paper_E1_precise` and `paper_E2_precise`, for the respective equilibria) and reran the original experiments with disturbances of $0.01$ until they converged to check the new equilibria are really reached (supplementary plots `paper_1_conv` and `paper_2_conv`, for the respective equilibria). Both experiments were showing expected results.

## Searching for the interesting behaviour of the model

After we had checked that our implementation of the model is correct, we tried to detect whether is the model capable of simulating more complex behaviour than just cenverging to given equilibria. Specifically, the authors of the paper mentioned two interesting phenomena as the motivation for model development: drug desentization and anaphylactic shock.

First, we tried to find initial conditions that would lead to the anaphylactic shock. As the first equilibrium was, in some sense, trivial and unrealistic (with zero concentration of all kinds of T cells), we focused on the second equilibrium.

We tried to increase the number of APCs and IL-6 ten times, leading to results shown in `lot_APCs` and `lot_Is`, respectivelly.

# References {-}

1. Abdullah, Rawan, Irina Badralexi, and Andrei Halanay. "Stability Analysis in a Mathematical Model for Allergic Reactions." Axioms 13.2 (2024): 102.
