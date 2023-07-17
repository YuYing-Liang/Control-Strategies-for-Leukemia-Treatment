## Description of algorithm and how to run code
all files with extension `.mlx` are main files that hold the main code of the Adaptive Model Predictive Control (AMPC) algorithm and generate important figures.

### Main files

`main.mlx` has the AMPC code for **fixed** on-drug and off-drug periods. Simulation parameters can be changed in the "Initialization" section (cycle length, initial states, number of cycles etc.). Run all sections prior to the "Plots" section to generate figures for neutrophil concentration over time, drug dose over time and parameter estimates over each cycle.

`comparison_individual.mlx` has AMPC code for **variable** on-drug and off-drug periods. This file also compares clinical data from Jost et al. to results using AMPC with the same cycle lengths as the clinical treatment. Run all sections, the final section generates the same figures as `main.mlx` with an addition of the clinical data. In this file you can only generate plots for one patient at a time.

`comparison.mlx` has AMPC code that generates comparison plots for all patients from Jost et al.'s clinical data. It also performs a high-level error analysis on the output of AMPC and clinical treatment strategies to the desired neturophil concentration.

### Helper Functions
`calc_params.m` contains the "adaptive" part of the AMPC algorithm. Here, we use all previous measurements to estimate model parameters.

`calc_control.m` contains the "model predictive" part of the AMPC algorithm. Using the estimated model parameters, we look ahead a certain number of cycles to determine the best amount of drug to give the patient for the next cycle.

`calc_control_comparison.m` this is the same as `calc_control.m` but with the ability to handle dynamic/changing cycle lengths.

`get_cycle_day.m` retrieves the day in the cycle given the cycle length and the day in treatment.

`g_nonlin.m` is the observation function g(x) in y=g(x). Outputs the measurement y given x.

`jost_discrete.m` is the system function f(x) in x'=f(x) which uses Jost et al.'s neutrophil model.

`jost_noisy.m` is the same as above but includes input noise.
