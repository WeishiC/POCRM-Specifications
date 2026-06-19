## Design specifications for the POCRM on two-dimensional grids
This repository contains code used for the paper _Design Specification of Partial Ordering Continual Reassessment Method Based on Consistency Conditions_ by Weishi Chen and Pavel Mozgunov.

## Functions
- `ListOrder.R` contains code to list all possible orderings on any given two-dimensional grids.
- `POCRM.R` contains code to simulate the operating characteristics of the two-stage Partial Ordering Continual Reassessment Method (Wages et al, 2011)
- `BLRM.R` contains code to simulate the operating characteristics of the Bayesian Logistic Regression Model (Riviere et al, 2014)
- `SFD.R` contains code to simulate the operating characteristics of the Surface-Free Design (Mozgunov et al, 2020)

The folder `RData` contains the scenarios used for the simulation studies.

## Dependencies
This project is written in R. Key R packages used (as seen in the code) include:
- dfcrm
- pocrm
- BOIN
- rjags
- nnet
- combinat
- tidyverse

## References
- Wages, N.A., Conaway, M.R. and O'Quigley, J. *Dose-finding design for multi-drug combinations*. Clinical Trials 2011; 8: 380-389.
- Mozgunov P, Gasparini M, and Jaki T. *A surface-free design for phase I dual-agent combination trials*. SMMR 2020; 29: 3093-109.
- Riviere M.K., Yuan Y., Dubois F., and Zohar S. *A Bayesian dose-finding design for drug combination clinical trials based on the logistic model*. PharmStats 2014; 13: 247-57.
- Lin R and Yin G. *Bayesian optimal interval design for dose finding in drug-combination trials*. SMMR 2017; 26(5): 2155-67.

##  License
MIT License. See `LICENSE` file for details.
