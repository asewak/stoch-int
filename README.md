# NorPEN 2021 IPW
This repository provides the functions for estimating causal effects of incremental propensity score interventions using IPW in a survival time setting.

## Prerequisites
* The data object is a `data.table` and contains the following variables:
    1. `id`: a unique individual identifier
    2. `t0`: a discrete time variable
    3. `rf`: a binary variable to specify the subset of the population for which the intervention is assigned to
    4. `A`: a binary variable indicating the initiation of treatment
    6. `Y`: a binary failure outcome
    7. `D`: a binary competing risk
* The treatment `A` is deterministically 1 for all future time intervals once treatment is initiated
* The individual has no records after the time interval in which the outcome `Y` occurs
* The outcome `Y` is deterministically 0 for all time intervals after the competing risk `D` has occurred
* `form`: the formula for the propensity score model e.g. `as.formula("A ~ L1 + L2 + t0")`
