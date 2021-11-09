# NorPEN 2021 IPW
This repository provides the functions for estimating causal effects of incremental propensity score interventions using inverse probability weighting (IPW) in a survival time setting.

## Prerequisites
* The input data object needs to be a `data.table` and contain the following variables:
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
* R package dependencies: `data.table`, `ggplot2`, `Epi`

## Usage
The main function `est_risk` estimates the cumulative incidence of outcome. This is the proportion of individuals who have failed at a specific interval after accounting for time-varying confounding. This has the following parameters:
* `dat`: the data set as a `data.table` object
* `mod_ps`: the propensity score model with `fitted.values` as one of the attributes (usually from `glm`)
* `dta`: the constant additive shift parameter in between 0 and 1
* `det`: `TRUE` if the intervention is deterministic
* `cuminc`: `TRUE` if the cumulative incidence is desired for all time points
* `comp`: `TRUE` if there is a competing risk outcome

For example
```r
est_risk(dat=dat, mod_ps=mod_ps, dta=0.4, det=F, cuminc=T, comp=F)
```

## What is in this repository?

Functions
* `est_stoch_int.R`: function to estimate the propensity score and to estimate the risk
* `est_boot.R`: nonparametric bootstrap for confidence intervals of estimates
* `plot_cuminc.R`: produces an adjusted cumulative incidence plot
* `plot_rrr.R`: produces a plot of the relative risk reduction varying by delta
* `example.R`: example usage of the estimation and plotting functions

Example data set
* `dat.RDS`: simulated example data set meeting the above prerequisites
* `boot_dta.RDS`: bootstrap output from 1000 resamples of the example data set
