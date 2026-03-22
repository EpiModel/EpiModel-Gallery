# CLAUDE.md

## Project Overview

EpiModel-Gallery is a collection of template examples for extending the EpiModel R package to model infectious disease dynamics over networks. Examples focus on stochastic network models using ERGMs (exponential random graph models). EpiModel itself is documented at https://github.com/statnet/EpiModel.

## Repository Structure

Each subdirectory is a self-contained example named `YYYY-MM-Description/` containing:
- `README.md` — explanation of the model, authors, and details
- `model.R` — main script for network estimation and epidemic simulation
- `module-fx.R` — custom module functions plugged into `control.net()`

Top-level files:
- `DESCRIPTION` — R package metadata (used for CI dependency management)
- `test.sh` — test runner that executes each example's `model.R`
- `EpiModel-Gallery.Rproj` — RStudio project config

### Gallery Examples

| Directory | Model |
|-----------|-------|
| `2018-08-AddingAnExposedState` | SEIR/SEIRS — adding exposed (latent) state to SIR |
| `2018-08-ObservedNetworkData` | Epidemics over observed (census) dynamic networks |
| `2018-08-SIwithVitalDynamics` | SI with aging, births, and deaths |
| `2018-08-TestAndTreatIntervention` | SIS with testing and treatment for bacterial STIs |
| `2018-09-CompetingStrains` | SIS with two competing pathogen strains |
| `2018-09-SocialDiffusion` | Information/behavior diffusion using SI framework |
| `2018-10-SEIRwithAONVax` | SEIR with all-or-nothing vaccination and vital dynamics |
| `2018-11-Syphilis` | Multi-stage syphilis model with diagnosis/treatment |
| `2018-12-SEIRSwithLeakyVax` | SEIRS with leaky vaccination and vital dynamics |
| `2019-03-HIV` | HIV with acute/chronic/AIDS stages and ART treatment |
| `2021-10-CostEffectivenessAnalysis` | SI with cost-effectiveness (costs, QALYs, discounting) |
| `2022-12-Multinets` | Multiple interacting networks with shared node set |

## Language and Dependencies

- **Language:** R (>= 4.0)
- **Primary dependency:** EpiModel (>= 2.4.0)
- **Other imports:** dplyr, networkDynamicData
- **Suggested:** testthat, dampack
- **Code style:** 2-space indentation, UTF-8 encoding

## Testing

Run all examples:
```bash
bash test.sh
```

This iterates through each subdirectory (excluding `renv/`), runs `Rscript <dir>/model.R` with an error handler, and reports pass/fail with timing.

### Unit Test Pattern

Every `model.R` must include these lines near the top (after loading EpiModel):
```r
# Standard Gallery unit test lines
rm(list = ls())
eval(parse(text = print(commandArgs(TRUE)[1])))
```

And use conditional parameters for interactive vs. CI runs:
```r
if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 800
} else {
  nsims <- 2
  ncores <- 2
  nsteps <- 200
}
```

## Contributing a New Example

1. Create a new subdirectory named `YYYY-MM-Description/`
2. Add three files: `README.md`, `model.R`, `module-fx.R`
3. Include the standard unit test lines in `model.R`
4. Ensure `bash test.sh` passes before submitting
