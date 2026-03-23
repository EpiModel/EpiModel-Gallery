# EpiModel-Gallery

[![Test Gallery Models](https://github.com/EpiModel/EpiModel-Gallery/workflows/Test-Gallery-Models/badge.svg?branch=main)](https://github.com/EpiModel/EpiModel-Gallery/actions)

Templates for extending the [EpiModel](https://github.com/statnet/EpiModel) platform to model infectious disease dynamics over networks. Each example demonstrates how to build custom modules for stochastic network models using exponential random graph models (ERGMs).

EpiModel provides built-in SIS/SIR models out of the box, but its module API supports arbitrarily complex epidemic systems. The learning curve for writing custom modules can be steep, so these examples teach by example -- from adding a single compartment to multi-stage disease models with interventions and cost-effectiveness analysis.


## Gallery Examples

| Example | Description |
|---------|-------------|
| [AddingAnExposedState](2018-08-AddingAnExposedState/) | SEIR/SEIRS: adding an exposed (latent) compartment, with optional waning immunity |
| [ObservedNetworkData](2018-08-ObservedNetworkData/) | Epidemics over observed (census) dynamic networks without ERGM estimation |
| [SIwithVitalDynamics](2018-08-SIwithVitalDynamics/) | SI with aging, births, deaths, and age-specific mortality |
| [TestAndTreatIntervention](2018-08-TestAndTreatIntervention/) | SIS with screening and antibiotic treatment for bacterial STIs |
| [CompetingStrains](2018-09-CompetingStrains/) | SIS with two pathogen strains differing in infectiousness and duration |
| [SocialDiffusion](2018-09-SocialDiffusion/) | SI framework repurposed for social diffusion with threshold and dose-response contagion |
| [SEIRwithAONVax](2018-10-SEIRwithAONVax/) | SEIR with all-or-nothing vaccination, vital dynamics, and herd immunity |
| [Syphilis](2018-11-Syphilis/) | Multi-stage syphilis with diagnosis, treatment, and recovery |
| [SEIRSwithLeakyVax](2018-12-SEIRSwithLeakyVax/) | SEIRS with leaky vaccination (reduced transmission probability) and vital dynamics |
| [HIV](2019-03-HIV/) | HIV with acute/chronic/AIDS stages and antiretroviral therapy (ART) |
| [CostEffectivenessAnalysis](2021-10-CostEffectivenessAnalysis/) | SI with cost-effectiveness analysis: costs, QALYs, discounting, and ICERs |
| [Multinets](2022-12-Multinets/) | Multilayer networks with cross-layer dependency (e.g., main vs. casual partnerships) |

Each example contains three files: `README.md` (model explanation), `model.R` (network estimation and simulation), and `module-fx.R` (custom module functions).


## Getting Started

### Prerequisites

- **R** >= 4.5
- **EpiModel** >= 2.6.0

Install EpiModel from CRAN:

```r
install.packages("EpiModel")
```

### Running an Example

Clone the repository, then run any example from the project root:

```r
source("2018-08-AddingAnExposedState/model.R")
```

Or from the command line:

```bash
Rscript 2018-08-AddingAnExposedState/model.R
```

### Running All Examples

```bash
bash test.sh
```

This runs each example's `model.R` and reports pass/fail with timing.


## How to Use

You may use and extend any of these examples in any way you want. We ask that if you publish a paper using EpiModel (with or without these gallery examples), you include a citation (see below).


## Contributing

Contributions of new gallery examples are welcome! To contribute:

1. Fork this repository on GitHub.
2. Create a new subdirectory named `YYYY-MM-Description/` containing three files:
   - `README.md` -- explanation of the model, authors, and details
   - `model.R` -- main script for network estimation and simulation
   - `module-fx.R` -- custom module functions plugged into `control.net()`
3. Include the standard unit test lines near the top of `model.R` (see existing examples).
4. Verify your example passes: `bash test.sh`
5. Submit a Pull Request. We will review, request changes if needed, and merge.

### Requesting New Examples

If you'd like to ask "how would you build a network model in EpiModel that does X?", file a [GitHub Issue](https://github.com/EpiModel/EpiModel-Gallery/issues) with a detailed description. We may generalize specific requests so the example is broadly useful.


## Citation

If using EpiModel for teaching or research, please include a citation:

> Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. *Journal of Statistical Software.* 2018; 84(8): 1-47. doi: 10.18637/jss.v084.i08


## License

[MIT](LICENSE)
