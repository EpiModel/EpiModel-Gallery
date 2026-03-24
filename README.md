# EpiModel Gallery

[![Test Gallery Models](https://github.com/EpiModel/EpiModel-Gallery/workflows/Test-Gallery-Models/badge.svg?branch=main)](https://github.com/EpiModel/EpiModel-Gallery/actions)

Templates for extending the [EpiModel](https://github.com/statnet/EpiModel) platform to model infectious disease dynamics over networks. Each example demonstrates how to build custom modules for stochastic network models using exponential random graph models (ERGMs).

EpiModel provides built-in SIS/SIR models out of the box, but its module API supports arbitrarily complex epidemic systems. The learning curve for writing custom modules can be steep, so these examples teach by example -- from adding a single compartment to multi-stage disease models with interventions and cost-effectiveness analysis.


## Gallery Examples

| Example | Description |
|---------|-------------|
| [SI with Vital Dynamics](examples/si-vital-dynamics/) | SI with aging, births, deaths, and age-specific mortality |
| [Adding an Exposed State](examples/seir-exposed-state/) | SEIR/SEIRS: adding an exposed (latent) compartment, with optional waning immunity |
| [Test and Treat](examples/sis-test-and-treat/) | SIS with screening and antibiotic treatment for bacterial STIs |
| [Competing Strains](examples/sis-competing-strains/) | SIS with two pathogen strains differing in infectiousness and duration |
| [Social Diffusion](examples/social-diffusion/) | SI framework repurposed for social diffusion with threshold and dose-response contagion |
| [Syphilis](examples/syphilis/) | Multi-stage syphilis with diagnosis, treatment, and recovery |
| [SEIR with AON Vaccination](examples/seir-aon-vaccination/) | SEIR with all-or-nothing vaccination, vital dynamics, and herd immunity |
| [SEIRS with Leaky Vaccination](examples/seirs-leaky-vaccination/) | SEIRS with leaky vaccination (reduced transmission probability) and vital dynamics |
| [HIV](examples/hiv/) | HIV with acute/chronic/AIDS stages and antiretroviral therapy (ART) |
| [Cost-Effectiveness](examples/cost-effectiveness/) | SI with cost-effectiveness analysis: costs, QALYs, discounting, and ICERs |
| [Observed Network Data](examples/observed-network-data/) | Epidemics over observed (census) dynamic networks without ERGM estimation |
| [Multinets](examples/multinets/) | Multilayer networks with cross-layer dependency (e.g., main vs. casual partnerships) |

Each example contains `model.R` (network estimation and simulation), `module-fx.R` (custom module functions), and `index.qmd` (annotated Quarto document).


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
source("examples/si-vital-dynamics/model.R")
```

Or from the command line:

```bash
Rscript examples/si-vital-dynamics/model.R
```

### Running All Examples

```bash
bash test.sh
```

This runs each example's `model.R` and reports pass/fail with timing.


## Contributing

Contributions of new gallery examples are welcome! To contribute:

1. Fork this repository on GitHub.
2. Create a new subdirectory under `examples/` containing:
   - `model.R` -- main script for network estimation and simulation
   - `module-fx.R` -- custom module functions plugged into `control.net()`
   - `index.qmd` -- annotated Quarto document
3. Include the standard unit test lines near the top of `model.R` (see existing examples).
4. Verify your example passes: `bash test.sh`
5. Submit a Pull Request.


## Citation

If using EpiModel for teaching or research, please include a citation:

> Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. *Journal of Statistical Software.* 2018; 84(8): 1-47. doi: 10.18637/jss.v084.i08


## License

[MIT](LICENSE)
