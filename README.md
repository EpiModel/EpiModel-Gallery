# EpiModel Gallery

[![Test Gallery Models](https://github.com/EpiModel/EpiModel-Gallery/workflows/Test-Gallery-Models/badge.svg?branch=main)](https://github.com/EpiModel/EpiModel-Gallery/actions)

### 📖 **[Browse the Gallery Website &rarr;](https://epimodel.github.io/EpiModel-Gallery/)**

A curated collection of tutorials for extending the [EpiModel](https://github.com/statnet/EpiModel) platform to model infectious disease dynamics over networks. Each example demonstrates how to build custom modules for stochastic network models using exponential random graph models (ERGMs), with fully annotated walkthroughs on the [Gallery website](https://epimodel.github.io/EpiModel-Gallery/).

EpiModel provides built-in SIS/SIR models out of the box, but its module API supports arbitrarily complex epidemic systems. The learning curve for writing custom modules can be steep, so these examples teach by example -- from adding a single compartment to multi-stage disease models with interventions and cost-effectiveness analysis. The examples are organized from beginner-friendly fundamentals through intermediate extensions to full applied disease models.


## Gallery Examples

**Fundamentals**

| Example | Description |
|---------|-------------|
| [SEIR/SEIRS: Adding an Exposed State](https://epimodel.github.io/EpiModel-Gallery/examples/seir-exposed-state/) | Adding an exposed (latent) compartment to SIR, with optional waning immunity |
| [SI with Vital Dynamics](https://epimodel.github.io/EpiModel-Gallery/examples/si-vital-dynamics/) | SI with aging, births, deaths, and age-specific mortality |

**Interventions**

| Example | Description |
|---------|-------------|
| [SEIR with AON Vaccination](https://epimodel.github.io/EpiModel-Gallery/examples/seir-aon-vaccination/) | SEIR with all-or-nothing vaccination, vital dynamics, and herd immunity |
| [SEIRS with Leaky Vaccination](https://epimodel.github.io/EpiModel-Gallery/examples/seirs-leaky-vaccination/) | SEIRS with leaky vaccination (reduced transmission probability) and vital dynamics |
| [Test and Treat](https://epimodel.github.io/EpiModel-Gallery/examples/sis-test-and-treat/) | SIS with screening and antibiotic treatment for bacterial STIs |

**Intermediate Extensions**

| Example | Description |
|---------|-------------|
| [Competing Strains](https://epimodel.github.io/EpiModel-Gallery/examples/sis-competing-strains/) | SIS with two pathogen strains differing in infectiousness and duration |
| [Social Diffusion](https://epimodel.github.io/EpiModel-Gallery/examples/social-diffusion/) | SI framework repurposed for social diffusion with threshold and dose-response contagion |

**Network Features**

| Example | Description |
|---------|-------------|
| [Observed Network Data](https://epimodel.github.io/EpiModel-Gallery/examples/observed-network-data/) | Epidemics over observed (census) dynamic networks without ERGM estimation |
| [Multilayer Networks](https://epimodel.github.io/EpiModel-Gallery/examples/multinets/) | Multilayer networks with cross-layer dependency (e.g., main vs. casual partnerships) |

**Full Disease Models**

| Example | Description |
|---------|-------------|
| [HIV](https://epimodel.github.io/EpiModel-Gallery/examples/hiv/) | HIV with acute/chronic/AIDS stages and antiretroviral therapy (ART) |
| [Syphilis](https://epimodel.github.io/EpiModel-Gallery/examples/syphilis/) | Multi-stage syphilis with diagnosis, treatment, and recovery |

**Advanced Extensions**

| Example | Description |
|---------|-------------|
| [Cost-Effectiveness Analysis](https://epimodel.github.io/EpiModel-Gallery/examples/cost-effectiveness/) | SI with cost-effectiveness analysis: costs, QALYs, discounting, and ICERs |

Each example contains `model.R` (network estimation and simulation), `module-fx.R` (custom module functions), and `index.qmd` (annotated tutorial on the website).


## Getting Started

The recommended way to learn is through the [Gallery website](https://epimodel.github.io/EpiModel-Gallery/), which presents each example as an annotated tutorial with model diagrams, code walkthroughs, and output plots. Start with the Fundamentals section and work through the examples in order.

To run the code yourself:

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
2. Create a new subdirectory under `examples/` named descriptively (e.g., `sir-vaccination/`) containing:
   - `model.R` -- main script for network estimation and simulation
   - `module-fx.R` -- custom module functions plugged into `control.net()`
   - `index.qmd` -- annotated Quarto tutorial (see existing examples for the format)
   - `thumbnail.png` -- thumbnail image for the gallery listing
3. Include the standard unit test lines near the top of `model.R` (see existing examples).
4. Verify your example passes: `bash test.sh`
5. Verify the website builds: `quarto render`
6. Submit a Pull Request.


## Citation

If using EpiModel for teaching or research, please include a citation:

> Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. *Journal of Statistical Software.* 2018; 84(8): 1-47. doi: 10.18637/jss.v084.i08


## License

[MIT](LICENSE)
