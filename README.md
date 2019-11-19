# EpiModel-Gallery
[![Build Status](https://travis-ci.org/EpiModel/EpiModel-Gallery.svg?branch=master)](https://travis-ci.org/EpiModel/EpiModel-Gallery)

This repository contains templates for extending the EpiModel research platform to address new and interesting research questions. EpiModel (http://epimodel.org) is an R package that provides tools for simulating mathematical models of infectious disease dynamics. Epidemic model classes include deterministic compartmental models, stochastic agent-based models, and stochastic network models. For now, the gallery examples here are focused on extensions for network models.

EpiModel has both several built-in models to run epidemic modeling simulations right out of the box, as well as a flexible application programming interface (API) for extending EpiModel to build out complex epidemic systems not included in those built-in models. For both built-in and extension models, EpiModel allows users to estimate and simulate statistical models for dynamic contact/partnership networks using exponential random graph models (ERGMs). But because the learning curve for building new EpiModel modules for the epidemic side of the things can be steep, we are providing templates here so you can learn by example.

## How to Use
You may use and extend any of the examples in any way you want! We do ask if you publish a paper using EpiModel with or without these gallery examples, you include a citation (see below). You can access all of the source code of this repository by forking the repository or downloading a ZIP file. Forking will be safer because any ongoing updates to the repository may be pulled after your initial repository download. Before running the examples, you must install `EpiModel` on the R statistical computing platform. Instructions are at our primary Github repository (https://github.com/statnet/EpiModel). If you have any questions, either file an issue or send us an email.

## Requests
Requests for new gallery examples are welcome! If you'd like to ask, "how would you build a network model in EpiModel that does X?", this is the place to do that. We ask that you file a [Github Issue](https://github.com/EpiModel/EpiModel-Gallery/issues) in this repository with a detailed description of the research question. Before filing the request, we ask that you read through the documentation for EpiModel (https://github.com/statnet/EpiModel), including our primary methods paper at the [Journal of Statistical Software](https://www.jstatsoft.org/article/view/v084i08). If your request is very specific, we may generalize it so that the Gallery example may be of broad use to the EpiModel community.

## Contributions
Contributions of gallery examples are encouraged! We ask that you contribute an example by:

1. Forking this repository on Github
2. Adding your example locally, putting it in a separate subfolder that is named in the convention (`YEAR-Month-BasicDescription`). The contribution should have three files: `README.md` that explains what's going on with the model, and also contain your name/affiliation and any other details you'd like to include; a `model.R` file that is the main script for estimating the network model and running the simulation; and a `module-fx.R` file that contains the new/edited module functions that are plugged into `control.net`). Within the `model.R` file, you should include the two `# Standard Gallery unit test lines` block that you see in all the existing examples in this Gallery (this sets up proper unit testing of your example on Travis CI).
3. Push your local changes up to your fork, and then do a Pull Request to the main `EpiModel-Gallery` repository. Your Pull Request should pass testing on Travis CI (that is, run without errors). We will then review, request any changes if necessary, and then merge into our repository. We'll then add your name to our list of contributors.

## Citations
If using EpiModel for teaching or research, please include a citation:

> Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. *Journal of Statistical Software.* 2018; 84(8): 1-47. doi: 10.18637/jss.v084.i08
