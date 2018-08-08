# EpiModel-Gallery

## About
This repository contains example templates to extend the EpiModel research platform to address novel research questions. EpiModel (http://epimodel.org) is an R package that provides tools for simulating mathematical models of infectious disease dynamics. Epidemic model classes include deterministic compartmental models, stochastic agent-based models, and stochastic network models. The gallery examples here are all focused on extending network models.

## How to Use
You may use and extend any of the examples in any way you want! We do ask if you publish a paper using EpiModel with or without these gallery examples, you include a citation (see below). You can access all of the source code of this repository either by cloning the repository or downloading a ZIP file. Both are available with the Green button at the top. Github cloning will be safer, of course, because our changes to the repository may be pulled after your initial clone. Before running the examples, you must install `EpiModel` and the R statistical computing platform. Instructions are at our primary Github repository (https://github.com/statnet/EpiModel). If you have any questions, either file an issue or send us an email.

## Requests
Requests for new gallery examples are welcome! If you'd like to ask, "how would you build a network model in EpiModel that does X?", this is the place to do? We ask that you file a [Github Issue](https://github.com/statnet/EpiModel-Gallery/issues) in this repository with a relatively detailed description of the problem. Before filing the request, we ask that you read through the documentation for EpiModel (https://github.com/statnet/EpiModel), including our primary methods paper at the [Journal of Statistical Software](https://www.jstatsoft.org/article/view/v084i08). After that ask away! We may take your specific request and generalize it a little bit so it may be of more broad use to others in the EpiModel community.

## Contributions
Contributions of gallery examples are encouraged! We ask that you contribute an example by:

1. Forking this repository on Github
2. Adding your example locally, putting it in a separate subfolder that is named in the convention (`YEAR-Month_DiseaseType-BasicDescription`). Ideally, the contribution should have three files: `README.md` that explains what's going on with the model, and also contain your name/affiliation and any other details you'd like to include; a `model.R` file that is the main script for estimating the network model and running the simulation; and a `module-fx.R` file that contains the new/edited module functions that are plugged into `control.net`).
3. Push your local changes up to your fork, and then do a Pull Request to the main `EpiModel-Gallery` repository. We will then review, request any changes if necessary, and then merge into our repository. We'll then add your name to our list of contributors.

## Citations
If using EpiModel for teaching or research, please include a citation:

> Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. *Journal of Statistical Software.* 2018; 84(8): 1-47. doi: 10.18637/jss.v084.i08
