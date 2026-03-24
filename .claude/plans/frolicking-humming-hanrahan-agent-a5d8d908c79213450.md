# EpiModel-Gallery Quarto Website Scaffold — Implementation Plan

## Overview

This plan scaffolds a Quarto-based website for the EpiModel-Gallery repository. It renames all 12 example directories to kebab-case under a new `examples/` parent, creates the full Quarto site configuration, converts the SI with Vital Dynamics example into a fully annotated proof-of-concept `index.qmd`, creates placeholder pages for the other 11 examples, and updates the test infrastructure and CI/CD pipelines.

Branch: `new-website` (currently identical to `main`)

---

## Step 0: Pre-Flight Checks

Before any changes, verify:
- On `new-website` branch: `git branch --show-current`
- Working tree is clean: `git status`
- All 12 example directories exist with expected files

---

## Step 1: Create the `examples/` Directory and Rename All 12 Examples

Create `examples/` at the repository root, then move each old directory into it with the new kebab-case name.

### Commands (one per directory)

```bash
mkdir examples
git mv 2018-08-AddingAnExposedState    examples/seir-exposed-state
git mv 2018-08-ObservedNetworkData     examples/observed-network-data
git mv 2018-08-SIwithVitalDynamics     examples/si-vital-dynamics
git mv 2018-08-TestAndTreatIntervention examples/sis-test-and-treat
git mv 2018-09-CompetingStrains        examples/sis-competing-strains
git mv 2018-09-SocialDiffusion         examples/social-diffusion
git mv 2018-10-SEIRwithAONVax          examples/seir-aon-vaccination
git mv 2018-11-Syphilis                examples/syphilis
git mv 2018-12-SEIRSwithLeakyVax       examples/seirs-leaky-vaccination
git mv 2019-03-HIV                     examples/hiv
git mv 2021-10-CostEffectivenessAnalysis examples/cost-effectiveness
git mv 2022-12-Multinets               examples/multinets
```

### Verification
- `ls examples/` should list exactly 12 directories
- Each directory should contain its original README.md, model.R, and (for 11 of 12) module-fx.R

---

## Step 2: Update `source()` Paths in All `model.R` Files

Every `model.R` uses a `source("OLD-DIR-NAME/module-fx.R")` call that assumes execution from the repo root. After the rename, these paths must change. The Multinets example has no `source()` call.

### Old → New Path Mapping

| File | Old source() path | New source() path |
|------|-------------------|-------------------|
| `examples/seir-exposed-state/model.R` | `source("2018-08-AddingAnExposedState/module-fx.R")` | `source("examples/seir-exposed-state/module-fx.R")` |
| `examples/observed-network-data/model.R` | `source("2018-08-ObservedNetworkData/module-fx.R")` | `source("examples/observed-network-data/module-fx.R")` |
| `examples/si-vital-dynamics/model.R` | `source("2018-08-SIwithVitalDynamics/module-fx.R")` | `source("examples/si-vital-dynamics/module-fx.R")` |
| `examples/sis-test-and-treat/model.R` | `source("2018-08-TestAndTreatIntervention/module-fx.R")` | `source("examples/sis-test-and-treat/module-fx.R")` |
| `examples/sis-competing-strains/model.R` | `source("2018-09-CompetingStrains/module-fx.R")` | `source("examples/sis-competing-strains/module-fx.R")` |
| `examples/social-diffusion/model.R` | `source("2018-09-SocialDiffusion/module-fx.R")` | `source("examples/social-diffusion/module-fx.R")` |
| `examples/seir-aon-vaccination/model.R` | `source("2018-10-SEIRwithAONVax/module-fx.R")` | `source("examples/seir-aon-vaccination/module-fx.R")` |
| `examples/syphilis/model.R` | `source("2018-11-Syphilis/module-fx.R")` | `source("examples/syphilis/module-fx.R")` |
| `examples/seirs-leaky-vaccination/model.R` | `source("2018-12-SEIRSwithLeakyVax/module-fx.R")` | `source("examples/seirs-leaky-vaccination/module-fx.R")` |
| `examples/hiv/model.R` | `source("2019-03-HIV/module-fx.R")` | `source("examples/hiv/module-fx.R")` |
| `examples/cost-effectiveness/model.R` | `source("2021-10-CostEffectivenessAnalysis/module-fx.R")` | `source("examples/cost-effectiveness/module-fx.R")` |
| `examples/multinets/model.R` | (no source() call) | (no change needed) |

### Implementation
Use `sed` or direct editing to replace the single `source()` line in each file.

---

## Step 3: Update Cross-References in README.md Files

Multiple README files use `../OLD-DIR-NAME` relative links to reference other examples. These must be updated to use `../NEW-DIR-NAME` since after the move, all examples share the `examples/` parent.

### Links to Update

| File | Old Link Target | New Link Target |
|------|----------------|-----------------|
| `examples/si-vital-dynamics/README.md` | `../2018-08-AddingAnExposedState` | `../seir-exposed-state` |
| `examples/si-vital-dynamics/README.md` | `../2019-03-HIV` | `../hiv` |
| `examples/seir-exposed-state/README.md` | `../2018-08-SIwithVitalDynamics` | `../si-vital-dynamics` |
| `examples/seir-exposed-state/README.md` | `../2018-10-SEIRwithAONVax` | `../seir-aon-vaccination` |
| `examples/seir-exposed-state/README.md` | `../2018-12-SEIRSwithLeakyVax` | `../seirs-leaky-vaccination` |
| `examples/seir-exposed-state/README.md` | `../2019-03-HIV` | `../hiv` |
| `examples/sis-competing-strains/README.md` | `../2018-08-SIwithVitalDynamics` | `../si-vital-dynamics` |
| `examples/sis-competing-strains/README.md` | `../2018-08-TestAndTreatIntervention` | `../sis-test-and-treat` |
| `examples/sis-test-and-treat/README.md` | `../2018-08-SIwithVitalDynamics` | `../si-vital-dynamics` |
| `examples/sis-test-and-treat/README.md` | `../2018-09-CompetingStrains` | `../sis-competing-strains` |
| `examples/sis-test-and-treat/README.md` | `../2018-11-Syphilis` | `../syphilis` |
| `examples/sis-test-and-treat/README.md` | `../2021-10-CostEffectivenessAnalysis` | `../cost-effectiveness` |
| `examples/observed-network-data/README.md` | `../2018-08-SIwithVitalDynamics` | `../si-vital-dynamics` |
| `examples/seir-aon-vaccination/README.md` | `../2018-12-SEIRSwithLeakyVax` | `../seirs-leaky-vaccination` |
| `examples/seir-aon-vaccination/README.md` | `../2018-08-SIwithVitalDynamics` | `../si-vital-dynamics` |
| `examples/seir-aon-vaccination/README.md` | `../2021-10-CostEffectivenessAnalysis` | `../cost-effectiveness` |
| `examples/seirs-leaky-vaccination/README.md` | `../2018-10-SEIRwithAONVax` | `../seir-aon-vaccination` |
| `examples/seirs-leaky-vaccination/README.md` | `../2018-08-SIwithVitalDynamics` | `../si-vital-dynamics` |

### Implementation
A batch `sed` operation across all `examples/*/README.md` files:
```bash
sed -i '' 's|\.\./2018-08-AddingAnExposedState|../seir-exposed-state|g' examples/*/README.md
sed -i '' 's|\.\./2018-08-ObservedNetworkData|../observed-network-data|g' examples/*/README.md
sed -i '' 's|\.\./2018-08-SIwithVitalDynamics|../si-vital-dynamics|g' examples/*/README.md
sed -i '' 's|\.\./2018-08-TestAndTreatIntervention|../sis-test-and-treat|g' examples/*/README.md
sed -i '' 's|\.\./2018-09-CompetingStrains|../sis-competing-strains|g' examples/*/README.md
sed -i '' 's|\.\./2018-09-SocialDiffusion|../social-diffusion|g' examples/*/README.md
sed -i '' 's|\.\./2018-10-SEIRwithAONVax|../seir-aon-vaccination|g' examples/*/README.md
sed -i '' 's|\.\./2018-11-Syphilis|../syphilis|g' examples/*/README.md
sed -i '' 's|\.\./2018-12-SEIRSwithLeakyVax|../seirs-leaky-vaccination|g' examples/*/README.md
sed -i '' 's|\.\./2019-03-HIV|../hiv|g' examples/*/README.md
sed -i '' 's|\.\./2021-10-CostEffectivenessAnalysis|../cost-effectiveness|g' examples/*/README.md
sed -i '' 's|\.\./2022-12-Multinets|../multinets|g' examples/*/README.md
```

---

## Step 4: Create `_quarto.yml`

Create `/Users/sjennes/git/EpiModel-Gallery/_quarto.yml` with the following content:

```yaml
project:
  type: website
  output-dir: _site

execute:
  freeze: auto

website:
  title: "EpiModel Gallery"
  description: "Template examples for extending EpiModel to model infectious disease dynamics over networks"
  repo-url: https://github.com/EpiModel/EpiModel-Gallery
  repo-actions: [source, issue]
  page-navigation: true

  navbar:
    background: primary
    search: true
    left:
      - text: "Home"
        href: index.qmd
      - text: "Examples"
        href: examples/si-vital-dynamics/index.qmd
      - text: "About"
        href: about.qmd
    right:
      - icon: github
        href: https://github.com/EpiModel/EpiModel-Gallery

  sidebar:
    - title: "Examples"
      style: docked
      search: true
      contents:
        - section: "Fundamentals"
          contents:
            - examples/si-vital-dynamics/index.qmd
            - examples/seir-exposed-state/index.qmd
            - examples/sis-test-and-treat/index.qmd
        - section: "Advanced Dynamics"
          contents:
            - examples/sis-competing-strains/index.qmd
            - examples/social-diffusion/index.qmd
            - examples/syphilis/index.qmd
        - section: "Vaccination"
          contents:
            - examples/seir-aon-vaccination/index.qmd
            - examples/seirs-leaky-vaccination/index.qmd
        - section: "Applied"
          contents:
            - examples/hiv/index.qmd
            - examples/cost-effectiveness/index.qmd
        - section: "Network Features"
          contents:
            - examples/observed-network-data/index.qmd
            - examples/multinets/index.qmd

format:
  html:
    theme: cosmo
    css: styles.css
    toc: true
    code-copy: true
    code-overflow: wrap
    highlight-style: github
```

### Design Decisions

- **`freeze: auto`**: Pre-rendered output is committed to `_freeze/` and only re-rendered when source changes. This is critical because R + EpiModel are not available in the Quarto publish CI step -- all R rendering happens locally and the frozen output is committed.
- **`cosmo` theme**: Clean, professional Bootstrap theme suitable for a technical gallery.
- **Navbar + Sidebar hybrid**: The navbar provides top-level site navigation (Home, Examples, About). Clicking "Examples" activates the sidebar which shows the tiered difficulty organization.
- **`repo-actions: [source, issue]`**: Each page gets "View source" and "Report issue" links to GitHub.
- **`page-navigation: true`**: Previous/next links at the bottom of each example page, following the sidebar order (which is the learning path order).

---

## Step 5: Create `index.qmd` (Landing Page)

Create `/Users/sjennes/git/EpiModel-Gallery/index.qmd`:

```yaml
---
title: "EpiModel Gallery"
subtitle: "Template examples for extending EpiModel to model infectious disease dynamics over networks"
listing:
  id: gallery-listing
  contents: examples/*/index.qmd
  type: grid
  grid-columns: 3
  image-height: 200px
  fields: [image, title, description, categories]
  sort: "order"
  categories: true
page-layout: full
---
```

Followed by Markdown body content:
- Brief introduction paragraph explaining what the Gallery is
- Reference to installing EpiModel
- The `:::{#gallery-listing}` div to place the listing
- A "Getting Started" section with quick instructions

### How the Listing Works

Each example's `index.qmd` frontmatter will contain `title`, `description`, `image` (pointing to `thumbnail.png`), `categories` (the tier name + model type tags), and a custom `order` field for sorting. The listing automatically generates cards from these metadata fields.

---

## Step 6: Create `about.qmd`

Create `/Users/sjennes/git/EpiModel-Gallery/about.qmd` with:
- Title "About the EpiModel Gallery"
- Sections on: what EpiModel is, how to use the Gallery, contributing, citation info
- Content adapted from the current `README.md`

---

## Step 7: Create `styles.css`

Create `/Users/sjennes/git/EpiModel-Gallery/styles.css` with minimal custom styles:
- Card hover effects for the gallery grid
- Category badge styling
- Any typography adjustments

This can be minimal at Phase 1 -- Quarto's cosmo theme provides a strong baseline.

---

## Step 8: Create `examples/_metadata.yml`

Create `/Users/sjennes/git/EpiModel-Gallery/examples/_metadata.yml` with shared defaults for all example pages:

```yaml
# Shared metadata for all example pages
execute:
  freeze: auto
  echo: true
  warning: false
  message: false

code-annotations: below

format:
  html:
    toc: true
    toc-depth: 3
    code-tools:
      source: true
    code-copy: true
    code-fold: false
```

### Design Decisions

- **`freeze: auto`**: Redundant with project-level setting but explicit for clarity; ensures no example re-renders in CI.
- **`echo: true`**: All code is shown (this is a teaching resource).
- **`warning: false` / `message: false`**: Suppress R startup messages and package warnings in rendered output.
- **`code-annotations: below`**: Enable Quarto's numbered code annotations for all examples. Annotations display as a numbered list below each code block.
- **`code-tools: source: true`**: Adds a "Code" button to toggle raw source view.

---

## Step 9: Create the Proof-of-Concept — `examples/si-vital-dynamics/index.qmd`

This is the most substantial file. It merges the content from `README.md`, `model.R`, and `module-fx.R` into a single annotated Quarto document.

### File: `examples/si-vital-dynamics/index.qmd`

#### Frontmatter

```yaml
---
title: "SI Model with Age-Specific Vital Dynamics"
description: "SI epidemic with aging, births, deaths, and age-specific mortality on a dynamic network with age-assortative mixing"
author: "Samuel M. Jenness"
date: "2018-08-01"
image: thumbnail.png
categories:
  - Fundamentals
  - SI
  - Vital Dynamics
order: 1

execute:
  freeze: auto

resources:
  - model.R
  - module-fx.R
---
```

#### Key Design Decisions for the Frontmatter

- **`image: thumbnail.png`**: Referenced by the gallery listing. Initially this file won't exist -- we'll generate it after the first render, or provide a placeholder.
- **`categories`**: Used for the gallery landing page category filters. "Fundamentals" is the tier; "SI" and "Vital Dynamics" are model-type tags.
- **`order: 1`**: Controls sort position in the gallery grid listing.
- **`resources: [model.R, module-fx.R]`**: Tells Quarto to copy these files to the output directory, making them downloadable.

#### Document Structure

The document body should follow this structure, interleaving narrative from the README with executable code from model.R and module-fx.R:

```
## Overview
(Adapted from README description section)

## Model Structure
### Disease Compartments
(Table from README)

### Flow Diagram
(Mermaid diagram from README -- Quarto renders mermaid natively)

### Vital Dynamics
(Narrative from README sections on aging, departures, arrivals)

## Module Functions

### Aging Module
```{r}
#| code-annotations: below
aging <- function(dat, at) {
  age <- get_attr(dat, "age")    # <1>
  age <- age + 1 / 52            # <2>
  dat <- set_attr(dat, "age", age)

  dat <- set_epi(dat, "meanAge", at, mean(age, na.rm = TRUE))  # <3>
  return(dat)
}
```
1. Retrieve current age attribute for all nodes
2. Increment by 1/52 years (one week per timestep)
3. Record population mean age as an epidemiological summary statistic

(Similar annotated blocks for dfunc and afunc)

### Download Standalone Scripts
::: {.callout-tip}
Download the standalone scripts: [model.R](model.R) | [module-fx.R](module-fx.R)
:::

## Network Setup
(Narrative + annotated code from sections 1-2 of model.R)

## Epidemic Parameters
(Narrative + annotated code from section 3 of model.R, including the source() call but using inline definition instead)

## Simulation: Baseline vs. Lethal Disease
(Narrative + code from sections 4-5 of model.R)

## Analysis & Results
(Narrative + plotting code from section 6 of model.R, with code-fold: true for plots)

## Parameters Reference
(Parameter tables from README)

## Module Execution Order
(From README)

## Next Steps
(From README, but with links updated to use Quarto site-relative paths)
```

#### Critical Implementation Notes

1. **The `source()` call**: In the rendered `.qmd`, the module functions are defined inline (with annotations), so there is no need for `source("examples/si-vital-dynamics/module-fx.R")`. The standalone `model.R` retains its `source()` call for direct execution.

2. **Code annotations syntax**: Each annotated line in R code blocks gets a `# <N>` comment at the end. An ordered list immediately after the code fence provides the annotation text. Example:
   ````
   ```{r}
   x <- 1 + 1        # <1>
   y <- sqrt(x)       # <2>
   ```
   1. Add one plus one
   2. Take the square root
   ````

3. **Mermaid diagrams**: Quarto renders mermaid natively with ` ```{mermaid} ` fenced blocks. The existing mermaid from the README can be used directly.

4. **Plot output**: Plots generated by the R code will be captured and displayed inline. For the initial freeze, run `quarto render examples/si-vital-dynamics/index.qmd` locally with R + EpiModel installed. This populates `_freeze/examples/si-vital-dynamics/` with the pre-rendered output.

5. **Thumbnail generation**: After the first render, take a screenshot of the prevalence comparison plot and save as `examples/si-vital-dynamics/thumbnail.png`. Alternatively, add R code at the end that uses `png()` to save a specific plot to `thumbnail.png` and exclude that chunk from the document output.

---

## Step 10: Create Placeholder `index.qmd` for the Other 11 Examples

For each of the 11 remaining examples, create an `index.qmd` with frontmatter only and a brief placeholder body. The content will be filled in during later phases.

### Template for Placeholder Pages

```yaml
---
title: "[Title from README.md H1]"
description: "[First sentence of README description]"
author: "[Author from README]"
date: "[Original date from directory name]"
image: thumbnail.png
categories:
  - [Tier name]
  - [Model type tags]
order: [N]
draft: true
execute:
  eval: false
---

::: {.callout-note}
This example is under construction. The standalone R scripts are available for download below.
:::

## Overview

[First paragraph from README.md]

## Download

- [model.R](model.R)
- [module-fx.R](module-fx.R)
```

### Specific Frontmatter for Each Placeholder

| Directory | title | categories | order |
|-----------|-------|-----------|-------|
| seir-exposed-state | SEIR/SEIRS: Adding an Exposed State | Fundamentals, SEIR, SEIRS | 2 |
| sis-test-and-treat | SIS with Test-and-Treat Intervention | Fundamentals, SIS, Intervention | 3 |
| sis-competing-strains | SIS with Competing Pathogen Strains | Advanced Dynamics, SIS, Multi-strain | 4 |
| social-diffusion | Social Diffusion on Networks | Advanced Dynamics, SI, Social Diffusion | 5 |
| syphilis | Multi-Stage Syphilis Model | Advanced Dynamics, Multi-stage, STI | 6 |
| seir-aon-vaccination | SEIR with All-or-Nothing Vaccination | Vaccination, SEIR, Vaccine | 7 |
| seirs-leaky-vaccination | SEIRS with Leaky Vaccination | Vaccination, SEIRS, Vaccine | 8 |
| hiv | HIV Transmission Model | Applied, HIV, Multi-stage, ART | 9 |
| cost-effectiveness | Cost-Effectiveness Analysis | Applied, SI, CEA | 10 |
| observed-network-data | Epidemics over Observed Network Data | Network Features, SI, Observed Data | 11 |
| multinets | Multiple Interacting Networks | Network Features, Multilayer | 12 |

### Note on `draft: true`
Setting `draft: true` in the frontmatter means these pages will be visible during local development (`quarto preview`) but excluded from the production build by default. This is ideal for Phase 1 where only the SI Vital Dynamics example is complete. When ready, remove `draft: true` from each example as it's converted.

**Alternative approach**: If you want all cards to appear in the gallery grid even as placeholders, omit `draft: true` and instead use `execute: eval: false` (already set) so no R code runs. The pages will render with just the narrative.

**Recommendation**: Omit `draft: true` so the gallery grid shows all 12 cards (demonstrating the full scope), but keep the "under construction" callout so visitors know the content isn't final. The `execute: eval: false` setting prevents any R code from running in placeholders.

---

## Step 11: Update `test.sh`

The current `test.sh` iterates over top-level directories. After the rename, examples live in `examples/`. The script must change its directory scanning.

### New `test.sh`

```bash
#!/bin/bash

set -e

d="$(ls -p examples/ | grep "/")"
dl="$(find examples/* -maxdepth 0 -type d | wc -l)"

echo " "
echo "#############################################"
echo "EpiModel Gallery Testing:" $dl "Directories"
echo "---------------------------------------------"

for i in $d; do
    dir="examples/$i"
    model="$dir/model.R"

    # Skip directories without model.R (shouldn't happen, but defensive)
    if [ ! -f "$model" ]; then
        continue
    fi

    echo -n "$dir ... "
    SECONDS=0
    if Rscript "$model" "options(error = function() q('no', 1, FALSE))" >& /dev/null
    then
        echo "OK" "($SECONDS seconds)"
    else
        test_failed=1
        echo "Failed"
    fi
    # Clean up generated PDFs
    rm -f "$dir"/*.pdf
done

rm -f *.pdf

echo "#############################################"
echo " "

if [ -z ${test_failed+x} ]
then
    echo "All Tests Succcessful"
    exit 0
else
    echo "Some Tests Failed"
    exit 1
fi
```

### Key Changes
- Scans `examples/` instead of the repo root
- Uses `examples/$i` prefix for all paths
- Removes the `renv/` skip (no longer needed since we only scan `examples/`)
- Cleans up PDFs in the correct subdirectory

---

## Step 12: Update `.gitignore`

Add Quarto-specific entries to `.gitignore`:

```
# Quarto
_site/
.quarto/

# Keep _freeze/ tracked (DO NOT add _freeze/ to .gitignore)
```

Also remove `*.html` from `.gitignore` since Quarto site build goes to `_site/` (already excluded) and we may want HTML files in examples. Actually, `*.html` in the gitignore may conflict with Quarto -- but since output goes to `_site/` and that's separately ignored, it's fine. However, remove the specific `2021-10-CostEffectivenessAnalysis/README.html` file that is currently tracked.

### Files to add to .gitignore
```
# Quarto output
/_site/
/.quarto/
```

### Note on `_freeze/`
The `_freeze/` directory MUST be committed to the repository. It contains the pre-rendered R output that the CI publish step uses without needing R installed. Do NOT add `_freeze/` to `.gitignore`.

---

## Step 13: Create `.github/workflows/publish.yml`

This workflow publishes the Quarto website to GitHub Pages. Because R rendering is done locally (via freeze), the CI step only needs Quarto, not R.

### File: `.github/workflows/publish.yml`

```yaml
on:
  workflow_dispatch:
  push:
    branches: main

name: Publish Website

jobs:
  build-deploy:
    runs-on: ubuntu-latest
    permissions:
      contents: write
    steps:
      - name: Check out repository
        uses: actions/checkout@v4

      - name: Set up Quarto
        uses: quarto-dev/quarto-actions/setup@v2

      - name: Render and Publish
        uses: quarto-dev/quarto-actions/publish@v2
        with:
          target: gh-pages
        env:
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
```

### Design Decisions

- **No R setup needed in CI**: Because `freeze: auto` is used and `_freeze/` is committed, Quarto skips R execution during the publish step. It only assembles the HTML from the frozen output.
- **Trigger on `main` only**: The website publishes when changes land on `main` (via PR merge).
- **`workflow_dispatch`**: Allows manual triggering for re-publishes.
- **Separate from `tests.yml`**: The test workflow validates R scripts; the publish workflow builds the website. They serve different purposes and have different dependencies.

### Pre-Requisite
Before the first publish, run `quarto publish gh-pages` locally once to create the `_publish.yml` file. Also ensure the GitHub repository settings have Pages configured to deploy from the `gh-pages` branch.

---

## Step 14: Update `.github/workflows/tests.yml`

The existing test workflow runs `bash test.sh`. Since we updated `test.sh` in Step 11, the workflow file itself needs minimal changes. However, we should verify it still works.

### Changes Needed
None to the workflow file itself -- it just runs `bash test.sh`, which we updated. But verify that the R dependency installation still works (it uses `setup-r-dependencies` which reads `DESCRIPTION`).

### One concern
The `DESCRIPTION` file lists the package dependencies. This should remain at the repo root and continue to work with `r-lib/actions/setup-r-dependencies@v2`. No changes needed.

---

## Step 15: Update the Top-Level `README.md`

Update the repository README to:
1. Add a link to the live website
2. Update the example table with new directory paths
3. Update the "Running an Example" instructions to use `examples/` paths
4. Keep the Contributing, Citation, and License sections

### Key Changes

```markdown
# EpiModel-Gallery

**[View the Gallery Website](https://epimodel.github.io/EpiModel-Gallery/)**

...

### Running an Example

```r
source("examples/si-vital-dynamics/model.R")
```

Or from the command line:

```bash
Rscript examples/si-vital-dynamics/model.R
```
```

Also update the Gallery Examples table:

| Example | Description |
|---------|-------------|
| [si-vital-dynamics](examples/si-vital-dynamics/) | SI with aging, births, deaths, and age-specific mortality |
| ... | ... |

---

## Step 16: Create a Placeholder `thumbnail.png`

For the proof-of-concept example (`si-vital-dynamics`), generate a thumbnail after the first local render. For the other 11 examples, we need placeholder images.

### Options
1. **Generate from R**: Add a code chunk in the proof-of-concept `index.qmd` that saves a plot as `thumbnail.png` (set `include: false` so it doesn't show in the document).
2. **Placeholder SVG/PNG**: Create a simple placeholder image for examples that aren't fully rendered yet.
3. **Skip for now**: The gallery listing will show cards without images. Quarto handles missing images gracefully.

**Recommendation**: For Phase 1, create thumbnails only for `si-vital-dynamics` (generate from the prevalence plot). For the other 11, omit the `image` field from the frontmatter until they are fully rendered. The gallery cards will still display title, description, and categories.

---

## Step 17: Delete the Legacy `README.html`

The file `examples/cost-effectiveness/README.html` (formerly `2021-10-CostEffectivenessAnalysis/README.html`) is a pre-rendered HTML file from the old setup. Delete it:

```bash
git rm examples/cost-effectiveness/README.html
```

---

## Step 18: Initial Local Render and Freeze

After all files are created, perform the initial render:

```bash
# Render the full site (requires R + EpiModel installed locally)
quarto render

# Or render just the proof-of-concept to populate _freeze/
quarto render examples/si-vital-dynamics/index.qmd
```

This creates:
- `_freeze/examples/si-vital-dynamics/index/` — frozen R output (commit this)
- `_site/` — built website (do not commit, it's in .gitignore)

Then commit `_freeze/` to the repository.

---

## Commit Strategy

### Commit 1: Directory renames and path updates
- `git mv` all 12 directories
- Update `source()` paths in all `model.R` files
- Update cross-references in all `README.md` files
- Update `test.sh`
- Update top-level `README.md`
- Delete `examples/cost-effectiveness/README.html`

### Commit 2: Quarto site scaffold
- `_quarto.yml`
- `index.qmd`
- `about.qmd`
- `styles.css`
- `examples/_metadata.yml`
- `.gitignore` updates
- 11 placeholder `index.qmd` files

### Commit 3: Proof-of-concept example
- `examples/si-vital-dynamics/index.qmd` (fully annotated)
- `examples/si-vital-dynamics/thumbnail.png` (generated)

### Commit 4: Frozen output + CI
- `_freeze/` directory (pre-rendered output)
- `.github/workflows/publish.yml`

### Commit 5 (optional): `_publish.yml`
- Generated by `quarto publish gh-pages` on first local publish

---

## File Inventory — All Files to Create or Modify

### New Files (13 files)
| File | Purpose |
|------|---------|
| `_quarto.yml` | Site configuration |
| `index.qmd` | Landing page with gallery grid |
| `about.qmd` | About page |
| `styles.css` | Custom CSS |
| `examples/_metadata.yml` | Shared example defaults |
| `examples/si-vital-dynamics/index.qmd` | Full proof-of-concept |
| `examples/seir-exposed-state/index.qmd` | Placeholder |
| `examples/sis-test-and-treat/index.qmd` | Placeholder |
| `examples/sis-competing-strains/index.qmd` | Placeholder |
| `examples/social-diffusion/index.qmd` | Placeholder |
| `examples/syphilis/index.qmd` | Placeholder |
| `examples/seir-aon-vaccination/index.qmd` | Placeholder |
| `examples/seirs-leaky-vaccination/index.qmd` | Placeholder |
| `examples/hiv/index.qmd` | Placeholder |
| `examples/cost-effectiveness/index.qmd` | Placeholder |
| `examples/observed-network-data/index.qmd` | Placeholder |
| `examples/multinets/index.qmd` | Placeholder |
| `.github/workflows/publish.yml` | GitHub Pages deployment |
| `examples/si-vital-dynamics/thumbnail.png` | Gallery card image |

### Modified Files (16 files)
| File | Change |
|------|--------|
| `test.sh` | Scan `examples/` instead of repo root |
| `README.md` | Update paths, add website link |
| `.gitignore` | Add `_site/`, `.quarto/` |
| `examples/seir-exposed-state/model.R` | Update `source()` path |
| `examples/observed-network-data/model.R` | Update `source()` path |
| `examples/si-vital-dynamics/model.R` | Update `source()` path |
| `examples/sis-test-and-treat/model.R` | Update `source()` path |
| `examples/sis-competing-strains/model.R` | Update `source()` path |
| `examples/social-diffusion/model.R` | Update `source()` path |
| `examples/seir-aon-vaccination/model.R` | Update `source()` path |
| `examples/syphilis/model.R` | Update `source()` path |
| `examples/seirs-leaky-vaccination/model.R` | Update `source()` path |
| `examples/hiv/model.R` | Update `source()` path |
| `examples/cost-effectiveness/model.R` | Update `source()` path |
| Multiple `README.md` files | Update cross-reference links |

### Deleted Files (1 file)
| File | Reason |
|------|--------|
| `examples/cost-effectiveness/README.html` | Legacy pre-rendered file |

### Moved/Renamed (12 directories)
All 12 example directories from repo root into `examples/` with new names.

### Generated (committed, not hand-authored)
| File | Purpose |
|------|---------|
| `_freeze/examples/si-vital-dynamics/index/` | Pre-rendered R output |

---

## Potential Issues and Mitigations

### 1. test.sh CWD Assumption
The current `model.R` files use `source("DIRNAME/module-fx.R")` which assumes execution from the repo root. The new paths (`source("examples/DIRNAME/module-fx.R")`) maintain this assumption. The `test.sh` script runs `Rscript` from the repo root, so this continues to work. **No issue.**

### 2. Quarto execute-dir
By default, Quarto's `execute-dir` is `file` (code executes with the `.qmd` file's directory as the working directory). However, `model.R` scripts expect the repo root as the working directory. For the `index.qmd` proof-of-concept, the module functions are defined inline, so `source()` is not needed. But if someone renders a `.qmd` that calls `source("examples/si-vital-dynamics/module-fx.R")`, it would fail because the CWD is the example directory.

**Mitigation**: In `_quarto.yml`, set `execute-dir: project` so all R code executes from the repo root. This makes the `source()` paths work consistently.

```yaml
project:
  type: website
  output-dir: _site
  execute-dir: project
```

### 3. Mermaid Rendering
Quarto renders mermaid diagrams natively. The existing README mermaid blocks use HTML-style `<b>` tags and `<br/>` inside mermaid nodes. Verify these render correctly in Quarto's mermaid implementation. If not, simplify the node labels.

### 4. Missing Thumbnails
The gallery listing will show cards without images for the 11 placeholder examples. This is acceptable for Phase 1. The `image: thumbnail.png` field should be omitted from placeholder frontmatter (not pointing to a nonexistent file).

### 5. Placeholder Pages and Listing
If placeholder pages have `draft: true`, they won't appear in the listing. Since we want all 12 cards in the gallery, we should NOT use `draft: true`. Instead, placeholders should have `execute: eval: false` and an "under construction" callout.

### 6. The EpiModel-Gallery.Rproj File
This file is in `.gitignore` already. No changes needed.

---

## Summary of Implementation Order

1. **Rename directories** (git mv) + update source/cross-ref paths
2. **Create Quarto config** (_quarto.yml, _metadata.yml, styles.css)
3. **Create landing and about pages** (index.qmd, about.qmd)
4. **Create 11 placeholder index.qmd files**
5. **Create the proof-of-concept** (si-vital-dynamics/index.qmd)
6. **Update test.sh** and verify `bash test.sh` passes
7. **Update .gitignore**
8. **Update top-level README.md**
9. **Delete legacy files** (README.html)
10. **Local render** to populate _freeze/
11. **Create publish.yml** workflow
12. **Commit and push**
