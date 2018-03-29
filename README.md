
gClinBiomarker
==============


Overview
--------

gClinBiomarker is an R package that contains functions to perform baseline and longitutinal biomarker analyses


Two documents are provided: - A user vignette demonstrates how gClinBiomarker package may help in your biomarker analysis - An exapmle use case document which contains more detailed use cases

pdf version of these two documents can be found  at inst/doc

The Rmarkdown templates to run the automatic workflows can be found here[https://github.com/lengning/gClinBiomarker-templates]


Getting Started
---------------

### Installation

To install this package from R, use `install_github()` function from the **devtools** package

In R, type:

``` r
## install.packages("devtools")
library(devtools)
install_github("lengning/gClinBiomarker")
```


### Getting the Documentation

Use the command `vignette(package = 'gClinBiomarker')` to view a list of avaialble vignettes for the `gClinBiomarker` package. Then use the command `vignette(<vignett_name>)` to view the relevant documentation.

gClinBiomarker Analysis and Workflows
-------------------------------------

### Supported endpoint and biomarker types

#### Endpoints

-   Time-to-event
-   Binary
-   Continuous
-   Longitudinal

#### Biomarkers

-   Continuous
-   Ordinal
-   Categorical
-   Longitudinal

#### Analysis

-   2-arm study: predictive and prognostic biomarker analysis
-   1-arm study: prognostic biomarker analysis

### Workflow

#### Step 1:

Is your biomarker population your full patient population? Take a look at these functions:

`SummaryVars()`, `CompareKM()`, `PlotRspBar()`

#### Step 2:

Is your biomarker response a dynamic range? Does it have a skewed distribution? Is it correlated with clinical variables?

`PlotProperty()`, `PlotTabForestMulti()`

#### Step 3:

Is your biomarker response associated with a clinical outcome? Is it prognostic or predictive? Is there an optimal biomarker cutoff with a consistent trend?

`PlotTabForestBiomarker()`, `PlotSTEPP()`

#### Step 4:

Estimate the clinical benefit within a biomarker subgroup.

`PlotKM()`, `PlotLong()`, `PlotRspBar()`, `CoxTab()`, `LogRankTab()`
