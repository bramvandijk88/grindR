# grindR

Phase plane analysis and parameter estimation for ODE models. A wrapper around
deSolve, rootSolve and FME providing five functions: `run()` (numerical
integration), `plane()` (phase plane and nullclines), `newton()` (steady
states), `continue()` (bifurcation diagrams) and `fit()` (parameter
estimation).

## Installation

grindR is published on [R-universe](https://bramvandijk88.r-universe.dev).
Install it together with its dependencies (from CRAN) using:

```r
install.packages("grindR",
  repos = c("https://bramvandijk88.r-universe.dev",
            "https://cloud.r-project.org"))
```

Both repositories are needed: grindR comes from R-universe, while its
dependencies (deSolve, rootSolve, FME) come from CRAN.

## Usage

```r
library(grindR)

model <- function(t, state, parms)
  with(as.list(c(state, parms)), list(c(r * N * (1 - N / K))))
p <- c(r = 0.5, K = 10)
s <- c(N = 0.1)

run(tmax = 20)               # integrate and plot the time course
plane(xmax = 12, ymax = 1)   # phase plane with nullclines
```

## Credits

Original author: Rob J. de Boer (Utrecht University). Maintained by Bram van
Dijk. See <https://tbb.bio.uu.nl/rdb/grind.html>.
