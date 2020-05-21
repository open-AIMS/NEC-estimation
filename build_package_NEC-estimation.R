

# This script includes the basic steps required to create and
#    build the jagsNEC package

library(devtools)
library(roxygen2)
library(knitr)
library(R.rsp)
library(digest)


devtools::document()

use_package("R2jags")

build()


