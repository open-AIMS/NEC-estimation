

# This script includes the basic steps required to create and
#    build the jagsNEC package

library(devtools)
library(roxygen2)
library(knitr)
library(R.rsp)
library(digest)

styler::style_pkg(filetype = c("R", "Rmd"))
#lintr::lint_package()
#devtools::test()

devtools::document()

rmarkdown::render("README.Rmd", output_format = "md_document")


use_package("R2jags")
#pkgdown::build_site()
#
build()

#devtools::check()

