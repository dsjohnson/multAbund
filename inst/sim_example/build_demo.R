setwd("~/research/projects/r_packages/multAbund/inst/sim_example")
knitr::purl("simulation_demo.Rmd")
knitr::knit("simulation_demo.Rmd")
rmarkdown::render("simulation_demo.Rmd")