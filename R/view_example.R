#' @title View simulation example
#' @description Users can view an .html file that contains an example with simulated abundance data.
#' The demo is similar to a vigette but it takes much too long to build. 
#' @author Devin S. Johnson
#' @export
#' 
view_example = function(){
  lib_loc = system.file(package="multAbund")
  browseURL(paste0('file://', file.path(lib_loc, "sim_example/simulation_demo.html")))
  rcode = paste0(system.file(package="multAbund"), "/sim_example/simulation_demo.R")
  message("R code available in the file:\n")
  message(paste0(system.file(package="multAbund"), "/sim_example/simulation_demo.R\n"))
} 