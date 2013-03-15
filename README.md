phylodesign
===========

build this package
R CMD build pkg: builds the package (generates an archive called foo.tar.gz)
R CMD build pkg --compact-vignettes --resave-data: same, making the package as small as possible
R CMD check phylodesign_1.0-0.gz --as-cran: runs package quality checks similar to the checks made by CRAN; must be passed without error/warning before 


install this R package:
R CMD INSTALL phylodesign_1.0-0.tar.gz

to view the package content in R:
require(phylodesign)
library(help=phylodesign)