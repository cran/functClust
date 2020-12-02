## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  out.width = "100%"
)

## ----pressure, echo=FALSE, fig.cap="Method for clustering functionally redundant components: (a) Component clustering, (b) Identification of assembly motifs, (c) Classification of sub-systems by assembly motif, (d) Computation of convergence criterion, (e) Iterative ajustment of component clustering", out.width = '100%'----
knitr::include_graphics("a.Figure.2.png")

## ----format_data---------------------------------------------------------
library(functClust)

# production of biomass in 2004, 2005 and 2006 in Biodiversity II experiment
data(CedarCreek.2004.2006.dat)
dim(CedarCreek.2004.2006.dat)
colnames(CedarCreek.2004.2006.dat)
rownames(CedarCreek.2004.2006.dat)

# the data for each component assemblage
line <- 10
CedarCreek.2004.2006.dat[line, ]

## ----foptions_1----------------------------------------------------------
getOption("verbose")
# to follow the computations
options(verbose = TRUE)
# to deactivate the option
options(verbose = FALSE)

## ----install-------------------------------------------------------------
#install.packages("functClust")

