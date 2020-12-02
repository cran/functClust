## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  out.width = "100%"
)

## ----simplest_use_fclust, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("The simplest use of fclust")----
library(functClust)
nbElt <- 16
dat.2004 <- CedarCreek.2004.2006.dat[ , c(1, 1 + 1:nbElt, 1 + nbElt + 1)]
res <- fclust(dat.2004, nbElt)

## ----simplest_use_plotfclust, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("The simplest use of fclust_plot")----
fclust_plot(res, main = "BioDiv2 2004")

## ----simplest_use_writefclust, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("The simplest use of fclust_plot")----
# fclust_write(res, filename = "myRecord")
# res <- fclust_read(filename = "myRecord")

## ----foptions_2----------------------------------------------------------
getOption("verbose")
# to follow the computations
options(verbose = TRUE)
# to deactivate the option
options(verbose = FALSE)

