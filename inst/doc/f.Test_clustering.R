## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  out.width = "100%"
)

## ----init_testfclust, echo=FALSE, results="hide"-------------------------
library(functClust)
fres <- CedarCreek.2004.2006.res

## ----optplottestcomp, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Raw tree on left, Tree with components sorted by decreasing effect on right")----
fclust_plot(fres = CedarCreek.2004.2006.res, opt.tree = "prd")

ftest_plot(fres = CedarCreek.2004.2006.res, 
           rtest = CedarCreek.2004.2006.test.components,
           main = "BioDIV2", 
           opt.var = "comp", opt.crit = "Jaccard", opt.comp = "sorted.tree")

## ----optplottestass, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Raw tree on left, Tree with components sorted by decreasing effect on right")----
ftest_plot(fres = CedarCreek.2004.2006.res, 
           rtest = CedarCreek.2004.2006.test.assemblages, 
           main = "BioDIV2", opt.var = "assemblages")

## ----foptions_6----------------------------------------------------------
getOption("verbose")
# to follow the computations
options(verbose = TRUE)
# to deactivate the option
options(verbose = FALSE)

