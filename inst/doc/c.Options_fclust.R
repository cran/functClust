## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  out.width = "100%"
)

## ----init_opt_fclust, echo=FALSE, results="hide"-------------------------
library(functClust)
nbElt <- 16
dat.2004 <- CedarCreek.2004.2006.dat[ , c(1:(nbElt + 2))]

## ----opt_method, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.method determines the method of clustering: on left 'divisive', on right 'agglomerative'")----
res <- fclust(dat.2004, nbElt, opt.method = "divisive")
fclust_plot(res, opt.tree = list("prd"))

res <- fclust(dat.2004, nbElt, opt.method = "agglomerative")
fclust_plot(res, opt.tree = list("prd"))

## ----opt_apriori, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.method='apriori' forces the hierarchical tree to include a given component partition")----
apriori <- c(1,3,2,4,1,3,3,2,1,2,1,4,2,3,4,4) 
apriori <- c("F","C3","L","C4","F","C3","C3","L","F","L","F","C4","L","C3","C4","C4")

res <- fclust(dat.2004, nbElt, opt.method = "apriori", affectElt = apriori)
fclust_plot(res, opt.tree = list("prd", "leg", cols = apriori))

## ----opt_model, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.model determines the model for predicting assemblage performance")----
res <- fclust(dat.2004, nbElt, opt.mod = "bymot")
fclust_plot(res, opt.tree = list("prd"), opt.perf = list("prd", "aov", pvalue = 0.01))

res <- fclust(dat.2004, nbElt, opt.mod = "byelt")
fclust_plot(res, opt.tree = list("prd"), opt.perf = list("prd", "aov", pvalue = 0.01))

## ----opt_mean, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.mean determines the formula to use in averaging")----
res <- fclust(dat.2004, nbElt, opt.mean = "amean")
fclust_plot(res, opt.tree = list("prd"), opt.perf = list("prd", "aov", pvalue = 0.01))

res <- fclust(dat.2004, nbElt, opt.mean = "gmean")
fclust_plot(res, opt.tree = list("prd"), opt.perf = list("prd", "aov", pvalue = 0.01))

## ----foptions_3----------------------------------------------------------
getOption("verbose")
# to follow the computations
options(verbose = TRUE)
# to deactivate the option
options(verbose = FALSE)

