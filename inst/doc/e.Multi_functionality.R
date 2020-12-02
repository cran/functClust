## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  out.width = "100%"
)

## ----init_plotfclust, echo=FALSE, results="hide"-------------------------
library(functClust)
nbElt <- 16
dat <- CedarCreek.2004.2006.dat

## ----init_opt_weight, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Weight = c(1,1,1) on left, = c(2,3,1) on right tree")----
res <- fclust(dat, nbElt, weight = c(1, 1, 1), opt.mod = "byelt", opt.mean = "gmean")
fclust_plot(res, main = "BioDiv2", opt.tree = list("prd"))

res2 <- fclust(dat, nbElt, weight = c(2, 3, 1), opt.mod = "byelt", opt.mean = "gmean")
fclust_plot(res2, main = "BioDiv2", opt.tree = list("prd"))

## ----plot_opt_xpr_perf, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.perf = list('xpr') plots the performances over time")----
fclust_plot(res, main = "BioDiv2", opt.perf = list("multi"))

## ----plot_opt_xpr_ass, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.ass = list(who = c('20', '138')) plots the performances over time of assemblages '20' and '138' ")----
fclust_plot(res, main = "BioDiv2", opt.ass = list(who = c("20", "138")))

## ----plot_opt_xpr_motif, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.motif = list('obs', 'seq') plots performances over time sorted by assembly motif")----
fclust_plot(res, main = "BioDiv2", opt.motif = list("prd", "hor", "multi", "aov", pvalue = 0.01))

## ----plot_opt_xpr_comp, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.comp = list('tree') plots performances over time sorted by component clusters as in hierarchical tree")----
fclust_plot(res, main = "BioDiv2", opt.comp = list("tree", "multi", "aov", pvalue = 0.01))

