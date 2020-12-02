## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  out.width = "100%"
)

## ----init_plotfclust, echo=FALSE, results="hide"-------------------------
library(functClust)
nbElt <- 16
dat.2004 <- CedarCreek.2004.2006.dat[ , c(1, 1 + 1:nbElt, 1 + nbElt + 1)]

## ----plot_opt_tree, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.tree manages the plot of primary and secondary trees")----
apriori <- c("F","C3","L","C4","F","C3","C3","L","F","L","F","C4","L","C3","C4","C4")
res <- fclust(dat.2004, nbElt, opt.mod = "byelt", opt.mean = "gmean")

fclust_plot(res, main = "BioDiv2", opt.tree = list("prd"))
fclust_plot(res, main = "BioDiv2", opt.tree = list("prd", cols = apriori))

## ----plot_opt_perf_stats, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.perf=list('stats_xx') manages the plot of statistics")----
fclust_plot(res, main = "BioDiv2", opt.perf = list("stats_II"))

## ----plot_opt_perf_prd, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.perf=list('prd') on left, =list('pub') on right")----
fclust_plot(res, main = "BioDiv2", nbcl = 3, opt.perf = list("prd", "aov", pvalue = 0.01, "pub"))

## ----plot_opt_perf_missing, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.perf=list('missing') on left, =list('calprd') on right")----
fclust_plot(res, main = "BioDiv2", nbcl = 4, opt.perf = list("missing"))
fclust_plot(res, main = "BioDiv2", nbcl = 4, opt.perf = list("calprd", "aov", pvalue = 0.01))

## ----plot_opt_perf_seq, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.perf=list('prd', 'seq') plots performances predicted from1 to nbOpt component clusters")----
fclust_plot(res, main = "BioDiv2", opt.perf = list("prd", "seq", "aov", pvalue = 0.01))

## ----plot_opt_perf_ass, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.perf=list('prd', 'ass') adds the names of assemblages on plot")----
fclust_plot(res, main = "BioDiv2", opt.perf = list("prd", "aov", pvalue = 0.01))
fclust_plot(res, main = "BioDiv2", opt.perf = list("prd", "ass"))

## ----plot_opt_motif, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.motif=list('obs', 'seq') plots performances sorted by assembly motif")----
fclust_plot(res, main = "BioDiv2", opt.motif = list("prd", "hor", "seq", "aov", pvalue = 0.01))

## ----plot_opt_comp, results="hide", out.width="48%", fig.align="default", fig.show="hold", fig.cap=c("Opt.comp=list('tree', 'leg') plots set of performances that contain a given component, sorted as the hierarchical tree of components")----
fclust_plot(res, main = "BioDiv2", opt.comp = list("tree", "leg", "aov", pvalue = 0.01))

