#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                             FEXPORT.R
#
# Exported functions                                                       ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include
#' fclust.R
#' plot_fclust.R
#' test_fclust.R
#' boot_fclust.R
#'
NULL


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# fclust functions                                                         ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Read a functional clustering for one or several performances
#'
#' @description Read the files resulting from a functional clustering
#'   and saved in text format in 6 different files
#'   by using function \code{fclust_write()}.
#'
#' @usage
#' fclust_read(filename = "")
#'
#' @inheritParams fclust_write
#'
#' @details The results are saved in 5 different files.
#'
#'  \itemize{
#'  \item \code{"filename.options.csv"}{: contains \code{nbElt}, \code{nbAss},
#'    \code{nbOpt}, \code{"opt.method"}, \code{"opt.mean"},
#'    \code{"opt.model"}.}
#'  \item \code{"filename.inputs.csv"}{: contains \code{fobs} and
#'   \code{names(fobs)}, \code{xpr} and \code{names(xpr)}.}
#'  \item \code{"filename.trees.csv"}{: contains the hierarchical tree
#'    \code{tree$aff} and \code{tree$cor}.}
#'  \item \code{"filename.matrices.csv"}{: contains the matrices
#'    \code{mCal}, \code{mPrd}, \code{mMotifs}, \code{tCal}, \code{tPrd},
#'     and \code{tNbcl}.}
#'  \item \code{"filename.stats.csv"}{: contains both statistical matrices
#'  \code{mStats} and \code{tStats}.}
#'  }
#'
#'  If only a file does not exist or is corrupted, the function is stopped.
#'
#'
#' @examples
#'
#' # save "res" in the files "myRecord.*" then read them again.
#'
#' res <- CedarCreek.2004.res
#' filename <- tempfile(pattern = "myRecord", tmpdir = tempdir())
#'
#' fclust_write(res, filename)
#' res <- fclust_read(filename)
#'
#' all.equal(res, CedarCreek.2004.res)
#'
#'
#' @seealso
#' \code{\link{fclust}}: make a functional clustering,\cr
#' \code{\link{fclust_plot}}: plot the results of a functional clustering,\cr
#' \code{\link{fclust_write}}: save the results of a functional clustering,\cr
#' \code{\link{fclust_read}}: read the results of a functional clustering.\cr
#'
#' @return The result of the functional clustering
#'   recorded in the files \code{"filename.*.csv"}.
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fclust_read <- function(filename = "") {

  if (nchar(filename) == 0) stop("'filename' cannot be an empty string")

  tmp <- c("nbElt", "nbAss", "nbXpr",
           "opt.method", "opt.mean", "opt.model",
           "opt.jack", "jack",
           "opt.na", "opt.repeat", "affectElt",
           "fobs", "mOccur", "xpr",
           "tree.I", "tree.II", "nbOpt",
           "mCal", "mPrd", "mMotifs", "mStats",
           "tCal", "tPrd", "tNbcl", "tStats")

  fres <- vector(mode = "list", length = length(tmp))
  names(fres) <- tmp


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Read options
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  tmp <- read_foptions(filename)

  fres$nbElt      <- tmp$nbElt
  fres$nbAss      <- tmp$nbAss
  fres$nbXpr      <- tmp$nbXpr
  fres$opt.method <- tmp$opt.method
  fres$opt.mean   <- tmp$opt.mean
  fres$opt.model  <- tmp$opt.model
  fres$opt.jack   <- tmp$opt.jack
  fres$jack       <- tmp$jack
  fres$opt.na     <- tmp$opt.na
  fres$opt.repeat <- tmp$opt.repeat
  fres$affectElt  <- tmp$affectElt


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Read inputs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  tmp <- read_finputs(filename, fres$nbElt)

  fres$fobs   <- tmp$fobs
  fres$mOccur <- tmp$mOccur
  fres$xpr    <- tmp$xpr


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Read trees
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  tmp <- read_ftrees(filename, fres$nbElt)

  fres$tree.I  <- tmp$tree.I
  fres$tree.II <- tmp$tree.II
  fres$nbOpt   <- tmp$nbOpt


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Read matrices and stats
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  tmp1 <- read_fmatrices(filename, fres$nbElt)
  tmp2 <- read_fstats(filename, fres$nbElt)

  fres$mCal    <- tmp1$mCal
  fres$mPrd    <- tmp1$mPrd
  fres$mMotifs <- tmp1$mMotifs
  fres$mStats  <- tmp2$mStats

  fres$tCal    <- tmp1$tCal
  fres$tPrd    <- tmp1$tPrd
  fres$tNbcl   <- tmp1$tNbcl
  fres$tStats  <- tmp2$tStats


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Outputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  attributes(fres)$class <- "list"

  return(fres)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Record a functional clustering for one or several performances
#'
#' @description Write the results of a functional clustering
#'   in text format in 6 different files.
#'
#' @usage
#' fclust_write(fres = NULL, filename = "")
#'
#' @param fres a list containing predictions of assembly performances
#'       and statistics computed by using a species clustering tree.
#'       The list is generated
#'       by the function \code{validate_ftree},
#'       also called by the function \code{fclust}.
#'
#' @param filename a string, used as radical for the 6 file names.
#'
#' @details The results are splitted in 5 different files.
#'
#'  \itemize{
#'  \item \code{"filename.options.csv"}{: contains \code{nbElt}, \code{nbAss},
#'    \code{nbXpr}, \code{"opt.method"}, \code{"opt.mean"},
#'    \code{"opt.model"}, \code{"opt.jack"},\code{"jack"},
#'    \code{"opt.na"}, \code{"opt.repeat"} and \code{"affectElt"}.}
#'  \item \code{"filename.inputs.csv"}{: contains \code{fobs}, \code{xpr}
#'     and \code{mOccur}.}
#'  \item \code{"filename.trees.csv"}{: contains
#'    the optimum number of functional clusters \code{nbOpt},
#'    and the hierarchical trees \code{tree.I} and \code{tree.II}.}
#'  \item \code{"filename.matrices.csv"}{: contains the matrices
#'    \code{mCal}, \code{mPrd}, \code{mMotifs}, \code{tCal}, \code{tPrd},
#'     and \code{tNbcl}.}
#'  \item \code{"filename.stats.csv"}{: contains both statistical matrices
#'  \code{mStats} and \code{tStats}.}
#'  }
#'
#' @examples
#'
#' # save "res" in the files "myRecord.*".
#'
#' res <-  CedarCreek.2004.res
#' filename <- tempfile(pattern = "myRecord", tmpdir = tempdir())
#'
#' fclust_write(res, filename)
#'
#'
#' @importFrom utils write.table
#'
#' @seealso
#' \code{\link{fclust}}: make a functional clustering,\cr
#' \code{\link{fclust_plot}}: plot the results of a functional clustering,\cr
#' \code{\link{fclust_write}}: save the results of a functional clustering,\cr
#' \code{\link{fclust_read}}: read the results of a functional clustering.\cr
#'
#' @return Nothing. It is a procedure.
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fclust_write <- function(fres = NULL, filename = "") {

  if (nchar(filename) == 0) stop("'filename' cannot be an empty string")
  if (length(fres) == 0) stop("'fres' cannot be NULL")

  # Inputs
  # nbElt, nbAss, nbXpr, "opt.method", "opt.mean", "opt.model": 1 word
  # "opt.na", "opt.repeat", affectElt : nbElt words


  tmp <- cbind(fres$nbElt, fres$nbAss, fres$nbXpr,
               fres$opt.method, fres$opt.mean, fres$opt.model,
               fres$opt.jack, fres$jack[1], fres$jack[2],
               fres$opt.na, fres$opt.repeat,
               matrix(fres$affectElt, nrow = 1))
  colnames(tmp) <- c("nbElt", "nbAss", "nbXpr",
                     "opt.method", "opt.mean", "opt.model",
                     "opt.jack", "jack.1", "jack.2",
                     "opt.na", "opt.repeat", names(fres$affectElt))
  utils::write.table(x = tmp,
                     file = paste(filename, "options", "csv", sep = "."),
                     append = FALSE, col.names = TRUE, row.names = FALSE,
                     sep = ",")

  # inputs
  # names(fobs), fobs, names(xpr), xpr, mOccur : nbAss words

  tmp <- cbind(names(fres$fobs), fres$fobs,
               names(fres$xpr), fres$xpr, fres$mOccur)
  colnames(tmp) <- c("names.fobs", "fobs", "names.xpr", "xpr",
                     colnames(fres$mOccur))
  utils::write.table(x = tmp,
                     file = paste(filename, "inputs", "csv", sep = "."),
                     append = FALSE, col.names = TRUE, row.names = FALSE,
                     sep = ",")


  # Trees
  # nbOpt, tree$aff, tree$cor : nbElt x nbElt words

  tmp <- cbind(fres$nbOpt, seq_len(fres$nbElt),
               fres$tree.I$aff, fres$tree.I$cor,
               fres$tree.II$aff, fres$tree.II$cor)
  colnames(tmp) <- c("nbOpt", "nbClu",
                     paste0(colnames(fres$tree.I$aff), ".I"), "R2.I",
                     paste0(colnames(fres$tree.II$aff), ".II"), "R2.II")
  utils::write.table(x = tmp,
                     file = paste(filename, "trees", "csv", sep = "."),
                     append = FALSE, col.names = TRUE, row.names = FALSE,
                     sep = ",")


  # Matrices
  # "mCal", "mPrd", "mMotifs", "tCal", "tPrd", "tNbcl" : nbClu x nbAss

  tmp <- cbind(c(rep("mCal", fres$nbElt),    rep("mPrd", fres$nbElt),
                 rep("mMotifs", fres$nbElt), rep("tCal", fres$nbElt),
                 rep("tPrd", fres$nbElt),    rep("tNbcl", fres$nbElt)),
               c(rep(seq_len(fres$nbElt), 6)),
               rbind(fres$mCal, fres$mPrd, fres$mMotifs,
                     fres$tCal, fres$tPrd, fres$tNbcl))
  colnames(tmp) <- c("mat", "nbClu", names(fres$fobs))
  utils::write.table(x = tmp,
                     file = paste(filename, "matrices", "csv", sep = "."),
                     append = FALSE, col.names = TRUE, row.names = FALSE,
                     sep = ",")


  # Stats
  # "mStats", "tStats" : nbClu x nbStat (=7)

  tmp <- cbind(c(rep("mStats", fres$nbElt), rep("tStats", fres$nbElt)),
               rep(seq_len(fres$nbElt), 2),
               rbind(fres$mStats, fres$tStats) )
  colnames(tmp)[c(1, 2)] <- c("mat", "nbClu")
  utils::write.table(x = tmp,
                     file = paste(filename, "stats", "csv", sep = "."),
                     append = FALSE, col.names = TRUE, row.names = FALSE,
                     sep = ",")
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot various graphs
#' of a functional clustering for one or several performances
#'
#' @description The function plots numerous useful graphs
#' for illustrating results and the ways by which they were obtained:
#' hierarchical trees of component clustering,
#' composition and mean performance of assembly motifs,
#' mean performance of assemblages containing a given components,
#' observed, simulated and predicted performances of assemblages
#' labelled by assembly motif,
#' performances of given assemblages...
#'
#' @usage
#' fclust_plot(fres, nbcl = 0, main = "",
#'             opt.tree  = NULL, opt.perf = NULL, opt.ass = NULL,
#'             opt.motif = NULL, opt.comp = NULL, opt.all = NULL )
#'
#' @param fres an object generated by the function \code{\link{fclust}}.
#'
#' @param nbcl an integer.
#' The integer indicates the number of component clusters
#' to take into account.
#' It can be lower than or equals to
#' the optimum number \code{fres$nbOpt} of component clusters.
#'
#' @param main a string, that is used
#' as the first, reference part of the title of each graph.
#'
#' @param opt.tree a list, that can include
#' \code{opt.tree = list("cal", "prd", cols, "zoom", window, "all")}.
#' This option list manages the plot of primary and secondary trees
#' of component clustering,
#' simplified or not, focussed on the main component clusters or not,
#' coloured by the user or not.
#' The item order in list is any. \cr
#'
#' \itemize{
#'  \item{\code{"cal"}}{  plots the primary tree of component clustering,
#'   from trunk until leaves.
#'   At trunk level, when all components are clustered
#'   into a large, trivial cluster,
#'   the coefficient of determination \code{R2} is low.
#'   At the leaves level, when each component is isolated in a singleton,
#'   the coefficient of determination is always equal to \code{1}.
#'   The primary tree is therefore necessarily over-fitted
#'   near the leaves level.
#'   The optimum number \code{fres$nbOpt} of component clusters
#'   is determined by the minimum \code{AICc}.
#'
#'   The blue dashed line indicates the level
#'   (optimum number \code{fres$nbOpt} of component clusters)
#'   where the tree must be optimally cut up.
#'   The red solid line indicates the value of tree efficiency \code{E}
#'   at the \code{nbcl}-level.
#'   The component clusters are named by lowercase letters,
#'   from left to right as \code{"a", "b", "c", ...}:
#'   the name and content of each component cluster
#'   is written on the following page.
#'   }
#'
#' \item{\code{"prd"}}{  plots the validated,
#'   secondary tree of component clustering,
#'   from trunk until validated leaves.
#'   Secondary tree is the primary tree cut
#'   at the level of the optimal number \code{nbOpt} of component clusters.
#'   \code{nbOpt} is determined
#'   by the first lowest value of \code{AIC} along the primary tree.
#'
#'   The red solid line indicates the value of tree efficiency \code{E}.
#'   \code{R2} and \code{E} are stored in \code{fres$tStats}.
#'   The component clusters are named by lowercase letters,
#'   from left to right as \code{"a", "b", "c", ...}:
#'   the name and content of each component cluster
#'   is written on the following page.
#'   }
#'
#' \item{\code{cols}}{  is a vector of colours, characters or integers,
#'   of same length as the number of components. This option specifies
#'   the colour of each component.
#'   The components labelled by the same integer
#'   have the same colour. If \code{cols} is not specified,
#'   the components that belong to a same cluster
#'   \emph{a posteriori} determined have the same colour.
#'   This option is useful when an \emph{a priori} clustering is known,
#'   to identify the components \emph{a priori} clustered
#'   into the \emph{a posteriori} clustering.
#'   }
#'
#' \item{\code{"zoom"}}{ if \code{"cal"} or \code{"prd"} is checked,
#'   this option allows
#'   to only plot the first, significant component clusters.
#'   The cluster on the far right (the cluster named by the last letter)
#'   is most often a large cluster, that includes many components
#'   of which the effects of assemblage performance are not significant.
#'   When the number of components is large, the tree is dense
#'   and the names of components are confusing.
#'   The option is useful to focus on the left, more signficant,
#'   part of the primary or secondary tree.
#'   If \code{"zoom"} is checked, \code{window} must be informed.
#'   If not, the function stops with an error message.
#'   Note that the large cluster, that includes many components,
#'   is always represented by at least one component.
#'   }
#'
#' \item{\code{window}}{  an integer, that
#'   specifies the number of components to plot.
#'   \code{window} must be informed when \code{"zoom"} is checked.
#'   If \code{window} is higher than the number of components, it is ignored.
#'   If \code{window} is lower than the number of significant components,
#'   it is ajusted in such a way that the large cluster,
#'   that includes many components,
#'   is at least represented by one component.
#'   }
#'
#' \item{\code{"all"}}{  plots all possible graphs.
#'   This option is equivalent to
#'   \code{opt.tree = list("cal", "prd", "zoom", window = 20)}.
#'   If the number of components is lower than 20,
#'   the option is equivalent to \code{opt.tree = list("cal", "prd")}.
#'   }
#' }
#'
#'
#' @param opt.perf a list, that can include
#' \code{opt.perf = list("stats_I", "stats_II",
#' "cal", "prd", "missing", "pub", "calprd",
#' "seq", "ass", "aov", pvalue, "all")}.
#' This option list manages the plot
#' of observed, modelled and predicted performances of assemblages,
#' and associated statistics. It also allows to plot performances of
#' some given, identified assemblages.
#' The item order in list is any. \cr
#'
#' \itemize{
#'  \item{\code{"stats_I", "stats_II":}}{  plot the statistics associated to
#'    fit of primary tree that best accounts
#'    for observed performances (\code{"stats_I"}),
#'    and of secondary tree that best predicts
#'    observed performances of assemblages (\code{"stats_II"}).
#'    Four graphs are plotted:
#'    1. coefficient of determination \code{R2}
#'    and efficiency \code{E} of models of component clustering
#'    (on y-axis) \emph{versus} the number of component clusters (on x-axis);
#'    2. the ratio of assemblage perfomances
#'    that cannot be predicted by cross-validation ("predicting ratio");
#'    3. and 4. the Akaike Information Criterion,
#'    corrected \code{AICc} or not \code{AIC} for small datasets.
#'    The green solid line indicates the first minimum of \code{AIC}
#'    that corresponds to the optimum number \code{nbOpt}
#'    of component clusters to consider.
#'    }
#'
#'  \item{\code{"cal", "prd":}}{  plot modelled performances
#'    \emph{versus} observed performances (\code{"cal"},
#'    or modelled and predicted by cross-validation performances
#'    \emph{versus} observed performances (\code{"prd"},
#'    for a number of component clusters increasing from \code{1}
#'    until the number of component clusters where efficiency \code{E}
#'    is maximum.
#'    Different symbols correspond to different assembly motifs.
#'    The prediction error induced by cross-validation is indicated
#'    by a short vertical line. \cr
#'
#'    The blue dashed lines are mean performances.
#'    The red solid line is 1:1 bissector line.
#'    The number of component clusters is indicated on graph left top.
#'    Predicting ratio and coefficient of determination \code{R2}
#'    of the clustering
#'    are indicated on graph right bottom.
#'    If \code{"prd"} is checked, efficiency \code{E}
#'    and \code{E/R2} ratio are added.
#'    If \code{"aov"} is checked, groups significantly different
#'    (at a p-value < \code{pvalue}) are indicated by differents letters
#'    on the right of graph.
#'    }
#'
#'  \item{\code{"missing":}}{  the option \code{"prd"} plot
#'    modelled and predicted by cross-validation performances
#'    \emph{versus} observed performances,
#'    using different symbols for different assembly motifs.
#'    The option \code{"missing"} plot the same data,
#'    but in using different symbols according to the clustering model
#'    used for predicting the performances of assemblages.
#'    This option allows to identify assemblages
#'    of which the performance cannot be predicted
#'    using the clustering model of the current level.
#'    The assemblages are plotted and named
#'    using the symbol corresponding to the level
#'    of the used clustering model. \cr
#'
#'    The blue dashed lines are mean performances.
#'    The red solid line is 1:1 bissector line.
#'    The number of component clusters is indicated on graph left top.
#'    Predicting ratio and coefficient of determination of the clustering
#'    are indicated on graph right bottom.
#'    If \code{"aov"} is checked, groups significantly different
#'    (at a p-value < \code{pvalue}) are indicated by differents letters
#'    on the right of graph.
#'    }
#'
#'  \item{\code{"pub":}}{  the option \code{"prd"} plot
#'    modelled and predicted by cross-validation performances
#'    \emph{versus} observed performances,
#'    using different symbols for different assembly motifs.
#'    The option \code{"pub"} plot the same data,
#'    but in using only one symbol.
#'    This option is useful for publication. \cr
#'
#'    The blue dashed lines are mean performances.
#'    The red solid line is 1:1 bissector line.
#'    The number of component clusters is indicated on graph left top.
#'    Predicting ratio and coefficient of determination of the clustering
#'    are indicated on graph right bottom.
#'    If \code{"aov"} is checked, groups significantly different
#'    (at a p-value < \code{pvalue}) are indicated by differents letters
#'    on the right of graph.
#'    }
#'
#'  \item{\code{"calprd":}}{  plot performances predicted by cross-validation
#'    \emph{versus} performances predicted by clustering model
#'    ("modelled performances"). This option is useful
#'    to identify which assembly motifs become difficult
#'    to predict by cross-validation. \cr
#'
#'    The blue dashed lines are mean performances.
#'    The red solid line is 1:1 bissector line.
#'    The number of component clusters is indicated on graph left top.
#'    Predicting ratio and coefficient of determination of the clustering
#'    are indicated on graph right bottom.
#'    If \code{"aov"} is checked, groups significantly different
#'    (at a p-value < \code{pvalue}) are indicated by differents letters
#'    on the right of graph. The letters are located
#'    at \code{mean(Fprd[motif == label])}.
#'    }
#'
#' \item{\code{"seq":}}{  plot performances of assembly motifs,
#'   from \code{1} to \code{nbMax} number of component clusters.
#'   Remember that number \code{m} of assembly motifs increases
#'   with the number \code{nbcl} of component clusters
#'   (\code{m = 2^nbcl - 1}). When the optimal number
#'   of component clusters is large,
#'   this option is useful to determine
#'   a number of component clusters lower
#'   than the optimal number of component clusters.
#'   Assembly motifs are named as the combinations of component clusters
#'   (see "opt.tree").
#'   }
#'
#'  \item{\code{"ass"}}{  plot the name of each assemblage
#'    close to its performance. This option can be used with
#'    the options \code{"cal"}, \code{"prd"}, \code{"pub"}
#'    and \code{"calprd"}. It must be used only
#'    if the number of assemblages is small.
#'    If the number of assemblages is large,
#'    the following option \code{"opt.ass"} is more convenient.
#'    }
#'
#'  \item{\code{"aov":}}{  does a variance analysis of assemblage performances
#'    by assembly motifs, and plot the result on the right of graphs.
#'    Different letters correspond to
#'    groups significantly different at a p-value < \code{pvalue}.
#'    If \code{"aov"} is checked, \code{pvalue} must be informed.
#'    If not, \code{pvalue = 0.001}.
#'   }
#'
#' \item{\code{pvalue:}}{  a probability used as threshold
#'   in the variance analysis. Then \code{pvalue} must be
#'   higher than \code{0} and lower than \code{1}.
#'   \code{pvalue} must be informed when \code{"aov"} is checked.
#'   Groups significantly different
#'   (at a p-value < \code{pvalue}) are then indicated by differents letters
#'   on the right of boxplots.
#'   }
#'
#' \item{\code{"all":}}{  plot all possible graphs.
#'   This option is equivalent to
#'   \code{opt.pref = list("cal", "prd", "pub", "calprd",
#'   "aov", pvalue = 0.001)}.
#'   }
#' }
#'
#'
#' @param opt.ass a list, that include
#' \code{opt.ass = list(sample, who)}.
#' This option plot modelled and predicted by cross-validation performances
#'  \emph{versus} observed performances,
#'  for a small sample of assemblages randomly drawn (\code{sample}),
#'  or for given, identified assemblages chosen by the user (\code{who}).
#'  The item order in list is any. \cr
#'
#' \itemize{
#'  \item{\code{sample:}}{  an integer. This integer
#'    specifies the number of assemblages
#'    to randomly drawn in the assemblage set,
#'    the plot as the option \code{opt.perf = list("prd")}.
#'    All chosen assemblages are plotted on a same graph.
#'    }
#'
#'  \item{\code{who:}}{  a list of assemblage names.
#'    The list contains the names of assemblages to plot.
#'    Each assemblage is plotted on a specific graph.
#'    This option is useful when ssemblage performances
#'    are observed over several experiments.
#'    }
#'  }
#'
#'
#' @param opt.motif a list, that can include
#' \code{opt.motif = list("obs", "cal", "prd", cols, "hor", "ver", "seq",
#' pvalue, "all")}.
#' This option list manages the plot of mean performances
#' of assembly motifs as boxplots,
#' observed, modelled or predicted by cross-validation,
#' horizontally or vertically,
#' sorted by increasing or decreasing mean values,
#' from \code{1} to \code{nbOpt} clusters of components.
#' The item order in list is any. \cr
#'
#' \itemize{
#'  \item{\code{"obs", "cal", "prd":}}{  plot the observed,
#'    modelled or predicted by cross-validation mean performances
#'    of assembly motifs as boxplots.
#'    Assembly motifs are named as the combinations of component clusters
#'    (see "opt.tree").
#'    The coloured squares are the mean performances of assembly motifs.
#'    Size (number of observed assemblages) of assembly motifs
#'    is indicated on the left of boxplots.
#'    The red dashed line is the mean performance of assembly motifs.
#'    If \code{"aov"} is checked, groups significantly different
#'    (at a p-value < \code{pvalue}) are indicated by differents letters
#'    on the right of boxplots.
#'    }
#'
#' \item{\code{"hor":}}{  plot boxplots as horizontal boxes:
#'   x-axis corresponds to assemblage performances,
#'   and y-axis corresponds to assembly motifs.
#'   It \code{"hor"} is not checked,
#'   boxplots are plotted as vertical boxes:
#'   x-axis corresponds to assembly motifs,
#'   and y-axis corresponds to assemblage performances.
#'   Option "ver" can also be used: "ver" = !"hor".
#'   }
#'
#' \item{\code{"seq":}}{  plot mean performances of assembly motifs,
#'   from \code{2} to \code{nbOpt} number of component clusters.
#'   Remember that number \code{m} of assembly motifs increases
#'   with the number \code{nbcl} of component clusters
#'   (\code{m = 2^nbcl - 1}). When the optimal number
#'   of component clusters is large,
#'   this option is useful to determine
#'   a number of component clusters lower
#'   than the optimal number of component clusters.
#'   Assembly motifs are named as the combinations of component clusters
#'   (see "opt.tree").
#'   }
#'
#' \item{\code{pvalue = value:}}{  a probability used as threshold
#'   in the variance analysis. Then \code{pvalue} must be
#'   higher than \code{0} and lower than \code{1}.
#'   \code{pvalue} must be informed when \code{"aov"} is checked.
#'   Groups significantly different
#'   (at a p-value < \code{pvalue}) are then indicated by differents letters
#'   on the right of boxplots.
#'   }
#'
#' \item{\code{"all":}}{  plot all possible graphs.
#'   This option is equivalent to
#'   \code{opt.motif = list("obs", "cal", "prd", "seq",
#'   "aov", pvalue = 0.001)}.
#'   }¶
#' }
#'
#'
#' @param opt.comp a list, that can include
#' \code{opt.comp = list("tree", "perf", "hor", "ver", cols,
#' pvalue, "zoom", window, "all")}.
#' This option list manages the plot as boxplot
#' of observed mean performances
#' of assemblages that contain a given component,
#' horizontally or vertically,
#' components sorted by increasing or decreasing mean values,
#' or components sorted like the clustering tree.
#' The item order in list is any. \cr
#'
#' \itemize{
#'  \item{\code{"tree", "perf":}}{  plot the observed mean performances
#'    of assemblages that contain a given component as boxplots.
#'    Each set of assemblages that contains a given component
#'    is named by the contained component.
#'    The coloured squares are the mean performances of assemblage sets.
#'    Size (number of observed assemblages) of assemblage sets
#'    is indicated on the left of boxplots.
#'    The red dashed line is the mean performance of assemblage sets.
#'    If \code{"aov"} is checked, groups significantly different
#'    (at a p-value < \code{pvalue}) are indicated by differents letters
#'    on the right of boxplots. \cr
#'
#'    If \code{"tree":} is checked, mean performances
#'    of assemblages that contain a given component
#'    are sorted like the clustering tree.
#'    If \code{"perf"} is checked, mean performances
#'    of assemblages that contain a given component
#'    are sorted by increasing mean performances.
#'    }
#'
#' \item{\code{"hor":}}{  plot boxplots as horizontal boxes:
#'   x-axis corresponds to assemblage performances,
#'   and y-axis corresponds to assemblage sets.
#'   It \code{"hor"} is not checked,
#'   boxplots are plotted as vertical boxes:
#'   x-axis corresponds to assemblage sets,
#'   and y-axis corresponds to assemblage performances.
#'   Option "ver" can also be used: "ver" = !"hor".
#'   }
#'
#' \item{\code{cols:}}{  is a vector of integers, of same length
#'   as the number of components. This option specifies
#'   the colour of each component.
#'   The components labelled by the same integer
#'   have the same colour. If \code{cols} is not specified,
#'   the components that belong to a same cluster
#'   \emph{a posteriori} determined have the same colour.
#'   This option is useful when an \emph{a priori} clustering is known,
#'   to identify the components \emph{a priori} clustered
#'   into the \emph{a posteriori} clustering.
#'   }
#'
#' \item{\code{pvalue = value:}}{  a probability used as threshold
#'   in the variance analysis. Then \code{pvalue} must be
#'   higher than \code{0} and lower than \code{1}.
#'   \code{pvalue} must be informed when \code{"aov"} is checked.
#'   Groups significantly different
#'   (at a p-value < \code{pvalue}) are then indicated by differents letters
#'   on the right of boxplots.
#'   }
#'
#' \item{\code{"all":}}{  plot all possible graphs.
#'   This option is equivalent to
#'   \code{opt.motif = list("tree", "aov", pvalue = 0.001,
#'   "zoom", window = 20)}.
#'   }
#' }
#'
#'
#' @param opt.all This option is equivalent to
#'         \code{opt.tree = "all", opt.comp = "all", opt.motif = "all",
#'         opt.perf = "all"}. This option is convenient to overview
#'         the different options of the function \code{fclust_plot}.
#'
#' @details If all the options are \code{NULL},
#' that is \code{opt.tree = NULL, opt.perf = NULL, opt.ass = NULL,
#' opt.motif = NULL, opt.comp = NULL, opt.all = NULL},
#' the function plot the main results, that are:
#' the secondary tree (\code{opt.tree = "prd"}),
#' assembly motifs as horizontal boxplots
#' (\code{opt.motif = list("obs", "hor")})), and
#' modelled and predicted by cross-validation mean performances
#' \emph{versus} observed performances (\code{opt.perf = "prd"}).
#'
#' @return Nothing. It is a procedure.
#'
#' @examples
#'
#' res <- CedarCreek.2004.res
#'
#' # plot the hierarchical tree of functionally redundant components
#' fclust_plot(res, main = "BioDiv2 2004", opt.tree = "prd")
#'
#' # plot AIC and AICc versus the number of clusters of components
#' layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
#' fclust_plot(res, main = "BioDiv2 2004", opt.perf = "stats_II")
#' layout(1)
#'
#' # plot the performances modelled and predicted versus observed performances
#' fclust_plot(res, main = "BioDiv2 2004", opt.perf = "prd")
#'
#' # plot the performances sorted by assembly motifs
#' layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
#' fclust_plot(res, main = "BioDiv2 2004",
#'             opt.motif = c("obs", "prd", "hor"))
#' layout(1)
#'
#'
#' @seealso
#' \code{\link{fclust}}: make a functional clustering,\cr
#' \code{\link{fclust_plot}}: plot the results of a functional clustering,\cr
#' \code{\link{fclust_write}}: save the results of a functional clustering,\cr
#' \code{\link{fclust_read}}: read the results of a functional clustering.\cr
#' \cr
#' \code{\link{plot_ftrees}} plot primary and secondary trees
#'  resulting from a functional clusterin, \cr
#' \code{\link{plot_fperf}} plot observed, modelled and predicted performances
#'  resulting from a functional clustering, \cr
#' \code{\link{plot_fass}} plot performances of some given assemblages, \cr
#' \code{\link{plot_fmotif}} plot as boxplot mean performances
#' of assemblages sorted by assembly motifs, \cr
#' \code{\link{plot_fcomp}} plot as boxplot mean performances
#' of assemblages containing a given component, \cr
#' \code{\link{fclust_plot}} plot all possible outputs
#' of a functional clustering.
#'
#' @references
#'
#' Jaillard, B., Richon, C., Deleporte, P., Loreau, M. and Violle, C. (2018)
#' \emph{An a posteriori species clustering
#' for quantifying the effects of species
#' interactions on ecosystem functioning}.
#' Methods in Ecology and Evolution, 9:704-715.
#' \url{https://doi.org/10.1111/2041-210X.12920}. \cr
#'
#' Jaillard, B., Deleporte, P., Loreau, M. and Violle, C. (2018)
#' \emph{A combinatorial analysis using observational data
#' identifies species that govern ecosystem functioning}.
#' PLoS ONE 13(8): e0201135.
#' \url{https://doi.org/10.1371/journal.pone.0201135}. \cr
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fclust_plot <- function(fres,
                        nbcl      = 0,
                        main      = "Title",
                        opt.tree  = NULL,
                        opt.perf  = NULL,
                        opt.ass   = NULL,
                        opt.motif = NULL,
                        opt.comp  = NULL,
                        opt.all   = NULL  )  {

  # plot by default
  if ( is.null(opt.tree) & is.null(opt.perf) &
       is.null(opt.ass)  & is.null(opt.motif) &
       is.null(opt.comp) & is.null(opt.all)  ) {
    opt.tree  <- list("prd", window = Inf, "leg")
    opt.perf  <- list("prd", "pub")
    opt.motif <- list("obs", "hor", "leg")
  }

  # plot all
  if (!is.null(opt.all))
    if ((opt.all == "all") | (opt.all == TRUE)) {
      opt.tree  <- list("all")
      opt.perf  <- list("all")
      opt.motif <- list("all")
      opt.comp  <- list("all")
    }

  # plots choisen by the user
  if (!is.null(opt.tree))  plot_ftrees(fres, nbcl, main, opt.tree)
  if (!is.null(opt.perf))  plot_fperf( fres, nbcl, main, opt.perf)
  if (!is.null(opt.ass))   plot_fass(  fres, nbcl, main, opt.ass)
  if (!is.null(opt.motif)) plot_fmotif(fres, nbcl, main, opt.motif)
  if (!is.null(opt.comp))  plot_fcomp( fres, nbcl, main, opt.comp)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Build a functional clustering for one or more performances
#'
#' @description Fit a primary tree of component clustering
#' to observed assemblage performances,
#' then prune the primary tree for its predicting ability and its parcimony,
#' finally retain a validated secondary tree
#' and the corresponding predictions, statistics and other informations.
#'
#' @usage
#' fclust(dat, nbElt,
#'        weight     = rep(1, dim(dat)[2] - nbElt - 1),
#'        opt.na     = FALSE,
#'        opt.repeat = FALSE,
#'        opt.method = "divisive",
#'        affectElt  = rep(1, nbElt),
#'        opt.mean   = "amean",
#'        opt.model  = "byelt",
#'        opt.jack   = FALSE,   jack = c(3,4) )
#'
#' @param dat a data.frame or matrix that brings together:
#' a vector of assemblage identity,
#' a matrix of occurrence of components within the system,
#' one or more vectors of observed performances.
#' Consequently, the data.frame or matrix dimensions are:
#' \code{dim(dat)[1]=} the number of observed assemblages,
#' \code{* dim(dat)[2]=} 1 + number of system components +
#' number of observed performances.
#' On a first line (colnames): assemblage identity,
#' a list of components identified by their names,
#' a list of performances identified by their names.
#' On following lines (a line by assemblage),
#' name of the assemblage (read as character),
#' a sequence of 0 (absence) and 1 (presence of component
#' within each assemblage)
#' (this is the matrix of occurrence of components within the system),
#' a sequence of numeric values for informed each observed performances
#' (this is the set of observed performances).
#'
#' @param nbElt an integer, that specifies the number of components
#' belonging to interactive system.
#' \code{nbElt} is used to know the dimension of matrix of occurrence.
#'
#' @param weight a vector of numerics,
#'  that specifies the weight of each performance.
#'  By default, each performance is equally weighted.
#'  If \code{weight} is informed, it must have the same length
#'  as the number of observed performances.
#'
#' @param opt.na a logical.
#' The records for each assemblage can have \code{NA}
#' in matrix of occurrence or in observed assemblage performances.
#' If \code{opt.na = FALSE} (by default), an error is returned.
#' If \code{opt.na = TRUE}, the records with \code{NA} are ignored.
#'
#' @param opt.repeat a logical.
#' in any case, the function looks for
#' different assemblages with identical elemental composition.
#' Messages indicate these identical assemblages.
#' If \code{opt.repeat = FALSE} (by default),
#' their performances are averaged.
#' If \code{opt.repeat = TRUE}, nothing is done,
#' and the data are processed as they are.
#'
#' @param opt.method a string that specifies the method to use.
#' \code{opt.method = c("divisive", "agglomerative", "apriori")}.
#' The three methods generate hierarchical trees.
#' Each tree is complete, running from a unique trunk
#' to as many leaves as components. \cr
#'
#' If \code{opt.method = "divisive"}, the components are clustered
#' by using a divisive method,
#' from the trivial cluster where all components are together,
#' towards the clustering where each component is a cluster.
#' This method gives the best result for several reasons,
#' exposed in detail in joined vignettes (see "The options of fclust").
#'
#' If \code{opt.method = "agglomerative"}, the components are clustered
#' by using an agglomerative method,
#' from the trivial clustering where each component is a cluster,
#' towards the cluster where all components are brought together
#' If all possible assemblages are not observed
#' (that is generally he case in practice),
#' the first clustering of few components can have no effect
#' on convergence criterion, indicing a non-optimum result.
#'
#' If \code{opt.method = "apriori"}, the user knows and gives
#' an "a priori" partitioning of the system components he is studying.
#' The partition is arbitrary, in any number of clusters of components,
#' but it must be specified (see following option \code{affectElt}).
#' The tree is then built:
#' \emph{(i)} by using \code{opt.method =  "divisive"}
#' from the defined component clustering towards as many leaves as components;
#' \emph{(ii)} by using \code{opt.method =  "agglomerative"}
#' from the  component clustering towards the trunk of tree.
#'
#' @param affectElt a vector of characters or integers,
#' as long as the number of components \code{nbElt},
#' that indicates the labels of different functional clusters
#' to which each component belongs.
#' Each functional cluster is labelled as a character or an integer, and
#' each component must be identified by its name in \code{names(affectElt)}.
#' The number of functional clusters defined in \code{affectElt}
#' determines an \emph{a priori} level of component clustering
#' (\code{level <- length(unique(affectElt))}).\cr
#'
#' If \code{affectElt = NULL} (by default),
#' the option \code{opt.method} must be specified.
#' If \code{affectElt} is specified,
#' the option \code{opt.method} switchs to \code{apriori}.
#'
#' @param opt.mean a character, equals to \code{"amean"} or \code{"gmean"}.
#' If \code{opt.mean = "amean"},
#' means are computed using an arithmetic formula,
#' if \code{opt.mean = "gmean"},
#' mean are computed using a geometric formula.
#'
#' @param opt.model a character equals to \code{"bymot"} or \code{"byelt"}.
#' If \code{opt.model = "bymot"},
#' the  modelled performances are means
#' of performances of assemblages
#' that share a same assembly motif
#' by including all assemblages that belong to a same assembly motif. \cr
#'
#' If \code{opt.model = "byelt"},
#' the modelled performances are the average
#' of mean performances of assemblages
#' that share a same assembly motif
#' and that contain the same components
#' as the assemblage to predict.
#' This procedure corresponds to a linear model within each assembly motif
#' based on the component occurrence in each assemblage.
#' If no assemblage contains component belonging to assemblage to predict,
#' performance is the mean performance of all assemblages
#' as in \code{opt.model = "bymot"}.
#'
#' @param opt.jack a logical,
#' that switchs towards cross-validation method.
#'
#' If \code{opt.jack = FALSE} (by default), a Leave-One-Out method is used:
#' predicted performances are computed
#' as the mean of performances of assemblages
#' that share a same assembly motif,
#' experiment by experiment,
#' except the only assemblage to predict. \cr
#'
#' If \code{opt.jack = TRUE}, a jackknife method is used:
#' the set of assemblages belonging to a same assembly motif is divided
#' into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#' Predicted performances of each subset of \code{jack[1]} assemblages
#' are computed, experiment by experiment,
#' by using the other (\code{jack[2] - 1}) subsets of assemblages.
#' If the total number of assemblages belonging
#' to the assembly motif is lower than \code{jack[1]*jack[2]},
#' predictions are computed by Leave-One-Out method.
#'
#' @param jack an integer vector of length \code{2}.
#' The vector specifies the parameters for jackknife method.
#' The first integer \code{jack[1]} specifies the size of subset,
#' the second integer \code{jack[2]} specifies the number of subsets.
#'
#' @details see Vignette "The options of fclust".
#'
#' @return Return a list containing the primary tree of component clustering,
#' predictions of assembly performances
#' and statistics computed by using the primary and secondary trees
#' of component clustering.
#'
#'  Recall of inputs:
#'  \itemize{
#'  \item \code{nbElt, nbAss, nbXpr}:
#'  the number of components that belong to the interactive system,
#'  the number of assemblages and the number of performances observed,
#'  respectively.
#'
#'  \item \code{opt.method, opt.mean, opt.model, opt.jack, jack, opt.na,
#'  opt.repeat, affectElt}: the options used
#'  for computing the resulting clustering trees,
#'  respectively.
#'
#'  \item \code{fobs, mOccur, xpr}:
#'  the vector or matrix of observed performances of assemblages,
#'  the binary matrix of occurrence of components, and
#'  the vector of weight of different performances,
#'  respectively.
#'  }
#'
#'  Primary and secondary, fitted and validated trees,
#'  of component clustering and associated statistics:
#'  \itemize{
#'  \item \code{tree.I, tree.II, nbOpt}:
#'  the primary tree of component clustering,
#'  the validated secondary tree of component clustering,
#'  and the optimum number of functional clusters,
#'  respectively.
#'  A tree is a list of a square-matrix of dimensions
#'  \code{nbLev * nbElt} (with \code{nbLev = nbElt}),
#'  and of a vector of coefficient of determination (of length \code{nbLev}).
#'
#'  \item \code{mCal, mPrd, tCal, tPrd}:
#'  the numeric matrix of modelled values,
#'  and of values predicted by cross-validation,
#'  using the primary tree (\code{mCal} and (\code{mPrd})
#'  or the secondary tree (\code{tCal} and (\code{tPrd}), respectively.
#'  All matrices have the same dimension \code{nbLev * nbAss}.
#'  \code{rownames} contains the number of component clusters,
#'  that is from \code{1} to \code{nbElt} clusters.
#'  \code{colnames} contains the names of assemblages.
#'
#'  \item \code{mMotifs, tNbcl}: the matrix
#'  of affectation of assemblages to different assembly motifs,
#'  coded as integers, and the matrices of the last tree levels
#'  used for predicting assemblage performances.
#'  All matrices have the same dimension \code{nbLev * nbAss}.
#'  \code{rownames} contains the number of component clusters,
#'  that is from \code{1} to \code{nbElt} clusters.
#'  \code{colnames} contains the names of assemblages.
#'
#'  \item \code{mStats, tStats}: the matrices of associated statistics.
#'  \code{rownames} contains the number of component clusters,
#'  that is from {1} to {nbElt} clusters.
#'  \code{colnames = c("missing", "R2cal", "R2prd", "AIC", "AICc")}.
#'  }
#'
#' @examples
#'
#' # Enable the comments
#' oldOption <- getOption("verbose")
#' if (!oldOption) options(verbose = TRUE)
#'
#' nbElt <- 16   # number of components
#' # index = Identity, Occurrence of components, a Performance
#' index <- c(1, 1 + 1:nbElt, 1 + nbElt + 1)
#' dat.2004 <- CedarCreek.2004.2006.dat[ , index]
#' res <- fclust(dat.2004, nbElt)
#' names(res)
#' res$tree.II
#'
#' options(verbose = oldOption)
#'
#'
#' @seealso
#' \code{\link{fclust}}: build a functional clustering,\cr
#' \code{\link{fclust_plot}}: plot the results of a functional clustering,\cr
#' \code{\link{fclust_write}}: save the results of a functional clustering,\cr
#' \code{\link{fclust_read}}: read the results of a functional clustering.\cr
#'
#' @references
#'
#' Jaillard, B., Richon, C., Deleporte, P., Loreau, M. and Violle, C. (2018)
#' \emph{An a posteriori species clustering
#' for quantifying the effects of species
#' interactions on ecosystem functioning}.
#' Methods in Ecology and Evolution, 9:704-715.
#' \url{https://doi.org/10.1111/2041-210X.12920}. \cr
#'
#' Jaillard, B., Deleporte, P., Loreau, M. and Violle, C. (2018)
#' \emph{A combinatorial analysis using observational data
#' identifies species that govern ecosystem functioning}.
#' PLoS ONE 13(8): e0201135.
#' \url{https://doi.org/10.1371/journal.pone.0201135}.
#'
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fclust <- function(dat, nbElt,
                   weight     = rep(1, dim(dat)[2] - nbElt - 1),
                   opt.na     = FALSE,
                   opt.repeat = FALSE,
                   opt.method = "divisive",
                   affectElt  = rep(1, nbElt),
                   opt.mean   = "amean",
                   opt.model  = "byelt",
                   opt.jack   = FALSE,   jack = c(3,4)   )  {

  dat <- format_fclust(dat, nbElt, weight, opt.na, opt.repeat)

  tree.I <- fit_ftree(dat$fobs, dat$mOccur, dat$xpr,
                      affectElt, opt.method, opt.mean, opt.model)

  res <- validate_ftree(tree.I,
                        dat$fobs, dat$mOccur, dat$xpr,
                        opt.method, opt.mean, opt.model,
                        opt.jack, jack)

  fres <- make_fclust(res, opt.na, opt.repeat, affectElt)

  return(fres)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  ftest functions                                                         ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Read the significance of different variables
#'  of a functional clustering
#'
#' @description Read a file of results obtained by
#' a test of significance of functional clustering.
#'
#' @usage
#' ftest_read(filename,
#'            opt.var = c("components", "assemblages", "performances") )
#'
#' @inheritParams ftest_write
#'
#' @details The functions
#'   \code{ftest_components}, \code{ftest_assemblages},
#'   \code{ftest_performances}, \code{fboot_assemblages}
#'   and \code{fboot_performances}.
#'   generate a list containing a matrix by clustering index
#'   ("Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'   "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'   "Sokal_Sneath1" and "Sokal_Sneath2" index).
#'   Only their dimensions change according the used functions.
#'   Consequently, a same function is used for recording and reading
#'   the results of both the test-functions.
#'
#' @examples
#'
#' # save "res" in the files "myRecord.*" then read it again.
#'
#' filename <- tempfile(pattern = "myRecord", tmpdir = tempdir())
#'
#' ftest_write(fres  = CedarCreek.2004.2006.res,
#'             rtest = CedarCreek.2004.2006.test.components,
#'             filename = filename,
#'             opt.var  = "components")
#'
#' rtest <- ftest_read(filename = filename, opt.var = "components")
#'
#' all.equal(rtest, CedarCreek.2004.2006.test.components)
#'
#'
#' @importFrom utils read.table
#'
#' @return a list of matrices,
#' each containing the results for a clustering index.
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ftest_read <- function(filename,
                    opt.var = c("components", "assemblages", "performances")) {

  # Matrices
  #  "R2", "E", "Czekanowski_Dice" "Folkes_Mallows"   "Jaccard"  "Kulczynski"
  #  "Precision"        "Rand"  "Recall"   "Rogers_Tanimoto"  "Russel_Rao"
  #  "Sokal_Sneath1"    "Sokal_Sneath2"    : nbClu x (nbElt or nbAss)

  all.var <- c("components", "assemblages", "performances")
  if (is.na(pmatch(opt.var, all.var)))
    stop("'opt.var' should be ", list_in_quote(all.var))
  opt.var <- match.arg(opt.var, all.var)

  file <- paste(filename, "test", opt.var, "csv", sep = ".")
  if (!file.exists(file)) stop("file '", file, "' does not exist.")

  tmp <- utils::read.table(file = file, header = TRUE, sep = ",",
                           check.names = FALSE, stringsAsFactors = FALSE)

  # build a list

  lnames <- unique(tmp[ , 1])
  index  <- c(1:2)

  rtest <- vector(mode = "list", length = length(lnames))
  names(rtest) <- lnames
  for (i in seq_along(lnames)) {
    rtest[[i]] <- as.matrix(tmp[tmp$mat == lnames[i], -index])
    rownames(rtest[[i]]) <- seq_len(dim(rtest[[i]])[1])
    #    colnames(rtest[[i]]) <- seq_len(dim(rtest[[i]])[2])
  }

  return(rtest)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Record the significance of different variables
#'  of a functional clustering
#'
#' @description Record in a file the results of
#' a test of significance of functional clustering.
#'
#' @usage
#' ftest_write(fres, rtest, filename,
#'             opt.var = c("components", "assemblages", "performances") )
#'
#' @inheritParams ftest
#'
#' @param rtest a list of matrices, each matrix containing
#'  the results for a given clustering index.
#'  The object \code{rtest} is generated by the function \code{ftest}.
#'
#' @param filename a string, used as radical for naming the file
#'   \code{"filename.components.csv"}, \code{"filemane.assemblages.csv"}
#'   or \code{"filemane.performances.csv"}
#'   according to the dimensions of matrices.
#'
#' @param opt.var a string, specifying the last part of the file-name.
#'  \code{opt} can only be equal to
#'  \code{"components"}, \code{"assemblages"} or \code{"performances"}.
#'
#' @details The functions \code{ftest},
#'   \code{ftest_components}, \code{ftest_assemblages} and
#'   \code{ftest_performances}
#'   generate a list containing a matrix by clustering index
#'   ("Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'   "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'   "Sokal_Sneath1" and "Sokal_Sneath2" index).
#'   Only their dimensions change according the used functions.
#'   Consequently, a same function is used for recording and reading
#'   the results of both the test-functions.
#'
#' @examples
#'
#' # save "rtest" in the file "myRecord.*".
#'
#' ftest_write(fres  = CedarCreek.2004.2006.res,
#'             rtest = CedarCreek.2004.2006.test.components,
#'             filename = tempfile(pattern = "myRecord", tmpdir = tempdir()),
#'             opt.var  = "components")
#'
#' @importFrom utils write.table
#'
#' @return Nothing. It is a procedure.
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ftest_write <- function(fres, rtest, filename,
                    opt.var = c("components", "assemblages", "performances")) {


  # Matrices
  #  "R2", "E", "Czekanowski_Dice" "Folkes_Mallows"   "Jaccard"  "Kulczynski"
  #  "Precision"        "Rand"  "Recall"   "Rogers_Tanimoto"  "Russel_Rao"
  #  "Sokal_Sneath1"    "Sokal_Sneath2"    : nbClu x (nbElt or nbAss)

  all.var <- c("components", "assemblages", "performances")
  if (is.na(pmatch(opt.var, all.var)))
    stop("'opt.var' should be ", list_in_quote(all.var))
  opt.var <- match.arg(opt.var, all.var)


  # all matrices are saved in a same file

  col1 <- NULL
  for (i in seq_along(names(rtest)))
    col1 <- c(col1, rep(names(rtest)[i], fres$nbElt))

  col2 <- rep(seq_len(fres$nbElt), length(names(rtest)))

  col3 <- NULL
  for (i in seq_along(names(rtest))) col3 <- rbind(col3, rtest[[i]])

  tmp <- cbind(col1, col2, col3)
  ttt <- switch(opt.var,
                components   = colnames(fres$mOccur),
                assemblages  = unique(rownames(fres$mOccur)),
                performances = unique(names(fres$xpr)) )
  colnames(tmp) <- c("mat", "nbClu", ttt)
  rownames(tmp) <- as.character(seq_len(length(col1)))

  utils::write.table(x = tmp,
                     file = paste(filename, "test", opt.var, "csv", sep = "."),
                     append = FALSE, col.names = TRUE, row.names = FALSE,
                     sep = ",")
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot the significance of different variables
#'  of a functional clustering
#'
#' @description Different plots are built according to the tested variable.
#'
#' @usage
#' ftest_plot(fres, rtest,
#'            main     = "Title",
#'            opt.var  = c("components", "assemblages", "performances"),
#'            opt.crit = "Jaccard",
#'            opt.comp = NULL, opt.ass = NULL, opt.perf = NULL)
#'
#' @inheritParams ftest
#'
#' @param rtest a list of matrices,
#'   each containing the results for a clustering index.
#'   \code{rtest} is an object generated by the function \code{ftest}.
#'
#' @param main a string, that is used
#' as the first, reference part of the title of each graph.
#'
#' @param opt.crit a list of strings,
#'  indicating the clustering indices to plot.
#'  The indices can be:
#'   "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'   "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'   "Sokal_Sneath1" or "Sokal_Sneath2".
#'   For more informations, see the notice of R-package \code{clusterCrit}.
#'
#' @param opt.comp a list, that can include
#' \code{opt.comp} = \code{list("all.together", "fgroups.together",
#' "comps.together",
#' "fgroups.byfg", "comps.byfg", "sorted.tree", "sorted.leg", "all")}.
#' This option list manages the plot
#' of results obtained using the function \code{ftest}
#' with \code{opt.var = "components"}.
#' The item order in list is any.
#'
#' \itemize{
#'  \item{\code{"all.together", "fgroups.together", "comps.together"}}
#'    {  plot \emph{(i)} the general mean index;
#'    \emph{(ii)} the mean indices for each functional group on a same plot;
#'    and \emph{(iii)} the mean index for each components
#'    on a same plot,
#'    when removing one after one each component from the dataset.
#'    This allows to evaluate the raw robustness of functional clustering
#'    to perturbation of dataset,
#'    and the weight of each cluster on functional clustering.
#'    }
#'
#' \item{\code{"fgroups.byfg", "comps.byfg"}}{  plot
#'    \emph{(i)} mean component clusters,
#'    functional group by functional group;
#'    \emph{(ii)} a graph by component, functional group by functional group;
#'    This allows to evaluate the weight of each component
#'     on functional clustering.
#'    }
#'
#' \item{\code{"sorted.tree", "sorted.leg"}}{  plot
#'   \emph{(i)} the hierarchical tree of components,
#'   with components decreasingly sorted according to their weight
#'   on functional clustering within each functional group;
#'   \emph{(ii)} the names of component decreasingly sorted
#'   according to their weight on functional clustering
#'   within each functional group.
#'   }
#'
#' \item{\code{"all"}}{  plot all possible graphs.
#'   This option is equivalent to
#'   \code{opt.comp} = \code{list("all.together", "fgroups.together",
#'   "comps.together", "fgroups.byfg", "comps.byfg",
#'   "sorted.tree", "sorted.leg")}.
#'   }
#' }
#'
#' @param opt.ass a list, that can include
#' \code{opt.ass} = \code{list("all.together", "motifs.together",
#' "assemblages.together",
#' "motifs.bymot", "assemblages.bymot", "sorted.leg", "all")}.
#' This option list manages the plot
#' of results obtained using the function \code{ftest}
#' with \code{opt.var = "assemblages"}.
#' The item order in list is any.
#'
#' \itemize{
#'  \item{\code{"all.together", "motifs.together", "assemblages.together"}}
#'    {  plot \emph{(i)} the general mean index;
#'    \emph{(ii)} the mean indices for each assembly motif on a same plot;
#'    and \emph{(iii)} the mean index for each assemblages
#'    on a same plot,
#'    when removing one after one each assemblage from the dataset.
#'    This allows to evaluate the raw robustness of functional clustering
#'    to perturbation of dataset,
#'    and the weight of each assemblage on functional clustering.
#'    }
#'
#' \item{\code{"motifs.bymot", "assemblages.bymot"}}{  plot
#'    \emph{(i)} mean assembly motifs,
#'    assembly motif by assembly motif;
#'    \emph{(ii)} a graph by removed assemblage,
#'    assembly motif by assembly motif;
#'    This allows to evaluate the weight of each assemblage
#'     on functional clustering.
#'    }
#'
#' \item{\code{"sorted.leg"}}{  plot
#'   the names of assemblages decreasingly sorted
#'   according to their weight on functional clustering.
#'   }
#'
#' \item{\code{"all"}}{  plot all possible graphs.
#'   This option is equivalent to
#'   \code{opt.ass} = \code{list("all.together", "motifs.together",
#'  "assemblages.together",
#'   "motifs.bymot", "assemblages.bymot", "sorted.leg")}.}
#'   }
#'
#' @param opt.perf a list, that can include
#' a list, that can include
#' \code{opt.comp} = \code{list("all.together", "performances.together",
#'  "sorted.leg")}.
#' This option list manages the plot
#' of results obtained using the function \code{ftest}
#' with \code{opt.var = "performances"}.
#' The item order in list is any.
#'
#' \itemize{
#'  \item{\code{"all.together", "performances.together"}}
#'    {  plot \emph{(i)} the general mean index;
#'    \emph{(ii)} the mean indices for each removed performance on a same plot,
#'    when removing one after one each performance from the dataset.
#'    This allows to evaluate the raw robustness of functional clustering
#'    to perturbation of dataset,
#'    and the weight of each performance on functional clustering.
#'    }
#'
#' \item{\code{"sorted.leg"}}{  plot
#'   the names of performances decreasingly sorted
#'   according to their weight on functional clustering.
#'   }
#'
#' \item{\code{"all"}}{  plot all possible graphs.
#'   This option is equivalent to
#'   \code{opt.comp} = \code{list("all.together", "performances.together",
#'  "sorted.leg")}.
#'   }
#' }
#'
#' @details The trees obtained by leaving out each element
#' are compared to the reference tree obtained
#' with all element of the variables
#' using different criteria of clustering:
#' "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'  For more informations, see the notice of R-package \code{clusterCrit}.
#'
#' @examples
#'
#' # Plot the hierachical tree of components
#' layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
#' fclust_plot(fres = CedarCreek.2004.2006.res, main  = "BioDIV2",
#'             opt.tree = "prd")
#'
#' # Plot the significance of each component within each components cluster
#' ftest_plot(fres  = CedarCreek.2004.2006.res,
#'            rtest = CedarCreek.2004.2006.test.components,
#'            main  = "BioDIV2",
#'            opt.var = c("components"), opt.crit = "Jaccard")
#' layout(1)
#'
#'
#' @references
#' Package "clusterCrit": Clustering Indices,
#'   by Bernard Desgraupes (University of Paris Ouest - Lab Modal'X)
#'
#' @return Nothing. It is a procedure.
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ftest_plot <- function(fres, rtest,
                       main     = "Title",
                      opt.var = c("components", "assemblages", "performances"),
                       opt.crit = "Jaccard",
                       opt.comp = NULL, opt.ass = NULL, opt.perf = NULL  ) {


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Check options
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  if (is.null(opt.var)) {

    if (is.null(opt.comp) & is.null(opt.ass) & is.null(opt.perf))
      stop("'opt.var' or other options should be non-null")

    if (!is.null(opt.comp)) {
      option <- "components"
    } else {
      if (!is.null(opt.ass)) {
        option <- "assemblages"
      } else {
        option <- "performances"
      }
    }

  } else {

    all    <- c("components", "assemblages", "performances")
    option <- check_foption(opt.var, all)
  }


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Switch to the function ftest_plot_components
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  if (option == "components") {

    if (is.null(opt.comp))
      opt.comp <- list("all.together", "sorted.tree", "sorted.leg")
    ftest_plot_components(fres, rtest, main, opt.crit, opt.comp)
  }


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Switch to the function ftest_plot_assemblages
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  if (option == "assemblages") {

    if (is.null(opt.ass)) opt.ass <- list("all.together", "sorted.leg")
    ftest_plot_assemblages(fres, rtest, main, opt.crit, opt.ass)
  }


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Switc to the function ftest_plot_performances
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  if (option == "performances") {

    if (is.null(opt.perf)) opt.perf <- list("all.together", "sorted.leg")
    ftest_plot_performances(fres, rtest, main, opt.crit, opt.perf)
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test the significance of different variables
#'  of a functional clustering
#'
#' @description The function allows to test the relative significance of
#'  each component, of each assemblage and of each performance
#'  on the result of the functional clustering.
#'  The method is based on removing one after the other
#'  each component, assemblage or performance,
#'  then evaluating the effect of these deletions
#'  on the functional clustering. Each new functional clustering
#'  is compared with the functional clustering obtained
#'  with the whole dataset. The process is time-consuming.
#'
#' @usage
#' ftest(fres,
#'       opt.var = c("components", "assemblages", "performances"),
#'       opt.nbMax = fres$nbOpt, opt.R2 = FALSE, opt.plot = FALSE )
#'
#' @param fres an object resulting from a functional clustering
#'  obtained with the whole dataset using the function \code{\link{fclust}}.
#'
#' @param opt.var a string, that indicates the variable to test.
#'   The option can be \code{"components"}, \code{"assemblages"}
#'   or \code{"performances"}.
#'
#' @param opt.nbMax a logical. If \code{opt.plot = TRUE},
#' at each test, the tree resulting from removing
#'  each component, assemblage or performance is plotted.
#'
#' @param opt.R2 a logical. If \code{opt.R2 = TRUE},
#' the primary tree is validated
#' and the vectors of coefficient of determination (\code{R^2})
#' and efficiency (\code{E}) are computed.
#'
#' @param opt.plot a logical. If \code{opt.plot = TRUE},
#' at each test, the tree resulting from removing
#' each component, assemblage or performance is plotted.
#'
#' @details None.
#'
#' @return a list of matrices, each matrix containing
#'  the results for a given clustering index.
#'
#' @examples
#'
#' # Enable the comments
#' oldOption <- getOption("verbose")
#' if (!oldOption) options(verbose = TRUE)
#' layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
#'
#' \dontshow{
#' # Test the significance of annual biomass production
#' ftest(fres = CedarCreek.2004.2006.res,
#'       opt.var = c("performance"), opt.nbMax = 1)
#' }
#'
#' \donttest{
#' # Test the significance of annual biomass production
#' test.perf <- ftest(fres = CedarCreek.2004.2006.res,
#'                    opt.var = c("performance"), opt.plot = TRUE)
#'
#' # Test the significance of each component within each component cluster
#' test.comp <- ftest(fres = CedarCreek.2004.res,
#'                    opt.var = c("components"), opt.plot = TRUE)
#' }
#'
#' layout(1)
#' options(verbose = oldOption)
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ftest <- function(fres,
                  opt.var = c("components", "assemblages", "performances"),
                  opt.nbMax = fres$nbOpt, opt.R2 = FALSE, opt.plot = FALSE ) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Check options
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  option <- check_foption(opt.var,
                          c("components", "assemblages", "performances"))

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Switch to appropriated function
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  rtest <-
    switch(option,
           components  = ftest_components(  fres, opt.nbMax, opt.R2, opt.plot),
           assemblages = ftest_assemblages( fres, opt.nbMax, opt.R2, opt.plot),
           performances = ftest_performances(fres, opt.nbMax, opt.R2, opt.plot)
    )

  return(rtest)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# fboot functions                                                          ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Read the robustness of a functional clustering
#'  evaluated by bootstrapping from 1 to (all-1) observations
#'
#' @description Read a file of results obtained by
#' a test of significance of functional clustering.
#'
#' @usage
#' fboot_read(filename,
#'            opt.var = c("assemblages", "performances"))
#'
#' @inheritParams fboot
#'
#' @details The function \code{fboot} ,
#'   generate a list of lists, each one containing a matrix by clustering index
#'   ("Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'   "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'   "Sokal_Sneath1" and "Sokal_Sneath2" index).
#'   Only their dimensions change according the used functions.
#'   Consequently, a same function is used for recording and reading
#'   the results of both the test-functions.
#'
#' @return a list of list of matrices,
#'  identical to this resulting from the function \code{fboot}.
#'
#' @examples
#'
#' # save "rtest" in the file "myRecord.*".
#'
#' filename <- tempfile(pattern = "myRecord", tmpdir = tempdir())
#'
#' fboot_write(fres = CedarCreek.2004.2006.res,
#'             lboot = CedarCreek.2004.2006.boot.performances,
#'             filename = filename,
#'             opt.var  = "performances")
#'
#' lboot <- fboot_read(filename = filename, opt.var = "performances")
#'
#' all.equal(lboot, CedarCreek.2004.2006.boot.performances)
#'
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fboot_read <- function(filename,
                       opt.var = c("assemblages", "performances")) {

  # evaluate the number of files

  i <- 1
  while (file.exists(
    paste(filename, "boot", opt.var, i, "csv", sep = "."))) i <- i + 1
  j <- i - 1


  # build a list and read all other matrices saved in a same file

  lboot <- vector(mode = "list", length = j)
  for (i in seq_len(j)) {
    file       <- paste(filename, "boot", opt.var, i, "csv", sep = ".")
    lboot[[i]] <- fboot_read_one_point(filename = file)
  }

  return(lboot)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Record the robustness of a functional clustering
#'  evaluated by bootstrapping from 1 to (all-1) observations
#'
#' @description Write a file of results obtained by
#' a test of significance of functional clustering.
#'
#' @usage
#' fboot_write(fres, lboot, filename,
#'             opt.var = c("assemblages", "performances"))
#'
#' @inheritParams fboot
#'
#' @param lboot a list of list of matrices,
#'  generated by the function \code{fboot}.
#'
#' @details The function \code{fboot} ,
#'   generate a list of lists, each one containing a matrix by clustering index
#'   ("Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'   "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'   "Sokal_Sneath1" and "Sokal_Sneath2" index).
#'   Only their dimensions change according the used functions.
#'   Consequently, a same function is used for recording and reading
#'   the results of both the test-functions.
#'
#' @return Nothing. It is a procedure.
#'
#' @examples
#'
#' # save "rtest" in the file "myRecord.*".
#'
#' fboot_write(fres = CedarCreek.2004.2006.res,
#'             lboot = CedarCreek.2004.2006.boot.performances,
#'             filename = tempfile(pattern = "myRecord", tmpdir = tempdir()),
#'             opt.var = "performances")
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fboot_write <- function(fres, lboot, filename,
                        opt.var = c("assemblages", "performances")) {


  all.var <- c("assemblages", "performances")
  if (is.na(pmatch(opt.var, all.var)))
    stop("'opt.var' should be ", list_in_quote(all.var))
  opt.var <- match.arg(opt.var, all.var)


  i <- 1
  for (i in seq_along(lboot)) {
    file <- paste(filename, "boot", opt.var, i, "csv", sep = ".")
    fboot_write_one_point(fres, lboot[[i]], filename = file)
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot the robustness of a functional clustering
#'  evaluated by bootstrapping from 1 to (all-1) observations
#'
#' @description Evaluate by bootstrapping
#' the robustness of a functional clustering
#' to perturbations of data. The perturbed data can be
#' the number of assemblages taken into account,
#' or the number of performances taken into account.
#'
#' @usage
#' fboot_plot(fres, lboot, main = "", opt.crit = "Jaccard",
#'            opt.var  = c("assemblages", "performances"))
#'
#'
#' @param fres an object resulting from a functional clustering
#'  obtained using the function \code{\link{fclust}}.
#'
#' @param lboot an object resulting from a functional clustering
#'  obtained using the function \code{\link{fclust}}.
#'
#' @param main a string, used as main title for all plots.
#'
#' @param opt.crit an object resulting from a functional clustering
#'  obtained using the function \code{\link{fclust}}.
#'
#' @param opt.var an object resulting from a functional clustering
#'  obtained using the function \code{\link{fclust}}.
#'
#' @details The trees obtained by bootstrapping of performances to omit
#'  are compared to the reference tree obtained with all components
#'  using different criteria :
#'  "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'  For more informations, see the notice of R-package \code{clusterCrit}.
#'
#' @references
#' Package "clusterCrit": Clustering Indices,
#'   by Bernard Desgraupes (University of Paris Ouest - Lab Modal'X)
#'
#' @importFrom stats median
#'
#' @importFrom graphics boxplot plot
#'
#' @return a list of lists, each containing a matrix by clustering index.
#'
#' @examples
#'
#' # Plot the significance of each component within each components cluster
#'
#' layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
#' fboot_plot(fres  = CedarCreek.2004.2006.res,
#'            lboot = CedarCreek.2004.2006.boot.assemblages,
#'            main  = "BioDIV2",
#'            opt.var = "assemblages", opt.crit = "Jaccard")
#' layout(1)
#'
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fboot_plot <- function(fres, lboot, main = "",
                       opt.crit = "Jaccard",
                       opt.var  = c("assemblages", "performances")) {


  opt.var <- check_foption(opt.var,
                           c("assemblages", "performances"))
  if (opt.var == "assemblages")  nbItem <- fres$nbAss
  if (opt.var == "performances") nbItem <- fres$nbXpr

  rboot   <- lboot[[1]]
  setCrit <- names(rboot)[names(rboot) %in% opt.crit]

  nrows   <- dim(rboot$Jaccard)[1]
  ncols   <- dim(rboot$Jaccard)[2]

  cols    <- extend_vector(fcolours(), nbItem)
  pchs    <- extend_vector(fsymbols(), nbItem)

  matPlot <- matrix(NA, nrow = ncols, ncol = (nbItem - 1))

  crit <- 1
  for (crit in seq_along(setCrit)) {

    main  <- paste0(setCrit[crit], " index")

    for (j in 1:fres$nbOpt) {

      for (i in 1:(nbItem - 1)) {
        mat <- matrix(as.numeric(unlist(lboot[[i]][setCrit[crit]])),
                       ncol = ncols, nrow = nrows, byrow = FALSE)
        matPlot[1:ncols, i] <- mat[j, ]
      }
      matPlot[matPlot == 0] <- NA

      boxplot(matPlot,
              main = paste(main, paste0("nbcl=", j), sep = " / "),
              ylim = c(0, 1), xlim = c(1, (nbItem - 1)),
              ylab = "Jaccard index",
              xlab = paste("number of", opt.var, "removed", sep = " "))
    }


    for (j in 1:fres$nbOpt) {

      for (i in 1:(nbItem - 1)) {
        mat <- matrix(as.numeric(unlist(lboot[[i]][setCrit[crit]])),
                      ncol = ncols, nrow = nrows, byrow = FALSE)
        matPlot[1:ncols, i] <- mat[j, ]
      }
      matPlot[matPlot == 0] <- NA

      plot(y = apply(matPlot, 2, median, na.rm = TRUE),
           ylab = "median Jaccard index",
           x = 1:(nbItem - 1),
           xlab = paste("number of", opt.var, "removed", sep = " "),
           main = paste(main, paste0("nbcl=", j), sep = " / "),
           ylim = c(0, 1), xlim = c(1, (nbItem - 1)),
           cex = 4, pch = fsymbols()[j], type = "b", bg = "white")
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Evaluate the robustness of a functional clustering
#'  by bootstrapping from 1 to (all-1) observations
#'
#' @description Evaluate by bootstrapping
#' the robustness of a functional clustering
#' to perturbations of data. The perturbed data can be
#' the number of assemblages taken into account,
#' or the number of performances taken into account.
#'
#' @usage
#' fboot(fres,
#'       opt.var = c("assemblages", "performances"), nbIter = 1,
#'       opt.nbMax = fres$nbOpt, opt.R2 = FALSE, opt.plot = FALSE,
#'       filename = "" )
#'
#' @param fres an object resulting from a functional clustering
#'  obtained using the function \code{\link{fclust}}.
#'
#' @param opt.var a string, that indicates the variable to test.
#' The option can be \code{"assemblages"} or \code{"performances"}.
#'
#' @param nbIter an integer, that indicates the number of random drawing to do.
#'
#' @param opt.nbMax a logical. If \code{opt.plot = TRUE},
#' the trees resulting from leaving out each performance is plotted.
#'
#' @param opt.R2 a logical. If \code{opt.R2 = TRUE},
#' the primary tree is validated
#' and the vectors of coefficient of determination (\code{R^2})
#' and efficiency (\code{E}) are computed.
#'
#' @param opt.plot a logical. If \code{opt.plot = TRUE},
#' the primary trees resulting from leaving out each performance are plotted.
#' If \code{opt.R2 = TRUE},
#' the secondary trees resulting from leaving out each performance are plotted.
#'
#' @param filename a string, used as radical for naming the file
#'   \code{"filename.components.csv"}, \code{"filemane.assemblages.csv"}
#'   or \code{"filemane.performances.csv"}
#'   according to the dimensions of matrices.
#'
#' @details The trees obtained by bootstrapping of performances to omit
#'  are compared to the reference tree obtained with all components
#'  using different criteria :
#'  "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'  For more informations, see the notice of R-package \code{clusterCrit}.
#'
#' @references
#' Package "clusterCrit": Clustering Indices,
#'   by Bernard Desgraupes (University of Paris Ouest - Lab Modal'X)
#'
#' @importFrom clusterCrit getCriteriaNames extCriteria
#'
#' @importFrom graphics text
#'
#' @return a list of lists, each containing a matrix by clustering index.
#'
#' @examples
#'
#' # Enable the comments
#' oldOption <- getOption("verbose")
#' if (!oldOption) options(verbose = TRUE)
#' layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
#'
#' \dontshow{
#' nbIter <- 1
#' fboot(fres = CedarCreek.2004.2006.res,
#'       opt.var   = "performances",
#'       nbIter    = nbIter,
#'       opt.nbMax = 1)
#' }
#'
#' \donttest{
#' test.boot <- fboot(fres = CedarCreek.2004.2006.res,
#'                    opt.var  = "performances",
#'                    nbIter   = 4,
#'                    opt.plot = TRUE)
#' }
#'
#' layout(1)
#' options(verbose = oldOption)
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fboot <- function(fres,
                  opt.var = c("assemblages", "performances"), nbIter = 1,
                  opt.nbMax = fres$nbOpt, opt.R2 = FALSE, opt.plot = FALSE,
                  filename = "") {


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Check options
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  opt.var <- check_foption(opt.var, c("assemblages", "performances"))

  if (nbIter < 1) stop("'nbIter' should be higher than 1")

  if (opt.var == "assemblages")  dimMat <- fres$nbAss - 1
  if (opt.var == "performances") dimMat <- fres$nbXpr - 1


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Computation
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  rboot <- vector(mode = "list", length = dimMat)
  for (i in seq_len(dimMat)) {
    rboot[[i]] <- fboot_one_point(fres, opt.var, nbIter,
                                  rm.number = i, seed = i,
                                  opt.plot = opt.plot, opt.nbMax = opt.nbMax)

    if (filename != "") {
      file <- paste(filename, "boot", opt.var, i, "csv", sep = ".")
      fboot_write_one_point(fres, rboot[[i]], filename = file)
    }
  }

  return(rboot)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
