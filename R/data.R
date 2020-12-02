
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Data of Cedar Creek experiment for the years 2004, 2005 and 2006
#'
#' @description Data from Cedar Creek Ecosystem Science Reserve.
#' The data used here as an example come from the \emph{e120} experiment,
#' more commonly known as \emph{BioDIV} experiment.
#' We have selected the "unsorted biomass" of the harvests
#' of years 2004, 2005 and 2006, and the plots planted without any oak.
#'
#' @source \url{https://www.cedarcreek.umn.edu/research/data/dataset?ple120}
#'
#' @format A data.frame. Each row corresponds to an \emph{species assemblage},
#'   that is here a plot.
#'   The different columns of the data.frame are:
#' \describe{
#'   \item{Plot:}{ numero of the plot. Each assemblage must have
#'         a specific elemental composition,
#'         different from this of other assemblages.
#'         This first column is the \emph{identity} of assemblages.}
#'   \item{Achmi, Agrsm, Amocan, Andge, Asctu, Elyca, Koecr, Lesca, Liaas,
#'         Luppe, Monfi, Panvi, Petpu, Poapr, Schsc, Sornu:}{ abbreviated
#'         names of species.
#'         The set of species makes up an ecosystem.
#'         Each species is informed by \code{1} if it belongs to the assemblage,
#'         and by \code{0} if it does not belong to the assemblage.
#'         The whole is the \emph{matrix of occurrence}
#'         of different \emph{components} that belong to, and thus make up,
#'         the \emph{interactive system} under consideration.}
#'   \item{y2004, y2005, y2006:}{ "unsorted biomass" of plots,
#'         \emph{i.e.} species assemblages, harvested in August 2004,
#'         July 2005 and June 2006. Each harvest is a
#'         \emph{collective, systemic performance}
#'         of species assemblages. The performances can be treated separately,
#'         or collectively.}
#' }
#'
#' @details The Cedar Creek Ecosystem Science Reserve work was supported
#'  by grants from the US National Science Foundation Long-Term Ecological
#'  Research Program (LTER) including DEB-0620652 and DEB-1234162.
#'  Further support was provided by the Cedar Creek Ecosystem Science
#'  Reserve and the University of Minnesota.
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"CedarCreek.2004.2006.dat"



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Functional clustering of species used in Cedar Creek experiment
#' for the year 2004
#'
#' @description Cedar Creek results obtained for the biomass production
#' in the year 2004 only.
#'
#' @format A hierarchical tree of species clustering:
#' \describe{
#' \item{tree$aff}{a square-matrix of dimensions species number x
#' species number}
#' \item{tree$cor}{a vector of coefficient of determination}
#' }
#'
#' @details None.
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"CedarCreek.2004.res"



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Functional clustering of species used in Cedar Creek experiment
#' for the years 2004, 2005 and 2006
#'
#' @description Cedar Creek results obtained for the biomass production
#' for the years 2004, 2005 and 2006. Each year is equally weighted.
#'
#' @format A hierarchical tree of species clustering:
#' \describe{
#' \item{tree$aff}{a square-matrix of dimensions species number x
#' species number}
#' \item{tree$cor}{a vector of coefficient of determination}
#' }
#'
#' @details None.
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"CedarCreek.2004.2006.res"


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test of significance of components (species)
#' used in Cedar Creek experiment for the years 2004, 2005 and 2006
#'
#' @description Each component is successively removed
#' from the data set, a functional clustering of remaining species is built,
#' and the effect of this species remove on remaining species clustering
#' is evaluated by comparing the resulting species clustering
#' with the species clustering obtained with all species.
#' The perturbation is measured by using different clustering indices.
#'
#' @format A list of matrices, each matrix containing
#'  the results for a given clustering index. The indices are
#'  "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'
#' @details None.
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"CedarCreek.2004.2006.test.components"



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test of significance of species assemblages (plots)
#' used in Cedar Creek experiment for the years 2004, 2005 and 2006
#'
#' @description Each species assemblage is successively removed
#' from the data set, a functional clustering of species is built,
#' and the effect of this assemblage remove on species clustering
#' is evaluated by comparing the resulting species clustering
#' with the species clustering obtained with all species assemblages.
#' The perturbation is measured by using different clustering indices.
#'
#' @format A list of matrices, each matrix containing
#'  the results for a given clustering index. The indices are
#'  "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'
#' @details None.
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"CedarCreek.2004.2006.test.assemblages"



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test of significance of different performances (yearly biomass)
#' used in Cedar Creek experiment for the years 2004, 2005 and 2006
#'
#' @description Each performance is successively removed
#' from the data set, a functional clustering of species is built,
#' and the effect of this performance remove on species clustering
#' is evaluated by comparing the resulting species clustering
#' with the species clustering obtained with all systemic performances.
#' The perturbation is measured by using different clustering indices.
#'
#' @format A list of matrices, each matrix containing
#'  the results for a given clustering index. The indices are
#'  "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'
#' @details None.
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"CedarCreek.2004.2006.test.performances"



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Evaluate by bootstrapping the robustness of species clustering
#' to the number of assemblages (plots) used in the cluster analysis.
#'
#' @description Species assemblages are randomly removed from the dataset,
#' from 1 to m-1 assemblages (with m = total number of assemblages),
#' a functional clustering of species is built,
#' and the effect of this assemblage remove on species clustering
#' is evaluated by comparing the resulting species clustering
#' with the species clustering obtained with all observed species assemblages.
#' The perturbation is measured by using different clustering indices.
#'
#' @format A list of matrices, each matrix containing
#'  the results for a given clustering index. The indices are
#'  "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'
#' @details None.
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"CedarCreek.2004.2006.boot.assemblages"



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Evaluate by bootstrapping the robustness of species clustering
#' to the number of performances (annual biomass) used in the cluster analysis.
#'
#' @description Performances are randomly removed from the dataset,
#' from 1 to n-1 performances (with n = total number of performances),
#' a functional clustering of species is built,
#' and the effect of this performance remove on species clustering
#' is evaluated by comparing the resulting species clustering
#' with the species clustering obtained with all systemic performances.
#' The perturbation is measured by using different clustering indices.
#'
#' Here, the total number of performances is 3 only,
#' then the interest of this process is very limited.
#'
#' @format A list of matrices, each matrix containing
#'  the results for a given clustering index. The indices are
#'  "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'
#' @details None.
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"CedarCreek.2004.2006.boot.performances"



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' End of file                                                             ####
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
