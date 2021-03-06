% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fclust.R
\name{plot_fcomp}
\alias{plot_fcomp}
\title{Plot mean performances
of assemblages containing a given component}
\usage{
plot_fcomp(fres, nbcl = 0, main = "Title", opt.comp = NULL )
}
\arguments{
\item{fres}{an object generated by the function \code{\link{fclust}}.}

\item{nbcl}{an integer.
The integer indicates the number of component clusters
to take into account.
It can be lower than or equals to
the optimum number \code{fres$nbOpt} of component clusters.}

\item{main}{a string, that is used
as the first, reference part of the title of each graph.}

\item{opt.comp}{a list, that can include
\code{opt.comp = list("tree", "perf", "hor", "ver", cols,
pvalue, "zoom", window, "all")}.
This option list manages the plot as boxplot
of observed mean performances
of assemblages that contain a given component,
horizontally or vertically,
components sorted by increasing or decreasing mean values,
or components sorted like the clustering tree.
The item order in list is any. \cr

\itemize{
 \item{\code{"tree", "perf":}}{  plot the observed mean performances
   of assemblages that contain a given component as boxplots.
   Each set of assemblages that contains a given component
   is named by the contained component.
   The coloured squares are the mean performances of assemblage sets.
   Size (number of observed assemblages) of assemblage sets
   is indicated on the left of boxplots.
   The red dashed line is the mean performance of assemblage sets.
   If \code{"aov"} is checked, groups significantly different
   (at a p-value < \code{pvalue}) are indicated by differents letters
   on the right of boxplots. \cr

   If \code{"tree":} is checked, mean performances
   of assemblages that contain a given component
   are sorted like the clustering tree.
   If \code{"perf"} is checked, mean performances
   of assemblages that contain a given component
   are sorted by increasing mean performances.
   }

\item{\code{"hor":}}{  plot boxplots as horizontal boxes:
  x-axis corresponds to assemblage performances,
  and y-axis corresponds to assemblage sets.
  It \code{"hor"} is not checked,
  boxplots are plotted as vertical boxes:
  x-axis corresponds to assemblage sets,
  and y-axis corresponds to assemblage performances.
  Option "ver" can also be used: "ver" = !"hor".
  }

\item{\code{cols:}}{  is a vector of integers, of same length
  as the number of components. This option specifies
  the colour of each component.
  The components labelled by the same integer
  have the same colour. If \code{cols} is not specified,
  the components that belong to a same cluster
  \emph{a posteriori} determined have the same colour.
  This option is useful when an \emph{a priori} clustering is known,
  to identify the components \emph{a priori} clustered
  into the \emph{a posteriori} clustering.
  }

\item{\code{pvalue = value:}}{  a probability used as threshold
  in the variance analysis. Then \code{pvalue} must be
  higher than \code{0} and lower than \code{1}.
  \code{pvalue} must be informed when \code{"aov"} is checked.
  Groups significantly different
  (at a p-value < \code{pvalue}) are then indicated by differents letters
  on the right of boxplots.
  }

\item{\code{"all":}}{  plot all possible graphs.
  This option is equivalent to
  \code{opt.motif = list("tree", "aov", pvalue = 0.001,
  "zoom", window = 20)}.
  }
}}
}
\value{
Nothing. It is a procedure.
}
\description{
The function plots, as vertical or horizontal boxplots,
 the mean performances
 of assemblages containing a given component.
}
\details{
None.
}
\seealso{
\code{\link{plot_ftrees}} plot primary and secondary trees
 resulting from functional clustering \cr
\code{\link{plot_fperf}} plot observed, modelled and predicted performances
 resulting from functional clustering \cr
\code{\link{plot_fass}} plot performances of some given assemblages \cr
\code{\link{plot_fmotif}} plot as boxplot mean performances
of assemblages sorted by assembly motifs \cr
\code{\link{plot_fcomp}} plot as boxplot mean performances
of assemblages containing a given component \cr
\code{\link{fclust_plot}} plot all possible outputs
of functional clustering
}
\keyword{internal}
