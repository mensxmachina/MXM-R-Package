\name{mammpc.output-class}
\Rdversion{1.1}
\docType{class}
\alias{mammpc.output-class}
\alias{mammpc.output}
\alias{plot,mammpc.output-method}
\alias{mammpc.output-method}
\alias{plot,mammpc.output,ANY-method}

\title{Class \code{"mammpc.output"}}
\description{
mammpc. output object class.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mammpc.output", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{selectedVars}:}{Object of class \code{"numeric"}}
    \item{\code{selectedVarsOrder}:}{Object of class \code{"numeric"}}
    \item{\code{hashObject}:}{Object of class \code{"list"}}
    \item{\code{pvalues}:}{Object of class \code{"numeric"}}
    \item{\code{stats}:}{Object of class \code{"numeric"}}
    \item{\code{univ}:}{Object of class \code{"list"}}
    \item{\code{max_k}:}{Object of class \code{"numeric"}}
    \item{\code{threshold}:}{Object of class \code{"numeric"}}
    \item{\code{test}:}{Object of class \code{"character"}}
    \item{\code{runtime}:}{Object of class \code{"proc_time"}}
}
}
\section{Methods}{
  \describe{
    \item{plot}{\code{plot(x = "mammpc.output", mode = "all")}: Generic function for plotting the generated pvalues of the MMPCoutput object. Argument mode = "all" for plotting all the pvalues or mode="partial" for partial plotting the first 500 pvalues }
	 }
}
\author{
Giorgos Athineou <athineou@csd.uoc.gr>
}

\seealso{
\code{\link{ma.mmpc}, \link{ma.ses}}
}
\examples{
showClass("mammpc.output")
}
\keyword{meta-analytic MMPC output}
