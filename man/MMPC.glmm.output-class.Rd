\name{MMPC.glmm.output-class}
\Rdversion{1.1}
\docType{class}
\alias{MMPC.glmm.output-class}
\alias{MMPC.glmm.output}
\alias{plot,MMPC.glmm.output-method}
\alias{MMPC.glmm.output-method}
\alias{plot,MMPC.glmm.output,ANY-method}

\title{Class \code{"MMPC.glmm.output"}}

\description{
MMPC.glmm output object class.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MMPC.glmm.output", ...)}.
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
    \item{\code{n.tests}:}{Object of class \code{"numeric"}}
    \item{\code{runtime}:}{Object of class \code{"proc_time"}}
    \item{\code{test}:}{Object of class \code{"character"}}
    \item{\code{slope}:}{Object of class \code{"logical"}}  }
}

\section{Methods}{
  \describe{
    \item{plot}{\code{plot(x = "MMPC.glmm.output", mode = "all")}: Generic function for plotting the generated pvalues of the MMPC.glmm.output object. Argument mode = "all" for plotting all the pvalues or mode="partial" for partial plotting the first 500 pvalues }
	 }
}

\author{
Giorgos Athineou <athineou@csd.uoc.gr>
}

\seealso{
\code{\link{MMPC.glmm}, \link{SES.glmm}}
}

\examples{
showClass("MMPC.glmm.output")
}

\keyword{MMPC.glmm output}
