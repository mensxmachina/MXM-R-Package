\name{MMPCoutput-class}
\Rdversion{1.1}
\docType{class}
\alias{MMPCoutput-class}
\alias{MMPCoutput}
\alias{plot,MMPCoutput-method}
\alias{MMPCoutput-method}
\alias{plot,MMPCoutput,ANY-method}


\title{Class \code{"MMPCoutput"}}
\description{
MMPC output object class.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MMPCoutput", ...)}.
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
}
}
\section{Methods}{
  \describe{
    \item{plot}{\code{plot(x = "MMPCoutput", mode = "all")}: Generic function for plotting the generated pvalues of the MMPCoutput object. Argument mode = "all" for plotting all the pvalues or mode="partial" for partial plotting the first 500 pvalues }
	 }
}
\author{
Giorgos Athineou <athineou@csd.uoc.gr>
}

\seealso{
\code{\link{MMPC}, \link{SES}}
}
\examples{
showClass("MMPCoutput")
}
\keyword{MMPC output}
