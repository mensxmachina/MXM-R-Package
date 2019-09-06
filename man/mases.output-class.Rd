\name{mases.output-class}
\Rdversion{1.1}
\docType{class}
\alias{mases.output-class}
\alias{mases.output}
\alias{plot,mases.output-method}
\alias{mases.output-method}
\alias{plot,mases.output,ANY-method}

\title{Class \code{"mases.output"}}
\description{
Meta analytic SES output object class.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mases.output", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{selectedVars}:}{Object of class \code{"numeric"}}
    \item{\code{selectedVarsOrder}:}{Object of class \code{"numeric"}}
    \item{\code{queues}:}{Object of class \code{"list"}}
    \item{\code{signatures}:}{Object of class \code{"matrix"}}
    \item{\code{hashObject}:}{Object of class \code{"list"}}
    \item{\code{pvalues}:}{Object of class \code{"numeric"}}
    \item{\code{stats}:}{Object of class \code{"numeric"}}
    \item{\code{univ}:}{Object of class \code{"list"}}
    \item{\code{max_k}:}{Object of class \code{"numeric"}}
    \item{\code{threshold}:}{Object of class \code{"numeric"}}
    \item{\code{runtime}:}{Object of class \code{"proc_time"}}
    \item{\code{test}:}{Object of class \code{"character"}}
}
}
\section{Methods}{
  \describe{
    \item{plot}{\code{plot(x = "mases.output", mode = "all")}: Generic function for plotting the generated pvalues of the mases.output object. Argument mode = "all" for plotting all the pvalues or mode="partial" for partial plotting the first 500 pvalues }
	 }
}
\author{
Giorgos Athineou <athineou@csd.uoc.gr>
}

\seealso{
\code{\link{ma.ses}, \link{ma.mmpc}}
}
\examples{
showClass("mases.output")
}
\keyword{meta analytic ses output}
