\name{SESoutput-class}
\Rdversion{1.1}
\docType{class}
\alias{SESoutput-class}
\alias{SESoutput}
\alias{plot,SESoutput-method}
\alias{SESoutput-method}
\alias{plot,SESoutput,ANY-method}


\title{Class \code{"SESoutput"}}
\description{
SES output object class.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SESoutput", ...)}.
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
    \item{\code{n.tests}:}{Object of class \code{"numeric"}}
    \item{\code{runtime}:}{Object of class \code{"proc_time"}}
    \item{\code{test}:}{Object of class \code{"character"}}
}
}
\section{Methods}{
  \describe{
    \item{plot}{\code{plot(x = "SESoutput", mode = "all")}: Generic function for plotting the generated pvalues of the SESoutput object. Argument mode = "all" for plotting all the pvalues or mode="partial" for partial plotting the first 500 pvalues }
	 }
}
\author{
Giorgos Athineou <athineou@csd.uoc.gr>
}

\seealso{
\code{\link{SES}}
}
\examples{
showClass("SESoutput")
}
\keyword{SES output}
