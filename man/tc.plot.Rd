\name{Plot of longitudinal data}
\alias{tc.plot}
\title{
Plot of longitudinal data
}

\description{
Plot of longitudinal data.
}

\usage{
tc.plot(target, tp, id, type = "l", ylab = "Values", xlab = "Time points",
       col = 2, lwd = 1, lty = 2, pch = 1, main = "Spaghetti plot")
}

\arguments{
\item{target}{
A numerical vector with the longitudinal data. 
}
\item{tp}{
The time points. It can either be a vector with length either equal to the number of time points or equal to the legnth of the target.
}
\item{id}{
A numerical vector specifying the subjects. It can either be a vector with length either equal to the number of subjects or equal to the legnth of the target.
}
\item{type}{
This is a graphical parameter. You can have lines "l" everywhere or lines with points at each time point "p". 
}
\item{ylab}{
This is a graphical parameter. The label on the y axis.
}
\item{xlab}{
This is a graphical parameter. The label on the x axis.
}
\item{col}{
This is a graphical parameter. The color of the lines. 
}
\item{lwd}{
This is a graphical parameter. The thickness of the lines.
}
\item{lty}{
This is a graphical parameter. The type of line, e.g. dashed, dotted, etc.
}
\item{pch}{
This is a graphical parameter. If the type is "b", then you can specify if you want different signs, for example circles, crosses, diamonds etc. 
}
\item{main}{
This is a graphical parameter. The title of the graph. 
}
}

\details{
The data must be longitudinal (the same subject measured multiple times at different time points) and for one variable only. For the graphical parameters see \code{\link{plot}} or \code{\link{par}}.
}

\value{
A plot with the longitudinal data over time.
}

\references{
Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). 
}


\author{
Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Michail Tsagris <mtsagris@uoc.gr>
}

\seealso{
\code{\link{testIndGLMMReg}, \link{SES.glmm} }
}

\examples{
y <- rnorm(70)
Days <- rep(1:7, each = 10)
id <- rep(1:10, 7)
tc.plot(y, Days, id)
tc.plot(y, Days, id, type = "b")
}

\keyword{ Longitudinal data }
\keyword{ Graphical visualisation }