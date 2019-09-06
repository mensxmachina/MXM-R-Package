\name{Read big data or a big.matrix object}
\alias{read.big.data}
\title{
Forward Backward Early Dropping selection regression
}

\description{
Forward Backward Early Dropping selection regression.
}

\usage{
read.big.data(path, sep = ",", header = FALSE)
}

\arguments{
\item{path}{
The path where the big.matrix object is.  
}
\item{sep}{
A field delimiter, for example " " (tab separated) or "," (comma separated). 
}
\item{header}{
If there are column names, then this should be TRUE.
}
}

\details{
The data (matrix) which will be read and compressed into a big.matrix object must be of type "numeric". We tested it and it works with "integer" as well. But, in general, bear in mind that only matrices will be read. We have not tested with data.frame for example. Woever, in the help page of "bigmemory" this is mentioned: Any non-numeric entry will be ignored and replaced with NA, so reading something that traditionally would be a data.frame won't cause an error. A warning is issued. In all cases, the object size is alwasy 696 bytes!
}

\value{
A big.matrix object.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>. For more information see the "bigmemory" package.
}

\seealso{
\code{ \link{big.gomp} }
}

\examples{
\dontrun{
dataset <- matrix( runif(100 * 50, 1, 100), ncol = 50 )
read.csv(data, "dataset.csv", sep = ",")
a <- read.big.data("dataset.csv")
}
}

