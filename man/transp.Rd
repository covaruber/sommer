\name{transp}
\alias{transp}
\title{Creating color with transparency}

\description{
This function takes a color and returns the same with a certain alpha grade transparency.
}
\usage{
transp(col, alpha=0.5)
}

\arguments{
  \item{col}{Color to be used for transparency}

\item{alpha}{Grade of transparency desired}
}
\details{
No major details.
}
\value{
If arguments are correctly specified the function returns:
\describe{
\item{$res}{ A new color with certain grade of transparency}
}
}

\references{
Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
transp("red", alpha=0.5)
}
\seealso{The core functions of the package \code{\link{mmer}} }
