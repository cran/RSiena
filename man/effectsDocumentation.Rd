\name{effectsDocumentation}
\alias{effectsDocumentation}
\title{
Function to create a table of documentation of effect names, short names etc.
}
\description{
  Produces a table of the shortnames and other information for effects,
  either in html or latex.
}
\usage{
effectsDocumentation(type = "html", display = type == "html",
filename = "effects")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
	type of output required. Valid options are "html" or "latex".
}
\item{display}{
  Boolean: should the output be displayed after creation. Only
  applicable to html output.
}
\item{filename}{
  Basename of the files to be written.
}
}
\details{
  The \code{allEffects} object is reformatted into a table, either latex
  or html.
}
\value{
  Nothing returned. Output files are created.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\examples{
\dontrun{effectsDocumentation()}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
