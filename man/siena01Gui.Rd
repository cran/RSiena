\name{siena01Gui}
\alias{siena01Gui}
\title{User interface }
\description{
Gui to allow entry of the data for a siena model fit.
}
\usage{
siena01Gui(getDocumentation=FALSE)
}
\arguments{
  \item{getDocumentation}{Flag to allow documentation of internal
  functions, not for use by users.}
}
\details{

  This function provides a graphical user interface for fitting Siena
  models.
  Note that this is less flexible and has fewer possibilities than
  creating data objects from matrices and
  vectors using \code{\link{sienaDependent}}, \code{\link{coCovar}} etc.,
  and finally \code{\link{sienaDataCreate}}.

  It can be run from within an R session, but is also called
  directly by \command{sienascript} (Linux or Mac) (the Windows
  interface via \command{siena.exe} has now been removed).  It allows
  entry of details of the data files required, either using the gui or
  by loading a session file. The \command{Apply} button causes a call to
  \code{\link{sienaDataCreateFromSession}} followed by a display of the
  \code{\link{sienaModelOptions}} screen.

  The required format for the column entries is described on the help
  page  for\cr
  \code{\link{sienaDataCreateFromSession}}, as this function can also be
  called directly.

  The entries for the table can be loaded from a file by using the
  buttons \code{Load new session from file} or \code{Continue session
	from file}. The former will create a new report file and produce a
  descriptive report. The latter will use an existing report
  file and omit the descriptive report.

 Alternatively, use \code{Add} and \code{Remove} buttons to enter the
 file names, and adjust the other columns to describe your data (see
 help page for \code{\link{sienaDataCreateFromSession}}).

 The \code{Save to file} button will save the entries in the table to a
 session file.

 The \code{Clear} button will empty the table.

 The \code{Apply} button will prompt to save the session, then create
 the data objects and display the \code{\link{sienaModelOptions}}
 screen.

 Exit by using the menu \code{File/Quit} or by closing the Window.
}
\value{
  None, although various objects made will still be in the directory if you
  are using this within an R session.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{    \code{\link{sienaDataCreateFromSession}},
    \code{\link{sienaDataCreate}}.}
\examples{
\dontrun{siena01Gui()}
}
\keyword{ misc}
