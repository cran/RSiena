\name{sienaDataCreateFromSession}
\alias{sienaDataCreateFromSession}
\title{Creates a Siena data object from a Siena session file}
\description{
Reads in a Siena session from file or \code{siena01Gui()}
and creates a Siena data or group object.}
\usage{
sienaDataCreateFromSession(filename = NULL, session = NULL,
    modelName = "Siena", edited = NULL, files = NULL)
}
\arguments{
  \item{filename}{Input session file}
  \item{session}{Input session (from \command{siena01Gui}) }
  \item{modelName}{Character string of project name }
  \item{edited}{Boolean, indicates whether a file has been edited and
	therefore should not be re-read. Used internally by
	\code{siena01Gui}.}
  \item{files}{List of data files, used internally by \code{siena01Gui}}
  }
\details{
  Allows creation of data objects of class "Siena" direct from data
  files rather than from the various Siena network and covariate
  objects. Is always called by \code{siena01Gui} but can also be used
  directly.

  The columns of the gui
  screen should have the format described below. If a session file is used
  for input, it should have columns with exactly the same names
  and in exactly the same order as those below with a row of column
  headings and no row numbers.
  \describe{
	\item{Group}{Used to identify the groups when using the
	  multi-group option described in the Manual. Must not
	  contain embedded blanks, and should be identical for all rows which
	  relate to the same group.}
	\item{Name}{Network files or dyadic covariates should use the same name
	  for each file of the set. Other files should have unique names, a list of
	  space separated ones for constant covariates.}
	\item{File Name}{in \code{siena01Gui}, usually entered by using a
	  file selection box, after clicking \command{Add}.}
	\item{Format}{Only relevant for networks or dyadic covariates.  Can be
	  matrix, a single Pajek network (.net) or a Siena network file
	  (and edgelist with three or four columns: from, to, value,
	  wave (optional)).
	  Not tested for dyadic covariates yet!}
	\item{Period(s)}{Only relevant for networks and dyadic covariates. All
	  other files cover all the relevant periods. Indicates the order of
	  the network and dyadic covariate files. Should range from 1 to
	  \emph{n} within each group. Enter multiple integers with spaces
	  between for Siena network multi-wave files. Use the value 1 or blank
	  for other files
	  which cover multiple periods.}
	\item{ActorSet}{If you have more than one set of nodes, use this column to
	  indicate which is relevant to each file. Should not contain
	  embedded blanks.}
	\item{Type}{Indicate here what type of data the file contains. Options
	  are "network", "behavior", "constant covariate",
	  "changing covariate",
	  "constant dyadic covariate", "changing dyadic covariate",
	  "exogenous event".
	}
	\item{Selected}{Yes or No. Only files with Yes will be included in
	  the model.}
	\item{Missing Values}{Enter any values which indicate missingness, with
	  spaces between different entries.}
	\item{Nonzero Codes}{Enter any values which indicate ties, with spaces
	  between different entries.}
	\item{NbrOfActors}{For Siena network files, enter the number of
  actors here.}
  }
  If using a file for input, it should be of one of the following types:
  \tabular{ll}{
	Extension \tab Type\cr
	\code{.csv}\tab Comma separated\cr
	\code{.dat} or \code{.prn} \tab Space delimited\cr
	\code{.txt} \tab Tab delimited
  }
  Network and covariate files should be text files with a row for each node. The
  numbers should be separated by spaces or tabs.
  Exogenous events should be specified by a file with a row for each
  node. Each row should be consist of a set of pairs of numbers which
  indicate the periods during which the corresponding actor was present. e.g.
\preformatted{
1 3
1.5 3
1 1.4 2.3 3
2.4 3
}
  would describe a network with 4 nodes, and 3 observations. Actor 1 is present
  all the time, actor 2 joins at time 1.5, actor 3 leaves and time 1.4 then
  rejoins at time 2.3, actor 4 joins at time 2.4. All intervals are treated as
  closed.
}
\value{
  A list with the following components:
	\item{OK}{Boolean, TRUE indicating success}
	\item{mydata}{A Siena data or group object, of class \code{\link{siena}} or
	\code{\link{sienaGroup}}}
	\item{myeff}{Effects object associated with \code{mydata}}
}
\references{See \url{http://www.stats.ox.ac.uk/~snidjers/siena/ }}
\author{Ruth Ripley }
\seealso{ \code{\link{sienaDataCreate}}, \code{\link{siena}}}
%\examples{
%}
\keyword{classes }
