\name{sienaModelOptions}
\alias{sienaModelOptions}
\title{Function to allow entry of model options}
\description{
  Displays a Gui with model options, and allows editing of effects plus
  running of Siena07
}
\details{
  Called from the \code{Apply} function in \code{\link{siena01Gui}}. An
  internal function of \code{\link{siena01Gui}}.

  Various parameters can be set on the upper part of the screen:
  \describe{
	\item{Estimation Method}{0 for Unconditional fitting, 1 for
	  Conditional. If there are multiple dependent variables, a list will be
	  displayed from which to choose.}
	\item{Standard starting value}{If checked, the estimation will ignore
	  the initial values of the parameters in the effects object, and use
	  the default ones.}
	\item{Specify random seed}{If you wish your run to be repeatable, check
	  this box and then choose any integer as the seed.}
	\item{Number of processors}{If checked, a box will appear for you to
	  select the number of processors to be used. All processes will run on
	  the same machine.}
	\item{Initial value of gain parameter}{A parameter to control the
	  Robbins-Monro algorithm. Will be multiplied by the number of
	  processors before use.}
	\item{Number of phase 2 subphases}{Default 4. To omit phase 2, set this
	  to 0.}
	\item{Derivative method}{0 for finite differences, 1 for score
	  function. Default 1.}
	\item{Number of phase 3 iterations}{Default 1000.}
  }
  If you wish to restrict the degree of the simulations, enter the value
  in the table on the bottom left.

  Desired effects can be selected by using the botton \code{Edit
	effects}. Change the \code{Include} column to a 1 to select, 0 to
  deselect.

  Initial values can be specified in the initialValues column.

  If it is desired to fix a parameter, set the \code{fix} column to 1.

  To request a test, set both the \code{test} and \code{fix} columns to 1
  and specify the value against which to test in the initialValue
  column.

  Some effects have parameter values: these can be specified in the
  \code{parm} column.

  Check the included effects by using the \code{Show included effects}
  button.

  The model can be fitted by using the \code{Estimate} button.

  The data objects can be saved to an R data set using \code{Save to
	file}.

  The results object can be saved to an R data set using \code{Save results}.

  The \code{Display Results} button is a toggle and should display or
  remove the display of the results file.

  \code{Exit Model Options} allows you to return to the previous screen.
}
\value{
None}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{\code{\link{siena01Gui}}, \code{\link{siena07}}}
\keyword{ misc }
