\name{sienaGroupCreate}
\alias{sienaGroupCreate}
\alias{sienaGroup}
\title{Function to group together several Siena data objects }
\description{Creates an object of class "sienaGroup" from a list of
  Siena data objects.
}
\usage{
sienaGroupCreate(objlist, singleOK = FALSE)
}
\arguments{
  \item{objlist}{ List of objects of class "siena" }
  \item{singleOK}{ Boolean: is it OK to only have one object?}
}
\details{

This function creates a Siena group object from several Siena data
objects, all of which use the same networks, covariates and actor sets
(although possibly different subsets of the actor set in each
object). It can be used as data input to \command{siena07} for the
multigroup option.  Also used internally for convenience with a single
Siena data object.}

\value{
  \item{... }{List containing the input objects, with attributes:}
  \item{netnames}{names of the dependent variables in each set}
  \item{symmetric}{vector of booleans, one for each dependent
	variable. TRUE if all occurrences of the network are symmetric.}
  \item{structural}{vector of booleans, indicating whether structurally
	fixed values occur in this network}
  \item{allUpOnly}{vector of booleans, indicating whether changes are
	all upwards in all the occurrences of this network}
  \item{allDownOnly}{similar to previous, but for downward changes}
  \item{anyUpOnly}{vector of booleans, indicating whether changes are
	all upwards in any of the occurrences of this network}
  \item{anyDownOnly}{similar to previous, but for downward changes}
  \item{types}{vector of network types of the dependent variables}
  \item{observations}{Total number of periods to process}
  \item{periodNos}{Sequence of numbers of periods which are not skipped
	in multigroup processing}
  \item{netnodeSets}{list of names of the node sets corresponding to the
	dependent variables}
  \item{cCovars}{names of the constant covariates, if any}
  \item{vCovars}{names of the changing  covariates, if any}
  \item{dycCovars}{names of the constant dyadic covariates, if any}
  \item{dyvCovars}{names of the changing dyadic covariates, if any}
  \item{ccnodeSets}{list of the names of the node sets corresponding to
	the constant covariates}
  \item{cvnodeSets}{list of the names of the node sets corresponding to
	the changing covariates}
  \item{dycnodeSets}{list of the names of the node sets corresponding to
	the constant dyadic covariates}
  \item{dyvcnodeSets}{list of the names of the node sets corresponding to
	the changing dyadic covariates}
  \item{compositionChange}{boolean: any composition change at all?}
  \item{exooptions}{named vector of composition change options for the
	node sets}
  \item{names}{Either from the input objects or "Data1", "Data2" etc}
  \item{class}{"sienaGroup" inheriting from "siena"}
  \item{balmean}{vector of means for balance calculations}
  \item{bRange}{vector of difference between maximum and minimum values
	for behavior variables, NA for other dependent variables}
  \item{behRange}{matrix of  maximum and minimum values
	for behavior variables, NA for other dependent variables}
  \item{bSim}{vector of similarity means for behavior variables, NA for
	other dependent variables}
  \item{bPoszvar}{vector of booleans indicating positive variance for
	behavior variables. NA for other dependent variables}
  \item{bMoreThan2}{vector of booleans indicating whether the behavior
	variables take more than 2 distinct values}
  \item{cCovarPoszvar}{vector of booleans indicating positive variance for
	constant covariates}
  \item{cCovarMoreThan2}{vector of booleans indicating whether the
	constant covariates take more than 2 distinct values}
  \item{cCovarRange}{vector of difference between maximum and minimum
	values for constant covariates}
  \item{cCovarRange2}{matrix of maximum and minimum values for constant
	covariates}
  \item{cCovarSim}{vector of similarity means for constant covariates}
  \item{cCovarMean}{vector of means for constant covariates}
  \item{vCovarRange}{vector of difference between maximum and minimum
	values for changing covariates}
  \item{vCovarSim}{vector of similarity means for changing covariates}
  \item{vCovarMoreThan2}{vector of booleans indicating whether the
	changing covariates take more than 2 distinct values}
  \item{vCovarPoszvar}{vector of booleans indicating positive variance for
	changing covariates}
  \item{vCovarMean}{vector of means for changing covariates}
  \item{dycCovarMean}{vector of means for constant dyadic covariates}
  \item{dycCovarRange}{vector of ranges for constant dyadic covariates}
  \item{dycCovarRange2}{matrix of maximum and minimum values
	for constant dyadic covariates}
  \item{dyvCovarRange}{vector of ranges for changing dyadic covariates}
  \item{dyvCovarMean}{vector of means for changing dyadic covariates}
  \item{anyMissing}{vector of booleans, one for each dependent variable,
	indicating the presence of any missing values}
  \item{netRanges}{matrix of maximum and minimum values for dependent
	networks, NA for behavior variables}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/ }}
\author{ Ruth Ripley }
\seealso{\code{\link{sienaDataCreate}} }
\examples{
## not very useful, but will work
mynet1 <- sienaNet(array(c(s501,s502),dim=c(50,50,2)))
mydata <- sienaDataCreate(mynet1)
mygroup <- sienaGroupCreate(list(mydata, mydata))
}
\keyword{ classes }
