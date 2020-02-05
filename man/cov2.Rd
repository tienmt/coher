\name{cov2}
\alias{cov2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fast covariance estimation for selected covariates
}
\description{
fast covariance estimation for selected covariates
}
\usage{
cov2(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 a matrix
}
}
\value{
covariance matrix is returned
}
\author{
The Tien Mai
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
crossprod(scale(x, TRUE, FALSE))/(NROW(x))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
