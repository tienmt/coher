\name{cor.scale}
\alias{cor.scale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
marginal correlation
}
\description{
calculating marginal correlation
}
\usage{
cor.scale(snp, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{
 a SNP covariate
}
  \item{y}{
phenotype
}
}

\value{
marginal correlation
}

\author{
The Tien Mai
}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (snp, y)
{
    a = snp - mean(snp)
    b = y - mean(y)
    abs(a \%*\% b)/sqrt(sum(a^2) * sum(b^2))
  }
}
