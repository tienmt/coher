\name{checkBinaryTrait}
\alias{checkBinaryTrait}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
checkBinaryTrait
}
\description{
checkBinaryTrait
}
\usage{
checkBinaryTrait(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
phenotype
}
}

\value{
binary or continuous
}

\author{
The Tien Mai
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y)
{
    if (!is.numeric(y))
        stop("Only numeric phenotypes are accepted.")
    if (length(levels(factor(y))) < 3)
        print("binomial")
    else print("continuous")
  }
}
