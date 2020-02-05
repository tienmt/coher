\name{coher}
\alias{coher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
coheritability estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
coher(Y, X, a.quantile = 0.25, alpha.enet = 1, para = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{a.quantile}{
%%     ~~Describe \code{a.quantile} here~~
}
  \item{alpha.enet}{
%%     ~~Describe \code{alpha.enet} here~~
}
  \item{para}{
%%     ~~Describe \code{para} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, X, a.quantile = 0.25, alpha.enet = 1, para = TRUE)
{
    library(parallel)
    require(doMC)
    registerDoMC(cores = detectCores() - 1)
    library(glmnet)
    q = length(Y)
    p = ncol(X)
    bhat = list()
    nameSNP = colnames(X)
    nameSAMPLe = rownames(X)
    selected.snps = c()
    for (j in 1:q) {
        matSNP = X[names(Y[[j]]), ]
        if (a.quantile > 0) {
            sam.cor <- apply(matSNP, 2, function(snp) cor.scale(snp,
                Y[[j]]))
            sam.cor[is.na(sam.cor)] = 0
            xsele = matSNP[, sam.cor > quantile(sam.cor, a.quantile)]
            sk = Matrix(xsele, sparse = T)
        }
        else if (a.quantile == 0) {
            sk = Matrix(X, sparse = T)
        }
        if (checkBinaryTrait(Y[[j]]) == "binomial") {
            enet <- cv.glmnet(x = sk, y = Y[[j]], parallel = para,
                alpha = alpha.enet, family = "binomial")
        }
        else {
            enet <- cv.glmnet(x = sk, y = Y[[j]], parallel = para,
                alpha = alpha.enet)
        }
        tam = coef(enet, s = "lambda.min")[-1]
        names(tam) = colnames(xsele)
        bhat[[j]] = rep(0, p)
        names(bhat[[j]]) = nameSNP
        bhat[[j]][names(tam)] <- tam
        selected.snps = c(selected.snps, names(tam[tam != 0]))
    }
    selected.snps = unique(selected.snps)
    covariance.selesnps = cov2(X[, selected.snps])
    mat.coher = matrix(0, nr = q, nc = q, dimnames = list(names(Y),
        names(Y)))
    for (i in 1:(q - 1)) {
        for (j in (i + 1):q) {
            mat.coher[i, j] = bhat[[i]][selected.snps] \%*\% covariance.selesnps \%*\%
                bhat[[j]][selected.snps]
            mat.coher[j, i] = bhat[[i]][selected.snps] \%*\% covariance.selesnps \%*\%
                bhat[[j]][selected.snps]/sqrt(bhat[[i]][selected.snps] \%*\%
                covariance.selesnps \%*\% bhat[[i]][selected.snps] *
                bhat[[j]][selected.snps] \%*\% covariance.selesnps \%*\%
                  bhat[[j]][selected.snps])
        }
    }
    return(coher.matrix = mat.coher, coefficients.B = bhat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
