
# scan for marginal correlation to remove some covariates
cor.scale <- function(snp,y){
  a = snp -mean(snp)
  b = y - mean(y)
  abs(a%*%b) / sqrt(sum(a^2)*sum(b^2))
}
cor.scale <- compiler::cmpfun(cor.scale)

# fast function to calculate covariance matrix
cov2 <- function(x)crossprod( scale(x , TRUE ,FALSE ) )/(NROW(x) )
cov2 <- compiler::cmpfun(cov2)

# checking binary phenos
checkBinaryTrait <- function(y){
  if (!is.numeric(y)) stop("Only numeric phenotypes are accepted.")
  if(length(levels(factor(y) ) ) < 3) print("binomial")
  else print("continuous")
}
checkBinaryTrait <- compiler::cmpfun(checkBinaryTrait)

#' @title coheritability estimation
#'
#' @description coheritability estimation
#'
#' @param Y,
#' @param X,
#' @param a.quantile = 0.25,
#' @param alpha.enet = 1,
#' @param para = TRUE
#'
#' @return
#'
#' @examples
#'
#' @export coher
coher = function(Y,X,
                 a.quantile = 0.25,
                 alpha.enet = 1,
                 para = TRUE){
  library(parallel)
  require(doMC)
  registerDoMC(cores=detectCores()-1)
  library(glmnet)

  q = length(Y)
  p = ncol(X)
  bhat = list()
  nameSNP = colnames(X)
  nameSAMPLe = rownames(X)
  selected.snps = c()
  for(j in 1:q){
    matSNP = X[names(Y[[j]]),]
    if(a.quantile > 0){
      sam.cor <- apply(matSNP, 2,function(snp) cor.scale(snp, Y[[j]]) )
      sam.cor[is.na(sam.cor)] = 0
      xsele = matSNP[, sam.cor> quantile(sam.cor, a.quantile)]
      sk = Matrix(xsele, sparse=T)
    }else if(a.quantile ==0){
      sk = Matrix(X, sparse=T)
    }
    if(checkBinaryTrait(Y[[j]]) == 'binomial' ){
      enet <- cv.glmnet(x= sk, y= Y[[j]] ,
                        parallel = para,
                        alpha = alpha.enet,
                        family='binomial')
    }else{
      enet <- cv.glmnet(x= sk, y= Y[[j]] ,
                        parallel = para,
                        alpha = alpha.enet)
    }
    tam = coef(enet, s='lambda.min')[-1]
    names(tam) = colnames(xsele)
    bhat[[j]] = rep(0, p)
    names(bhat[[j]]) = nameSNP
    bhat[[j]][names(tam)] <- tam
    selected.snps = c(selected.snps,names(tam[tam!=0]) )
  }
  selected.snps = unique(selected.snps )

  # colculate the covariance matrix
  covariance.selesnps = cov2(X[,selected.snps] )

  #colculate the cohers
  mat.coher = matrix(0,nr=q,nc=q, dimnames = list(names(Y),names(Y)) )
  for (i in 1:(q-1)){
    for (j in (i+1):q){
      mat.coher[i, j] = bhat[[i]][selected.snps] %*% covariance.selesnps %*% bhat[[j]][selected.snps]
      mat.coher[j,i] = bhat[[i]][selected.snps] %*% covariance.selesnps %*% bhat[[j]][selected.snps]/ sqrt(bhat[[i]][selected.snps] %*% covariance.selesnps %*% bhat[[i]][selected.snps] * bhat[[j]][selected.snps] %*% covariance.selesnps %*% bhat[[j]][selected.snps]  )
    }
  }
  return(coher.matrix = mat.coher, coefficients.B = bhat)
}
coher <- compiler::cmpfun(coher)


