#' coher
#'
#' @description coheritability estimation
#'
#' @param Y a list of vectors indicating the values for the phenotypes of intererest for each of the genomes
#' @param X binary SNP/unitig presence/absence matrix
#' @param a.quantile (default = 0.25)
#' @param alpha.enet (default = 1)
#' @param para (default = TRUE)
#'
#' @return a list with the inferred co-heritability matrix and the coefficients of the elastic-net model for each phenotype
#'
#' @examples
#' data('coher_example')
#' result <- coher(coher_example$Y, coher_example$X)
#'
#'
#' @export
coher = function(Y, X,
                 a.quantile = 0.25,
                 alpha.enet = 1,
                 para = TRUE){

  doMC::registerDoMC(cores=detectCores()-1)

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
      sk = Matrix::Matrix(xsele, sparse=TRUE)
    }else if(a.quantile ==0){
      sk = Matrix::Matrix(X, sparse=TRUE)
    }
    if(checkBinaryTrait(Y[[j]]) == 'binomial' ){
      enet <- glmnet::cv.glmnet(x= sk, y= Y[[j]] ,
                        parallel = para,
                        alpha = alpha.enet,
                        family='binomial')
    }else{
      enet <- glmnet::cv.glmnet(x= sk, y= Y[[j]] ,
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

  # calculate the covariance matrix
  covariance.selesnps = cov2(X[,selected.snps] )

  #calculate the cohers
  mat.coher = matrix(0,nr=q,nc=q, dimnames = list(names(Y),names(Y)) )
  for (i in 1:(q-1)){
    for (j in (i+1):q){
      mat.coher[i, j] = bhat[[i]][selected.snps] %*% covariance.selesnps %*% bhat[[j]][selected.snps]
      mat.coher[j,i] = bhat[[i]][selected.snps] %*% covariance.selesnps %*% bhat[[j]][selected.snps]/ sqrt(bhat[[i]][selected.snps] %*% covariance.selesnps %*% bhat[[i]][selected.snps] * bhat[[j]][selected.snps] %*% covariance.selesnps %*% bhat[[j]][selected.snps]  )
    }
  }
  return(list(coher.matrix = mat.coher, coefficients.B = bhat))
}

# scan for marginal correlation to remove some covariates
cor.scale <- function(snp,y){
  a = snp -mean(snp)
  b = y - mean(y)
  abs(a%*%b) / sqrt(sum(a^2)*sum(b^2))
}

# fast function to calculate covariance matrix
cov2 <- function(x) Matrix::crossprod( scale(x , TRUE ,FALSE ) )/(NROW(x) )

# checking binary phenos
checkBinaryTrait <- function(y){
  if (!is.numeric(y)) stop("Only numeric phenotypes are accepted.")
  if(length(levels(factor(y) ) ) < 3) print("binomial")
  else print("continuous")
}
