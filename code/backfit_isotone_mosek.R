
## IN:
##   X  n--by--p
##

## 
##

source("isotone1d_mosek.R")

backfit_isotone_mosek <- function(y, X, tforder){

    n = nrow(X)
    p = ncol(X)

    ordX = matrix(0, n, p)
    rankX = matrix(0, n, p)
    
    for (j in 1:p){
        ordX[, j] = order(X[, j])
        rankX[, j] = rank(X[, j])
    }

    pattern = rep(0, p)
    tol = 1e-4
    resid = y
    obj = sum(resid^2)/n
    obj_drop = Inf
    cur_fs = matrix(0, n, p)
    
    while (obj_drop > tol){
        
        for (j in 1:p){

            resid = resid + cur_fs[, j]
            
            curfit1 = isotone1d_mosek( resid[ordX[, j]], c(0), tforder)
            curfit2 = isotone1d_mosek( -resid[ordX[, j]], c(0), tforder)
            
            fhat1 = curfit1$fhats[rankX[, j], 1]
            fhat2 = -curfit2$fhats[rankX[, j], 1]
            
            err1 = sum(( resid - fhat1 )^2)
            err2 = sum(( resid - fhat2 )^2)

            if (err1 < err2){
                cur_fs[, j] = fhat1
                pattern[j] = 1
            } else {
                cur_fs[, j] = fhat2
                pattern[j] = -1
            }

            resid = resid - cur_fs[, j]
        }


        new_obj = min(obj, sum(resid^2)/n)
        obj_drop = obj - new_obj
        obj = new_obj
        print(new_obj)
        print(pattern)
    }

    return(list(fhats = cur_fs, pattern=pattern))
}

##
n = 200
p = 10
X = matrix(rnorm(n*p), n, p)
Sigma = toeplitz(0.6^(0:(p-1)))
X = X %*% Sigma^(1/2)

betastar = rnorm(p)

y = X %*% betastar

myfit = backfit_isotone_mosek(y, X, 0)
myfit2 = backfit_tf_mosek(y, X, 0.01, 0)

0
sign(betastar)
myfit$pattern
myfit2$pattern
