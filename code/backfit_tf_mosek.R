
## IN:
##   X  n--by--p
##

## 
##

source("tf_mosek.R")

backfit_tf_mosek <- function(y, X, lambda, tforder){

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

        betanorm = 0
        for (j in 1:p){

            resid = resid + cur_fs[, j]
            
            curfit = tf_mosek( resid[ordX[, j]], c(lambda), tforder)
            
            fhat = curfit$fhats[rankX[, j], 1]
            betahat = curfit$betahats[, 1]

            if (sum(betahat) > 0){
                pattern[j] = 1
            } else {
                pattern[j] = -1
            }

            cur_fs[, j] = fhat
            resid = resid - cur_fs[, j]
            betanorm = betanorm + sum(abs(betahat))
        }


        new_obj = min(obj, sum(resid^2)/n + lambda*betanorm)
        obj_drop = obj - new_obj
        obj = new_obj
        print(new_obj)
        print(pattern)
    }

    return(list(fhats = cur_fs, pattern=pattern))
}

##
n = 500
p = 20
X = matrix(rnorm(n*p), n, p)
betastar = rnorm(p)

y = X %*% betastar

myfit = backfit_tf_mosek(y, X, 0.001, 0)
