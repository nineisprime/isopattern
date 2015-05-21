## IN:
##   X  (n--by--p)
##   newx (nnew--by--p)
##   y  (n)
##   lambdas (nlambda) vector

## OUT:
##   pattern (p) 1/-1 vector
##   ynew (nnew) vector

## divides (X,y) into 4 folds for cross-validation
## minimize predictive error to select best lambda

source("../eval_additive_tf.R")
source("../backfit_tf_admm.R")

cv_pattern <- function(X, newx, y, lambdas, opt){
    n = nrow(X)
    p = ncol(X)
    nlambda = length(lambdas)
    nnew = nrow(newx)

    num_fold = 4
    nfold = floor(n/num_fold)

    errs = rep(0, nlambda)
    for (ii in 1:num_fold){
        validixs = (1 + (ii-1)*nfold):(ii*nfold)
        tmp = 1:n
        trainixs = tmp[!(tmp %in% validixs)]

        trainx = X[trainixs, ]
        trainy = y[trainixs]
        validx = X[validixs, ]
        validy = y[validixs]

        fitobj = backfit_tf_admm(trainy, trainx, lambdas, 0)
        ## fitobj$fhats (n--by--p--by--nlambda)
        ## fitobj$pattern (p--by--nlambda)

        for (il in 1:nlambda){
            evalobj = eval_additive_tf(fitobj$fhats[, , il], trainx, validx)
            errs[il] = errs[il] + sum( (evalobj$ynew - validy)^2 )/sum(validy^2)
        }
    }

    print(errs)
    
    
    ## select best lambda to use
    best_il = 0
    if (opt == "min"){
        best_il = which(errs == min(errs))
    } else {
        tmp = which(errs < min(errs) + sd(errs))
        best_il = tmp[length(tmp)]
    }

    print(best_il)
    
    fitobj = backfit_tf_admm(y, X, c(lambdas[best_il]), 0)
    evalobj = eval_additive_tf(fitobj$fhats[, , 1], X, newx)

    return(list(pattern = fitobj$pattern[, 1], ynew = evalobj$ynew))
}


###
n = 100
p = 3
lambdas = c(0.001, 0.01, 0.05, 0.1, 1)
X = matrix(rnorm(n*p), n, p)
betastar = runif(p) - 1/2
y = X %*% betastar + rnorm(n)

newx = matrix(rnorm(n*p), n, p)
newy = newx %*% betastar 

res = cv_pattern(X, newx, y, lambdas, opt="min")

sum((newy - res$ynew)^2)/sum(newy^2)
