## IN:
##   covariance, n, p

## OUT:
##   pattern recovery accuracy

source("monotone_functions.R")
source("cv_pattern.R")

run_trial <- function(covariance, n, p){

    ntest = 1000
    ntotal = n + ntest
    
    Z = matrix(rnorm(ntotal*p), nrow=ntotal, ncol=p)
    X = Z %*% chol(covariance)
    
    out = additive_monotone_fn(X)
    y = out$y
    true_pattern = out$pattern
    
    SNR = 4
    y = SNR*y/sd(y)
    y = y + rnorm(ntotal)
    
    Xtest = X[(n+1):ntotal, ]
    ytest = y[(n+1):ntotal]
    X = X[1:n, ]
    y = y[1:n]

    prederrs = c(0, 0)
    paterrs = c(0, 0, 0)
    
    
    lambdas = c(1e-5, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2,
                1e-1, 5e-1, 1)
    res = cv_pattern(X, Xtest, y, lambdas, opt="min")
    yhat = res$ynew
    patternhat = res$pattern

    prederrs[1] = sum( (ytest - yhat)^2)/sum(ytest^2)
    paterrs[1] = sum(abs(patternhat - true_pattern))/2

    ## test linear model
    lin_beta = solve(t(X) %*% X, t(X) %*% y)
    yhat = Xtest %*% lin_beta
    patternhat = sign(lin_beta)

    prederrs[2] = sum( (ytest - yhat)^2)/sum(ytest^2)
    paterrs[2] = sum(abs(patternhat - true_pattern))/2


    ## test marginal baseline
    patternhat = rep(0, p)
    for (j in 1:p){
        ordx = order(X[, j])
        yord = y[ordx]
        if (yord[n] - yord[1] > 0)
            patternhat[j] = 1
        else
            patternhat[j] = -1
    }
    paterrs[3] = sum(abs(patternhat - true_pattern))/2
    
    return(list(prederrs = prederrs, paterrs = paterrs))
}

##
covariance = matrix(c(1, .2, .7,
                      .2, 1, .1,
                      .7, .1, 1), 3, 3)
p = 3
n = 120
res = run_trial(covariance, n, p)
res


##
p = 10
n = 70
tmp = matrix(rnorm(p*p), p, p)
covar = t(tmp) %*% tmp
D = diag(diag(covar)^(-1/2))
covar = D %*% covar %*% D
covar = (covar + matrix(1, p, p))/2
res = run_trial(covar, n, p)
res

##
p = 50
n = 100
covar = toeplitz(0.6^(0:(p-1)))
res = run_trial(covar, n, p)
res

