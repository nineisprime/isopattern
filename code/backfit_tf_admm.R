
## IN:
##   X  n--by--p
##   lambdas  p-list of lambda

## OUT:
##   fhats -- (n--by--p--by--nlambda)
##   pattern -- (p--by--nlambda)
##

#source("tf_mosek.R")
library(glmgen)

backfit_tf_admm <- function(y, X, lambdas, tforder){

    n = nrow(X)
    p = ncol(X)
    nlambda = length(lambdas)
    
    ordX = matrix(0, n, p)
    rankX = matrix(0, n, p)
    
    for (j in 1:p){
        ordX[, j] = order(X[, j])
        rankX[, j] = rank(X[, j])
    }

    pattern = matrix(0, p, nlambda)
    tol = 1e-4
    MAX_ITER = 1000
    
    resid = matrix(y, n, 1) %*% matrix(1, 1, nlambda); # (n--by--nlambda)
        
    objs = rep(sum(resid^2)/n, nlambda)
    obj_drops = rep(Inf, nlambda)
    cur_fs = array(0, dim=c(n,p,nlambda))
    #betahats = matrix(0, n-1, p)
    
    it = 0
    while ( sum(obj_drops > tol) || it < 40){
        it = it + 1;
        if (it > MAX_ITER)
            break
        
        betanorms = rep(0, nlambda)
        for (j in 1:p){
            fhats = matrix(0, n, nlambda)
            resid = resid + cur_fs[, j, ]

            for (il in 1:nlambda){
                mylambda = lambdas[il]*n*(1 + exp(-0.25*it + 5))
                
                curfit = trendfilter(y=resid[ordX[, j], il], x=1:n,
                                 lambda=mylambda, k=tforder,
                                 family="gaussian", verbose=0)

                fhats[, il] = 1
                fhats[, il] = predict(curfit, lambda=mylambda)
                                       
            }

            
            betahats = fhats[2:n, ] - fhats[1:(n-1), ];
                                        #betahat is (n-1 --by-- nlambda)

            fhats = fhats[rankX[, j], ] #fhats is (n--by--nlambda)

            if (nlambda == 1)
                betahats = matrix(betahats, n-1, 1)

            pattern[j, ] = 1
            zero_ixs = apply(abs(betahats), 2, sum) < 1e-6
            pattern[j, zero_ixs] = 0
            neg_ixs = apply(betahats, 2, sum) < -1e-6
            pattern[j, neg_ixs] = -1

            cur_fs[, j, ] = fhats
            resid = resid - cur_fs[, j, ]
            
            betanorms = betanorms + apply(abs(betahats), 2, sum)
        }

        new_objs = pmin(objs, apply(resid^2 /n, 2, sum) + lambdas * betanorms)

        #print(apply(resid^2 / n, 2, sum))
        obj_drops = objs - new_objs
        objs = new_objs
        print(objs)
    }

    return(list(fhats = cur_fs, pattern=pattern))
}

##
n = 50
p = 4
X = matrix(rnorm(n*p), n, p)
betastar = rnorm(p)

y = X %*% betastar

lambdas = c(0.01)
myfit = backfit_tf_admm(y, X, lambdas, 0)
