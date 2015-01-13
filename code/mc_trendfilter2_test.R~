#Create the matrix mat which plays the role of
#design matrix X in the Lasso formulation.
isomat2 = function(n){

    X = matrix(0, n, n)
    X[lower.tri(X, diag=TRUE)] = 1

    X = X %*% X

    X = (1/n) * X
    X[, 1] = 1
    return(X)
}

isomat1 = function(n){

    X = matrix(0, n, n)
    X[lower.tri(X, diag=TRUE)] = 1
    return(X)
}


library(glmnet)
library(doMC)

registerDoMC(4)

niter = 40

n_ls = c(300, 600, 900, 1200, 1800)
s0_ls = c(3)
lambda_ls = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1, 1.5)

results = array(0, c(length(s0_ls), length(n_ls), length(lambda_ls)))

for (ii in 1:length(s0_ls)){
    for (jj in 1:length(n_ls)){
        for (kk in 1:length(lambda_ls)){

        risks = foreach (it=1:niter) %dopar% {
    
            print(c(it, ii, jj, kk));
    
            n = n_ls[jj]
            ##n = 400
            X = isomat1(n)
            s0 = s0_ls[ii]
            ##s0 = 3
            c = 1
            sigma = 0.2
            betastar = rep(0, n);
            delta = floor(n/s0)
            ixs = (0:s0)*delta
            ixs[1] = 1
            ixs = ixs[-(s0+1)]
            ##
            betastar[ixs] = 1
            betastar[ixs] = betastar[ixs] * ((-1)^(1:s0))
            ##betastar[ixs] = sign(runif(n=s0)-0.5)
            betastar = c*betastar/sum(abs(betastar))
            ##betastar[ixs] = 0
            epsilon = rnorm(n, 0, sigma)

            #const = 0.02
            const = 0.2
            lambda = const*sigma*(log(n)/n)^(1/2)
            fstar = X %*% betastar
            fstar = fstar - mean(fstar)
            y = fstar + epsilon
            ##
            ##

            glmfit = glmnet(x=X, y=y, lambda=lambda, standardize=FALSE,
                family="gaussian", thresh=1e-8, maxit=200000)
            betahat = as.vector(glmfit$beta)
            
            fhat1 = X %*% betahat
            fhat1 = fhat1 - mean(fhat1)

            (100/n)*sum( (fstar - fhat1)^2)
        }

        results[ii,jj,kk] = mean(unlist(risks))
        save(results, file="Res.RData")
        
    }

    }
}



## n=600
## p=200
## sigma=0.5
## X = matrix(rnorm(n*p), n, p)
## betastar = rep(0, p)
## betastar[1] = 1
## betastar = as.matrix(betastar)
## y = X %*% betastar + rnorm(n, sigma)
## glmfit = glmnet(x=X, y=y, lambda=sigma*sqrt(log(p)/n))
## betahat = glmfit$beta
## norm(betahat - betastar)
## norm(betastar)
## norm(X %*% (betahat - betastar))
