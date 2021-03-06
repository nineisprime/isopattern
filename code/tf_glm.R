##
## tf_glm
##
## IN: y vector
##     lambda
##     tforder integer

## OUT: fhat
##      X
##      betahat


tf_glm <- function(y, lambdas, tforder){

    n = length(y);
    X = xmat(n, tforder);

    glmfit = glmnet(x=X, y=y, lambda=lambdas, standardize=FALSE,
            family="gaussian", thresh=1e-8, maxit=200000)

    betahat = as.matrix(glmfit$beta);
    
    fhat = X %*% betahat;

    return(list(fhat = fhat, betahat = betahat, x = X))
}



xmat <- function(n, tforder) {

    X = matrix(0, n, n-1);
    X[lower.tri(X, diag=FALSE)] = 1;

    X = X - matrix(1/n, n, n) %*% X 
    
    if (tforder == 0){
        return(X)
    }

    i = 0;
    while (i < tforder) {
        i = i + 1;
        Xcur = matrix(0, n-i, n-i-1);
        Xcur[lower.tri(Xcur, diag=FALSE)] = 1;
        X = X %*% Xcur;
    }

    return(X)
}
    

## Example Usage:
library(glmnet)
n = 200;
##y = rnorm(n);

X = xmat(n, 0);
betastar = rep(0, n-1);
betastar[floor(n/2)] = 1;
y = X %*% betastar + 0.2*rnorm(n);

lambdas = c(0.5, 0.1, 0.05, 0.01);

ptm = proc.time()
res = tf_glm(y, lambdas, 1);
proc.time() - ptm

