##
## tf_glm
##
## IN: y vector
##     lambda
##     tforder integer

## OUT: fhat
##      X
##      betahat

isotone1d_mosek <- function(y, lambdas, tforder){


    ## variable ordering:
    ## f (n), s (n-tforder-1), t(n), t'(1), t''(1)
    n = length(y);
    ns = n-tforder-1;
    
    D = create_D(n, tforder);

    fhats = matrix(0, n, length(lambdas));
    betahats = matrix(0, ns, length(lambdas));

    j = 0;
    for (lambda in lambdas){
        j = j+1;
        
        obj = c(rep(0, n), lambda*rep(1,ns), rep(0,n), 1/n, 0);

        # regularization inequalities
        A = cBind(D, -Diagonal(ns), Matrix(0, ns, n+2));
        #A = rBind(A, cBind(-D, -Diagonal(ns), Matrix(0, ns, n+2)));
        buc = rep(0, ns);
        blc = rep(0, ns);

        ## objective equality
        A = rBind(A, cBind(Diagonal(n), Matrix(0,n,ns), Diagonal(n), Matrix(0,n,2)));
        blc = c(blc, y);
        buc = c(buc, y);

        A = rBind(A, cBind(Matrix(1,1,n), Matrix(0, 1, n+ns+2)));
        blc = c(blc, 0);
        buc = c(buc, 0);

        bux = c(rep(Inf, 2*n+ns+1), 1)
        blx = c( rep(-Inf, n), rep(0, ns), rep(-Inf, n), 0, 1);

        cqo = list(sense="min");
        cqo$A = A;
        cqo$bc = rbind(blc = blc, buc = buc);
        cqo$bx = rbind(blx = blx, bux = bux);
        cqo$cones = cbind(
            list("RQUAD", c(2*n + ns + 1, 2*n + ns + 2, seq(n + ns + 1, 2*n + ns)))
            );
        cqo$c = obj

        res = mosek(cqo, list(verbose=0));

        fhats[, j] = res$sol$itr$xx[1:n];

        betahats[, j] = as.vector(D %*% fhats[, j]);
        
        
     }
    return(list(fhats = fhats, betahats = betahats))
    
}


create_D <- function(n, tforder){

    i = 0;

    Dfinal = diag(1,n,n);
    while (i <= tforder){

        D = matrix(0, n-1-i, n-i);
        diag(D[, -1]) = 1
        diag(D[, -(n-i)]) = -1

        Dfinal = D %*% Dfinal;
        i = i + 1;
    }

    Dfinal = Matrix(Dfinal, sparse=TRUE)
    return(Dfinal)
}


## Example Usage:
library(Rmosek)
n = 200;
##y = rnorm(n);

X = xmat(n, 0);
betastar = rep(0, n-1);
betastar[floor(n/2)] = 1;
y = X %*% betastar + 0.2*rnorm(n);

#lambdas = c(0.5, 0.1, 0.05, 0.01);
lambdas = c(0);

ptm = proc.time()
res = isotone1d_mosek(y, lambdas, 1);
proc.time() - ptm
