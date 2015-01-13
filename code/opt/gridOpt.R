library(Matrix)

## "D" is m--by--m density
## "Y" is m--by--m function values

gridopt <- function(Y, D, lambda){

    m = nrow(Y)
    maxiter = 30

    ## store current fits
    z = matrix(0, nrow=m, ncol=2)
    obj = Inf

    progress = rep(0, maxiter)
    for (it in 1:maxiter) {

        zmat2 = matrix(1, nrow=m, ncol=1) %*% t(z[, 2])
        wt1 = apply(D, 1, sum)
        resfor1 = apply( (Y - zmat2)*D, 1, sum)

        opt1res = mono_mosek(resfor1, m, wt1, lambda)
        z[, 1] = opt1res

        #post.val = sum( - 2*sum(resfor1*z[,1]) + sum(wt1*z[,1]^2) +
        #     lambda*(z[m,1] - z[1,1]) )
        
        zmat1 = matrix(1, nrow=m, ncol=1) %*% t(z[, 1])
        resfor2 = apply( (t(Y) - zmat1)*t(D), 1, sum)
        wt2 = apply(D, 2, sum)       

        opt2res = mono_mosek(resfor2, m, wt2, lambda)
        z[, 2] = opt2res

        zmat2 = matrix(1, nrow=m, ncol=1) %*% t(z[, 2])
        
        new_obj = sum( D*(Y - t(zmat1) - zmat2)^2 ) +
                    lambda*( z[m,1] - z[1,1] + z[m,2] - z[1,2] )
        
        
        progress[it] = new_obj
        
        if ((obj - new_obj) < 1e-20)
            break

        obj = new_obj
    }
    print(it)
    print(progress)
    flush.console()
    return(z)
}

## i = 1,...,m
## min_f  sum_i -2*yi*f(xi) + wti*f(xi)^2 + lambda(f(xm) - f(x1))
##   s.t. f(xi) is increasing

last <- function(vec){
    return(vec[length(vec)])
}


mono_mosek <- function(resfor1, m, wt1, lambda){

    f.index = 1:m
    r.index = last(f.index) + 1:m

    t.index = last(r.index) + 1
    s.index = t.index + 1

    num.vars = s.index

    iso.fit <- list(sense="min")

    vec = rep(0, num.vars)
    vec[t.index] = 1
    vec[last(f.index)] = lambda
    vec[f.index[1]] = -lambda

    iso.fit$c <- vec
    
    A1 <- Matrix(0, nrow=m, ncol=num.vars)

    
    A1[f.index, f.index] = diag(wt1)
    A1[f.index, r.index] = -diag(sqrt(wt1))
        
    A2 = Matrix(0, nrow=m-1, ncol=num.vars)
    A2[1:(m-1), f.index] = bandSparse(m-1, m, c(0,1), list(rep(-1,m-1), rep(1,m-1)))

    iso.fit$A <- rBind(A1,A2)

    iso.fit$bc <- rbind(
        blc = c(resfor1, rep(0, m-1)),
        buc = c(resfor1, rep(Inf, m-1))
        )

    iso.fit$bx <- rbind(
        blx = c(rep(-Inf, num.vars-2), 0, 1/2),
        bux = c(rep(Inf, num.vars-2), Inf, 1/2)
        )

    iso.fit$cones <- cbind(
        list("RQUAD", c(t.index, s.index, r.index))
        )

    opts = list()
    opts$verbose = 1
    result <- mosek(iso.fit, opts)

    fhat = result$sol$itr$xx[f.index]
    return(fhat)
}
        

smoothDensity <- function(D){

    m = nrow(D)

    Dnew = D
    
    for (i in 2:(m-1)){
        for (j in 2:(m-1)){
            
            Dnew[i,j] = 0.4*D[i,j] +
                  0.6*0.25*(D[i-1,j] + D[i,j-1] + D[i+1,j] + D[i,j+1])
            
        }
    }
    return(Dnew)
}
