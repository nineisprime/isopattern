## IN:
##   X  (n--by--p) vector
##
## OUT:
##   y  (n) vector, centered
##  pattern (p) vector, random 0,1


additive_monotone_fn <- function(X){
    n = nrow(X)
    p = ncol(X)

    tmp = rnorm(p)
    pattern = rep(0, p)
    pattern[tmp < 0] = -1
    pattern[tmp > 0] = 1

    y = rep(0, n)
    fully = matrix(0, n, p)
    cases = rep(0,p)
    for (j in 1:p){
        case = sample(1:4, 1)
        cases[j] = case
        if (case == 1){
            fully[, j] = exp2fun(X[, j]) * pattern[j]
            #fully[, j] = stepfun(X[, j]) * pattern[j]            
            y = y + fully[, j]
            
        } else if (case == 2){
            fully[, j] = expfun(X[, j]) * pattern[j]
            #fully[, j] = twostepfun(X[, j]) * pattern[j]
            y = y + fully[, j]
            
        } else if (case == 3){
            fully[, j] = stepfun(X[, j]) * pattern[j]
            y = y + fully[, j]
            
        } else if (case == 4){
            fully[, j] = twostepfun(X[, j]) * pattern[j]
            y = y + fully[, j]
        }
    }
    print(cases)
    return(list(y=y, pattern=pattern, fully = fully))
}

linfun <- function(x){
    return(x - mean(x))
}


exp2fun <- function(x){
   y = - exp(-x)
   y = y - mean(y)
   y = 3*y/max(abs(y))
   return(y)
}


expfun <- function(x){
    y = exp(x)
    y = y - mean(y)
    y = 3*y/max(abs(y))
    #y = y - mean(y)
    return(y)
}


stepfun <- function(x){
    n = length(x)
    y = rep(0, n)
    ordx = order(x)
    
    ixs = seq(floor(n/5), n)
    y[ordx[ixs]] = 1

    return(y - mean(y))
}

twostepfun <- function(x){
    n = length(x)
    y = rep(0, n)
    ordx = order(x)
    
    ixs1 = seq(1, floor(n/2))
    ixs2 = seq(floor(4*n/5), n)
    y[ordx[ixs1]] = -1
    y[ordx[ixs2]] = 1

    return(y - mean(y))
}


##
n = 20
p = 10
X = matrix(rnorm(n*p), n, p)
res = additive_monotone_fn(X)
