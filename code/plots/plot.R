
##

#nls = c(100, 200, 300, 400, 500, 600)
#p = 100

p = 10
nls = c(30, 40)

ntrial = 2

paterr_mat = matrix(0, 3, ntrial)
prederr_mat = matrix(0, 2, ntrial)

for (it in 1:ntrial){

    tmp = matrix(rnorm(p*p), p, p)
    covar = t(tmp) %*% tmp
    D = diag(diag(covar)^(-1/2))
    covar = D %*% covar %*% D
    covar = (2*covar + matrix(1,p,p))/3
    
    for (ii in 1:length(nls)){
        n = nls[ii]        
        out = run_trial(covar, n, p)

        paterr_mat[, it] = out$paterrs
        prederr_mat[, it] = out$prederrs
    }
    
}
