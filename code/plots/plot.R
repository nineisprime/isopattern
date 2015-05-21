library(glmgen)
source("run_trial.R")
##

nls = c(100, 200, 300, 400, 500, 600)
p = 100

#p = 10
#nls = c(30, 40)

ntrial = 15

paterr_mat = array(0, c(3, length(nls), ntrial))
prederr_mat = array(0, c(2, length(nls), ntrial))

for (it in 1:ntrial){

    tmp = matrix(rnorm(p*p), p, p)
    covar = t(tmp) %*% tmp
    D = diag(diag(covar)^(-1/2))
    covar = D %*% covar %*% D
    covar = (2*covar + matrix(1,p,p))/3
    
    for (ii in 1:length(nls)){
        n = nls[ii]        
        out = run_trial(covar, n, p)

        paterr_mat[, ii, it] = out$paterrs
        prederr_mat[, ii, it] = out$prederrs
    }
    save.image(file="tmp.RData")
}

paterr_ave = apply(paterr_mat, c(1,2), mean)
prederr_ave = apply(prederr_mat, c(1,2), mean)

paterr_sd = apply(paterr_mat, c(1,2), sd)
prederr_sd = apply(prederr_mat, c(1,2), sd)

save("paterr_mat", "prederr_mat", "paterr_ave", "prederr_ave", "nls", "p", "paterr_sd", "prederr_sd", file="nls_experiment.RData")


