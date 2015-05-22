library(glmgen)
source("run_trial.R")
##

pls = c(40, 80, 120, 160, 200)
n = 400

#p = 10
#nls = c(30, 40)

ntrial = 15

paterr_mat = array(0, c(3, length(pls), ntrial))
prederr_mat = array(0, c(2, length(pls), ntrial))

for (it in 1:ntrial){
    
    for (ii in 1:length(pls)){
        print("new run.")
        print(c(it, ii))

        p = pls[ii]    
        tmp = matrix(rnorm(p*p), p, p)
        covar = t(tmp) %*% tmp
        D = diag(diag(covar)^(-1/2))
        covar = D %*% covar %*% D
        covar = (2*covar + matrix(1,p,p))/3        
        
    
        out = run_trial(covar, n, p)

        paterr_mat[, ii, it] = out$paterrs
        prederr_mat[, ii, it] = out$prederrs
    }
    save.image(file="tmp2.RData")
}

paterr_ave = apply(paterr_mat, c(1,2), mean)
prederr_ave = apply(prederr_mat, c(1,2), mean)

paterr_sd = apply(paterr_mat, c(1,2), sd)
prederr_sd = apply(prederr_mat, c(1,2), sd)

save("paterr_mat", "prederr_mat", "paterr_ave", "prederr_ave", "pls", "n", "paterr_sd", "prederr_sd", file="pls_experiment.RData")


