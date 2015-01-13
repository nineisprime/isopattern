
library(glmnet)
library(MASS)

d = 3

n = 200
success = 0

Sigma = matrix( c(1, 0.85, -.95,
                  0.85, 1, -.95,
                -.95,  -.95,  1), 3, 3)

Sigma = matrix( c(1, .99, .3,
                  .99, 1, .3,
                  .3,  .3,  1), 3, 3)

ntrial = 300
M = matrix(0, n-1, n-1)
isomatX = isomat(n)
for (it in 1:ntrial){

Z = mvrnorm(n=n, mu=rep(0,d), Sigma=Sigma)

P_ls = list()
allX = c()
#isomatX = isomat(n)

for (j in 1:d){
    P = matrix(rep(0, n*n), n, n)

    for (i in 1:n) {
        P[order(Z[, j])[i], i] = 1
    }

    P_ls[[j]] = P
    allX = cbind(allX, P %*% isomatX)
}

M = M - t(isomatX) %*% t(P_ls[[1]]) %*% P_ls[[2]] %*% isomatX;
#M = M + t(P_ls[[1]]) %*% P_ls[[2]] 

}

#M1 = -isomatX
M = M/ntrial

j = 100
plot(M[,j])
#lines(1:n, M1[,j])
