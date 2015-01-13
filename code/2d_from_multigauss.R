#Create the matrix mat which plays the role of
#design matrix X in the Lasso formulation.
isomat = function(n){
m1 = matrix(0,n,(n-1))
for (i in 1:(n-1)){
  for (j in i:(n-1)){
    m1[i,j] = -1
  }
}

m2 = matrix(0,n,(n-1))
for (j in 1:(n-1)){
  m2[,j] = rep(j,n)/n
}

mat = m1 + m2
return(mat)
}




library(glmnet)
library(MASS)

n = 2000
success = 0

corr = -0.9
Sigma = matrix(c(1,corr,corr,1), 2,2)

ntrial = 10

#for (it in 1:ntrial){
#print(it)
##
Z = mvrnorm(n=n, mu=rep(0,2), Sigma=Sigma)
P1 = matrix(rep(0, n*n), n, n)
P2 = matrix(rep(0, n*n), n, n)
for (i in 1:n) {
    P1[rank(Z[ ,1])[i], i] = 1
    P2[i, rank(Z[, 2])[i]] = 1
}

X1 = isomat(n)
X2 = P1 %*% P2 %*% X1
X = cbind(X1,X2)



f <- list()
f[[1]] = function(x){
  #return(2*x)
  return(10 * x^3)
}
f[[2]] = function(w){
  ans = 4*w
  return(ans)
}

Zunif1 = sapply(Z[,1], FUN=pnorm) - 0.5
Zunif2 = sapply(Z[,2], FUN=pnorm) - 0.5

f1 = sapply(Zunif1, FUN=f[[1]])
f2 = sapply(Zunif2, FUN=f[[2]])

f1 = f1[ order(Z[,1]) ]
f2 = f2[ order(Z[,2]) ]

delta1star = f1[2:n] - f1[1:(n-1)]
delta2star = f2[2:n] - f2[1:(n-1)]
deltastar = c(delta1star, delta2star)


keym = round(solve( t(X1) %*% X1) %*% t(X2) %*% X1)

epsilon = rep(0, n)
y = X %*% deltastar 


lambda = 0.01
#lambda = 0.1/n^(1/5)
#lambda = 0.05
glmfit = glmnet(x=X, y=y, lambda=lambda, standardize=FALSE, family="gaussian", thresh=1e-7, maxit=5000, intercept=FALSE)

deltahat = as.vector(glmfit$beta)
deltahat1 = deltahat[1:(n-1)]
deltahat2 = deltahat[n:(2*n-2)]

f1hat = X1 %*% deltahat1;
f2hat = X1 %*% deltahat2;

par(mfrow=c(1,2))

x = sort(Zunif1)
y = y - rep(mean(y), n)
plot(x,y, main="bleh")
lines(x, f1, lwd = 2, col = "green")
lines(x, f1hat, lwd = 2, col = "blue")

x = sort(Zunif2)
y = y - rep(mean(y), n)
plot(x,y, main="bleh")
lines(x, f2, lwd = 2, col = "green")
lines(x, f2hat, lwd = 2, col = "blue")

