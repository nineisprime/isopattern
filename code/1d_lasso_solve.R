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
n = 900
X = isomat(n);
deltastar = 1*rep(1,n-1)/(n-1);

epsilon = rep(0, n);

lambda = 0.2/sqrt(n)

y = X %*% deltastar + epsilon;
glmfit = glmnet(x=X, y=y, lambda=lambda, standardize=FALSE, family="gaussian", thresh=1e-15, maxit=200000, intercept=FALSE)
xtx = t(X) %*% X;

deltahat = as.vector(glmfit$beta)
