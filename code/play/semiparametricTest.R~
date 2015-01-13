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

d = 3

n = 1500
success = 0

Sigma = matrix( c(1, 0.85, -.95,
                  0.85, 1, -.95,
                -.95,  -.95,  1), 3, 3)

Sigma = matrix( c(1, 0.8, -0.8,
                  0.8, 1, -0.8,
                  -0.8,  -0.8,  1), 3, 3)

#Sigma = matrix( c(1, -0.9, -0.9, 1), 2, 2)

ntrial = 10

#for (it in 1:ntrial){
#print(it)
##
Z = mvrnorm(n=n, mu=rep(0,d), Sigma=Sigma)


P_ls = list()
allX = c()
isomatX = isomat(n)


#slits = c(n/10, n/2, floor(n/8))
slits = c(n/2, floor(n/2), floor(n/4))

for (j in 1:d){
    P = matrix(rep(0, n*n), n, n)
    for (i in 1:n) {
        #P[i, order(Z[,j])[i]] = 1
        P[order(Z[, j])[i], i] = 1
    }
    P_ls[[j]] = P
    allX = cbind(allX, P %*% isomatX[, slits])
}



#######
## Define f* and delta*

f <- list()
f[[1]] = function(x){
  #return(2*x)
  return(x^3)
}
f[[2]] = function(w){
  ans = 4*w
  return(ans)
}


#Zunif = matrix(0, n, d)
#fvals = matrix(0, n, d)
#deltastar = matrix(0, n-1, d)

deltastar = matrix(0, d, d)
diag(deltastar) = 1
#deltastar[slits[1], 1] = 1
#deltastar[slits[2], 2] = 1
#deltastar[(9/10)*n, 3] = 1
#deltastar[slits[3], 3] = 1


Zunif = matrix(0, n, d)
for (j in 1:d){
    Zunif[, j] = sapply(Z[, j], FUN=pnorm) - 0.5
#    fvals[, j] = sapply(Zunif[, j], FUN=f[[1]])
#    fvals[, j] = fvals[order(Z[, j]), j]
#    deltastar[, j] = fvals[2:n, j] - fvals[1:(n-1), j]
}

#deltastar[, 1] = matrix(0, n-1, 1)
#keym = round(solve( t(X1) %*% X1) %*% t(X2) %*% X1)

epsilon = rep(0, n)
y = allX %*% matrix(deltastar, d*d, 1)


#######
### Solve!
##

lambda = 0.05
#lambda = 0.1/n^(1/5)
#lambda = 0.01
glmfit = glmnet(x=allX, y=y, lambda=lambda, standardize=FALSE, family="gaussian", thresh=1e-8, maxit=5000, intercept=FALSE)
deltahat = matrix(as.vector(glmfit$beta), d, d)





f1hat = isomatX %*% deltahat[, 1]
f2hat = isomatX %*% deltahat[, 2]

f3hat = isomatX %*% deltahat[, 3]

f1 = isomatX %*% deltastar[, 1]
f2 = isomatX %*% deltastar[, 2]

f3 = isomatX %*% deltastar[, 3]

#par(mfrow=c(1,d))

x = sort(Zunif[, 1])
y = y - rep(mean(y), n)
plot(x,y, main="bleh")
lines(x, f1, lwd = 2, col = "green")
lines(x, f1hat, lwd = 2, col = "blue")

x = sort(Zunif[, 2])
y = y - rep(mean(y), n)
plot(x,y, main="bleh")
lines(x, f2, lwd = 2, col = "green")
lines(x, f2hat, lwd = 2, col = "blue")

x = sort(Zunif[, 3])
y = y - rep(mean(y), n)
plot(x,y, main="bleh")
lines(x, f3, lwd = 2, col = "green")
lines(x, f3hat, lwd = 2, col = "blue")


yy3 = t(isomatX) %*% (f3hat - f3)
yy2 = t(isomatX) %*% t(P_ls[[3]]) %*% P_ls[[2]] %*% (f2hat - f2)
yy1 = t(isomatX) %*% t(P_ls[[3]]) %*% P_ls[[1]] %*% (f1hat - f1)
#yy1 = t(isomatX) %*% (f1hat - f1)

plot(1:(n-1), yy3, col="green", ylim=c(-70, 50))
lines(1:(n-1), yy2, col="blue")
lines(1:(n-1), yy1, col="black")
lines(1:(n-1), yy1+yy2+yy3, col="red")


yy3 = t(isomatX) %*% (f1hat - f1)
yy2 = t(isomatX) %*% t(P_ls[[1]]) %*% P_ls[[2]] %*% (f2hat - f2)
yy1 = t(isomatX) %*% t(P_ls[[1]]) %*% P_ls[[3]] %*% (f3hat - f3)


plot(1:(n-1), yy3, col="green", ylim=c(-70, 50))
lines(1:(n-1), yy2, col="blue")
lines(1:(n-1), yy1, col="black")
lines(1:(n-1), yy2+yy3+yy1, col="red")



yy3 = t(isomatX) %*% (f2hat - f2)
yy2 = t(isomatX) %*% t(P_ls[[2]]) %*% P_ls[[1]] %*% (f1hat - f1)
#yy1 = t(isomatX) %*% (f1hat - f1)

plot(1:(n-1), yy3, col="green", ylim=c(-70, 50))
lines(1:(n-1), yy2, col="blue")
lines(1:(n-1), yy3+yy2, col="red")
#lines(1:(n-1), yy1, col="red")
