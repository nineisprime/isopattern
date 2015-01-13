library(Rmosek)


m = 100
D = matrix(2*abs(rnorm(m*m)), m, m)
#D = matrix(1,m,m)
diagmat = diag(rep(1,m))
diagmat = diagmat[m:1,]
D = D + diagmat*2
D = smoothDensity(D)
#D = smoothDensity(D)

D = D/sum(D)

D1 = rowSums(D)
D2 = rowSums(t(D))

isomat = function(n) {
    
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

X = isomat(m)

betastar1 = runif(m-1)
for (i in 1:(m-1)) {
    cut = runif(1)*0.8
    betastar1[i][betastar1[i] > cut] = 0
}
betastar1 = betastar1/sum(betastar1)
#betastar1 = c(rep(0,12),1,rep(0,16))

fstar1 = X %*% betastar1 
fstar1 = fstar1 - sum(fstar1 * D1)


betastar2 = runif(m-1)
for (i in 1:(m-1)) {
    cut = runif(1)*0.2
    betastar2[i][betastar2[i] > cut] = 0
}
betastar2 = betastar2/sum(betastar2)
#betastar2 = c(rep(0,15),1,rep(0,13))

fstar2 = X %*% betastar2
fstar2 = fstar2 - sum(fstar2 * D2)

Y = matrix(fstar1, m, 1) %*% matrix(1, 1, m) + 
    matrix(1, m, 1) %*% matrix(fstar2, 1, m)

fhats = gridopt(Y, D, 0.01)


par(mfrow=c(1,2))

plot(1:m, fstar1, type="l", col="red")
lines(1:m, fhats[,1])


plot(1:m, fstar2, type="l", col="red")
lines(1:m, fhats[,2])
