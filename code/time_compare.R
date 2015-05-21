

library(Rmosek)
n = 2000;
##y = rnorm(n);

tforder = 1;
X = xmat(n, tforder);

betastar = rep(0, n-tforder-1);
betastar[floor(n/2)] = 2;
ixs = floor(n/5):floor(n/4)
betastar[ixs] = -0.05


y = (1/n) * X %*% betastar + rnorm(n);

lambdas = c(0.5, 0.1, 0.05, 0.01);
#lambdas = c(0.1);

ptm = proc.time()
res = tf_mosek(y, lambdas, 1);
proc.time() - ptm


plot(y)
lines(res$fhats[,2])


ptm = proc.time()
res2 = trendfilter(y, 1:n, k=1, lambda=n*lambdas)
proc.time() - ptm

yy = predict(res2, x.new=1:n, lambda=n*lambdas[2])
plot(y)
lines(yy)

