n = 100
## generate P, X
# Import Rmosek
if (!require("Rmosek")) {
  stop ("Rmosek not installed.")
}

######################################
#Generates the n by n-1 Isotonic matrix.
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

################################ 
X = isomat(n)


################################
#Generate permutation P#########
new.ord = sample(1:n, size=n, replace=FALSE, prob=rep(1,n))
#new.ord = sample(1:n, size=n, replace=FALSE, prob=exp( 5*(1:n)/n))
#new.ord = n+1 - 1:n
#new.ord = c(10,5,2,8,9,3,4,6,7,1)
new.ord = as.integer(new.ord)
sigma = invPerm(new.ord)
PX = X[new.ord, ]

# third variable
new.ord2 = sample(1:n, size=n)
P2X = X[new.ord2, ]


################################
beta1star = runif(n-1)
beta1star = beta1star/sum(beta1star)


beta2star = runif(n-1)
beta2star[beta2star > 0.05] = 0
beta2star[floor(n/2)] = 0.1
beta2star = beta2star/sum(beta2star)

beta3star = rep(0, n-1)
beta3star[floor(n/2)] = 1
beta3star[floor(n/3)] = 1
beta3star = beta3star/sum(beta3star)

#beta3star = rep(0, n-1)

p = 3


################################
y = X %*% beta1star + PX %*% beta2star + P2X %*% beta3star
## rmosek call 1 L1 minimization
## # Set up the program
 
  
beta1pos.index <- seq(1,(n-1))
beta1neg.index <- seq(1,(n-1)) + (n-1)

beta2pos.index <- seq(1,(n-1)) + 2*(n-1)
beta2neg.index <- seq(1,(n-1)) + 3*(n-1)

beta3pos.index <- seq(1,(n-1)) + 4*(n-1)
beta3neg.index <- seq(1,(n-1)) + 5*(n-1)

num.vars = 2*p*(n-1)
##############################################
noiseless.pattern <- list(sense = "min")
noiseless.pattern$c <- rep(1,num.vars)
##############################################
# Affine constraint 1: auxiliary variables [no cost]
# r = y - X beta
A1 <- Matrix(0, nrow = n, ncol = num.vars)
A1 <- cbind(X,-X,PX,-PX,P2X,-P2X)
A1 <- Matrix(A1)
###############################################
noiseless.pattern$A <- A1
noiseless.pattern$bc <- rbind(
  blc = t(y),
  buc = t(y)
)

#constraints on the program variables
noiseless.pattern$bx <- rbind(
  blx = c(rep(0, num.vars)),
  bux = c(rep(Inf, num.vars))
)
##################################
r <- mosek(noiseless.pattern)
sol1 <- r$sol$itr$xx
obj1 <- sum(sol1)


beta1.index <- seq(1,(n-1))
beta2.index <- seq(1,(n-1)) + (n-1)
beta3.index <- seq(1,(n-1)) + 2*(n-1)

num.vars = p*(n-1)
##############################################
noiseless.pattern <- list(sense = "min")
noiseless.pattern$c <- rep(1, num.vars)
##############################################
# Affine constraint 1: auxiliary variables [no cost]
# r = y - X beta
A1 <- Matrix(0, nrow = n, ncol = num.vars)
A1 <- cbind(X,PX,P2X)
A1 <- Matrix(A1)
###############################################
noiseless.pattern$A <- A1
noiseless.pattern$bc <- rbind(
  blc = t(y),
  buc = t(y)
)
#constraints on the program variables
noiseless.pattern$bx <- rbind(
  blx = c(rep(0, num.vars)),
  bux = c(rep(Inf, num.vars))
)
##################################
r <- mosek(noiseless.pattern)
sol2 <- r$sol$itr$xx
obj2 <- sum(sol2)
#####################################
print(c(obj1, obj2, obj1 - obj2))
#print(sol1[beta1neg.index] > 1e-5)
#print(sol1[beta2neg.index] > 1e-5)

beta1a = sol1[beta1pos.index] - sol1[beta1neg.index]
beta1b = sol2[beta1.index]
beta2a = sol1[beta2pos.index] - sol1[beta2neg.index]
beta2b = sol2[beta2.index]

beta3a = sol1[beta3pos.index] - sol1[beta3neg.index]
beta3b = sol2[beta3.index]

ya = X %*% beta1a + PX %*% beta2a + P2X %*% beta3a
yb =  X %*% beta1b + PX %*% beta2b + P2X %*% beta3b

#print(round(cbind(beta1a, beta2a, beta1b, beta2b), 3))
