# f4: sample outside the reference database
source("tools.r")

K = 4 # The real world has 4 populations, but the reference database only contains the first 3
N=100 # This is the number of diploids per population in the training data
M=500
set.seed(2)
iterations = 4

cols = rainbow(3)
p = matrix(0, nrow = K, ncol = M) # p[k,m] is the allele frequency in population k for locus m
for(k in 1:max(K)) {
  for(m in 1:M) {
    p[k,m] = rbeta(1,2,2)
  }
}

pdf(file = "f4.pdf", width=14, height=4, family="Times", onefile=FALSE)
par(ann=TRUE, cex=1.4, mar = c(0,0,0,0))
plot(c(0.5, 14), c(2,6), type="n", axes=FALSE, main="", xlab="", ylab="")

# sample 1-4
plotq(x=2, y=5, c(1/3, 1/3, 1/3), height=.5, width=3.5, cols=cols, var=NULL) 
text(1, 5.25, "true q", cex=1.5)
for(i in 1:iterations) {
  q = c(1/3, 1/3, 1/3,0)
  x = 0 * 1:M
  loc = q %*% p # these are the success probabilities
  for(m in 1:M) {
    x[m] = rbinom(1, 2, loc[m])
  }
  p.loc = p[1:3,]
  K.loc = 3
  # effects of finite N
  p.hat = p.loc
  for(k in 1:K.loc) {
    for(m in 1:M) {
      p.hat[k,m] = rbinom(1, 2*N, p.loc[k,m])/2/N
    }
  }
  q.hat = get.admixtureproportions(x, p.hat, tol=1e-6, verbose = TRUE)
  var.q.hat = varq(p.hat, q.hat, N=rep(N, K.loc), x=x, lim="MN")
  loc = eigen(var.q.hat)
  var = 2*sqrt(loc$values[-K.loc])*t(loc$vectors[,-K.loc])

#  text(1, length(K)+1-.75*i, paste("K=", K[i]))
  plotq(x=2, y=iterations - .75*i + 1, q.hat, height=.5, width=3.5, cols=cols, var=var) 
  text(1, iterations - .75*i + 1.25, paste("Ind.", i), cex=1.5)
  cat("Log-likelihood: ", get.loglik.admixture(x, p.hat, q.hat, sum=TRUE), "\n")
}

# sample 5-8
rect(6, 5, 9.5, 5.5)
text(10, 6, "outside reference", cex=1.5)
lines(c(9,7.5), c(5.75,5.25), lty=2)
lines(c(11,12.5), c(5.75,5.25), lty=2)

for(i in 1:iterations) {
  q = c(0,0,0,1)
  x = 0 * 1:M
  loc = q %*% p # these are the success probabilities
  for(m in 1:M) {
    x[m] = rbinom(1, 2, loc[m])
  }
  p.loc = p[1:3,]
  K.loc = 3
  # effects of finite N
  p.hat = p.loc
  for(k in 1:K.loc) {
    for(m in 1:M) {
      p.hat[k,m] = rbinom(1, 2*N, p.loc[k,m])/2/N
    }
  }
  q.hat = get.admixtureproportions(x, p.hat, tol=1e-6, verbose = TRUE)
  var.q.hat = varq(p.hat, q.hat, N=rep(N, K.loc), x=x, lim="MN")
  loc = eigen(var.q.hat)
  var = 2*sqrt(loc$values[-K.loc])*t(loc$vectors[,-K.loc])

  plotq(x=6, y=iterations - .75*i + 1, q.hat, height=.5, width=3.5, cols=cols, var=var)
  cat("Log-likelihood: ", get.loglik.admixture(x, p.hat, q.hat, sum=TRUE), "\n")

}

# sample 9-12
rect(10, 5, 13.5, 5.5)
rect(10, 5, 11.75, 5.5, col = cols[1], border = NA)

for(i in 1:iterations) {
  q = c(0.5,0,0,.5)
  x = 0 * 1:M
  loc = q %*% p # these are the success probabilities
  for(m in 1:M) {
    x[m] = rbinom(1, 2, loc[m])
  }
  p.loc = p[1:3,]
  K.loc = 3
  # effects of finite N
  p.hat = p.loc
  for(k in 1:K.loc) {
    for(m in 1:M) {
      p.hat[k,m] = rbinom(1, 2*N, p.loc[k,m])/2/N
    }
  }
  q.hat = get.admixtureproportions(x, p.hat, tol=1e-6, verbose = TRUE)
  var.q.hat = varq(p.hat, q.hat, N=rep(N, K.loc), x=x, lim="MN")
  loc = eigen(var.q.hat)
  var = 2*sqrt(loc$values[-K.loc])*t(loc$vectors[,-K.loc])

#  text(1, length(K)+1-.75*i, paste("K=", K[i]))
  plotq(x=10, y=iterations - .75*i + 1, q.hat, height=.5, width=3.5, cols=cols, var=var) 
  cat("Log-likelihood: ", get.loglik.admixture(x, p.hat, q.hat, sum=TRUE) , "\n")

}

# mean difference of log-likelihoods between left and middle: 90.70412
# mean difference of log-likelihoods between left and right:  28.73852

dev.off()

