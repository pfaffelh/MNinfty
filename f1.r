# f1: trial plot to display variance in a barplot
library("colorspace") 
source("tools.r")

q = c(.4,.3,0.3) # This is the real q
K=length(q)
set.seed(121)
iterations=1000

comp = TRUE

pdf(file = "f1.pdf", width=14, height=3, family="Times", onefile=FALSE)
par(ann=TRUE, cex=1.4, mar = c(0,0,0,0))
plot(c(0.5, 13.5), c(0,3), type="n", axes=FALSE, main="", xlab="", ylab="")

# round1

M=500
N=c(10, 10, 10) # These are diploid individuals

if(comp) {
  p = matrix(0, nrow = K, ncol = M) # p[k,m] is the allele frequency in population k for locus m
  for(m in 1:M) {
     p[1,m] = rbeta(1,5,5)
     p[2,m] = rbeta(1,5,5)
     p[3,m] = rbeta(1,5,5)
  }
  rownames(p) = c("A", "B", "C")

  # Compute FST between A and B
  FST = 0*(1:M)
  for(m in 1:M) {
    FST[m] = 1 - mean(p[1:2,m]*(1-p[1:2,m])) / (mean(p[1:2,m])*(1-mean(p[1:2,m])))
  }
  cat("The average FST between A and B is", mean(FST), ".\n")

  # create training.sample
  training.samLoc = c(rep("A", N[1]), rep("B", N[2]), rep("C", N[3]))
  training.sample = matrix(0, nrow = sum(N), ncol=M)
  for(i in 1:nrow(training.sample)) {
    for(m in 1:ncol(training.sample)) {
      training.sample[i,m] = rbinom(1, 2, p[training.samLoc[i],m])
    }
  }
  p.hat = getfreqs(training.sample, training.samLoc)

  # produce sample and estimates
  x = 0 * 1:M
  loc = q %*% p # these are the success probabilities
  for(m in 1:M) {
    x[m] = rbinom(1, 2, loc[m])
  }
  q.hat = get.admixtureproportions(x, p.hat, tol=1e-6, verbose = TRUE)
  var.q.hat = varq(p.hat,q.hat,N,x=x,lim="MN")
  loc = eigen(var.q.hat)
  var = 2*sqrt(loc$values[-K])*t(loc$vectors[,-K])

  # get bootstrap estimates
  
  q.hat.loc=matrix(0, ncol=K, nrow=iterations)
  for(i in 1:iterations) {
    loc = bootstrap.sample(x, training.sample, training.samLoc, lim="MN") 
    x.loc = loc$test.sample
    p.hat.loc = getfreqs(loc$training.sample, loc$training.samLoc)
    q.hat.loc[i,] = get.admixtureproportions(x.loc, p.hat.loc, tol=1e-6, verbose = TRUE)
  }

  # get covariance matrix from bootstrap
  var.q.hat.boot = matrix(0, nrow=K, ncol=K)
  for(k in 1:K) {
    for(l in 1:K) {
      var.q.hat.boot[k,l] = sum(q.hat.loc[,k]*q.hat.loc[,l])/iterations - mean(q.hat.loc[,k])*mean(q.hat.loc[,l])
    }
  }
  loc = eigen(var.q.hat.boot)
  var.boot = 2*sqrt(loc$values[-K])*t(loc$vectors[,-K])
  save(q.hat, var, var.boot, file = "f1.1")
}
load("f1.1")

text(3.5, 2.75, "M=500, N=30", cex=1.5)
text(1, .75, "bootstrap", cex=1.5)
text(1, 1.5, "Thm 1", cex=1.5)
text(1, 2.25, "true q", cex=1.5)
plotq(x=2, y=2, q, height=.5, width=3.5, cols=rainbow(3), var=NULL) 
plotq(x=2, y=1.25, q.hat, height=.5, width=3.5, cols=rainbow(3), var=var) 
plotq(x=2, y=.5, q.hat, height=.5, width=3.5, cols=rainbow(3), var=var.boot) 


# round2

M=5000
N=c(10, 10, 10) # These are diploid individuals

if(comp) {
  p = matrix(0, nrow = K, ncol = M) # p[k,m] is the allele frequency in population k for locus m
  for(m in 1:M) {
     p[1,m] = rbeta(1,5,5)
     p[2,m] = rbeta(1,5,5)
     p[3,m] = rbeta(1,5,5)
  }
  rownames(p) = c("A", "B", "C")

  # create training.sample
  training.samLoc = c(rep("A", N[1]), rep("B", N[2]), rep("C", N[3]))
  training.sample = matrix(0, nrow = sum(N), ncol=M)
  for(i in 1:nrow(training.sample)) {
    for(m in 1:ncol(training.sample)) {
      training.sample[i,m] = rbinom(1, 2, p[training.samLoc[i],m])
    }
  }
  p.hat = getfreqs(training.sample, training.samLoc)

  # produce sample and estimates
  x = 0 * 1:M
  loc = q %*% p # these are the success probabilities
  for(m in 1:M) {
    x[m] = rbinom(1, 2, loc[m])
  }
  q.hat = get.admixtureproportions(x, p.hat, tol=1e-6, verbose = TRUE)
  var.q.hat = varq(p.hat,q.hat,N,x=x,lim="MN")
  loc = eigen(var.q.hat)
  var = 2*sqrt(loc$values[-K])*t(loc$vectors[,-K])

  # get bootstrap estimates
  q.hat.loc=matrix(0, ncol=K, nrow=iterations)
  for(i in 1:iterations) {
    loc = bootstrap.sample(x, training.sample, training.samLoc, lim="MN") 
    x.loc = loc$test.sample
    p.hat.loc = getfreqs(loc$training.sample, loc$training.samLoc)
    q.hat.loc[i,] = get.admixtureproportions(x.loc, p.hat.loc, tol=1e-6, verbose = TRUE)
  }

  # get covariance matrix from bootstrap
  var.q.hat.boot = matrix(0, nrow=K, ncol=K)
  for(k in 1:K) {
    for(l in 1:K) {
      var.q.hat.boot[k,l] = sum(q.hat.loc[,k]*q.hat.loc[,l])/iterations - mean(q.hat.loc[,k])*mean(q.hat.loc[,l])
    }
  }
  loc = eigen(var.q.hat.boot)
  var.boot = 2*sqrt(loc$values[-K])*t(loc$vectors[,-K])
  save(q.hat, var, var.boot, file = "f1.2")
}
load("f1.2")

text(7.5, 2.75, "M=5000, N=30", cex=1.5)
#text(1, .75, "bootstrap")
#text(1, 1.5, "Thm 1")
#text(1, 2.25, "true q")
plotq(x=6, y=2, q, height=.5, width=3.5, cols=rainbow(3), var=NULL) 
plotq(x=6, y=1.25, q.hat, height=.5, width=3.5, cols=rainbow(3), var=var) 
plotq(x=6, y=.5, q.hat, height=.5, width=3.5, cols=rainbow(3), var=var.boot) 


# round3

M=500
N=c(1000, 1000, 1000) # These are diploid individuals

if(comp) {
  p = matrix(0, nrow = K, ncol = M) # p[k,m] is the allele frequency in population k for locus m
  for(m in 1:M) {
     p[1,m] = rbeta(1,5,5)
     p[2,m] = rbeta(1,5,5)
     p[3,m] = rbeta(1,5,5)
  }
  rownames(p) = c("A", "B", "C")

  # create training.sample
  training.samLoc = c(rep("A", N[1]), rep("B", N[2]), rep("C", N[3]))
  training.sample = matrix(0, nrow = sum(N), ncol=M)
  for(i in 1:nrow(training.sample)) {
    for(m in 1:ncol(training.sample)) {
      training.sample[i,m] = rbinom(1, 2, p[training.samLoc[i],m])
    }
  }
  p.hat = getfreqs(training.sample, training.samLoc)

  # produce sample and estimates
  x = 0 * 1:M
  loc = q %*% p # these are the success probabilities
  for(m in 1:M) {
    x[m] = rbinom(1, 2, loc[m])
  }
  q.hat = get.admixtureproportions(x, p.hat, tol=1e-6, verbose = TRUE)
  var.q.hat = varq(p.hat,q.hat,N,x=x,lim="MN")
  loc = eigen(var.q.hat)
  var = 2*sqrt(loc$values[-K])*t(loc$vectors[,-K])

  # get bootstrap estimates
  q.hat.loc=matrix(0, ncol=K, nrow=iterations)
  for(i in 1:iterations) {
    loc = bootstrap.sample(x, training.sample, training.samLoc, lim="MN") 
    x.loc = loc$test.sample
    p.hat.loc = getfreqs(loc$training.sample, loc$training.samLoc)
    q.hat.loc[i,] = get.admixtureproportions(x.loc, p.hat.loc, tol=1e-6, verbose = TRUE)
  }

  # get covariance matrix from bootstrap
  var.q.hat.boot = matrix(0, nrow=K, ncol=K)
  for(k in 1:K) {
    for(l in 1:K) {
      var.q.hat.boot[k,l] = sum(q.hat.loc[,k]*q.hat.loc[,l])/iterations - mean(q.hat.loc[,k])*mean(q.hat.loc[,l])
    }
  }
  loc = eigen(var.q.hat.boot)
  var.boot = 2*sqrt(loc$values[-K])*t(loc$vectors[,-K])
  save(q.hat, var, var.boot, file = "f1.3")
}
load("f1.3")

text(11.5, 2.75, "M=500, N=3000", cex=1.5)
#text(1, .75, "bootstrap")
#text(1, 1.5, "Thm 1")
#text(1, 2.25, "true q")
plotq(x=10, y=2, q, height=.5, width=3.5, cols=rainbow(3), var=NULL) 
plotq(x=10, y=1.25, q.hat, height=.5, width=3.5, cols=rainbow(3), var=var) 
plotq(x=10, y=.5, q.hat, height=.5, width=3.5, cols=rainbow(3), var=var.boot) 

dev.off()

