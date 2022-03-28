# F3: explanation plot for variance bars
library("colorspace") 
source("tools.r")

q = c(.4,.3,0.3) # This is the real q
K=length(q)
set.seed(121)

pdf(file = "f3.pdf", width=7, height=7, family="Times", onefile=FALSE)
par(ann=TRUE, cex=1.4, mar = c(0,0,0,0))
plot(c(1, 4), c(0,5), type="n", axes=FALSE, main="", xlab="", ylab="")

M=500
N=c(10, 10, 10) # These are diploid individuals

comp=TRUE

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

  save(q.hat, var, file = "f3")
}
load("f3")

text(2.5, 3.8, "means that the confidence region is within", cex=1.5)
plotq(x=1, y=4.4, q.hat, height=.5, width=3, cols=rainbow(3), var=var) 
plotq(x=1, y=.5, q.hat-var[2,], height=.5, width=3, cols=rainbow(3), var=NULL) 
plotq(x=1, y=1.25, q.hat+var[2,], height=.5, width=3, cols=rainbow(3), var=NULL) 
plotq(x=1, y=2, q.hat-var[1,], height=.5, width=3, cols=rainbow(3), var=NULL) 
plotq(x=1, y=2.75, q.hat+var[1,], height=.5, width=3, cols=rainbow(3), var=NULL) 
for(i in 1:2) {
  lines(1+3*c(cumsum(q.hat+var[1,])[i],cumsum(q.hat+var[1,])[i]), c(4.9, 2.75), lty=3)
  lines(1+3*c(cumsum(q.hat-var[1,])[i],cumsum(q.hat-var[1,])[i]), c(4.9, 2), lty=3)
  lines(1+3*c(cumsum(q.hat+var[2,])[i],cumsum(q.hat+var[2,])[i]), c(4.9, 1.25), lty=3)
  lines(1+3*c(cumsum(q.hat-var[2,])[i],cumsum(q.hat-var[2,])[i]), c(4.9, .5), lty=3)
}

dev.off()

