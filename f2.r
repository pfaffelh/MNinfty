# f2: effects of a large number of reference populations
source("tools.r")

q = c(.5,.5) # This is the real q, add 0s 
K=c(2, 3, 5, 10)
cols = rainbow(max(K))

N=10 # This is the number of diploids per population in the training data
M=500

set.seed(11)

p = matrix(0, nrow = max(K), ncol = M) # p[k,m] is the allele frequency in population k for locus m
for(k in 1:max(K)) {
  for(m in 1:M) {
    p[k,m] = rbeta(1,2,2)
  }
}

# create sample
x = 0 * 1:M
loc = q %*% p[1:2,] # these are the success probabilities
for(m in 1:M) {
  x[m] = rbinom(1, 2, loc[m])
}

pdf(file = "f2.pdf", width=7, height=2.5, family="Times", onefile=FALSE)
par(ann=TRUE, cex=1.4, mar = c(0,0,0,0))
plot(c(1.4, 5.5), c(2.2,4.5), type="n", axes=FALSE, main="", xlab="", ylab="")

# iterate all numbers of reference populations
for(i in 1:length(K)) {
  q.loc = c(q, rep(0, K[i]-2))
  p.loc = p[1:K[i],]
  p.hat = p.loc
  for(k in 1:K[i]) {
    for(m in 1:M) {
      p.hat[k,m] = rbinom(1, 2*N, p.loc[k,m])/2/N
    }
  }
  q.hat = get.admixtureproportions(x, p.hat, tol=1e-6, verbose = TRUE)
  var.q.hat = varq(p.hat, q.hat, N=rep(N, K[i]), x=x, lim="MN")
  loc = eigen(var.q.hat)
  var = 2*sqrt(loc$values[-K[i]])*t(loc$vectors[,-K[i]])

  text(1.5, length(K)+1-.6*i, paste("K=", K[i]), cex=1.5)
  plotq(x=2, y=length(K)+.75-.6*i, q.hat, height=.4, width=3, cols=cols[1:K[i]], var=var) 
}

dev.off()

