# F3: explanation plot for variance bars, and triangle plots
library("colorspace") 
source("tools.r")

q = c(.4,.3,0.3) # This is the real q
K=length(q)
set.seed(121)

pdf(file = "f3.pdf", width=12, height=7, family="Times", onefile=FALSE)
par(ann=TRUE, cex=1.4, mar = c(0,0,0,0))
plot(c(0, 9), c(-.5,5), type="n", axes=FALSE, main="", xlab="", ylab="")

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

cols=rainbow(3)

text(2.5, 3.8, "means that the confidence region is within", cex=1.5)
plotq(x=1, y=4.4, q.hat, height=.5, width=3, cols=cols, var=var) 
plotq(x=1, y=.5, q.hat-var[2,], height=.5, width=3, cols=cols, var=NULL) 
plotq(x=1, y=1.25, q.hat+var[2,], height=.5, width=3, cols=cols, var=NULL) 
plotq(x=1, y=2, q.hat-var[1,], height=.5, width=3, cols=cols, var=NULL) 
plotq(x=1, y=2.75, q.hat+var[1,], height=.5, width=3, cols=cols, var=NULL) 
for(i in 1:2) {
  lines(1+3*c(cumsum(q.hat+var[1,])[i],cumsum(q.hat+var[1,])[i]), c(4.9, 2.75), lty=3)
  lines(1+3*c(cumsum(q.hat-var[1,])[i],cumsum(q.hat-var[1,])[i]), c(4.9, 2), lty=3)
  lines(1+3*c(cumsum(q.hat+var[2,])[i],cumsum(q.hat+var[2,])[i]), c(4.9, 1.25), lty=3)
  lines(1+3*c(cumsum(q.hat-var[2,])[i],cumsum(q.hat-var[2,])[i]), c(4.9, .5), lty=3)
}

#dev.off()

#pdf(file = "f3a.pdf", width=7, height=8, family="Times", onefile=FALSE)
#par(ann=TRUE, cex=1.4, mar = c(0,0,0,0))
#plot(c(-1.65, 1.3), c(-1.9, 2.5), type="n", axes=FALSE, main="", xlab="", ylab="")

xoffset = 7
yoffset = 2
lines(xoffset+c(0, 1), yoffset + c(sqrt(3)-1/2, -1/2))
lines(xoffset+c(0, -1), yoffset + c(sqrt(3)-1/2, -1/2))
lines(xoffset+c(-1, 1), yoffset + c(-1/2, -1/2))

points(xoffset-0, yoffset + 1.3, col=cols[1], cex=4, pch=20)
points(xoffset-1, yoffset -.5, col=cols[2], cex=4, pch=20)
points(xoffset+1, yoffset -.5, col=cols[3], cex=4, pch=20)
#text(0, 1.4, "A", col="white", cex=2)
#text(-1.2, -1/2, "B", col="white", cex=2)
#text(1.2, -1/2, "C", cex=2)

co = coords(q.hat) 
co.true = coords(q) 

points(xoffset+co[1], yoffset + co[2], pch=20, col = cols[1], cex=2)
#points(co.true[1], co.true[2], pch=20, col = cols[1], cex=2)

co1 = coords(q.hat)
co2 = coords(q.hat+var[1,])
co3 = coords(q.hat+var[2,])
co4 = coords(q.hat-var[1,])
co5 = coords(q.hat-var[2,])
  
arrows(xoffset+co1[1], yoffset + co1[2], xoffset+co2[1], yoffset + co2[2], length=0.15, cex=3, lwd=2)
arrows(xoffset+co1[1], yoffset + co1[2], xoffset+co3[1], yoffset + co3[2], length=0.15, cex=3, lwd=2)
arrows(xoffset+co1[1], yoffset + co1[2], xoffset+co4[1], yoffset + co4[2], col="black", length=0.15, cex=3, lwd=2)
arrows(xoffset+co1[1], yoffset + co1[2], xoffset+co5[1], yoffset + co5[2], length=0.15, col="black", cex=3, lwd=2)

library("plotrix")

draw.ellipse(xoffset+co1[1], yoffset + co1[2], 1.5*sqrt(sum(co2-co1)^2), .82*sqrt(sum(co3-co1)^2), angle = -15)
text(7, 1, "confidence region", cex=1.5)
lines(c(7.1, 7.1), c(1.2, 1.9), lty=2)

lines(c(1, 4), c(0, 0))
lines(c(1, 1), c(-.1, 0.1))
lines(c(4, 4), c(-.1, 0.1))
text(1, -.3, "0")
text(4, -.3, "1")

lines(c(4, xoffset+co2[1]), c(3, yoffset + co2[2]), lty=3)
lines(c(4, xoffset+co3[1]), c(1.5, yoffset + co3[2]), lty=3)
lines(c(4, xoffset+co4[1]), c(2.25, yoffset + co4[2]), lty=3)
lines(c(4, xoffset+co5[1]), c(.75, yoffset + co5[2]), lty=3)

dev.off()