# f7: display effects of finite M, N
source("tools.r")
library("colorspace") 

q = c(.4,.3,0.3) # This is the real q
K=length(q)
set.seed(12)

M=c(10, 20, 50, 100, 200, 500)
N=c(10, 20, 50, 100, 200, 500, 1000, 2000) # These are diploid individuals

p = matrix(0, nrow = K, ncol = max(M)) # p[k,m] is the allele frequency in population k for locus m
for(m in 1:max(M)) {
   p[1,m] = rbeta(1,5,5)
   p[2,m] = rbeta(1,5,5)
   p[3,m] = rbeta(1,5,5)
}
rownames(p) = c("A", "B", "C")

# produce sample and estimates
x = 0 * 1:max(M)
loc = q %*% p # these are the success probabilities
  for(m in 1:max(M)) {
  x[m] = rbinom(1, 2, loc[m])
}

# compute total variance 
M.vec = N.vec = res.vec = NULL
for(i in 1:length(M)) {
  for(j in 1:length(N)) {
    M.vec = c(M.vec, as.character(M[i]))
    N.vec = c(N.vec, as.character(N[j]))
    loc = varq(p[,1:M[i]], q, rep(N[j],3), x=x[1:M[i]], lim="MN")
    res.vec = c(res.vec, sqrt(loc[1,1] + loc[2,2] + loc[3,3]))
#    res.vec = c(res.vec, i+j)
  }
}

var = data.frame(x=N.vec, y=M.vec, z=res.vec)
library(ggplot2)
ggplot(data = var, aes(x=x, y=y, fill = z)) + geom_tile() + xlab("N") + ylab("M") +
 scale_fill_gradient2(low = "white", high = "blue",
  limit = c(0,max(var$z)), space = "Lab", name=expression(sqrt(paste(V[q], "[||", hat(q), -q, , "||]")))) + theme_minimal() + ylim(c("10", "20", "50", "100", "200", "500")) + xlim(c("10", "20", "50", "100", "200", "500", "1000", "2000")) + theme(text = element_text(family = "Times"))+ theme(axis.text.x = element_text(angle = 90))

ggsave("f7a.pdf", width=4, height=4)

# compute total variance 
M.vec = N.vec = res.vec = NULL
for(i in 1:length(M)) {
  for(j in 1:length(N)) {
    M.vec = c(M.vec, as.character(M[i]))
    N.vec = c(N.vec, as.character(N[j]))
    loc1 = varq(p[,1:M[i]], q, rep(N[j],3), x=x[1:M[i]], lim="MN")
    loc2 = varq(p[,1:M[i]], q, rep(N[j],3), x=x[1:M[i]], lim="N")
    res.vec = c(res.vec, (sqrt(loc2[1,1] + loc2[2,2] + loc2[3,3]))/sqrt(loc1[1,1] + loc1[2,2] + loc1[3,3]))
#    res.vec = c(res.vec, i+j)
  }
}

var = data.frame(x=N.vec, y=M.vec, z=res.vec)
library(ggplot2)
ggplot(data = var, aes(x=x, y=y, fill = z)) + geom_tile() + xlab("N") + ylab("M") +
 scale_fill_gradient2(low = "white", high = "blue",
  limit = c(0,max(var$z)), space = "Lab", name="contribution \n of finite N") + theme_minimal() + ylim(c("10", "20", "50", "100", "200", "500")) + xlim(c("10", "20", "50", "100", "200", "500", "1000", "2000")) + theme(text = element_text(family = "Times"))+ theme(axis.text.x = element_text(angle = 90))

ggsave("f7b.pdf", width=5, height=4)

