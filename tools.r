# Throughout, we assume bi-allelic data

bootstrap.sample <- function(test.sample, training.sample, training.samLoc, lim="MN") {
  locM = 1:ncol(training.sample)
  locN = 1:nrow(training.sample)
  if(lim=="M" || lim=="MN") locM = sample(ncol(training.sample), ncol(training.sample), replace = TRUE)
  if(lim=="N" || lim=="MN") locN = sample(nrow(training.sample), nrow(training.sample), replace = TRUE)
  list(test.sample=test.sample[locM], training.sample = training.sample[locN,locM], training.samLoc = training.samLoc[locN])
}

# compute frequencies from sample along samLoc
getfreqs = function(sample, samLoc){
  t(simplify2array(by(sample,samLoc,colMeans,na.rm=TRUE))) / 2
}

# sample has nss columns and ninds rows; entries are 0, 1, 2
# freqs has npops rows and nss columns
# countsAncestral is a vector of length npops
# value is a matrix admixture.proportions with nsam.mixed rows and npops columns;
# admixture.proportion[i,k] gives the proportion of the genome of individual i coming from population k

# Here, x is a vector

get.admixtureproportions<-function(x, p, tol=1e-6, verbose = FALSE) {
  p = p[,!is.na(x)]
  x = x[!is.na(x)]
  
  K = nrow(p)
  M = ncol(p)

  
  # initialize admixture proportions
  res = as.vector(rdirichlet(1,K))
  #  res = rep(1/K, K)
  names(res) = rownames(p)

  j=1
  err = 1
  while(err>tol) {
    j = j+1
    loc = fun2(res, p, x)
    err = sum(abs(res - loc))
    res = loc
  }
  if(verbose) cat("\t Iterations: ", j, "\n")
  res
}

fun2<-function(q, p, loc.x) {
  K = length(q)
  M = length(loc.x)  
  E = matrix(0, nrow=K, ncol = M)

  # admixture.proportions = admixture.proportions / sum(admixture.proportions)

  loc = q %*% p # has M cols
  loc[loc==0] = 1e-16
  loc[loc==1] = 1-1e-16
  for(k in 1:K) {
    E[k,] = (loc.x * p[k,] / loc + (2 - loc.x) * (1-p[k,])/(1 - loc))
  }
  res = rowSums(E)/M * q / 2
  res/(sum(res))
}

get.loglik.admixture<-function(loc.x, p, q, sum=TRUE) {
  loc = t(p) %*% q
  if(sum) res = sum(log(choose(2, loc.x) * loc^loc.x * (1-loc)^(2-loc.x)))
  else res = log(choose(2, loc.x) * loc^loc.x * (1-loc)^(2-loc.x))
  res
}

rdirichlet<-function(n, classes){
  res = matrix(rexp(classes*n), ncol=classes, nrow=n)
  res / rowSums(res)
}

bound<-function(x, lower=0, upper=1) {
  res = x
  for(i in 1:length(x)) {
    res[i] = min(max(x[i],lower), upper)
  }
  res
}

# var is a matrix containing the eigenvalues to be drawn
plotq<-function(x, y, ia, height, width, cols, var=NULL) {
  npops = length(ia)
  rect(x, y, x+width*ia[1], y+height, col = cols[1], border = NA)
  for(i in 2:npops) {
    rect(x+width*cumsum(ia)[i-1], y, x+width*cumsum(ia)[i], y+height, col = cols[i], border = NA)
  }
  if(length(var)) {
    dy = height / (nrow(var)+1)
    for(i in 1:nrow(var)) {
      if(var[i,1]>0) var[i,]=-var[i,]
      loc1 = ia + var[i,]
      loc1 = loc1/sum(loc1)
      loc2 = ia - var[i,]
      loc2 = loc2/sum(loc2)
      eps = rnorm(1, mean=0, sd = 0.1*dy)
      if(abs(loc1[1]-loc2[1])>0.0001) {
        lines(bound(c(x+width*ia[1], x+width*loc2[1]), lower=x, upper=x+width), c(y+height-i*dy+eps, y+height-i*dy+eps), lwd=2, col="darkgrey")
        lines(bound(c(x+width*ia[1], x+width*loc1[1]), lower=x, upper=x+width), c(y+height-i*dy+eps, y+height-i*dy+eps), lwd=2, col="black")
	lines(bound(c(x+width*loc1[1], x+width*loc1[1]), lower=x, upper=x+width), c(y+height-i*dy+eps-.2*dy, y+height-i*dy+eps+.2*dy), lwd=2, col="black")
	lines(bound(c(x+width*loc2[1], x+width*loc2[1]), lower=x, upper=x+width), c(y+height-i*dy+eps-.2*dy, y+height-i*dy+eps+.2*dy), lwd=2, col="darkgrey")
      }	
      for(j in 2:(npops-1)) {
        eps = rnorm(1, mean=0, sd = 0.1*dy)
        if(abs(cumsum(loc1)[j]-cumsum(loc2)[j])>0.0001) {
          lines(bound(c(x+width*cumsum(ia)[j], x+width*cumsum(loc2)[j]), lower=x, upper=x+width), c(y+height-i*dy+eps, y+height-i*dy+eps), lwd=2, col="darkgrey")
          lines(bound(c(x+width*cumsum(ia)[j], x+width*cumsum(loc1)[j]), lower=x, upper=x+width), c(y+height-i*dy+eps, y+height-i*dy+eps), lwd=2, col="black")
          lines(bound(c(x+width*cumsum(loc1)[j], x+width*cumsum(loc1)[j]), lower=x, upper=x+width), c(y+height-i*dy+eps-.2*dy, y+height-i*dy+eps+.2*dy), lwd=2, col="black")
          lines(bound(c(x+width*cumsum(loc2)[j], x+width*cumsum(loc2)[j]), lower=x, upper=x+width), c(y+height-i*dy+eps-.2*dy, y+height-i*dy+eps+.2*dy), lwd=2, col="darkgrey")
	}
      }
    }
  }
  rect(x, y, x+width, y+height)
}

# Now we come to all matrices from the manuscript

DeltaTilde<-function(p,q,x=NULL) { # p and q are vectors of length K, x \in {0,1,2}
  K = length(q)
  loc = (q %*% p)[1,]
  res = 0
  if(!is.null(x)) {
    if(x > 0) res = res + x * (diag(rep(1,K))/loc - (p %*% t(q)) / loc^2)
    if(x < 2) res = res - (2-x) * (diag(rep(1,K))/(1-loc) - ((1-p) %*% t(q)) / (1-loc)^2)
  } else {
    res = 2 * (((1-p) %*% t(q)) / (1-loc) - (p %*% t(q)) / loc)
  }
  res
}

Pi<-function(p,r) { # p and r are vectors of length K
  diag(p*(1-p)/r)
}

Gamma<-function(p,q,r,x=NULL){ # here, p is a matrix, and x is a vector
  K = length(q)
  res = matrix(0, nrow=K, ncol=K)
  if(!is.null(x)) {
    p = p[,!is.na(x)]
    x = x[!is.na(x)]
    M = ncol(p)
    for(i in 1:M) {
      loc = DeltaTilde(p[,i], q, x[i]) %*% Pi(p[,i], r) %*% t(DeltaTilde(p[,i], q, x[i]))
      if(is.nan(loc[1,1])) cat("NaN produced for ", colnames(p)[i], "\n")
      res = res + loc
    }
  } else {
    M = ncol(p)
    for(i in 1:M) {
      loc = DeltaTilde(p[,i], q) %*% Pi(p[,i], r) %*% t(DeltaTilde(p[,i], q))
      if(is.nan(loc[1,1])) cat("NaN produced for ", colnames(p)[i], "\n")
      res = res + loc
    }

  }
  res/M/2
}

# if x==NULL, return Sigma from the manuscript, otherwise, return Sigma_x

Sigma<-function(p, q, x=NULL) {
  K = length(q)
  res = matrix(0, nrow=K, ncol=K)
  if(!is.null(x)) {
    p = p[,!is.na(x)]
    x = x[!is.na(x)]
    M = ncol(p)
    for(i in 1:M) {
      loc = (q %*% p[,i])[1,1]
      if(x[i] > 0) res = res + x[i] * p[,i]%*%t(p[,i]) / loc^2
      if(x[i] < 2) res = res + (2-x[i]) * (1-p[,i])%*%t(1-p[,i]) / (1-loc)^2 
    }
  } else {
    M = ncol(p)
    for(i in 1:M) {
      loc = (q %*% p[,i])[1,1]
      if(loc>0) res = res + 2 * p[,i]%*%t(p[,i]) / loc
      if(loc<1) res = res + 2 * (1-p[,i])%*%t(1-p[,i]) / (1-loc)
    }
  }
  solve(res/M/2)
}

# N is a vector of length K for the size of the reference database
varq<-function(p,q,N,x=NULL,lim="MN"){ 
  if(lim!="M" && lim!="N" && lim!="MN") {
    cat("lim must be M, N or MN.\n")
    exit()
  }
  r = N/sum(N)
  N = sum(N)
  K = nrow(p)
  M = ncol(p)
  res = matrix(0, ncol=K, nrow=K)
  if(lim == "M" || lim == "MN") {
    res = res + (Sigma(p, q, x) - q %*% t(q))/2/M
  }
  if(lim == "N" || lim == "MN") {
    res = res + Sigma(p, q, x) %*% Gamma(p,q,r,x) %*% Sigma(p, q, x)/4/N/M
  }
  res
}
