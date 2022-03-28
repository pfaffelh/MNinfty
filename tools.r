# Throughout, we assume bi-allelic data

require("parallel")
require("vcfR")

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

# Here, sample is a vector

get.admixtureproportions.multi<-function(x, p, tol=1e-6, verbose = FALSE) {
  res = NULL
  for(i in 1:nrow(x)) {
    if(verbose) cat("Method: admixture \t Name: ", rownames(x)[i])
    res = rbind(res, get.admixtureproportions(x[i,], p, tol=tol, verbose=verbose))
  }
  colnames(res) = rownames(p)
  rownames(res) = rownames(x)
  res
}

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

coords<-function(vec)  {
   x = vec[1]
   y = vec[2]
   z = vec[3]
   x*c(0,sqrt(3)-1/2) + y*c(-1,-1/2) + z*c(1,-1/2)
}

getData<-function(filenameDataVCFGZ, skip=0, windowSize=-1, inds, diploid = TRUE) {
  if(diploid) {
    trans = c(0,1,1,2, NA)
    names(trans)=c("0|0", "0|1", "1|0", "1|1", "NA")
  } else {
    trans = c(0,1, NA)
    names(trans)=c("0", "1", "NA")
  }
  samLoc = inds[,2]
  
  a = read.vcfR(filenameDataVCFGZ, skip = skip, nrows = windowSize)
  genotypes = a@gt[,-1]
  row.names(genotypes) = a@fix[,"ID"]
  sample = matrix(0, ncol = nrow(genotypes), nrow = nrow(inds))
  rownames(sample) = as.character(inds[,1])
  colnames(sample) = as.character(a@fix[,"ID"])
  for(k in 1:nrow(inds)) {
    sample[k,] = as.numeric(trans[genotypes[,as.character(inds[k,1])]])
  }
  snps = colnames(sample) 
  for(i in c(which(snps == ""),which(is.na(snps)))) {
    colnames(sample)[i] = snps[i] = paste(filenameDataVCFGZ, "_", skip + i, sep="") # sometimes SNPs have no identifier
  }
  ref.allele = a@fix[,"REF"]
  alt.allele = a@fix[,"ALT"]
  names(ref.allele) = a@fix[,"ID"]
  
  res = list(sample = sample, samLoc = samLoc, ref.allele = ref.allele, alt.allele = alt.allele)
  res
}

rdirichlet<-function(n, classes){
  res = matrix(rexp(classes*n), ncol=classes, nrow=n)
  res / rowSums(res)
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
        lines(c(x+width*ia[1], x+width*loc2[1]), c(y+height-i*dy+eps, y+height-i*dy+eps), lwd=2, col="darkgrey")
        lines(c(x+width*ia[1], x+width*loc1[1]), c(y+height-i*dy+eps, y+height-i*dy+eps), lwd=2, col="black")
	lines(c(x+width*loc1[1], x+width*loc1[1]), c(y+height-i*dy+eps-.2*dy, y+height-i*dy+eps+.2*dy), lwd=2, col="black")
	lines(c(x+width*loc2[1], x+width*loc2[1]), c(y+height-i*dy+eps-.2*dy, y+height-i*dy+eps+.2*dy), lwd=2, col="darkgrey")
      }	
      for(j in 2:(npops-1)) {
        eps = rnorm(1, mean=0, sd = 0.1*dy)
        if(abs(cumsum(loc1)[j]-cumsum(loc2)[j])>0.0001) {
          lines(c(x+width*cumsum(ia)[j], x+width*cumsum(loc2)[j]), c(y+height-i*dy+eps, y+height-i*dy+eps), lwd=2, col="darkgrey")
          lines(c(x+width*cumsum(ia)[j], x+width*cumsum(loc1)[j]), c(y+height-i*dy+eps, y+height-i*dy+eps), lwd=2, col="black")
          lines(c(x+width*cumsum(loc1)[j], x+width*cumsum(loc1)[j]), c(y+height-i*dy+eps-.2*dy, y+height-i*dy+eps+.2*dy), lwd=2, col="black")
          lines(c(x+width*cumsum(loc2)[j], x+width*cumsum(loc2)[j]), c(y+height-i*dy+eps-.2*dy, y+height-i*dy+eps+.2*dy), lwd=2, col="darkgrey")
	}
      }
    }
  }
  rect(x, y, x+width, y+height)
}

plotq.multi.old<-function(x, y, ia, height, width, cols, alpha=NULL) { # here, ia has multiple rows and we use transparency
  if(is.null(alpha)) alpha = rep(1/nrow(ia), nrow(ia))
  cols.rgb = col2rgb(cols)/255
  for(i in 1:nrow(ia)) {
    cols.loc = cols
    for(j in 1:length(cols)) {
      cols.loc[j] = rgb(cols.rgb[1,j], cols.rgb[2,j], cols.rgb[3,j], alpha[i])
    }
    plotq(x, y, ia[i,], height, width, cols.loc)
  }
}

# here, ia has multiple rows and we use the first
plotq.multi<-function(x, y, ia, height, width, cols, steps = 1000, weights=NULL) { 
  if(is.null(weights)) weights = rep(1/nrow(ia), nrow(ia))
  cols.rgb = col2rgb(cols)/255
  npops = length(cols)
  res = NULL
  ia.sum = round(t(apply(ia, 1, cumsum)),3)
  for(i in 1:steps) {
    loc = t(apply(ia.sum, 1, function(x) min(npops,1+sum(x < i/steps))))
#    cat("i = ", i, ", loc = ", loc, "\n")
    res = rbind(res,t(round(cols.rgb[,loc] %*% weights,3)))
  }

  for(i in 1:steps) {
    cols.loc = rgb(res[i,1], res[i,2], res[i,3])
    rect(x+width*(i-1)/steps, y, x+width*i/steps, y+height, col = cols.loc, border = NA)
  }
  rect(x, y, x+width, y+height)
}




# if x==NULL, return Sigma from the manuscript, otherwise, return Sigma_x

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
