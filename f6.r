source("tools.r")

# Define the relevant populations in the training sample
# "AFR" "AMR" "EAS" "EUR" "MEA" "OCE" "SAS"
cols=c("red", "white", "green", "blue", "lightblue", "orange", "yellow")
names(cols) = c("AFR", "AMR", "EAS", "EUR", "MEA", "OCE", "SAS")

# Abbreviate superpops in all samples
trans = c("AFR", "EUR", "SAS", "EAS", "MEA", "OCE", "AMR", "SIB", "NEU", "EEU", "CAS", "ROM", "NEA", "SEA")
names(trans) = c("African", "European", "South Asian", "East Asian", "Middle East", "Oceanian", "American", "Siberian", "Northen European", "Eastern European", "Central Asian", "Roma", "Northeast Asian", "Southeast Asian")


# load the test sample from one sheet of the Snipper-xls-file.
a = read.csv("Verogen ForenSeq DNA Signature Kit 56-AIM SGDP Samples.csv", header = FALSE)
test.sample = as.matrix(a[2:nrow(a), 5:60])
rownames(test.sample) = as.vector(a[2:nrow(a),4])
for(i in 1:ncol(test.sample)) colnames(test.sample)[i] = as.vector((a[1,5:60])[1,i])
test.sample[test.sample == "NN"] = NA # some fields contain NN, which means NA
test.samLoc.full = as.vector(a[2:nrow(a),2])
test.samLoc.pop = as.vector(a[2:nrow(a),3])
test.samLoc = trans[test.samLoc.full]
names(test.samLoc) = names(test.samLoc.full)

# now load the training data from the other sheet of the Snipper-xls-file
a = read.csv("Verogen ForenSeq DNA Signature Kit 56-AIM Training Set Reference.csv", header = FALSE)
training.sample = as.matrix(a[6:nrow(a), 4:(ncol(a)-1)])
rownames(training.sample) = as.vector(a[6:nrow(a),3])
for(i in 1:ncol(training.sample)) colnames(training.sample)[i] = as.vector((a[1,4:(ncol(a)-1)])[1,i])
training.sample[training.sample == "NN"] = NA # some fields contain NN, which means NA
training.samLoc.full = as.vector(a[6:nrow(a),2])
training.samLoc = trans[training.samLoc.full]
names(training.samLoc) = names(training.samLoc.full)

# check if all SNPs are bi-allelic (in the full sample)
full.sample = rbind(training.sample, test.sample)
full.samLoc = c(training.samLoc, test.samLoc)

for(AIM in colnames(full.sample)) {
  dat = as.vector(full.sample[,AIM])
  dat = dat[!is.na(dat)]
  dat = unlist(strsplit(dat, ","))
  dat = unlist(strsplit(dat, ""))
  if(length(unique(dat))!=2) {
    cat("Warning: ", AIM, " is not bi-allelic. Present alleles are ")
    for(i in unique(dat)) cat(i, " ")
    cat("\n")
  }  
}

# Determine reference alleles (by first occurrence in the training set)
ref.allele = training.sample[1,]
for(j in 1:ncol(training.sample)) {
  loc = unlist(strsplit(training.sample[1,j],""))
  ref.allele[j] = loc[1]
}

# translate training.sample and test.sample to number of occurrences of reference allele 
training.sample.x = matrix(0, ncol=ncol(training.sample), nrow=nrow(training.sample))
colnames(training.sample.x)=colnames(training.sample)
rownames(training.sample.x)=rownames(training.sample)
test.sample.x = matrix(0, ncol=ncol(test.sample), nrow=nrow(test.sample))
colnames(test.sample.x)=colnames(test.sample)
rownames(test.sample.x)=rownames(test.sample)

for(AIM in colnames(training.sample)) {
  for(i in 1:nrow(training.sample)) {
    if(!is.na(training.sample[i, AIM])) {
      loc = as.character(training.sample[i, AIM])
      loc = unlist(strsplit(loc, ","))
      loc = unlist(strsplit(loc, ""))
      training.sample.x[i, AIM] = as.numeric(sum(loc == ref.allele[AIM]))
    }
  }
  for(i in 1:nrow(test.sample)) {
    if(!is.na(test.sample[i, AIM])) {
      loc = as.character(test.sample[i, AIM])
      loc = unlist(strsplit(loc, ","))
      loc = unlist(strsplit(loc, ""))
      test.sample.x[i, AIM] = sum(loc == ref.allele[AIM])
      test.sample.x[i, AIM] = as.numeric(test.sample.x[i, AIM])
    }
  }
}

# generate a matrix for distinguishing two populations

M = ncol(training.sample.x)
p = getfreqs(training.sample.x, training.samLoc)
K = 7

pops = sort(unique(training.samLoc))

x=NULL
y=NULL
res=NULL
for(k in 1:(K-1)) {
  for(l in (k+1):K) {
    x = c(x, pops[k])
    y = c(y, pops[l])
    # set q to mixture of i and k
    t = q = 0*(1:K)
    q[k]=q[l] = 0.5
    t[k]=1/sqrt(2)
    t[l]=-1/sqrt(2)
#    res = c(res, (t %*% (Sigma(p,q) - q %*% t(q)) %*% t)[1,1]/2/M)
#    res[k,l] = (Sigma(p,q) - q %*% t(q))[k,l]
    res = c(res, (Sigma(p,q) - q %*% t(q))[k,l]/2/M)
  }
}

covariances = data.frame(x=x, y=y, cov=res)
library(ggplot2)
ggplot(data = covariances, aes(x=x, y=y, fill = cov)) + geom_tile() + xlab("") + ylab("") + scale_fill_gradient2(low = "blue", high = "white", limit = c(min(res), max(res)), space = "Lab", name="-Covariance") + theme_minimal() + coord_fixed() + theme(text = element_text(family = "Times")) + ylim(rev(unique(covariances$y))) 

ggsave("f6.pdf", width=4, height=4)

