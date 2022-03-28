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

identifier = c("LP6005592-DNA_H03", "LP6005619-DNA_C01", "LP6005442-DNA_D02", "SS6004477", "LP6005443-DNA_A03", "LP6005443-DNA_H03")
#identifier = identifier[1:3]
p = getfreqs(training.sample.x, training.samLoc)

K = 7
cols = rainbow(K)
pdf(file = "f5.pdf", width=12, height=5, family="Times", onefile=FALSE)
par(ann=TRUE, cex=1.4, mar = c(0,0,0,0))
plot(c(1.4, 13.5), c(-2,5), type="n", axes=FALSE, main="", xlab="", ylab="")

trans = c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)")
x.loc = test.sample.x # [identifier,]
sd = 2

identifier=c(17, 38, 4, 101, 89, 111)

#for(i in 1:nrow(x.loc)) {
for(i in 1:length(identifier)) {
  eps=0.01
  N = table(training.samLoc)
  x = x.loc[identifier[i],]
  q.hat = get.admixtureproportions(x, p, tol=1e-6, verbose = TRUE)
  N.loc = N[q.hat>eps]
  p.loc = p[q.hat>eps,]
  cols.loc = cols[q.hat>eps]
  q.hat = q.hat[q.hat>eps]
  q.hat = q.hat/sum(q.hat)
  K = length(q.hat)
  if(K>1) {
    var.q.hat = varq(p.loc, q.hat, N=N.loc, x=x, lim="MN")
    loc = eigen(var.q.hat)
    var = sd*sqrt(loc$values[-K])*t(loc$vectors[,-K])

    text(1.5, 4.2 -.7*i, trans[i])
    plotq(x=5, y=4-.7*i, q.hat, height=.4, width=3.5, cols=cols.loc, var=var)
    text(3, 4.2 - .7*i, rownames(test.sample)[identifier[i]])
    text(10, 4.2 - .7*i, paste(test.samLoc.full[identifier[i]], ", ", test.samLoc.pop[identifier[i]], sep=""))
  }
}

text(3, 4.5, "Identifier")
text(7, 4.5, "IA-estimate")
text(10, 4.5, "population")

pops = rownames(p)
for(i in 1:7) {
  text(1.15 +1.3*i, -1, pops[i], cex=1.3)
  rect(.8+1.3*i, -1.6, 1.5+1.3*i, -1.3, col=cols[i])
}

dev.off()



# Take home: EUR, MEA and SAS are hard to distinguish. Hence, variance
# often goes into these directions.


