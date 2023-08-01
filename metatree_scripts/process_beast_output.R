# make pinniped supertree time trees 

setwd("Pinniped_metatree/TimeTree/BEAST_run2/")
x <- read.table("pinniped_combined.log", header=1)
trees <- read.nexus("pinniped_combined.trees")
write.tree(trees[[which.max(x$posterior)]], file="pinniped_map.tre")


library(treespace)


KCdist05 <- multiDist(trees, lambda=0.5)

d<-matrix(data=NA, nrow=length(trees), ncol=length(trees))
d[lower.tri(d, diag = F)] <- KCdist05
d<-t(d)
d[lower.tri(d)] <- t(d)[(lower.tri(d))]
diag(d) <- rep(0,length(trees))
median.tree05<-trees[[which.min(colSums(d))]] #378
write.tree(median.tree05, "median_tree.tre")

