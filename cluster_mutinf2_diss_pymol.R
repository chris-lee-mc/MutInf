library(marray)
library(cluster)
prec <- read.table("3FJQ_cart_hotatoms.reslist-nsims6-structs3332-bin20_bootstrap_avg_dist_stdev_matrix.txt")
nrow(prec)
ncol(prec)
blah = hclust(dist(prec,method="euclidean"),method="ward")
blah2 = cutree(blah, k=11)
blah3=as.table(blah2)
blah4=data.frame(blah3)
blah4
write.table(blah4, "output_clusters.txt", quote=FALSE, sep='\t')

#END
