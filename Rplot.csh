#!/bin/tcsh

set mymatrix=`ls *bootstrap_avg_mutinf_res_sum_0diag.txt | tail -n 1`
set myprefix = `echo $mymatrix | sed -e 's/bootstrap_avg_mutinf_res_sum_0diag.txt/mutinf/'`
cat > makeplot.R <<EOF
 source("http://bioconductor.org/biocLite.R")
 biocLite("marray")
 library(marray)
 S1M47_mutent <- read.table("${mymatrix}")
 pdf("${myprefix}_unclustered.pdf",width=48,height=48)
 heatmap(as.matrix(S1M47_mutent), col=maPalette(low="white",mid="blue",high="red"), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA, cexRow=1.05, cexCol=1.05, oldstyle=FALSE,symm = TRUE)
 dev.off()

pdf("${myprefix}_clustered.pdf",width=48,height=48)
heatmap(x=as.matrix((S1M47_mutent)),zlim=c(0,max(S1M47_mutent)*1.0), col = maPalette(low = "white", mid = "blue", high = "red", k =50),     add = FALSE, xaxs = "i", yaxs = "i",xaxt= "n", yaxt = "n", reorderfun = function(d,w) rev(reorder(d,w)),revC=TRUE, cexRow=0.8,  cexCol=0.8,      oldstyle = FALSE,symm=TRUE )
dev.off()

EOF


cat makeplot.R | R --no-save

#END



