#!/bin/tcsh

set mymatrix=${1}
set myprefix = `echo $mymatrix | sed -e 's/reslist.*bootstrap/mutinf/'`
cat > makeplot.R <<EOF
 library(marray)
 S1M47_mutent <- read.table("${mymatrix}")
 pdf("${myprefix}_hist.pdf",width=48,height=48)
 heatmap(as.matrix(S1M47_mutent), col=maPalette(low="white",mid="blue",high="red"), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA, cexRow=1.05, cexCol=1.05, oldstyle=FALSE,symm = TRUE)
 dev.off()


EOF


cat makeplot.R | R --no-save

#END



