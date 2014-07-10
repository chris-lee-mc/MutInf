#!/bin/bash

cat > doRplot_triplet_hist.R <<EOF
blah <- read.table("mutinfs_for_histograms.txt")

pdf("triplet_vs_12_plot.pdf")
plot(blah\$V1,blah\$V6,xlab="mutual information between dihedral 1 and 2", ylab="triplet mutual info")
dev.off()

pdf("ind_vs_12_plot.pdf")
plot(blah\$V1,blah\$V5,xlab="mutual information between dihedral 1 and 2", ylab="Triplet mutual info")
dev.off()

pdf("uncorrected_vs_12_plot.pdf")
dev.off()
pdf("uncorrected_vs_12_plot.pdf")
plot(blah\$V1,blah\$V4,xlab="mutual information between dihedral 1 and 2", ylab="Triplet mutual info")
dev.off()

pdf("mutinf_1_2_hist.pdf")
hist(blah\$V1, xlab="mutual info between dihedral 1 and 2")
dev.off()


pdf("mutinf_2_3_hist.pdf")
hist(blah\$V2, xlab="mutual info bewteen dihedral 2 and 3")
dev.off()

pdf("mutinf_1_3_hist.pdf")
hist(blah\$V3, xlab="mutual info between dihedral 1 and 3")
dev.off()


pdf("uncorrected_triplet_mutinf_hist.pdf")
hist(blah\$V4, xlab="Triplet mutual information", main="Uncorrected Triplet Mutual Information Values")
dev.off()

pdf("independent_triplet_mutinf_hist.pdf")
hist(blah\$V5, xlab="Triplet mutual information", main="Independent Triplet Mutual Information Values")
dev.off()


pdf("corrected_triplet_mutinf_hist.pdf")
hist(blah\$V6, xlab="Triplet mutual information", main="Corrected Triplet Mutual Information Values")
dev.off()


EOF

cat  doRplot_triplet_hist.R | R --no-save

#EOF

