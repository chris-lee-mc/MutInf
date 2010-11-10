#!/usr/bin/R --no-restore --no-save --no-readline
# Create a plot for a pair of dihedral angles

make_plot <- function(fn, contour_min, contour_max) {
	x <- as.matrix(read.table(fn))
	#print(paste("MATRIX: ", x))

	nx <- dim(x)[1]
	ny <- dim(x)[2]
	x2 <- x[2:nx,2:ny]
	
	pdf(paste(fn,'.pdf', sep=""), width=5, height=5, version = "1.4",pointsize = 8)

	filled.contour(x2, col=rainbow(20), nlevels=11,
					plot.axes = { axis(1, at=seq(0,1, length.out=nx-1), labels=x[1,2:ny]) 
        					      axis(2, at=seq(0,1, length.out=ny-1), labels=x[2:nx,1]) })

        dev.off()
}


args <-	commandArgs()
print(paste("ARGS: ", args))

fn1 <- args[5]

make_plot(fn1, 0, .1)
