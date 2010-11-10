#!/usr/bin/R --no-restore --no-save --no-readline
# Create 3 plots for a pair of dihedral angles: Pij, PiPj, and Pij/PiPj

make_plot <- function(fn, contour_min, contour_max) {
	x <- as.matrix(read.table(fn))
	#print(paste("MATRIX: ", x))

	nx <- dim(x)[1]
	ny <- dim(x)[2]
	x2 <- x[2:nx,2:ny]
	
	pdf(paste(fn,'.pdf', sep=""), width=5, height=5, version = "1.4",pointsize = 8)

	filled.contour(x2, zlim=c(0,.2), col=rainbow(20), levels=seq(contour_min, contour_max, length=11),
					plot.axes = { axis(1, at=seq(0,1, length.out=nx-1), labels=x[1,2:ny]) 
        					      axis(2, at=seq(0,1, length.out=ny-1), labels=x[2:nx,1]) })
#	axis(1, at=seq(0,1, length.out=nx-1), labels=x[1,2:ny])
#	axis(2, at=seq(0,1, length.out=ny-1), labels=x[2:nx,1])

    dev.off()
}


args <-	commandArgs()
print(paste("ARGS: ", args))

fn1 <- args[5]
#fn2 <- args[6]
#fn3 <- args[7]

make_plot(fn1, 0, .1)
#make_plot(fn2, 0, .1)
#make_plot(fn3, 0, 2)

