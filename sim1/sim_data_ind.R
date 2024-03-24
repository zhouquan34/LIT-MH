args <- commandArgs(TRUE)
seed <- as.numeric(args[1])
snr <- as.numeric(args[2])
set.seed(seed)
n <- as.numeric(args[3])
p <- as.numeric(args[4])
dir <- args[5]
m <- as.numeric(args[6])
s = 10
b0 = c(2,-3,2,2,-3,3,-2,3,-2,3)
beta = c(b0, rep(0, p - s)) * sqrt( log(p)/n ) * snr

X = matrix(rnorm(n*p), nrow=n)
X = scale(X, center=TRUE)
y = X %*% beta + rnorm(n)
X = round(X, 6)

pre = paste(dir, "/ind", args[1], "_", args[6], sep="")

true = cbind(1:10, beta[1:10])


write.table(X, file=paste(pre, ".mat", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=" ")
write.table(y, file=paste(pre, ".ph", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(true, file=paste(pre, ".beta", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

