args <- commandArgs(TRUE)
seed <- as.numeric(args[1])
sigma <- as.numeric(args[2])
n <- as.numeric(args[3])
p <- as.numeric(args[4])
s <- as.numeric(args[5])

set.seed(seed)

d = 20
C = diag(d)
for (i in 1:d){
	for (j in 1:d){
		C[i,j] = exp(-abs(i - j)/3)
	}
}
R = chol(C)

X = matrix(rnorm(n*p), nrow=n)  
for (i in 1:(p/d)){
	start = (i-1)*d + 1
	pick = start:(start + d - 1)
	X[,pick] = X[,pick] %*% R 
}
X = scale(X, center=TRUE)
X = round(X, 6)

sel = sample(1:p, s, replace=FALSE)
beta = rnorm(s, mean=0, sd=sigma)
y0 = X[,sel] %*% beta 
y = y0 + rnorm(n)
y = y - mean(y)
cat("pve =", var(y0)/var(y), "\n") 

yy = t(y) %*% y
kappa = 1
g = p - 1

marginal.ll.gprior <- function(s){
	if (length(s) == 0){
		ll = -0.5 * n * log( 1 + 1 / g ) 
		return(list("beta"=c(), "ll"=ll))
	}
	z = as.matrix(X[,s])
	zy = t(z) %*% y
	beta = solve(t(z) %*% z, zy)
	k = length(s)
	ll = -0.5 * n * log( (yy - t(zy) %*% beta + yy / g)/yy )  - k * kappa * log(p) - 0.5 * k * log(g + 1)
	return(list("beta"=beta, "ll"=ll))
}

model.ll <- function(s){
	s = s[order(s)]
	name = paste(s, collapse=",")
	res = marginal.ll.gprior(s)
	ll = as.numeric(res$ll)
	return(ll)
}

pre = paste("dat/cor", seed, "_", as.character(10*sigma), sep="")
true = cbind(sel, beta)

cat("True", model.ll(sel), "\n")
mbf = matrix(0, nrow=p, ncol=3)
for (k in 1:p){
	ifs = 0
	if (k %in% sel){ifs = 1}
	mbf[k, ] = c(k, ifs, model.ll(k) )
#	if (k<=10){cat(k, model.ll(k), "\n") }
}

write.table(mbf, file=paste(pre, "bf", sep="."), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(X, file=paste(pre, "mat", sep="."), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=" ")
write.table(y, file=paste(pre, "ph", sep="."), quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(true, file=paste(pre, "beta", sep="."), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

