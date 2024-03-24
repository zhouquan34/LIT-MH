args <- commandArgs(TRUE)
pre <- args[1]
dat.file = paste(pre, "mat", sep=".")
ph.file = paste(pre, "ph", sep=".")

X <- read.table(dat.file, header=TRUE)
y <- scan(ph.file)
y <- y - mean(y)
n = nrow(X)
p = ncol(X)
cat("n =", n, "; p =", p, "\n")
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
	if (length(s) > 1){
		s = s[order(s)]
	}
	name = paste(s, collapse=",")
	res = marginal.ll.gprior(s)
	ll = as.numeric(res$ll)
	return(ll)
}

m <- c()
ll = model.ll(m)

while (1){
	ll2 = numeric(p)
	cat("Model size = ", length(m), "; like = ", model.ll(m),  "\n")
	for (k in 1:p){
		if (k %in% m){
			ll2[k] = -1e10
		}else{
			m2 = c(m, k)
			ll2[k] = model.ll(m2)
		}
	}
	if (max(ll2) > ll){
		m = c(m, which.max(ll2))
		ll = max(ll2)
	}else{
		break 
	}
}


while (length(m) > 0){
	cat("Model size = ", length(m), "; like = ", model.ll(m),  "\n")
	size = length(m)
	ll2 = numeric(size)
	for (k in 1:size){
		m2 = m[-k]
		ll2[k] = model.ll(m2)
	}
	if (max(ll2) > ll){
		m = m[-which.max(ll2)]
		ll = max(ll2)
	}else{
		break 
	}			
}

write.table(as.matrix(m), file=paste(pre, ".im", sep=""), quote=F, row.names=F, col.names=F)

