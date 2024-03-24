library(LaplacesDemon, lib.loc="~/local/R/")

ess <- function(stat, life){		
	n = length(life)
	x = numeric(sum(life))
	now = 1
	for (i in 1:n){
		if (life[i] > 0){
			x[now:(now + life[i]-1)] = rep(stat[i], life[i])
			now = now + life[i]
		}
	}
	ess = ESS(x)
	if (ess < 1e-10){ess = 0}
	return(ess)
}

pre0 <- "dat/out/s"

f1 = "ess.log"
f2 = "ess.table"
for (sigma in 1:5){
for (ss in c("r", "i1", "i2", "sq")){
	cat("Parsing", sigma, ss, "\n")
	sum = 0
	for (j in 1:20){
		pre = paste(pre0, j, "_", sigma, "_", ss, sep="")
		path <- paste(pre, ".path.txt", sep="")
		d <- read.table(path, header=FALSE, skip=1)
		life <- d[,2]
		pve <- d[,6]
		beta.se <- d[,10]
		stat = c(ess(beta.se, life), ess(pve, life))
		cat(sigma, ss, j, stat, "\n", file=f1, append=TRUE)
		sum = stat + sum
	}
	cat(sigma, ss, sum/20, "\n")
	cat(ss, sigma, sum/20, "\n", file=f2, append=TRUE)
}
}

