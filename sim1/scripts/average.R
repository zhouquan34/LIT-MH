args <- commandArgs(TRUE)
group <- args[1]
id <- args[2]
pre0 <- paste("sim", group, "/out", id, "/s", sep="")

convert <- function(sel, iter, val){
	n = length(iter)
	now = 1
	ns = length(sel)
	v = numeric(ns)
	for (k in 1:n){
		while (now <= ns && sel[now] <= iter[k]){
			v[now] = val[k]
			now = now + 1
		}
	}
	return(v)
}


for (ss in c('r', 'i1', 'i2', 'sq')){
	pick = 0:2000
	if (ss == 'r'){
		pick = seq(0, 100000, by = 10)
	}
	np = length(pick)
	s0 = numeric(np) ## 10:beta.mse 
	s1 = numeric(np) ## 11:beta.mse.rb
	for (j in 1:100){
		pre = paste(pre0, j, "_", ss, sep="")
		path <- paste(pre, ".path.txt", sep="")

		d <- read.table(path, header=FALSE, skip=1)
		s0 = s0 + convert(pick, d[,1], d[,10])	
		s1 = s1 + convert(pick, d[,1], d[,11])	
	}
	ave = cbind(s0, s1)
	ave = ave / 100
	ave = cbind(pick, ave)
	file = paste("sim", group, "/", group, id, ss, ".bmse", sep="")
	write.table(ave, file=file, quote=F, row.names=F, col.names=F, sep="\t")
}


