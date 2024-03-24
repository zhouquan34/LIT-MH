args <- commandArgs(TRUE)

A0 = 1e5
A1 = 2e3
qs = c(0.025, 0.975)
null.a = -2.5e-7
null.b = -4.00002e-9

qstring <- function(x){
	m = round(median(x))
	if (m == A0){m = "$10^5+$"}
	qx = round(quantile(x, 0.95))
	if (qx == A0){qx = "$10^5+$"}
	if (qx == A1){qx = "2000+"}
	ss = paste(m, " (", qx, ")", sep="")
	return(ss)
}

for (sim in c('A', 'B', 'E', 'F')){
	null = 0
	if (sim == 'A' || sim == 'B'){
		null = null.a
	}else{
		null = null.b
	}

for (j in 1:3){
	file = paste('sim', sim, '/', 'out', j, '.txt', sep="")
	d <- read.table(file)
	true <- d[,2]	
	pick = c(3,5,7, 9)
	ns = length(pick)
	d0 = d[,pick]
	high = as.numeric(apply(d0, 1, max))
	stat = numeric(ns)
	for (s in 1:ns){
		fail = (d0[,s] < high)
		iter = d[,pick[s] + 1]
		if (s == 1){
			iter[fail] = A0
		}else{
			iter[fail] = A1
		}
		stat[s] = qstring(iter)
	}
	cat(sim, j, qstring(true - null), qstring(high - null), stat[1], stat[2], stat[3], stat[4], sep="\t")
	cat("\n")
}

}


