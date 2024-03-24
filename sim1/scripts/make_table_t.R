A0 = 1e5
A1 = 2e3
qs = c(0.025, 0.975)
null.a = -2.5e-7
null.b = -4.00002e-9

time.all <- read.table('time.txt')
g <- c('r', 'i1', 'i2', 'sq')
ns = 4

for (sim in c('A', 'B', 'E', 'F')){
	null = 0
	if (sim == 'A' || sim == 'B'){
		null = null.a
	}else{
		null = null.b
	}

	for (j in 1:3){
		file = paste('sim', sim, '/', 'out', j, '.txt', sep="")
		t = time.all[which(time.all[,1] == sim & time.all[,2] == j), 3:5]	
		sec = matrix(0, nrow=100, ncol=4)
		for (s in 1:ns){
			t2 = t[which(t[,2] == g[s]), c(1,3)]
			sec[,s] = t2[order(t2[,1]),2]
		}

		d <- read.table(file)
		true <- d[,2]	
		pick = c(3,5,7, 9)
		d0 = d[,pick]
		high = as.numeric(apply(d0, 1, max))
		stat = numeric(ns)
		stat.med = numeric(ns)
		for (s in 1:ns){
			fail = (d0[,s] < high)
			iter = d[,pick[s] + 1]
			total = 0
			if (s == 1){
				iter[fail] = A0
				total = A0
			}else{
				iter[fail] = A1
				total = A1
			}
			stat[s] = mean(iter * sec[,s]/ total)
			stat.med[s] = median(iter * sec[,s]/ total)
		}
		cat(sim, j, signif(stat.med,4), "\n")
	}
}


