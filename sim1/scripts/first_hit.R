args <- commandArgs(TRUE)
pre0 <- args[1]

for (j in 1:100){
	stat = c()
	true <- 0
	for (ss in c('r', 'i1', 'i2',  'sq')){
		pre = paste(pre0, j, "_", ss, sep="")
		log <- paste(pre, ".log.txt", sep="")
		path <- paste(pre, ".path.txt", sep="")

		a <- readLines(log)
		for (i in 1:length(a)){
			line = a[[i]]	
			ch = unlist(strsplit(line, "\\s+"))
			if (length(ch) < 2){next}
			if (ch[1] == "True" && ch[2] == "ll"){
				true = as.numeric(ch[4])
				break
			}
		}

		d <- read.table(path, header=FALSE, skip=1)
		m = max(d[,5])
		diff = m - true 
		k = min(which(d[,5] >= m - 1e-6)) - 1
		stat = c(stat, m, d[k, 1])
	}
	stat = c(j, true, stat)
	cat(stat, "\n")
}


