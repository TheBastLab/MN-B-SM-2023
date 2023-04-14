
arg = commandArgs(TRUE)
fname = arg[1]
outname = arg[2]
ssize = as.numeric(arg[3])

ne = 1000 #Population size (diploid)
seq.leng = 1e6 #Sequence size 1MB

#Read populations
pop.sites = as.matrix(read.table(paste(fname,".txt",sep=""),header=FALSE))
pop.pos = read.table(paste(fname,"_pos.txt",sep=""),header=FALSE)[,1]

samp = sample(ne,ssize)
samp.sites = pop.sites[c(samp,ne+samp),]

afq = apply(samp.sites,2,mean)

samp.sites = samp.sites[,(afq*(afq-1) != 0)]
samp.pos = pop.pos[(afq*(afq-1) != 0)]

write(samp.pos,paste(outname,"_pos.txt",sep=""),ncolumns=1)
write.table(samp.sites,paste(outname,".txt",sep=""),row.names=FALSE,col.names=FALSE)
