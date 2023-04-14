
arg = commandArgs(TRUE)
fname = arg[1]
driftline.file = arg[2]

ne = 1000 #Population size (diploid)
seq.leng = 1e6 #Sequence size 1MB
mu = 0.05 #Gives theta=200 for the entire sequence, 2e-4 per nucleotide
rec.rate = 0.01 #Recombination rate 0.01 across the entire sequence
gc.rate=0.002

#Read populations
pop.sites = as.matrix(read.table(paste(fname,".txt",sep=""),header=FALSE))
pop.pos = read.table(paste(fname,"_pos.txt",sep=""),header=FALSE)[,1]

#Read driftline
driftline = read.table(driftline.file,header=FALSE)[,1]
driftline = c(driftline,rep(1000,10000))

pop.cloner = rep(0,1000)
pop.cloner[1] = 1 #Initialize the cloner mark

#Evolve population 

for (gen in 1:10000) {

 par1 = sample(1:(2*ne),2*ne,replace=TRUE,prob=rep(1,2*ne)) #Parentage sampling: from haplotype to gamete. Because of neutrality, all haplotypes have same probability to be sampled.
 loci = length(pop.pos)
 
 sex.num = ne-driftline[gen]
 
 if (sex.num > 0) {
  temp = which(pop.cloner==0)
  par1 = sample(c(temp,temp+ne),2*sex.num,replace=TRUE) #Parentage sampling: from haplotype to gamete.
 }
 if (sex.num < ne) {
  temp = which(pop.cloner==1)
  if (length(temp)==1) {par2 = rep(temp,(ne-sex.num))}
  else{par2 = sample(temp,(ne-sex.num),replace=TRUE)} #Clonal sampling: from individual to individual.
 }
 
 
 newsites = matrix(0,2*ne,loci) #Empty matrix to put the sexual offspring first and clones second.
 #Recombinations
 ptnr = c(ne+(1:ne),1:ne)
 if (sex.num > 0) {
  rec.info = rbinom(2*sex.num,1,rec.rate)*runif(2*sex.num,0,1e6) #0 if no recombination; otherwise position of recombination
  rec.sites = (t(matrix(pop.pos,loci,2*sex.num))>=matrix(rec.info,2*sex.num,loci))
  newsites[1:sex.num,] = pop.sites[par1[1:sex.num],] * rec.sites[1:sex.num,] + pop.sites[ptnr[par1[1:sex.num]],] * (1-rec.sites[1:sex.num,]) #Here the new gametes have to be separated into two "discontinuous" blocks in the offspring matrix
  newsites[(ne+1):(ne+sex.num),] = pop.sites[par1[(sex.num+1):(2*sex.num)],] * rec.sites[(sex.num+1):(2*sex.num),] + pop.sites[ptnr[par1[(sex.num+1):(2*sex.num)]],] * (1-rec.sites[(sex.num+1):(2*sex.num),])
 }
 #Clones
 if (sex.num < ne) {
  newsites[(sex.num+1):ne,] = pop.sites[par2,]
  newsites[(ne+sex.num+1):(2*ne),] = pop.sites[(par2+ne),]
 } 
 
 if(loci == 1) {newsites = as.matrix(newsites)} #Make sure the sites is a matrix even if only 1 site
 
 #Mutations
 mut.inds = (runif(2*ne,0,1) < mu) #Which individuals have mutation? The total number is binomial, close to Poisson
 mut.ct = 0 #Counting number of individuals already mutated, to know which locus to mutate.
 if (sum(mut.inds) > 0) {
  newsites = cbind(newsites,matrix(0,2*ne,sum(mut.inds)))
  for (i in 1:(2*ne)) {
   if (mut.inds[i]) {
    mut.ct = mut.ct + 1
	newsites[i,(loci+mut.ct)] = 1
   }
  }
 }
 pop.pos = c(pop.pos,ceiling(runif(mut.ct,0,1e6)))
 
 #Clean up fixed/lost sites
 afq = apply(newsites,2,mean)
 if (sum(afq) != 0 & sum(afq>0) != 1) {
  newsites = newsites[,(afq*(afq-1) != 0)]
  pop.pos = pop.pos[(afq*(afq-1) != 0)]
 }
 
 pop.sites = newsites
 pop.cloner = c(rep(0,sex.num),rep(1,ne-sex.num))
 

 #Gene conversion
 loci = length(pop.pos)
 
 gc.num = rpois(1,gc.rate*ne)
 
 if (gc.num > 0) {
  gci.ind = sample(1:ne,gc.num,replace=TRUE)
  gci.leng = rgeom(gc.num,1/1000)
  gci.dir = sample(1:2,gc.num,replace=TRUE)
  for (i in 1:gc.num) {
   gci.loc = runif(1,0,(1e6-gci.leng[i]))
   temp = (pop.pos>gci.loc & pop.pos<=(gci.loc+gci.leng[i]))
   if (gci.dir[i] == 1) {pop.sites[gci.ind[i],temp] = pop.sites[(gci.ind[i]+ne),temp]}
   else {pop.sites[(gci.ind[i]+ne),temp] = pop.sites[gci.ind[i],temp]}
  }
 }


 if (gen %% 1000 == 0) {print(gen)} #Print progress per kgen


 if (gen %% 2000 == 0) {
  write(pop.pos[order(pop.pos)],paste(fname,"clonegen_",gen,"_pos.txt",sep=""),ncolumns=1)
  write.table(pop.sites[,order(pop.pos)],paste(fname,"clonegen_",gen,".txt",sep=""),row.names=FALSE,col.names=FALSE)
 }




}
