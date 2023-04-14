
arg = commandArgs(TRUE)
fname = arg[1]


ne = 1000 #Population size (diploid)
seq.leng = 1e6 #Sequence size 1MB
mu = 0.05 #Gives theta=200 for the entire sequence, 2e-4 per nucleotide
rec = 0.01 #Recombination rate 0.01 across the entire sequence

##### Part 1: produce sexually reproducing populations

pop.sites = matrix(rbinom(20*ne,1,0.5),2*ne,10)
pop.pos = c(1e5,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,1e6)
#Initialize the population with artificial SNPs. During the burn-in, these fake SNPs should all be fixed/lost.

for (gen in 1:16000) { #Burn-in process. Each run of the loop is a generation.
 #In my usual simulations, all these are done with a function that is invoked once per generation. This would make editing the script easier. But here we put everything sequentially for ease of understanding.
 
 par1 = sample(1:(2*ne),2*ne,replace=TRUE,prob=rep(1,2*ne)) #Parentage sampling: from haplotype to gamete. Because of neutrality, all haplotypes have same probability to be sampled.
 loci = length(pop.pos)
 #Recombinations
 rec.info = rbinom(2*ne,1,rec)*runif(2*ne,0,1e6) #0 if no recombination; otherwise position of recombination
 ptnr = c(ne+(1:ne),1:ne)
 rec.sites = (t(matrix(pop.pos,loci,2*ne))>=matrix(rec.info,2*ne,loci)) #Sites on the right side of the recombination point. For non-recombined gametes, all sites taken.
 newsites = pop.sites[par1,] * rec.sites + pop.sites[ptnr[par1],] * (1-rec.sites)
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

 if (gen %% 1000 == 0) {print(gen)} #Print progress per kgen
}

write(pop.pos[order(pop.pos)],paste(fname,"_pos.txt",sep=""),ncolumns=1)
write.table(pop.sites[,order(pop.pos)],paste(fname,".txt",sep=""),row.names=FALSE,col.names=FALSE)

