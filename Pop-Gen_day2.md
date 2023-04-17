# Population genetics day 2

## VCF and SNP matrix

VCF contains more information than we usually need. For most downstream studies, the important data is what genotype each individual has on each variable site. This is usually presented as a matrix.

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	T501	T502	T505	T506	T507
chr1	10405	.	T	A	119.50	.	AC=2;AF=0.250;AN=8;DP=56;ExcessHet=0.3218;FS=0.000;MLEAC=1;MLEAF=0.125;MQ=34.26;QD=29.11;SOR=1.179	GT:AD:DP:GQ:PGT:PID:PL:PS	0/0:21,0:21:57:.:.:0,57,855	0/0:17,0:17:48:.:.:0,48,720	0/0:14,0:14:42:.:.:0,42,436	./.:1,0:1:.:.:.:0,0,0	1|1:0,3:3:9:1|1:10405_T_A:135,9,0:10405
chr1	10430	.	A	T	330.61	.	AC=3;AF=0.375;AN=8;BaseQRankSum=-1.147e+00;DP=90;ExcessHet=1.0474;FS=8.080;MLEAC=3;MLEAF=0.375;MQ=45.77;MQRankSum=-3.970e+00;QD=11.02;ReadPosRankSum=0.291;SOR=1.137	GT:AD:DP:GQ:PGT:PID:PL:PS	0/0:23,0:23:60:.:.:0,60,900	0/0:24,0:24:23:.:.:0,23,746	0|1:17,7:25:99:0|1:10403_A_T:121,0,598:10403	./.:10,0:10:.:.:.:0,0,0	1|1:0,6:6:18:1|1:10405_T_A:226,18,0:10405
chr1	10491	.	G	A	473.27	.	AC=1;AF=0.100;AN=10;BaseQRankSum=-5.780e-01;DP=161;ExcessHet=3.0103;FS=1.412;MLEAC=1;MLEAF=0.100;MQ=49.89;MQRankSum=-1.590e-01;QD=15.78;ReadPosRankSum=0.656;SOR=1.082	GT:AD:DP:GQ:PL	0/0:35,0:35:99:0,99,1112	0/0:46,0:46:99:0,105,1397	0/1:15,15:31:99:483,0,504	0/0:27,0:27:78:0,78,1170	0/0:17,0:17:48:0,48,720

```

## Some basic but useful R commands

```
mat.1 = matrix(0,20,10) #Creating matrix (aka 2-dimensional array)
mat.1[3,4] = 10 #Member of matrix

mat.1 = matrix(rgeom(200,0.1),20,10)
apply(mat.1,1,sum) #Doing math operations on matrix, row by row (1) or column by column (2)

x = 101:200
y = sample(x,10,replace=TRUE) #Random sampling from a vector, with or without replacement

x = runif(20,0,1) #Uniformly distributed ("Random") betweeen 0 and 1
(x > 0.5) #Boolean vector, can be used to subsetting data
mat.1[(x>0.5),] #Only the rows that correspond to x entries larger than 0.5

plot (1:20, x) #Default plot with scattered dots
plot (1:20, x, type="l", xlim=c(0,40), ylim=c(-1,2), xlab = "X axis", ylab = "Y axis", main = "Test Plot") #Line plot with custom "field of view" and labels
lines(1:20, sort(runif(20,0,1)), col=2) #Add line to plot

x = 8
paste("Text_",x,"_text.txt",sep="") #Concatenate character strings, useful for file names

```

## R reading of SNP matrix

```
pop.sites = as.matrix(read.table("pop_file.txt",header=FALSE))
pop.pos = read.table("pop_file_pos.txt",header=FALSE)[,1]

```

## Calculating Basic Population Statistics

Here are the R commands used to calculate the relevant statistics. These will be added after class.

Theta_pi

Theta_w

Allele Frequency Spectrum

Fis (Hardy-Weinburg)


