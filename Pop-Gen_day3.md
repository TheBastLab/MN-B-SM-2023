# Add module

    module avail
    module add
    module list

# Data download 
    https://uni-koeln.sciebo.de/s/IDcUK28aB5B8Dwx
    https://uni-koeln.sciebo.de/s/ZqheSDmcXZwAj7e
    passwd atcg

# PCA
## What's PCA?
  Principal component analysis, or PCA, is a statistical procedure that allows you to summarize the information content in large data tables by means of a smaller set of “summary indices” that can be more easily visualized and analyzed.

## commands
    les $.vcf |awk '/^#/||$7=="PASS"||$7=="MQRankSum-12.5;ReadPosRankSum-8"||$7=="MQRankSum-12.5"||$7=="ReadPosRankSum-8"' > $.hardfiltered.vcf
    awk '$5 !~ /([[:alpha:]]|*)+,([[:alpha:]]|*)/{print}' $hardfiltered.vcf > $hardfiltered.bi.vcf
    vcftools --gzvcf 9samples_chr9.vcf.f.f.bi.gz --plink --out 9samples_chr9.vcf.f.f.bi.gz.plink
    plink --noweb --file $hardfiltered.bi.vcf --make-bed --out $hardfiltered.bi.vcf_bfile 
    plink --threads 16 --bfile $hardfiltered.bi.vcf_bfile --pca 3 --out $hardfiltered.bi.pca3.vcf_bfile
    
## plot
    R
    library(data.table)
    install.packages("tidyverse")
    library(tidyverse)
    install.packages("ggrepel")
    library(ggrepel) 
    re1a = fread("all_samples_mapped2Ger.hardfiltered.bi.plink.vcf_bfile_pca20.eigenval")
    re1b = fread("all_samples_mapped2Ger.hardfiltered.bi.plink.vcf_bfile_pca20.eigenvec")
    re1a$por = re1a$V1/sum(re1a$V1)*100
    pdf('pca12.pdf') 
    ggplot(re1b,aes(x = V3,y = V4,color=V1)) + geom_point(size=1.5) + xlab(paste0("PC1 (",round(re1a$por[1],2),"%)")) +ylab(paste0
    ("PC2 (",round(re1a$por[2],2),"%)"))+scale_color_manual(values = c("#953CCC", "#74BBF1", "#FFD500", "#A8F1AD"))+scale_x_contin
    uous(limits=c(-0.4, 0.4))+geom_text_repel(aes(label=V2))+labs(color="Country")+theme_bw()
    dev.off()
    

# Fst
## What's Fst?
  FST is the proportion of the total genetic variance contained in a subpopulation (the S subscript) relative to the total genetic variance (the T subscript). Values can range from 0 to 1. High FST implies a considerable degree of differentiation among populations.
  
## commands
    vcftools --vcf all_samples_mapped2Ger.hardfiltered.bi.vcf --weir-fst-pop de.txt --weir-fst-pop ru.txt  --out de-ru --fst-window-size 500000 --fst-window-step 50000
    R
    a<-read.table('de-ru.windowed.weir.fst',header=T)
    head(a)
      CHROM BIN_START BIN_END N_VARIANTS WEIGHTED_FST MEAN_FST
    1  chr9         1  500000      13973     0.241968 0.174812
    2  chr9     50001  550000      14418     0.222682 0.163322
    3  chr9    100001  600000      16315     0.209036 0.147515
    4  chr9    150001  650000      16859     0.204524 0.142917
    5  chr9    200001  700000      15382     0.204730 0.143377
    6  chr9    250001  750000      15061     0.202856 0.140516
    pdf('plt.pdf')
    plot((a$BIN_START-1)/50000+1,a$MEAN_FST,type='p')
    dev.off()

# Pi
## What's Pi?
Nucleotide diversity (often referred to using the symbol π) is the average pairwise difference between all possible pairs of individuals in your sample. It is a very intuitive and simple measure of genetic diversity, and is accurately estimated even with very few samples.
## commands
    vcftools --vcf ../all_samples_mapped2Ger.hardfiltered.bi.vcf --keep $i.txt --recode --recode-INFO-all --out pop_${i}
    vcftools --vcf pop_${i}.recode.vcf --out ${i}_pi_500kb --window-pi 500000 --window-pi-step 50000
    
# Tajima's D
## What's Tajima's D?
  Tajima's D is a population genetic test statistic created by and named after the Japanese researcher Fumio Tajima.[1] Tajima's D is computed as the difference between two measures of genetic diversity: the mean number of pairwise differences and the number of segregating sites, each scaled so that they are expected to be the same in a neutrally evolving population of constant size.
## commands
    vcftools --vcf ../all_samples_mapped2Ger.hardfiltered.bi.vcf --keep $i.txt --recode --recode-INFO-all --out pop_${i}
    vcftools --vcf pop_${i}.recode.vcf --out ${i}_TajimaD_500kb --TajimaD 500000
    
    
