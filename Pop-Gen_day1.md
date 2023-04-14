# Population genetics day 1

## Assembly statistics 

[assembly-stats](https://github.com/sanger-pathogens/assembly-stats)

Installation
```sh
conda install -c bioconda assembly-stats
```

Getting assembly statistics
```sh
assembly-stats Ppr.eifel.hap0.chr9.fasta
```
```sh
stats for Ppr.eifel.hap0.chr9.fasta
sum = 18765566, n = 1, ave = 18765566.00, largest = 18765566
N50 = 18765566, n = 1
N60 = 18765566, n = 1
N70 = 18765566, n = 1
N80 = 18765566, n = 1
N90 = 18765566, n = 1
N100 = 18765566, n = 1
N_count = 1530
Gaps = 6
``` 

## Data cleaning

[FastQC](https://github.com/s-andrews/FastQC)

Installation 
```sh
conda create -n fastqc_env 
conda activate fastqc_env
conda install -c bioconda fastqc

# and when you are done with it
conda deactivate
```

Running quality control
```sh
mkdir fastqc_raw
fastqc -o fastqc_raw reads1.fastq.gz reads2.fastq.gz
```

## Mapping

[BWA](https://github.com/lh3/bwa)
[SAMtools](https://github.com/samtools/samtools)

Installation
```sh
conda create -n bwa_env
conda activate bwa_env
conda install -c bioconda bwa
conda install -c bioconda samtools

# when you are done with it
conda deactivate
```

Mapping the reads to the reference assembly
```sh
bwa index Ppr.eifel.hap0.chr9.fasta
bwa mem -t 2 -o mapping.sam Ppr.eifel.hap0.chr9.fasta reads1.fastq.gz reads2.fastq.gz
```

Converting the SAM file (human readable) to a BAM file (binary, much lighter)
```sh
samtools view -@ 1 -b -o mapping.sam mapping.bam
```

Sorting the mapping fille
```sh
samtools sort -@ 1 -o mapping.sorted.bam mapping.bam
```

Getting mapping statistics
```sh
samtools flagstat -@ 1 mapping.bam 
```

## Generate VCF

[GATK](https://gatk.broadinstitute.org/hc/en-us)

```sh
gatk CreateSequenceDictionary -R Ppr.eifel.hap0.chr9.fasta -O Ppr.dict
```
  
    module avail
    module add
    module list
  
### call gvcf
    $gatk HaplotypeCaller -R $ref --emit-ref-confidence GVCF -I $bam2 -O $gvcf2
### merge them 
    $gatk CombineGVCFs -R $ref -V $gvcf1 -V $gvcf2 -O merged_gvcf/$mergedgvc
### detect SNPs
    $gatk GenotypeGVCFs -R $ref -V $mergedgvcf -O $vcf
### compress
    bgzip -f $vcf
    tabix -p vcf ${vcf}.gz    
### Select only SNPs from VCF 
    $gatk SelectVariants -select-type SNP -V ${vcf}.gz -O ${vcf}.snp.gz   
### filter SNPs by parameters
    $gatk VariantFiltration -V ${vcf}.snp.gz --filter-expression "QD <2.0 || MQ <40.0 || FS >60.0 || SOR >3.0 || ReadPosRankSum < -8.0" --filter-name "PASS" -O ${vcf}.snp.f.gz
  
    module purge
