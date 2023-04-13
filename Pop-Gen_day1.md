# Population genetics day 1

## Assembly statistics 

[assembly-stats](https://github.com/sanger-pathogens/assembly-stats)
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

```sh
mkdir fastqc_raw
fastqc -o fastqc_raw reads1.fastq.gz reads2.fastq.gz
```

## Mapping
bowtie2, bwa
samtools view, samtools sort

## Generate VCF
GATK
