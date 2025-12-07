# Final project - code - supplementary materials
# Rachel Lekanoff, F694, December 2025

## Read Counts of Leukoma staminea Transcriptomes Generated via de novo and Reference Genome Assembly Approaches
# annotated code and scripts used to analyze L. staminea RNAseq data
# UAF's RCS Chinook HPC cluster was used to run all analyses
# envirnoments are generally set up with one program per environment, but could be set up in one environment


Download the files from NCBI with known SRR number, using prefetch in sratoolkit.
sratoolkit must be configured, and then can be used to download SRR files from NCBI. convert RNAseq data from .SRR to .fastq.gz.

```bash
prefetch SRR18446124 # download SRR file
sbatch srr_fastq_gz.sh      # script to convert
```

srr_fastq_gz.sh 
```bash
#!/bin/bash

#SBATCH --partition=t1small
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=12
#SBATCH --time=24:00:00

export PATH=$PATH:/home/rmlekanoff/sratoolkit.3.2.1-ubuntu64/bin  # location of downloaded sratoolkit on Chinook

fastq-dump --split-files --gzip --outdir /home/rmlekanoff/final_proj /home/rmlekanoff/final_proj/srr/SRR18446124

# fastq-dump takes the raw RNAseq data in SRR format and separates them to forward and reverse reads and outputs those paired read files into the specified directory.
```

run_fastqc.sh to run quality checks on the raw fastq.gz files
remember to activate environment: 
conda activate fastqc_env
sbatch run_fastqc.sh
```bash
#!/bin/bash
#SBATCH --partition=t1small
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=12
#SBATCH --time=08:00:00
#SBATCH --dependency=afterok:483901 ssr_fastq_gz.sh

INPUT_DIR="/home/rmlekanoff/final_proj"
OUTPUT_DIR="/home/rmlekanoff/final_proj/fastqc_reports"

# Perform FastQC on each specified file line by line
fastqc "${INPUT_DIR}/SRR18446124_1.fastq.gz" -o "$OUTPUT_DIR"
fastqc "${INPUT_DIR}/SRR18446124_2.fastq.gz" -o "$OUTPUT_DIR"

# download .html files to local machine to view in your browser
```

Trimming the reads with BBDuk
BBDuk installation steps (https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide)
test install: /home/rmlekanoff/bbmap/stats.sh in=bbmap/resources/phix174_ill.ref.fa.gz

```bash
sbatch bbduk_trim.sh
```

```bash
#!/bin/bash

#SBATCH --partition=t1small
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=12
#SBATCH --time=24:00:00

export PATH=$PATH:/home/rmlekanoff/bbmap # location of bbduk program in Chinook

bbduk.sh in1=/home/rmlekanoff/final_proj/SRR18446124_1.fastq.gz in2=/home/rmlekanoff/final_proj/SRR18446124_2.fastq.gz out1=clean_SRR18446124_1.fastq.gz out2=clean_SRR18446124_2.fastq.gz ktrim=r k=23 mink=11 hdist=1 ref=/home/rmlekanoff/bbmap/resources/adapters.fa tbo tpe

# default common adapters are included in bbduk as adapters.fa, which is used here

```

after trimming, use fastqc on the clean SRR18446124 .fastq.gz files
```bash
sbatch clean_fastq.sh
```

```bash
#!/bin/bash

#SBATCH --partition=t1small
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=12
#SBATCH --time=24:00:00

INPUT_DIR="/home/rmlekanoff/final_proj"
OUTPUT_DIR="/home/rmlekanoff/final_proj/clean_fastqc_reports"

# Perform FastQC on each specified file line by line
fastqc "${INPUT_DIR}/clean_SRR18446124_1.fastq.gz" -o "$OUTPUT_DIR"
fastqc "${INPUT_DIR}/clean_SRR18446124_2.fastq.gz" -o "$OUTPUT_DIR"

```


de novo assembly of transcriptome with Trinity

```bash
conda activate trinity
sbatch trinity_clean.sh
```
 
```bash
#!/bin/bash
#SBATCH --partition=t1standard
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=28
#SBATCH --ntasks=84
#SBATCH --mail-user=rmlekanoff@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

Trinity --seqType fq --max_memory 150G --left /home/rmlekanoff/final_proj/clean_SRR18446124_1.fastq.gz --right /home/rmlekanoff/final_proj/clean_SRR18446124_2.fastq.gz --output /center1/FISH694/rmlekanoff/trinity_out_clean_24

# memory usage for Trinity is recommended at 1 Gb Ram per 100,000 reads. This dataset has ~150,000,000 reads. Earlier run attempts did not take this into account, meaning the individual steps in Trinity completed, but the final assembly of the de novo transcriptome did not happen. In this scenario, the logs of the Trinity run read that all individual steps were completely successfully (even with a happy litle :-) included). But the final assembly is not compiled and there is no indication of this except for missing the final 'Trinity.fasta' file needed for downstream analyses. Fortunately Trinity will pick up where it left off on the last run, so the analysis can continue after properly optimizing the script!
```

the reference genome needs to be indexed. Downlad the whole genome from NCBI. Run hisat2-build to generate indexed genome files.

```bash
conda activate hisat2
sbatch reference.sh
```

```bash
#!/bin/bash

#SBATCH --partition=t1small
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=12

hisat2-build /center1/FISH694/rmlekanoff/man-clam/ncbi_dataset/data/GCF_026571515.1/GCF_026571515.1_ASM2657151v2_genomic.fna /home/rmlekanoff/final_proj/index/manila-clam

# index the reference genome
```

assemble the Leukoma RNAseq reads to the indexed Manila clam reference genome. This still takes place in the hisat2 environment.

```bash
sbatch assemble.sh
```

```bash
#!/bin/bash

#SBATCH --partition=t1small
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=24

hisat2 -x /home/rmlekanoff/final_proj/index/manila-clam -1 /home/rmlekanoff/final_proj/clean_SRR18446124_1.fastq.gz -2 /home/rmlekanoff/final_proj/clean_SRR18446124_2.fastq.gz -S /center1/FISH694/rmlekanoff/hisat2_out.sam
 
```

index the trinity made de novo assembly with hisat2-build.

```bash
conda activate trinity
sbatch trinity_present.sh
```

```bash
#!/bin/bash

#SBATCH --partition=t1standard
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=rmlekanoff@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

hisat2-build -f /center1/FISH694/rmlekanoff/trinity_out_clean_24/Trinity.fasta /center1/FISH694/rmlekanoff/trinity_index/leukoma_rnaseq

# -f argument indicates the reference input is a fasta file.
```

align the original RNAseq reads to the indexed trinity de novo transcriptome
map_trinity.sh
```bash
#!/bin/bash

#SBATCH --partition=t1small
#SBATCH --ntasks=24
#SBATCH --tasks-per-node=24

hisat2 -x /center1/FISH694/rmlekanoff/trinity_index/leukoma_rnaseq -1 /home/rmlekanoff/final_proj/clean_SRR18446124_1.fastq.gz -2 /home/rmlekanoff/final_proj/clean_SRR18446124_2.fastq.gz -S /center1/FISH694/rmlekanoff/trinity_hisat2_out.sam

```

the output of hisat2 is in .sam format, which is not human readable or ready to use in downstream analyses. Convert the .sam to .bam using samtools. 

```bash
conda activate samtools_env
sbatch bamsam.sh
```

```bash
#!/bin/bash

#SBATCH --partition=t1small
#SBATCH --ntasks=12
#SBATCH --tasks-per-node=12

samtools view -b -S /center1/FISH694/rmlekanoff/hisat2_out.sam > /home/rmlekanoff/final_proj/hisat2_out.bam

samtools view -b -S /center1/FISH694/rmlekanoff/trinity_hisat2_out.sam > /home/rmlekanoff/final_proj/trinity_hisat2_out.bam

```

print stats from flagstat analysis on the Leukoma reads mapped to both the Manila clam reference genome and the de novo assembly and save the outputs to a txt file 
conda activate samtool_env
samstats.sh
```bash
#!/bin/bash

#SBATCH --partition=t1small
#SBATCH --ntasks=2
#SBATCH --tasks-per-node=12

samtools flagstat hisat2_out.bam > hisat2_out_stats.txt

samtools flagstat /center1/FISH694/rmlekanoff/trinity_hisat2_out.bam > trinity_hisat2_out_stats.txt

```

Some additional qc with qualimap on the Manila clam reference genome assembly.
No annotated genome for Leukoma, currently, so this step is only done on the Leukoma RNAseq data mapped to the Manila clam reference.

```bash
conda activate qualimap
sbatch qualimap.sh
```

```bash
#!/bin/bash

#SBATCH --partition=t1small
#SBATCH --ntasks=2
#SBATCH --tasks-per-node=24

qualimap rnaseq -bam /home/rmlekanoff/final_proj/hisat2_out.bam \
    -gtf /center1/FISH694/rmlekanoff/man-clam/ncbi_dataset/gff/genomic.gff \
    -outdir /center1/FISH694/rmlekanoff/qualimap_out \
    -outfile qualimap_report.pdf \
    -outformat "HTML" \
    -p non-strand-specific

```
