# RNAseq_QuantSeq
- Lexogen's QuantSeq (3'mRNA-Seq) analysis pipeline
- Command line for each processing steps are described in below. This script follows and modified from Lexogen's data analysis guide (https://www.lexogen.com/quantseq-data-analysis/).
- Please see the script for the nextflow.

## Install required software
```
conda install -c bioconda star subread multiqc samtools cutadapt
conda install nextflow fastqc bbmap
```

### Convert GTF to bed
```
cat gencode.vM25.primary_assembly.annotation.gtf | awk 'BEGIN{FS="\t";OFS="\t"}{if($1 ~ /^ch/){split($9,a,";");print $1,$4-1,$5,a[1],"0",$7,$4-1,$4-1,"255,0,0",1,$5-$4,$4-1}else{print $0}}' | sed 's/gene_id //g' | tr -d '"' > gencode.vM25.primary_assembly.annotation.bed
```

### Create STAR index
```
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./STARindex --genomeFastaFiles GRCm38.primary_assembly.genome.fa --sjdbGTFfile gencode.vM25.primary_assembly.annotation.gtf --sjdbOverhang 49 --outFileNamePrefix ./STARindex/GRCm38_50bp
```

---
## Workflow with nextflow script
Run pipeline with nextflow (https://www.nextflow.io/)

```
nextflow run nf_QuantSeq.sh --read '/data/fastq/*.fastq.gz' --genome '/data/genome/STARindex'
```

---
## Workflow in the command line (Based on Lexogen's analysis guide)
Current workflow consider 50bp single-end sequencing and gencode mouse genome.

1. Run FastQC for all samples
```
for sample in *.fastq.gz; do cat ${sample} | fastqc ${sample} -o qc -t 2; done
```

2. BBDuk - remove the adapter contamination, polyA read through, and low quality tails
```
for fileName in *.fastq.gz; do sample=${fileName%.*.*}; bbduk.sh -Xmx2g in=${fileName} out=cleaned/${sample}.fq.gz ref=/Users/kicheol/miniconda3/pkgs/bbmap-38.87-hf29c6f4_0/opt/bbmap-38.87-0/resources/polyA.fa.gz,/Users/kicheol/miniconda3/pkgs/bbmap-38.87-hf29c6f4_0/opt/bbmap-38.87-0/resources/truseq_rna.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 &>> bbduk_log.txt; done
```

3. STAR mapping
```
mkdir /data/star_out
for sample in sample1 sample2 sample3 ; do \
STAR --runThreadN 8 --genomeDir /data/gencode_M25/STARindex --readFilesCommand gunzip -c --readFilesIn /data/fastq/"$sample".fq.gz \
--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix /data/star_out/"$sample"/"$sample"_ ;\
done
```

4. Indexed bam files are necessary for many visualization and downstream analysis tools
```
cd star_out
for bamfile in /data/star_out/*/*_Aligned.sortedByCoord.out.bam ; do samtools index -@ 8 ${bamfile}; done
```

5. featureCount for gene counting
```
featureCounts -T 8 -s 1 -t exon -g gene_id -a /data/gencode_M25/gencode.vM25.primary_assembly.annotation.gtf -o gene_featureCounts_output.txt ./sample1/sample1_Aligned.sortedByCoord.out.bam ./sample2/sample2_Aligned.sortedByCoord.out.bam ./sample3/sample3_Aligned.sortedByCoord.out.bam
```

6. Further analysis (differential expression)


---
[Created by Kicheol Kim - Updated Jan 26. 2021]
