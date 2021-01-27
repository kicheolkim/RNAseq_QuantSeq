#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.read = "/data/fastq/*.fastq.gz"
params.gtf = "/data/genome/gencode.vM25.primary_assembly.annotation.gtf"
params.staridx = "/data/genome/STARindex"


process fastqc {
  tag "FASTQC for raw reads"

  input:
    tuple val(sample_id), file(reads)
  output:
    file("fastqc_${sample_id}_logs")

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -t ${task.cpus} -q ${reads} -f fastq -o fastqc_${sample_id}_logs
  """
}

process multiqc_raw {
  publishDir "$baseDir/qc_fastq_raw", mode:'copy'

  input:
  file('*')

  output:
  file("multiqc_*")

  script:
  """
  multiqc .
  """
}

process trimming {
  tag "bbduk - adaptor and polyA trimming"

  input:
    tuple val(sample_id), file(reads)
  output:
    tuple val(sample_id), file("cleaned-*.fq.gz")

  script:
  """
  bbduk.sh -Xmx2g in=${reads} out=cleaned-${sample_id}.fq.gz ref=/Users/miniconda3/pkgs/bbmap-38.87-hf29c6f4_0/opt/bbmap-38.87-0/resources/polyA.fa.gz,/Users/miniconda3/pkgs/bbmap-38.87-hf29c6f4_0/opt/bbmap-38.87-0/resources/truseq_rna.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20
  """
}

process mapping {
  tag "STAR - genome mapping"
  publishDir "$baseDir", mode:'copy', pattern:'*_Aligned.sortedByCoord.out.bam'
  publishDir "$baseDir", mode:'copy', pattern:'*_Aligned.sortedByCoord.out.bam.bai'
  publishDir "$baseDir", mode:'copy', pattern:'*_Log.final.out'

  input:
    tuple val(sample_id), file(reads)
  output:
    file("*_Aligned.sortedByCoord.out.bam")

  script:
  """
  STAR --runThreadN ${task.cpus} --genomeDir $params.staridx --readFilesCommand gunzip -c --readFilesIn ${reads} --outFileNamePrefix ${sample_id}_ --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate
  samtools index -@ ${task.cpus} ${sample_id}_Aligned.sortedByCoord.out.bam
  """
}

process counting {
  tag "featureCount - gene counting"
  publishDir "$baseDir", mode:'copy', pattern:"*"

  input:
    file(bam)
  output:
    file('*')

  script:
  """
  featureCounts -T ${task.cpus} -s 1 -t exon -g gene_id -a $params.gtf -o gene_featureCounts_output.txt ${bam}
  """
}


workflow {
  read_ch = channel
    .fromPath( params.read )
    .map { file -> tuple(file.simpleName, file) }

    read_ch
      fastqc(read_ch)
      multiqc_raw(fastqc.out.collect())

    read_ch
      trimming(read_ch)
      mapping(trimming.out)
      counting(mapping.out.collect())
}
