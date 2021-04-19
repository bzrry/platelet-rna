#! /bin/bash

set -e

if [ ! -f "hg38_Gencode_V28.bed" ]; then
    echo "Download and unzip hg38_Gencode_V28.bed.gz from sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_Gencode_V28.bed.gz/download"
    exit 1
fi

echo "Downloading references"
curl http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > hg38.fa
curl http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz | gunzip -c > hg38.ensGene.gtf

echo "Generating genome index"
STAR \
  --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir human_STARindex \
  --genomeFastaFiles hg38.fa \
  --sjdbGTFfile hg38.ensGene.gtf \
  --sjdbOverhang 100 \
  --genomeSAsparseD 2 \
  --genomeSAindexNbases 10 \
  --limitGenomeGenerateRAM 32212254720

echo "Processing samples:"
cat $1

cat $1 | while IFS=, read -r srr sample
do
  echo "Downloading $sample"
  fastq-dump --stdout $srr | gzip > fastqs/${sample}.fq.gz
done

cat $1 | while IFS=, read -r srr sample
do
  echo "Trimming $sample"
  trim_galore --illumina --phred33 --stringency 5 --cores 4 -o trimmed_fastqs fastqs/${sample}.fq.gz
done

cat $1 | while IFS=, read -r srr sample
do
  echo "Running fastqc on $sample"
  fastqc trimmed_fastqs/${sample}_trimmed.fq.gz --extract
done

cat $1 | while IFS=, read -r srr sample
do
  echo "Aligning $sample"
  STAR \
    --runMode alignReads \
    --runThreadN 7 \
    --genomeDir human_STARindex \
    --readFilesIn trimmed_fastqs/${sample}_trimmed.fq.gz \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix aligned/${sample}_trimmed. \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All
done

cat $1 | while IFS=, read -r srr sample
do
  echo "Indexing $sample"
  samtools index aligned/${sample}_trimmed.Aligned.sortedByCoord.out.bam
done

cat $1 | while IFS=, read -r srr sample
do
  echo "Running geneBody_coverage on $sample"
  geneBody_coverage.py -r hg38_Gencode_V28.bed -i aligned/${sample}.Aligned.sortedByCoord.out.bam -o rseqc/gbc.${sample}
done

cat $1 | while IFS=, read -r srr sample
do
  echo "Running read_distribution on $sample"
  read_distribution.py -r hg38_Gencode_V28.bed -i aligned/${sample}.Aligned.sortedByCoord.out.bam > rseqc/${sample}.read_distribution.txt
done

featureCounts \
  -T 7 \
  --minOverlap 10 \
  -a hg38.ensGene.gtf \
  -o feature-counts.txt \
  aligned/*.bam

multiqc .


