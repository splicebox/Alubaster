##  Software
##------------

# NOTE: oases is included with the distribution
SCRIPTS=/path/to/Alubaster/
SIM4DB=/path/to/sim4db/bin/directory
OASES=${SCRIPTS}/oases/oases                 # copy your oases software here
VELVETH=${SCRIPTS}/oases/velvet/tvelveth
VELVETG=${SCRIPTS}/oases/velvet/velvetg
TOPHAT=tophat2

##  Required reference files directory
##-------------------------------------
REFDIR=$SCRIPTS/reference

# Alu databases and consensus Alu sequences
ALUDB=$REFDIR/ALUdb.fasta
ALUFA=$REFDIR/ALU
KRAKENDB=$REFDIR/krakendb_alu                # if not provided with the software, you can find a link to a copy from the README.md file;
                                             # then 'gunzip' an d 'tar -xvf' the file

# genes, genomes and indices
BWT2IDX=/path/to/genome/dir/bowtie2prefix    # e.g., /data/genomes/hg38/hg38
GENOME=/path/to/genome/genomefasta           # e.g., /data/genome/hg38/hg38.fa
GENDEX=/path/to/genome/genomefasta.hdrs      # generate this with 'grep "^>" genome.fa > genome.hdrs' ; e.g., /data/genome/hg38/hg38.fa.hdrs
ANNOT=/path/to/gene/annotation/gencodegtf    # e.g., /data/genome/hg38/annotation/gencode.v22.annotation.gtf

# in the reference directory
TXPT2GENE=$REFDIR/gencode.vXX.Txpt2Gene      # generate this with 'perl ./get_annot_index.pl < gencode.v36.annotation.gtf > gencode.v36.Txpt2Gene'
EXONS=$REFDIR/gencode.vXX.mx.singl_exon.bed  # this file contains the exons in multi-exon transcripts, one per line; generate this with
                                             # 'gtf2bed gencode.vXX.annotation.gtf  |  awk '{ if ($10!=1) print $_; }' | bedtools bed12tobed6 > gencode.vXX.mx.singl_exon.bed'
                                             # or using your own scripts. Expected format:
                                             # chr1    11868   12227   ENST00000456328.2       0       +


