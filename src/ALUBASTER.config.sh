##  Software  
##------------

# Alubaster pipeline (this script) location 
SCRIPTS=/ccb/salz4-1/cbc_core/Florea_Burns_Jan2020/Pipeline4git/software      #/path/to/Alubaster/scripts/directory
# /scratch/groups/lflorea1/shared/sw/src/sim4db-r2008/Linux-amd64/bin
SIM4DB=/ccb/salz4-1/cbc_core/scripts/sim4db/sim4db-r2008/Linux-amd64/bin  # /path/to/sim4db/bin/directory
# /path/to/oases/executable   
#/scratch/groups/lflorea1/projects/Repeat/Alubaster/software/oases/oases
OASES=${SCRIPTS}/oases/oases
# /path/to/velveth/executable
VELVETH=${SCRIPTS}/oases/velvet/tvelveth
# /path/to/velvetg/executable
VELVETG=${SCRIPTS}/oases/velvet/velvetg
#/path/to/tophat/executable
TOPHAT=tophat2 

##  Required Reference files directory 
##-------------------------------------

GENDIR=/ccb/salz4-1/florea/repeat/Data             #/path/to/required/reference/files/forAlubaster

KRAKENDB=$GENDIR/krakendb_alu
ALUDB=$GENDIR/ALUdb.fasta
ALUFA=$GENDIR/ALU

BWT2IDX=/ccb/salz7-data/genomes/hg38/hg38c   #$GENDIR/hg38/hg38c
GENOME=/ccb/salz7-data/genomes/hg38/hg38c.fa     #GENOME=$GENDIR/hg38/hg38c.fa
GENDEX="/ccb/salz7-data/genomes/hg38/hg38c.fa.hdrs"                  #$GENDIR/hg38/hg38c.fa.hdrs
ANNOT=/home/florea/gencode.v22.annotation.gtf      #$GENDIR/gencode.v22.annotation.gtf

TXPT2GENE=$GENDIR/gencode.v22.Txpt2Gene
EXONS=$GENDIR/gencode.v22.mx.singl_exon.bed
GINDEX=$GENDIR/gencode.v22.Txpt2Gene2


