SET=$1  #  SRR1284895  # $1

BASE=/ccb/salz4-1/florea/repeat/
DATADIR=$BASE/Neuro/Data
KKDIR=$BASE/Neuro/Kraken
THDIR=$BASE/Neuro/Tophat
WORKDIR=$BASE/Pipeline


FASTQ1=$DATADIR/${SET}_1.fastq.gz
FASTQ2=$DATADIR/${SET}_2.fastq.gz

KRAKENDB=/ccb/salz4-1/cbc_core/Utils/krakendb_alu/

GENOME=/ccb/salz7-data/genomes/hg38/hg38c.fa
ANNOT=/home/florea/gencode.v22.annotation.gtf
TXPT2GENE=$BASE/Data/gencode.v22.Txpt2Gene
EXONS=$BASE/Data/gencode.v22.mx.singl_exon.bed

### Step 0. Run kraken

#mkdir -p $KKDIR
#mkdir -p $KKDIR/$SET

#kraken --db $KRAKENDB --threads 8 --min-hits 3 --quick --only-classified-output --classified-out $KKDIR/$SET/${SET}_1.Alu.fastq --output $KKDIR/$SET/${SET}_1.kraken_out --gzip_compressed $FASTQ1 2> $KKDIR/$SET/${SET}_1.kraken.summary &
#kraken --db $KRAKENDB --threads 8 --min-hits 3 --quick --only-classified-output --classified-out $KKDIR/$SET/${SET}_2.Alu.fastq --output $KKDIR/$SET/${SET}_2.kraken_out --gzip_compressed $FASTQ2 2> $KKDIR/$SET/${SET}_2.kraken.summary

#fastq2fasta.pl < $KKDIR/$SET/${SET}_1.Alu.fastq > $KKDIR/$SET/${SET}_1.Alu.fasta
#fastq2fasta.pl < $KKDIR/$SET/${SET}_2.Alu.fastq > $KKDIR/$SET/${SET}_2.Alu.fasta
 
#### Start the actual pipeline

mkdir -p $WORKDIR/$SET
cd $WORKDIR/$SET

## 1. Clean up the Alu read files further
##### NOTE: Make minid 75 or 70, or otherwise just accept any match; did this by removing -minid 80 in the sim4db call below
#nohup sim4db -gen $BASE/Data/ALUdb.fasta -cdna $KKDIR/$SET/${SET}_1.Alu.fasta -minlength 22 -aligns -o - | EM2cov.pl > ./${SET}_1.Alu.sim4db.stats &
#nohup sim4db -gen $BASE/Data/ALUdb.fasta -cdna $KKDIR/$SET/${SET}_2.Alu.fasta -minlength 22 -aligns -o - | EM2cov.pl > ./${SET}_2.Alu.sim4db.stats
#wait

## 2. Sort original file by readid
#samtools sort -n -o $THDIR/$SET/accepted_hits.sorted-by-name.bam $THDIR/$SET/accepted_hits.bam

## 3. Create the $SET.concordant.sam and $SET.nonconcordant.sam files (both sorted by name from here on)
#samtools view -h $THDIR/$SET/accepted_hits.sorted-by-name.bam | get_aligns.pl -k $KKDIR/$SET/$SET ./$SET -u  # added -u
#samtools view -bT $GENOME ./$SET.concordant.sam > ./$SET.concordant.bam
#samtools view -bT $GENOME ./$SET.nonconcordant.sam > ./$SET.nonconcordant.bam

## first must be sorted by coordinate, then index
#samtools sort -o ./$SET.concordant.so.bam ./$SET.concordant.bam ; samtools index ./$SET.concordant.so.bam
#samtools sort -o ./$SET.nonconcordant.so.bam ./$SET.nonconcordant.bam ; samtools index ./$SET.nonconcordant.so.bam

## 4. Find annotation anchor exons; used to be -abam 
#intersectBed -wo -bed -split -a ./$SET.nonconcordant.bam -b $EXONS > nonconcordant; cat nonconcordant | take_last.pl | sort -k 13 -k 14n -k 15n > ./$SET.nonconcordant.genes.bed
#cut -f13-16,19 ./$SET.nonconcordant.genes.bed | sort -u > ./$SET.nonconcordant.exons.bed
#cut -f13-16,19 ./$SET.nonconcordant.genes.bed | sort | uniq -c > ./$SET.nonconcordant.genes.counts

## 5. extract mates of mapped nonconcordant reads (convert to fasta)
#samtools view ./$SET.nonconcordant.bam | pull_read_list_from_bam.pl -o ./$SET.nonconcordant.mates.fastq -d ./$SET.nonconcordant.fastq --gzip \
#                                        $FASTQ1 $FASTQ2 --gzip -o ./$SET.nonconcordant.mates.fastq
#fastq2fasta.pl < ./$SET.nonconcordant.mates.fastq > ./$SET.nonconcordant.mates.fa
#fastq2fasta.pl < ./$SET.nonconcordant.fastq > ./$SET.nonconcordant.fa
#leaff -F ./$SET.nonconcordant.mates.fa
#leaff -F ./$SET.nonconcordant.fa
#
## Example then of retrieval by id:
## leaff -F NA12814.nonconcordant.mates.fa -s ERR188081.1393181

# 6. Run the search for signal (*** Might want to break the id by columns: chrom f t tid later)
#bed2groups.pl $ANNOT < ./$SET.nonconcordant.genes.bed | sort -k 1,1 -k 2n -k 3n | cut -f4-6 | nr_groups.pl > ./$SET.nonconcordant.genes.nr.sorted.groups


# Notes: the only reason for including the overlap file is to get the orientation of the match; maybe it could be obtained and recorded earlier, in bed2groups? Done.
#search_signal.pl ./$SET ./$SET.nonconcordant.genes.nr.sorted.groups ./$SET.nonconcordant.mates.fa $ANNOT $GENOME 2> $SET.log  # generates the files $SET.{SIGNAL,UNSPLICED,REGION}.sim4db  

#exit;

# It might be necessary to clean this up: cat SRR1284895sim.SIGNAL.sim4db | perl parse | perl parse1 | perl parse2 > signal 2> signal.log
em2stats_SIGNAL.pl < $SET.SIGNAL.sim4db > $SET.SIGNAL.stats
em2stats_REGION.pl < $SET.REGION.sim4db > $SET.REGION.stats
em2stats_UNSPLICED.pl < $SET.UNSPLICED.sim4db  > $SET.UNSPLICED.stats

# prepare for manual inspection and debug
sort_stats.pl < $SET.SIGNAL.stats > $SET.SIGNAL.sorted.stats
sort_stats.pl < $SET.REGION.stats > $SET.REGION.sorted.stats
sort_stats.pl < $SET.UNSPLICED.stats > $SET.UNSPLICED.sorted.stats

# apply filters
##echo "cov:0.8 pctid:93 "; perl call_REGION_new.pl 0.8 93 < REGION3.sorted.stats | cut -d ' ' -f4 | sort | uniq -c
call_REGION.pl 0.8 93 < $SET.REGION.sorted.stats > $SET.REGION.c0.8.p93.stats.nr
call_UNSPLICED.pl 0.8 93 < $SET.UNSPLICED.sorted.stats > $SET.UNSPLICED.c0.8.p93.stats.nr
call_SIGNAL.pl 0.8 80 10 < $SET.SIGNAL.sorted.stats > $SET.SIGNAL.c0.8.p80.r10.stats.nr

pick_read_hapmap.pl $TXPT2GENE $SET.SIGNAL.c0.8.p80.r10.stats.nr \
                            $SET.UNSPLICED.c0.8.p93.stats.nr \
                            $SET.REGION.c0.8.p93.stats.nr > $SET.Selection.c0.8.p80.r10.txt

pick_read_insert_filter_new_hapmap.pl -b 2 5 0 2 0.5 -debug -p $SET.u \
                                       $SET.nonconcordant.genes.bed \
                                       $SET.Selection.c0.8.p80.r10.txt \
                                       $SET.SIGNAL.c0.8.p80.r10.stats.nr \
                                       $SET.REGION.c0.8.p93.stats.nr \
                                       $SET.UNSPLICED.c0.8.p93.stats.nr \
                                       > $SET.u.all_selection.c0.8.p80.r10.txt

exit;

##perl pick_read.pl | tr ';' ' ' | cut -d ' ' -f1,3,5 | sort -k 1,1 -k 2,2 | uniq | sort -k 3 > Selection.short.txt
###cut -d ' ' -f1,3 Selection.short.txt | sort | uniq -c | sort -nrk1 | more

###
