set -o allexport   # export all vars defined from now on
SDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
. $SDIR/ALUBASTER.config.sh  #import variables from ALUBASTER.config.sh

#USAGE:
# pipeline_4git.sh SRR820468 /ccb/salz4-1/cbc_core/Florea_Burns_Jan2020/data/Tophat_cerebellum/ /ccb/salz3/florea/NewJulip/Cerebellum /ccb/salz4-1/cbc_core/Florea_Burns_Jan2020/Test4git/cereb &> test.log
SET=$1             # sample ID;
THDIR=$2         # Path to TOPHAT output directory containing accepted_hits.bam and unmapped.bam files for this sample;
FASTQDIR=$3    # Path to fastqfiles directory;
WORKDIR=$4      # Main dir for this DATASET analysis. For each sample a subdirectory will be created.

t=4                  # No of threaths used with kraken, samtools ..
# Fastq files (make sure their name is in this format or change it here: 
FASTQ1=${FASTQDIR}/${SET}_1.fastq.gz
FASTQ2=${FASTQDIR}/${SET}_2.fastq.gz
# Kraken output will be stored here:
KKDIR=${WORKDIR}/${SET}/Kraken

set +o allexport   
mkdir -p ${WORKDIR}
mkdir -p ${WORKDIR}/${SET}
mkdir -p ${KKDIR}

#------------------------

# date
# echo "=========================="
# echo "NEW RUN"
# echo "$WORKDIR : $SET"
# echo "=========================="
# 
# 
# ## 1. Run kraken
# 
# #comment next lines if just testing sim4db on already generated kraken file :
# #date
# echo "STEP1 - kraken"
# kraken --db $KRAKENDB --threads $t --min-hits 3 --quick --only-classified-output --classified-out ${KKDIR}/${SET}_1.Alu.fasta --output ${KKDIR}/${SET}_1.kraken_out $FASTQ1 2> ${KKDIR}/${SET}_1.kraken.summary 
# wait
# kraken --db $KRAKENDB --threads $t --min-hits 3 --quick --only-classified-output --classified-out ${KKDIR}/${SET}_2.Alu.fasta --output ${KKDIR}/${SET}_2.kraken_out $FASTQ2 2> ${KKDIR}/${SET}_2.kraken.summary
# wait
# echo -e "Done -step1 \n"

cd ${WORKDIR}/${SET}

# ## 2. Clean up the Alu read files further
# date
# echo "STEP2 - sim4db"
# 
# ${SIM4DB}/sim4db -gen $ALUDB -cdna ${KKDIR}/${SET}_1.Alu.fasta -minlength 22 -aligns -o - | ${SCRIPTS}/EM2cov_coords.pl > ./${SET}_1.Alu.sim4db.coords.stats
# wait
# echo "sim4db - R1 DONE"
# 
# ${SIM4DB}/sim4db -gen $ALUDB -cdna ${KKDIR}/${SET}_2.Alu.fasta -minlength 22 -aligns -o - | ${SCRIPTS}/EM2cov_coords.pl > ./${SET}_2.Alu.sim4db.coords.stats
# wait
# echo -e "sim4db - R2 DONE \n"
# 
# ## 3. Partial mappings: Extract and trim unmapped; re-align; generate corresponding NEW concordant and nonconcordant BAM files
# 
# date
# echo "STEP3 -- get_unmapped_readpairs"
# 
# # generates files unmapped_*.{NEW,orig}.fastq.gz in the working dir
# 
# ${SCRIPTS}/get_unmapped_readpairs.pl ${SET} ${KKDIR}/${SET} ${THDIR}/unmapped.bam
# wait 
# 
# if [[ ! -s "${SET}_unmapped_1.NEW.fastq" ]] && [[ ! -s "${SET}_unmapped_2.NEW.fastq" ]]; then echo "ERROR: ${SET}_unmapped_*.NEW.fastq files are EMPTY  --  Exiting pipeline "; exit; fi;
# 
# mkdir -p ./${SET}.unmapped_realign
# 
# #date
# echo "STEP3 -- tophat in ./${SET}.unmapped_realign directory"
# $TOPHAT -o ./${SET}.unmapped_realign -p $t --max-multihits 10 $BWT2IDX  ${SET}_unmapped_1.NEW.fastq ${SET}_unmapped_2.NEW.fastq
# wait

### Mainframe: the same for full alignments:
# echo
# ls -ld $KRAKENDB
# #ls -l ${KKDIR}/${SET}_1.Alu.fastq
# ls -l ${KKDIR}/${SET}_1.Alu.fasta
# ls -l ${SET}_1.Alu.sim4db.coords.stats
# ls -l ${KKDIR}/${SET}_*
# ls -ld ${SET}.unmapped_realign
# ls -l ${SET}_unmapped_1.NEW.fastq
# ls -l ${SET}_unmapped_2.NEW.fastq
echo -e "Done -step3 \n"


## 4. Sort original (and unmapped_realign) files by readid
date
echo "STEP4 - samtools sort"
 
samtools sort -n -@ $t -o ${THDIR}/accepted_hits.sorted-by-name.bam ${THDIR}/accepted_hits.bam
wait
samtools sort -n -@ $t -o ./${SET}.unmapped_realign/accepted_hits.sorted-by-name.bam ./${SET}.unmapped_realign/accepted_hits.bam
wait
echo -e "Done -step4"

## 5. Create the ${SET}.concordant.sam and ${SET}.nonconcordant.sam files (both sorted by name from here on)
##    Create separately for initial set of reads and then for the unmapped_realign reads, then concatenate
##   (because they are sorted by readid)

date
echo "STEP 5  -- get_aligns"
samtools view -h ${THDIR}/accepted_hits.sorted-by-name.bam | ${SCRIPTS}/get_aligns.pl -k $KKDIR/${SET} -d ./${SET} -o ./${SET} -u  # added -u
wait
date
echo "STEP 5  -- samtools view"
samtools view -bT $GENOME ./${SET}.concordant.sam > ./${SET}.concordant.bam
samtools view -bT $GENOME ./${SET}.nonconcordant.sam > ./${SET}.nonconcordant.bam
wait
date
echo "STEP 5  -- unmapped, get_aligns"
samtools view -h ./${SET}.unmapped_realign/accepted_hits.sorted-by-name.bam | ${SCRIPTS}/get_aligns.pl -k $KKDIR/${SET} -d ./${SET} -o ${SET}.unmapped_realign -u
wait
samtools view -bT $GENOME ./${SET}.unmapped_realign.nonconcordant.sam > ./${SET}.unmapped_realign.nonconcordant.bam
samtools view -bT $GENOME ./${SET}.unmapped_realign.concordant.sam > ./${SET}.unmapped_realign.concordant.bam
wait
date
echo "STEP 5  -- merge"
samtools merge ./tmp.concordant.bam ./${SET}.concordant.bam ./${SET}.unmapped_realign.concordant.bam
samtools merge ./tmp.nonconcordant.bam ./${SET}.nonconcordant.bam ./${SET}.unmapped_realign.nonconcordant.bam
wait
date
echo "STEP 5  -- move"
mv ./${SET}.concordant.bam ./${SET}.tmp_orig.concordant.bam
mv ./${SET}.nonconcordant.bam ./${SET}.tmp_orig.nonconcordant.bam
mv ./tmp.concordant.bam ./${SET}.concordant.bam
mv ./tmp.nonconcordant.bam ./${SET}.nonconcordant.bam
wait
echo -e "Done -step5 \n"

## 6. Find annotation anchor exons; used to be -abam 
date 
echo "STEP 6  -- intersectbed"
intersectBed -wo -bed -split -a ./${SET}.nonconcordant.bam -b $EXONS > nonconcordant; cat nonconcordant | ${SCRIPTS}/take_last.pl | sort -k 13 -k 14n -k 15n > ./${SET}.nonconcordant.genes.bed
wait
cut -f13-16,19 ./${SET}.nonconcordant.genes.bed | sort -u > ./${SET}.nonconcordant.exons.bed
cut -f13-16,19 ./${SET}.nonconcordant.genes.bed | sort | uniq -c > ./${SET}.nonconcordant.genes.counts
wait
echo -e "Done -step6 \n"

## 7. extract mates of mapped nonconcordant reads (convert to fasta)
date
echo "STEP 7  -- pull_read_list_from_bam"
samtools view ./${SET}.nonconcordant.bam | ${SCRIPTS}/pull_read_list_from_bam.pl -o ./${SET}.nonconcordant.mates.fastq -d ./${SET}.nonconcordant.fastq \
                                  $FASTQ1 $FASTQ2 --gzip #### -o ./${SET}.nonconcordant.mates.fastq
wait
date
echo "STEP 7  -- 2x fastq2fasta"
${SCRIPTS}/fastq2fasta.pl < ./${SET}.nonconcordant.mates.fastq > ./${SET}.nonconcordant.mates.fa
wait
${SCRIPTS}/fastq2fasta.pl < ./${SET}.nonconcordant.fastq > ./${SET}.nonconcordant.fa
wait
date
echo "STEP 7  -- leaff"
${SIM4DB}/leaff -F ./${SET}.nonconcordant.mates.fa
wait
${SIM4DB}/leaff -F ./${SET}.nonconcordant.fa
wait
echo -e "Done -step7 \n"

## 8. Run the search for signal (*** Might want to break the id by columns: chrom f t tid later)
date
echo "STEP 8  -- bed2groups"
${SCRIPTS}/bed2groups.pl $ANNOT < ./${SET}.nonconcordant.genes.bed | sort -k 1,1 -k 2n -k 3n | cut -f4-6 | ${SCRIPTS}/nr_groups.pl > ./${SET}.nonconcordant.genes.nr.sorted.groups
wait
date
echo "STEP 8  -- search_signal"
#echo "${SCRIPTS}/search_signal.pl ./${SET} ./${SET}.nonconcordant.genes.nr.sorted.groups ./${SET}.nonconcordant.mates.fa $ANNOT $GENOME $SIM4DB 2> ${SET}.log"
${SCRIPTS}/search_signal.pl ./${SET} ./${SET}.nonconcordant.genes.nr.sorted.groups ./${SET}.nonconcordant.mates.fa $ANNOT $GENOME $SIM4DB 2> ${SET}.log  # generates the files ${SET}.{SIGNAL,UNSPLICED,REGION}.sim4db
wait 
echo -e "Done -step8 \n"

# 9. Gather stats, and then apply the filters. 
#  NOTE: It might be necessary to clean this up
#  OLD: cat SRR1284895sim.SIGNAL.sim4db | parse | parse1 | perl parse2 > signal 2> signal.log
date
echo "STEP 9  -- parse_fix_aln; em2stats_SIGNAL"
cat ${SET}.SIGNAL.sim4db | ${SCRIPTS}/parse_fix_aln.pl | ${SCRIPTS}/em2stats_SIGNAL.pl > ${SET}.SIGNAL.stats
wait
date
echo "STEP 9  -- parse_fix_aln; em2stats_UNSPLICED"
cat ${SET}.UNSPLICED.sim4db | ${SCRIPTS}/parse_fix_aln.pl | ${SCRIPTS}/em2stats_UNSPLICED.pl > ${SET}.UNSPLICED.stats
wait
date
echo "STEP 9  -- parse_fix_aln; em2stats_REGION"
cat ${SET}.REGION.sim4db | ${SCRIPTS}/parse_fix_aln.pl | ${SCRIPTS}/em2stats_REGION.pl > ${SET}.REGION.stats
wait
date
echo "STEP 9  -- 3x sort_stats"
## prepare for manual inspection and debug
${SCRIPTS}/sort_stats.pl < ${SET}.SIGNAL.stats > ${SET}.SIGNAL.sorted.stats
wait
${SCRIPTS}/sort_stats.pl < ${SET}.REGION.stats > ${SET}.REGION.sorted.stats
wait
${SCRIPTS}/sort_stats.pl < ${SET}.UNSPLICED.stats > ${SET}.UNSPLICED.sorted.stats
wait
echo -e "Done -step9 \n"

## 10. Apply the two filters

#Filter 1:
### NOTE: Implemented the sim4db match to ALU requirement here, in call_SIGNAL, but could probably go into the search (to avoid spending time with alignments)
date
echo "STEP 10  -- call_SIGNAL"
${SCRIPTS}/call_SIGNAL.pl 0.8 80 10 < ${SET}.SIGNAL.sorted.stats > ${SET}.SIGNAL.c0.8.p80.r10.stats.nr
wait
date
echo "STEP 10  -- call_UNSPLICED"
${SCRIPTS}/call_UNSPLICED.pl 0.8 93 < ${SET}.UNSPLICED.sorted.stats > ${SET}.UNSPLICED.c0.8.p93.stats.nr
wait
date
echo "STEP 10  -- call_REGION"
${SCRIPTS}/call_REGION.pl 0.8 93 < ${SET}.REGION.sorted.stats > ${SET}.REGION.c0.8.p93.stats.nr
wait
# Filter 2:
date
echo "STEP 10  -- pick_read_hapmap"
${SCRIPTS}/pick_read_hapmap.pl $TXPT2GENE \
                    ${SET}.SIGNAL.c0.8.p80.r10.stats.nr \
                    ${SET}.UNSPLICED.c0.8.p93.stats.nr \
                    ${SET}.REGION.c0.8.p93.stats.nr > ${SET}.Selection.c0.8.p80.r10.txt

wait
date
echo "STEP 10  -- pick_read_insert_filter_assemble_new_hapmap"
${SCRIPTS}/pick_read_insert_filter_assemble_new_hapmap.pl -b 2 5 0 2 0.5 -debug -p ${SET}.u.asm \
                                       ${SET}.nonconcordant.genes.bed \
                                       ${SET}.Selection.c0.8.p80.r10.txt \
                                       ${SET}.SIGNAL.c0.8.p80.r10.stats.nr \
                                       ${SET}.REGION.c0.8.p93.stats.nr \
                                       ${SET}.UNSPLICED.c0.8.p93.stats.nr \
                                       > ${SET}.u.asm.all_selection.c0.8.p80.r10.txt
wait
echo -e "Done - step10 \n"

# 11. Run the assemblies
date
echo "STEP 11  -- call_assembly"
for i in 21 23 25 27
do
   ${SCRIPTS}/call_assembly.pl ${SET}.u.asm.all_selection.c0.8.p80.r10.txt \
                    ${SET}.nonconcordant.fa ${SET}.nonconcordant.mates.fa \
                    -gi $TXPT2GENE -ge $GENDEX -gen $GENOME -p ${SET}.u.asm.k${i} ${i} 2> ${SET}.u.asm.k${i}.log
done
wait
echo -e "Done step11 \n"

## 12. Screen for candidate reads

## Part 1: Identify (exon,Alu) junction reads:
##
## Start from all reads in the Selection file (YES), since we know they
## (the mates) have both a match on the exon and a match on the Alu; 
## Then, Filter those that appear both in the reads and the mates fastq file
## (i.e., remove those where the anchor read itself contains some repeat portion).
date
echo "STEP 12  -- screen_yes"
grep Yes ./${SET}.Selection.c0.8.p80.r10.txt |  cut -d ' ' -f1 > screen.yes
wait
echo "STEP 12  -- filter_R_fastq"
cat screen.yes | ${SCRIPTS}/filter_R_fastq.pl ./${SET}.nonconcordant.mates.fastq > screen.unique
wait
echo "STEP 12  -- em2stats_SIGNAL_screen" 
## Collect information about composition of mates from the SIGNAL sim4db file.
cat ./${SET}.SIGNAL.sim4db | ${SCRIPTS}/parse_fix_aln.pl | ${SCRIPTS}/em2stats_SIGNAL_screen.pl > ${SET}.SIGNAL.screen.stats
wait
date
echo "STEP 12  -- screen_SIGNAL ; add_gene_name"
## Screen SIGNALs to detect and prioritize potential candidates.
cat ${SET}.SIGNAL.screen.stats | ${SCRIPTS}/screen_SIGNAL.pl -s screen.unique | ${SCRIPTS}/add_gene_name.pl $TXPT2GENE > ${SET}.SIGNAL.unique.screen.stats
wait
### Part 2: Refine the alignment, starting from the SIGNAL.unique.screen.stats file
##
## Extract the fasta sequences of the reads, in the order given
date
echo "STEP 12 part 2; output: *.SIGNAL.screened.stats"
grep -v "::" ${SET}.SIGNAL.unique.screen.stats > ${SET}.SIGNAL.screened.stats
wait
## run run_parse.pl to process each line at a time, to avoid any interferences ...
date
echo "STEP 12  -- run_parse"
${SCRIPTS}/run_parse.pl ${SET} $GENOME $GENDEX $ALUFA ${SET}.nonconcordant.mates.fa ${SET}.SIGNAL.screened.stats ${SET}.screened.out 
wait

date
echo "ALL DONE !!!"
exit;

## then adjust manually; note that RAC and rAC are faulty; for future, to correct must check in the main (search_signal.pl) code
#cp ${SET}.screened.out ${SET}.screened.clean.out
#exit;
