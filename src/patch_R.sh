TISSUE=Lung_test
BASE=/scratch/groups/lflorea1/projects/Repeat/Alubaster
WORKDIR=$BASE/$TISSUE
SCRIPTS=$BASE/software

ALUBED=/scratch/groups/lflorea1/shared/genomes/forAlubaster/Alu_repeats_hg38.bed

for i in SRR1072367
do
  cd ${WORKDIR}/${i}
  ${SCRIPTS}/pick_read_insert_filter_R_assemble_new_hapmap.pl -b 2 5 0 2 0.5 \
                                       ./${i}.nonconcordant.genes.bed \
                                       ./${i}.Selection.c0.8.p80.r10.txt \
                                       ./${i}.SIGNAL.c0.8.p80.r10.stats.nr \
                                       ./${i}.REGION.c0.8.p93.stats.nr \
                                       ./${i}.UNSPLICED.c0.8.p93.stats.nr \
                                     > ./${i}.tmp_all_selection.c0.8.p80.r10.txt
  ${SCRIPTS}/parse_selection.pl < ./${i}.tmp_all_selection.c0.8.p80.r10.txt > ./${i}.tmp_selection.bed
  intersectBed -wo -bed -split -a ./${i}.tmp_selection.bed -b $ALUBED > ./${i}.selection_Alu_overlaps.bed
  $SCRIPTS/modify_selection.pl ./${i} < ./${i}.tmp_all_selection.c0.8.p80.r10.txt > ./${i}.u.asm.all_selection.c0.8.p80.r10.R.txt
  rm ./${i}.tmp_all_selection.c0.8.p80.r10.txt
done
