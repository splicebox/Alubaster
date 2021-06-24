## Alubaster
Detection of novel *Alu* exonization events from RNA-seq data.

Described in:

  Florea L, Payer L, Antonescu C, Yang G and Burns K. (2021) Detection of *Alu* exonization events in human frontal cortex from RNA-seq data, *Submitted*.

Supplementary data described in the article can be found [here](http://ccb.jhu.edu/software/Alubaster/data/).

```
Copyright (C) 2018-2021 and GNU GPL v3.0 by Liliana Florea, Lindsay Payer
```

### Synopsis
Alubaster identifies candidate gene loci by locating mapped reads (‘anchors’) whose *Alu*-containing mates could not be found in the genome. It then applies two types of filters to select a more accurate subset of loci. The first, a *signal filter*, identifies read evidence for an Alu exonization event, by searching the mate’s sequence against a concatenation of the neighboring exons’ and *Alu* sequences while concomitantly ruling out false positive matches that could have resulted from local or more distant *Alu* elements in the genome. The second, a *context filter*, evaluates the likelihood of an event based on the strength of the signal versus the local context, in particular repeat content and proportion of signal-to-‘context’ matches. 

### Installation
Alubaster is written in Perl and shell. To download the source code, clone the current GitHub repository:

```
git clone https://github.com/splicebox/Alubaster.git
```

#### Prerequisites
Required bioinformatics packages: [sim4db](https://sourceforge.net/projects/kmer), [oases](https://github.com/dzerbino/oases), [tophat2](https://github.com/infphilo/tophat) and [kraken](https://github.com/DerrickWood/kraken). A copy of the 'oases' program is included with this software. Follow the instructions for each program to install and compile, then update the paths in the file 'ALUBASTER.config.sh'.

#### Before you start
Before you start your first Alubaster run, you will need to prepare several reference data files:

1. Download the GENCODE gene annotation file (version of your choice) from [gencodegenes.org](http://gencodegenes.org) (GTF file). Save this file in a directory of your choice.
2. Generate a header file for the genome fasta file, ideally in the same directory as the genome:
```
grep "^>" genome.fa > genome.fa.hdrs
```

Next, perform the following steps in the reference directory:

3. Generate the annotation index:
```
 perl make_annot_index.pl  < gencode.vXX.annotation.gtf > gencode.vXX.annotation.Txpt2Gene
```
4. Generate the exons file, for instance by converting the GTF annotation file to multi-block BED, filtering out single exon transcripts, and then extracting the individual exons in a BED file:
```
gtf2bed gencode.vXX.annotation.gtf  |  awk '{ if ($10!=1) print $_; }' | bedtools bed12tobed6 > gencode.vXX.mx.singl_exon.bed
```
or you can use your own scripts. The format for the output file is:
```
# chrom    from    to   enstranscriptID     0      strand
chr1    11868   12227   ENST00000456328.2       0       +
```
5. You should find a copy of the kraken Alu sequence database in the ‘reference’ directory. If not already there, for instance because of file size restrictions, you can download a copy [here](https://ccb.jhu.edu/software/Alubaster/data/). Copy this in the 'reference' directory, or elsewhere, and extract the files, for instance with 'gunzip' and 'tar -xvf'.
6. Update the paths in the 'ALUBASTER.config.sh' file to point to your local files and software.


### Usage
```
  runAlubaster.sh <SampleName> <TophatDir> <FastqDir> <OutputDir>

  Required parameters:
  
  <SampleName>  Sample name as it apears in fastq files (i.e., sample name would be ABC if fastq files are ABC_{1,2}.fastq.gz)
  <TophatDir>   Path to a directory containing tophat output ('accepted_hits.bam' and 'unmapped.bam') for this sample
  <FastqDir>    Path to fastq files directory
  <OutputDir>   Path to a directory where the analysis is done.  If it does not exist, it will be created.
                For each sample, a subdirectory will be created and all the out files for this run will be written there;
                e.g., for sample ABC the results will be written to <OutputDir>/ABC/
```

### Terms of use
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received (LICENSE.txt) a copy of the GNU General Public License v3.0 along with this program; if not, you can obtain one from http://www.gnu.org/licenses/gpl.txt or by writing to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
