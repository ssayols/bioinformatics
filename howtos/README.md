## Table of Contents

* [Convert BAM to BigWig](#convert-bam-to-bigwig)
* [Convert BAM to BigWig with R](#convert-bam-to-bigwig-with-r)
* [Strand specific tracks](#strand-specific-tracks)
* [LogFC tracks](#logfc-tracks)
* [Convert GTF to BED](#convert-gtf-to-bed)
* [Deduplicate UMIs](#deduplicate-umis)
* [Downsample](#downsample)
   * [Fastq](#fastq)
   * [BAM](#bam)
   * [BAM (other)](#bam-other)
* [Merge BAM files](#merge-bam-files)
* [Merge BigWig tracks](#merge-bigwig-tracks)
* [Retrieve chromosome sizes](#retrieve-chromosome-sizes)
* [LiftOver coordinates between different assemblies](#liftover-coordinates-between-different-assemblies)
* [Submit jobs to random machines in the cluster](#submit-jobs-to-random-machines-in-the-cluster)
* [Count reads on repetitive regions](#count-reads-on-repetitive-regions)
* [Predict genes network with genemania](#predict-genes-network-with-genemania)
* [Annotation of principal isoforms from Biomart](#annotation-of-principal-isoforms-from-biomart)
* [Match a pattern in the genome](#match-a-pattern-in-the-genome)
* [Reverse complement a DNA sequence](#reverse-complement-a-dna-sequence)
* [Query the Uniprot Rest API](#query-the-uniprot-rest-api)
* [Calculate best primers using melting temperature](#calculate-best-primers-using-melting-temperature)
* [Broad Institute GSEA](#broad-institute-gsea)
* [Problematic and blacklisted regions from UCSC](#problematic-and-blacklisted-regions-from-UCSC)
* [Signal summary around regions of interest](#signal-summary-around-regions-of-interest)

## Convert BAM to BigWig
This script will loop over the BAM files in a directori, and convert them to BigWig files or visualization in Genome Browsers.

To queue it into LSF, the for loop has to be moved into a master script that will queue with bsub this script passing the BAM file to be processed as a parameter.

**Requirements**

* bedtools
* Picard tools
* bedGraphToBigWig from the UCSC Genome Browser tools
* The cromosome sizes, which can be retrieved as described [here](#retrieve-chromosome-sizes).

**Source**

```bash
#!/bin/bash
BEDTOOLS=/opt/bedtools-2.17.0/bin
PICARD=/opt/picard-tools-1.103
BEDGRAPHBW=/opt/ucsc/bedGraphToBigWig
GENOME=./mm9.chrom.sizes

for BAM in *.bam
do
        ## scale by Reads per 1000 Mapped Reads
        TOTAL_MAPPED=$(java -Xmx5000m -jar ${PICARD}/BamIndexStats.jar I=${BAM} | awk '{ SUM += $5 } END { print SUM }')
        echo "TOTAL MAPPED READS:" ${TOTAL_MAPPED} | tee scale_log.txt
        SCALE=$(echo "1000000/$TOTAL_MAPPED" | bc -l)
        echo SCALE: $SCALE / CALCULATING COVERAGE IN $BAM | tee scale_log.txt;
        ${BEDTOOLS}/genomeCoverageBed -bg -split -scale ${SCALE} -ibam $BAM -g $GENOME > ${BAM%_sorted.bam}_scaled.bedgraph
        echo CONVERT $BAM to BIGWIG;
        $BEDGRAPHBW ${BAM%_sorted.bam}_scaled.bedgraph $GENOME ${BAM%_sorted.bam}_scaled.bw
done
```

## Convert BAM to BigWig with R

This method doesn't include normalization to total coverage.

```R
rtracklayer::export.bw(con="sample.bw",               ## export to bigwig
  GenomicAlignments::coverage(                        ## convert to coverages
    GenomicAlignments::readGAlignments("sample.bam")  ## read in BAM file (use readGAlignmentPairs for paired-end files)
  )
)
```

## Strand specific tracks

The trick here's to use `deeptools`'s `bamCoverage` with `--filterRNAstrand` to get only forward/reverse reads into 2 separate tracks. 

```sh
module load deepTools/3.1.0_debian9

bamCoverage --filterRNAstrand forward --numberOfProcessors 16 --outFileFormat bigwig --binSize 1 --skipNonCoveredRegions --normalizeUsing CPM --bam polyAOligotex.bam -o polyAOligotex.forward.bw &
bamCoverage --filterRNAstrand reverse --numberOfProcessors 16 --outFileFormat bigwig --binSize 1 --skipNonCoveredRegions --normalizeUsing CPM --bam polyAOligotex.bam -o polyAOligotex.reverse.bw &

wait $(jobs -p)
```

## LogFC tracks

Using Bioconductor, simply tile the genome into a number of bins of a specific width, and calculate the log FC:

```R
bins <- tileGenome(seqinfo(BSgenome.Mmusculus.UCSC.mm10::Mmusculus), tilewidth=1e6, cut.last.tile=TRUE)
bins <- keepStandardChromosomes(bins, pruning.mode="coarse")
seqlevelsStyle(bins) <- "UCSC"

sub_igg_bin <- function(ip_file, igg_file) {
  ip_log2_igg_file <- sub(".bw$", ".1M_smoothed.log2_igg.bw", ip_file)

  # read tracks
  ip  <- keepStandardChromosomes(coverage(import.bw(file.path("tracks/filtered", ip_file )), weight="score"))
  igg <- keepStandardChromosomes(coverage(import.bw(file.path("tracks/filtered", igg_file)), weight="score"))
  
  # smooth
  ip  <- binnedAverage(bins, ip , "score")$score
  igg <- binnedAverage(bins, igg, "score")$score

  # calculate ratio adding 1e-6 as pseudocount
  bins$score <- log2(1e-6 + ip) - log2(1e-6 + igg)

  export.bw(bins, ip_log2_igg_file)
}
```

A version at 1-bp resolution (without bins), would be possible simple avoiding the call to `binnedAverage()` and calculating the logFC directly on the `coverage()` objects.

## Convert GTF to BED

See https://github.com/ExpressionAnalysis/ea-utils/blob/master/clipper/gtf2bed

```bash
gtf2bed.pl Saccharomyces_cerevisiae.R64-1-1.86.gtf > Saccharomyces_cerevisiae.R64-1-1.86.bed
```

## Deduplicate UMIs
Suposing we have the UMIs and the sequences demultiplexed in separate reads (fastq files).

**Merge.pl**

Read both files, and write output in the format:

<pre>
UMI  SEQ  Q1  Q2  +  QUAL  HEAD
</pre>

Were Q1=(sum of bases with quality>30) and Q2=(sum of all individual base qualities).

A perl approach would look like:

```perl
open FQ , "gunzip -c $ARGV[0] |";
open UMI, "gunzip -c $ARGV[1] |";
while(!eof(FQ) and !eof(UMI)) {
    @fq=split /\n/, <FQ>.<FQ>.<FQ>.<FQ>;
    @umi=split /\n/, <UMI>.<UMI>.<UMI>.<UMI>;
    $q1=eval join("+", map { (ord($_) - 33) > 30 || 0 } split //, $fq[3]);
    $q2=eval join("+", map { ord($_) } split //, $fq[3]);
    print join "\t", ($umi[1], $fq[1], $q1, $q2, $fq[0], $fq[2], $fq[3]), "\n";
}
```

**Split.pl**
Split the *reversely* sorted result of Merge.pl

```perl
my($key1, $key2) = ('', '');
while($l = <>) {
    chomp $l;
    @l = split /\t/, $l;
    if(!($key1 eq $l[0] && $key2 eq $l[1])) {
        print "$l[4]\n$l[1]\n$l[5]\n$l[6]\n";
        ($key1, $key2) = ($l[0], $l[1]);
    }
}
```

**Loop over the input files and run the tasks on LSF**

```perl
#!/bin/bash
hosts=(imbc1 imbc4 imbc5 imbc6)
for f in *R?.fastq.gz; do
    name=${f%.fastq.gz}
    echo "perl dedup_merge.pl $f ${name}.umi.fastq.gz | sort -k1,1 -k2,2 -k3,3nr -k4,4nr | perl dedup_split.pl | gzip > ${name}.dedup.fastq.gz" | bsub -J $name -W5:00 -o ${name}.out -e ${name}.err -m ${hosts[$((RANDOM % 2))]}
done 
```

## Downsample

### Fastq

The code works in parallel and downsamples to 10% of the original reads. It's not random, takes only 1 out of 10 reads from the fastq file.

**Source**
```bash
$ function PROCESS { zcat $1 | perl -e'while(my $l=<>.<>.<>.<>) { if($. % 40 == 0) { print "$l" }}' | gzip > ${1%.fastq.gz}.subset.fastq.gz; }
$ export -f PROCESS
$ find . -name "*.fastq.gz" | parallel -j 4 PROCESS
```

### BAM

This script will loop over the BAM files in a directory and downsample them according to the percentage described in the 3rd column of a csv file.

**Requirements**

* samtools (0.1.18 has a known bug when downsampling with >50% likelihood)

**Source**

```bash
#!/bin/bash

# Input file should look like (sample name,total reads,downsampling prob):
#
#36_HEK293T_REX_shdicer,1962219,10.19
#42_HEK293T_REX_parent,10198296,1.96
#14_HEK293T_shRandom,11939348,1.68

BAMS=..
PREF=Sample_richly_2014_03__
SUFF=.cut.bam

for l in $(cat downsample20M.csv); do
  # get the selection probability (downsampling factor)
  f=$(echo $l|cut -f1 -d,)
  x=$(echo $l|cut -f3 -d,)

  # get the samplename and call samtools to downsample
  if (( $(echo "$x >= 1.0" | bc -l) )); then
    echo "Sample $f cannot be downsampled due to not enough reads. Creating symbolic link instead"
    ln -s ${BAMS}/${PREF}${f}${SUFF} ${PREF}${f}.20M${SUFF}
  else
    xx=$(echo "$x * 100" |bc)
    echo "Sample $f being downsampled ${xx}%"
    echo "samtools view -s ${SEED}.${x} -bh ${BAMS}/${PREF}${f}${SUFF} > ${PREF}${f}.20M${SUFF}" | bsub -cwd $(pwd) -J ${f} -n1 -W2:00
  fi
done
```

## BAM (other)

Downsample to a predefined number of reads.

**Requirements**

* samtools (0.1.18 has a known bug when downsampling with >50% likelihood)

**Source**

```bash
#!/bin/bash
function downsample {
  bam=$1
  bamcoord=${bam%.sortedByName.out.markDups.bam}.sortedByCoord.out.markDups.bam
  counts=$2
  seed=42

  # Calculate the sampling factor based on the intended number of reads:
  prob=$(samtools idxstats $bamcoord | cut -f3 | awk -v SEED=$seed -v COUNT=$counts 'BEGIN {total=0} {total += $1} END {print SEED+COUNT/(total/2)}')
  if [[ $FACTOR > 1 ]]; then
    prob=1
  fi

  # and call samtools to downsample based on the factor calculated
  echo "downsampling $(basename $bam) to $counts with probability $prob"
  samtools view -bs ${prob} $bam > ${bam%.bam}.${counts}.bam
}

for f in ../rawdata/*.sortedByName.out.markDups.bam; do
  for x in 500000 1000000 5000000; do
    downsample $f $x
  done
done
```

## Merge BAM files

This script takes the BAM files from several replicates and merges them together.

The trick is to subsample to the smallest library in order to avoid overrepresentation of one library compared to the rest.

It consists of 3 steps:

* get minimum number of reads mapped
* downsample to the minimum
* merge the downsampled bams

**Source**
```bash
#!/bin/bash
SAMTOOLS=/opt/samtools/1.3/samtools
INPUTDIR=/projects/xxx/with_UMIs/mapped
SAMPLES=(sample_??_WT_0_2h sample_??_WT_24h sample_??_ERC)
SEED=666

cd $INPUTDIR

for SAMPLE in ${SAMPLES[@]}; do
    OUT=$(echo $SAMPLE | sed 's/??_//')
    echo "== $OUT =="

    # get minimum number of reads mapped
    REPS=() # will fill out later
    MAPPEDREADS=()  # will fill out later
    for REP in ${SAMPLE}*.bam; do
        M=$($SAMTOOLS flagstat $REP | grep mapped | cut -f1 -d" " | head -n1)
        MAPPEDREADS+=($M)
        REPS+=($REP)
        echo "$REP has $M reads"
    done
    MIN=$(echo ${MAPPEDREADS[@]} | tr " " "\n" | sort -n | head -n1)

    # downsample to the minimum
    REPSDOWNSAMPLED=()
    for i in $(seq 1 ${#REPS[@]}); do
        X=${REPS[$i-1]}
        REPSDOWNSAMPLED+=(${X%.bam}.subsampled.bam)
        if [ "$MIN" -eq "${MAPPEDREADS[$i-1]}" ]; then
            echo "not downsampling $X"
            $SAMTOOLS view -F4 -bh $X > ${REPSDOWNSAMPLED[$i-1]}
        else
            PROB=$(echo "$MIN / ${MAPPEDREADS[$i-1]}" | bc -l | cut -c2-3)
            echo "downsampling $X to ${PROB}% of reads"
            $SAMTOOLS view -s ${SEED}.${PROB} -F4 -bh $X > ${REPSDOWNSAMPLED[$i-1]}
        fi
    done

    # merge the downsampled bams
    echo "merging ${REPSDOWNSAMPLED[@]} into ${OUT}.bam"
    $SAMTOOLS merge ${OUT}.bam ${REPSDOWNSAMPLED[@]} && rm ${REPSDOWNSAMPLED[@]}
done
```

## Merge BigWig tracks

This script takes the bigwig tracks from several replicates and merges them together.

To queue it into LSF, the for loop has to be moved into a master script that will queue with bsub this script passing the BAM file to be processed as a parameter.

**Requirements**

* bedGraphToBigWig from the UCSC Genome Browser tools
* The cromosome sizes, which can be retrieved as described [here](#retrieve-chromosome-sizes).

**Source**

```bash
#!/bin/bash
UCSC=/opt/ucsc/latest/
INPUTDIR=/projects/xxx/with_UMIs/tracks
SAMPLES=(sample_01_??_WT_0_2h sample_01_??_WT_24h sample_01_??_ERC)
CHRSIZES=/projects/xxx/with_UMIs/chr.sizes  # got it with: samtools idxstats sample.bam | cut -f1-2 > ./chr.sizes

cd $INPUTDIR

for SAMPLE in ${SAMPLES[@]}; do
    echo $SAMPLE

    REP=$(ls ${SAMPLE}*.bw) # get list of replicates to merge
    N=$(echo $REP | wc -w)  # how many?
    NORM=$(echo "1 / $N" | bc -l)   # average each signal by the number of replicates
    OUT=$(echo $SAMPLE | sed 's/??_//')

    ${UCSC}/bigWigMerge $REP stdout \
        | awk -v NORM=$NORM '$4=NORM*$4' \
        | LC_COLLATE=C sort -k1,1 -k2,2n > ${OUT}.bed
    ${UCSC}/bedGraphToBigWig ${OUT}.bed $CHRSIZES ${OUT}.bw && rm ${OUT}.bed
done
```

## Retrieve chromosome sizes

Use the UCSC toolkit from the deps tree.

The chromosome sizes are often used to generate tracks from a bam file.

**Requirements**

* fetchChromSizes from the UCSC tools, or
* a mysql client, or
* curl, or
* R + a BSGenome package, or
* samtools and a bam file

**Source**

```bash
$ fetchChromSizes hg18 > hg18.chrom.sizes
```

Or:

```bash
$ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  > hg19.genome
```

Or:

```bash
$ curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
```

Or:

```bash
Rscript -e 'library(BSgenome.Hsapiens.UCSC.hg38); write.table(seqinfo(BSgenome.Hsapiens.UCSC.hg38), col.names=F, quote=F, sep="\t")' | cut -f1-2 > hg38.chrom.sizes
```

Or:

```
samtools idxstats sample.bam | cut -f1-2 > ./chr.sizes
```

## LiftOver coordinates between different assemblies

Basically, use the UCSC utils to convert the bigWig to bedGraph, lifting it over, and back to bigWig.

```bash
ml kentUtils

# or download them from UCSC's ftp
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v287/bigWigToBedGraph
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v287/liftOver
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v287/bedGraphToBigWig

# get the chain and chr sizes
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz | gzip -cd > hg19ToHg38.over.chain
wget https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes

# and liftover
bigWigToBedGraph sample1.hg19.bw sample1.hg19.bedGraph
liftOver -bedPlus=4 sample1.hg19.bedGraph hg19ToHg38.over.chain sample1.hg38.bedGraph unmapped
sort -k1,1 -k2,2n sample1.hg38.bedGraph > sample1.hg38.sorted.bedGraph
bedGraphToBigWig sample1.hg38.sorted.bedGraph hg19.chrom.sizes sample1.hg38.bw
```

There's on issue that prevents easily lifting over bigWigs, and that's overlapping ranges, which are not allowed. To overcome the issue, [do it in R](https://github.com/ssayols/bioinformatics/tree/master/Rsnippets#liftover-coordinates-between-different-assemblies) using the `rtracklayer` package.

## Submit jobs to random machines in the cluster

If you have high demanding IO jobs, it's advisable to distribute the load across different machines. Otherwise, the system could become unusable.

```bash
MACHINES=(imbc1 imbc4 imbc5)
echo "zcat x1.fastq.gz | fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 15 | gzip > x1.trimmed.fastq.gz" | \
  bsub -J fastx -W5:00 -app Reserve1G -n1 -m ${MACHINES[$((RANDOM % 3))]}
echo "zcat x2.fastq.gz | fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 15 | gzip > x2.trimmed.fastq.gz" | \
  bsub -J fastx -W5:00 -app Reserve1G -n1 -m ${MACHINES[$((RANDOM % 3))]}
echo "zcat x3.fastq.gz | fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 15 | gzip > x3.trimmed.fastq.gz" | \
  bsub -J fastx -W5:00 -app Reserve1G -n1 -m ${MACHINES[$((RANDOM % 3))]}
echo "zcat x4.fastq.gz | fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 15 | gzip > x4.trimmed.fastq.gz" | \
  bsub -J fastx -W5:00 -app Reserve1G -n1 -m ${MACHINES[$((RANDOM % 3))]}
echo "zcat x5.fastq.gz | fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 15 | gzip > x5.trimmed.fastq.gz" | \
  bsub -J fastx -W5:00 -app Reserve1G -n1 -m ${MACHINES[$((RANDOM % 3))]}
```

Please remind that's not the best strategy for high IO load. That should go in combination with:

* using the scratch disks
* LSF doesn't let you book for IO resources. Thus, limit the number of concurrent tasks

## Count reads on repetitive regions

Let's do it in two parts:

* Count the number of TTAGGG repetitions:

```sh
for f in *.fastq.gz; do
    echo $f
    counts=$(zcat $f | awk 'NR%4==2' | perl -aln -e'print $c = () = $F[0] =~ /(TTAGGG|CCCTAA)/gi' | sort | uniq -c | paste -s)
    echo "$f,$counts" >> ttaggg.csv
done
```

* Plot showing the fraction of TTAGGG reads:

```R
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(WriteXLS)

CWD <- "/fsimb/groups/imb-bioinfocf/projects/butter/imb_Butter_2014_03_U2OS_ChIP_51set/"
setwd(CWD)

pdf("./results/ttaggg.pdf")

# prepare data structure
ttaggg <- read.csv("./rawdata/ttaggg.csv", head=F)
files  <- ttaggg$V1
ttaggg <- lapply(ttaggg$V2, function(x) {
    x <- strsplit(gsub("^\\s+|\\s+$", "", unlist(strsplit(x, "\t"))), " ")
    positions <- as.integer(sapply(x, function(x) x[2]))
    counts    <- as.integer(sapply(x, function(x) x[1]))
    out <- numeric(9)   # has to match max(repeats) + 1; check before
    out[positions + 1] <- counts
    out
})

ttaggg <- do.call(rbind, ttaggg)
colnames(ttaggg) <- as.character(0:8)
rownames(ttaggg) <- gsub("\\.fastq\\.gz", "", files)

ttaggg.perc <- apply(ttaggg[, -1], 1, function(x) x / sum(x))    # fraction of 1+ repeats!!

out <- list(absolute=as.data.frame(ttaggg),
            fraction=as.data.frame(t(ttaggg.perc)))
WriteXLS("out", "./results/ttaggg.xls", row.names=TRUE)

# do plot
df <- melt(ttaggg.perc)
ggplot(df, aes(x=X2, y=value, fill=as.factor(X1))) +
    geom_bar(stat="identity") +
    scale_fill_brewer("repeats", palette="Greens") +
    xlab("") + ylab("fraction") + coord_flip() + ggtitle("U2OS ChIP") +
    theme_bw()

dev.off()
```

## Predict genes network with genemania

A bit convoluted to setup. First get the [command line tools jar](http://pages.genemania.org/command-line-tools/). Then, install the [Cytoscape app](http://apps.cytoscape.org/apps/genemania), open Cytoscape and download the network sources for the desired organism (I didn't find any other way to get them from genemania.org). Once downloaded, put them together with the binary (the jar file).

```R
GM_THREADS <- 4
GM_PATH <- "/home/sayolspu/src/tox_prediction/networks/genemania-cytoscape-plugin-2.18.jar"
GM_DATA_PATH <- "/home/sayolspu/src/tox_prediction/networks/genemania-data/gmdata-2017-07-13-core/"

gmn <- tapply(drugtargets$target_id, drugtargets$drug_id, function(x) {
  # temp file with the query:
  # S. Cerevisiae
  # CDC27 APC11 APC4  XRS2  RAD54 APC2  RAD52 RAD10 MRE11 APC5
  # preferred
  # 150
  # bp
  input   <- tempfile()
  output  <- tempfile()
  targets <- get.nwobj.genes(unique(x))$genesymbol
  writeLines(paste("H. Sapiens",
                   paste(targets, collapse="\t"),
                   "preferred",  # coexp (Co-expr), gi (Genetic int), pi (Physical int)
                   "150", # related-gene-limit
                   "bp",  # Networks are weighted in an attempt to reproduce GO BP co-annotation patterns
                   sep="\n"),
             input)

  # run command
  cmd <- paste("java -Xmx1800M",
               paste("-cp", GM_PATH, "org.genemania.plugin.apps.QueryRunner"),
               paste("--data", GM_DATA_PATH),
               "--out xml",
               "--threads", GM_THREADS,
               "--results", output,
               input)
  dir.create(output)
  try(system(cmd))

  # parse output
  f <- file.path(output, paste0(basename(input), "-results.report.xml"))
  if(file.exists(f)) {
    gm <- xmlToList(xmlParse(f))
    nw <- do.call(rbind, gm[["results"]][["interactions"]])
    nw <- as_adjacency_matrix(graph_from_edgelist(nw[, 1:2, drop=FALSE]), sparse=FALSE)
  } else {
    nw <- matrix(nrow=0, ncol=0)
  }

  nw
})

names(gmn) <- tapply(drugtargets$drug_name, drugtargets$drug_id, unique)
saveRDS(gmn, file="DILI_drugtargets_genemania_networks.RDS")
```

## Annotation of principal isoforms from Biomart

Efforts for curating transcripts available only in GRCH38 are
  * [Gencode basic flag](https://www.ensembl.org/info/genome/genebuild/transcript_quality_tags.html#basic): subset of the GENCODE transcript set, containing only 5' and 3' complete transcripts. 
  * [MANE Select](https://www.ncbi.nlm.nih.gov/refseq/MANE/): default transcript per human gene that is representative of biology, well-supported, expressed and highly-conserved
  * [APPRIS](http://appris-tools.org/) selects a single CDS variant for each gene as the 'PRINCIPAL' isoform based on the range of protein features

```R
appris <- read.delim("http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh37/appris_data.principal.txt", head=FALSE)
appris <- appris[appris$V5 == "PRINCIPAL:1", ]

mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org"))
regions  <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "chromosome_name",
                               "transcript_gencode_basic", "transcript_start", "transcript_end",
                               "strand", "transcription_start_site"),
                  filters="biotype", values="protein_coding", mart=mart)
                  
regions <- regions[regions$ensembl_transcript_id %in% appris$V3, ]  # main isoform
regions <- regions[!is.na(regions$transcript_gencode_basic), ]      # only complete transcripts (redundant with appris)
regions$strand <- ifelse(regions$strand[1] > 0, "+", "-")
regions <- makeGRangesFromDataFrame(regions, keep.extra.columns=TRUE)
```

## Match a pattern in the genome

Option 1: using a short read mapper like bowtie:

```sh
#!/bin/bash
REF=/fsimb/common/genomes/homo_sapiens/ensembl/grch37/canonical/index/bowtie1/genome
MOTIF=GCGATCGC

ml bowtie
bowtie $REF -v 0 -ca $MOTIF | cut -f3-4 | sort -u > asisi.txt
```

Option2: using Bioconductor

```R
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

asisi_motif <- "GCGATCGC"

asisi <- do.call(rbind, mclapply(1:24, function(i) {  # loop over chr 1:22+XY
  m <- matchPattern(asisi_motif, Hsapiens[[i]])
  data.frame(chr=seqnames(Hsapiens)[i], start=start(m), end=end(m))
}))

asisi <- makeGRangesFromDataFrame(asisi)
seqlevelsStyle(asisi) <- "Ensembl"
```

## Reverse complement a DNA sequence

in Perl:

```perl
sub revcomp {
  my ($dna) = @_;
  my $rc = reverse($dna);
  $rc =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
  return $rc;
}
```

in R:

```R
chartr("ACTG", "TGAC", paste(rev(unlist(strsplit(s, ""))), collapse=""))
```

or using Bioconductor for more complex [IUPAC strings](https://www.bioinformatics.org/sms/iupac.html):

```R
library(Biostrings)
x <- DNAString("ACGT-YN-")
reverseComplement(x)
```

## Query the Uniprot Rest API

Retrieving stuff using the [Rest API](https://www.uniprot.org/help/api_queries) provided by Uniprot:

```sh
wget -qO- "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cprotein_existence%2Cid%2Cprotein_name%2Cgene_primary%2Corganism_name%2Clength%2Csequence&format=tsv&query=%28%28organism_id%3A9606%29%29" | gunzip > scripts/iupred.protein_sequences.txt'
```

[Multiple fields](https://www.uniprot.org/help/return_fields) can be retrieved, also from [external references](https://www.uniprot.org/help/return_fields_databases).

See API URL generated in here: https://www.uniprot.org/uniprotkb?query=(organism_id:9606). The interesting thing of using this API call is that length of the results is not limited and not paginated, which complicates retrieving eg. whole human set of proteins.

## Calculate best primers using melting temperature

```R
##############################
##
## Do primers from a multifasta file
##
## 1-read multifasta file
## 2-loop, for every sequence:
##   2.1-take n=MINNUC from 5'. Take n=n+1 until TM(seq_n) > TM
##   2.2-take n=MINNUC from 3'. Take n=n+1 until TM(RC(seq_n)) > TM
##   2.3-output: name;primer_5';TM_5';primer_3';TM_3'
##
## TM=4*(C+G)+2*(A+T)
##
##############################
library("ShortRead")
FASTA <- "feat.fasta"   # input file with gene sequences
MINL  <- 15
TM    <- 58
OUT   <- "doprim.csv"   # output file with primers

## 1-read multifasta file
fasta <- readFasta(FASTA)

# calculate melting temperature
tm <- function(s) {

        nC <- length(gregexpr("C",s)[[1]])
        nG <- length(gregexpr("G",s)[[1]])
        nA <- length(gregexpr("A",s)[[1]])
        nT <- length(gregexpr("T",s)[[1]])

        return(4*(nC+nG)+2*(nA+nT))
}

# calculate reverse complimentary
RC <- function(s) {
        chartr("ACTG","TGAC",paste(rev(unlist(strsplit(s,""))),collapse=""))
}

## 2-loop, for every sequencia:
for(i in 1:length(fasta)) {

        s <- as.character(sread(fasta[i]))      # sequence
        x <- as.character(id(fasta[i]))         # id

        ##   2.1-take n=MINNUC from 5'. Take n=n+1 until TM(seq_n) > TM
        n  <- MINL
        t5 <- 0
        while(t5 < TM) {
                s5 <- substr(s,1,n)
                t5 <- tm(s5)
                n  <- n + 1
        }

        ##   2.2-take n=MINNUC from 3'. Take n=n+1 until TM(RC(seq_n)) > TM
        n  <- MINL
        t3 <- 0
        while(t3 < TM) {
                s3 <- RC(substr(s,nchar(s) - n,nchar(s)))       # idem, but Reverse Complimentary
                t3 <- tm(s3)
                n  <- n + 1
        }

        ##   2.3-output: name;primer_5';TM_5';primer_3';TM_3'
        cat(paste(x,s5,t5,nchar(s5),s3,t3,nchar(s3),s,sep=";"),file=OUT,fill=T,append=T)
}
```

## Broad Institute GSEA
Command line gene set enrichment analysis (GSEAPreranked test) using Broad's software and datasets.

```bash
#!/bin/bash

export PROJECT=/fsimb/groups/imb-bioinfocf/projects/sfb/kraemer/sfb_kraemer_2021_01_dzulko_pr130_RNAseq

[ -d ${PROJECT}/tmp ] || mkdir -p ${PROJECT}/tmp
[ -d ${PROJECT}/bin/gsea ] || mkdir -p ${PROJECT}/bin/gsea
[ -d ${PROJECT}/ref/gsea ] || mkdir -p ${PROJECT}/ref/gsea
[ -d ${PROJECT}/results/gsea ] || mkdir -p ${PROJECT}/results/gsea

# download gsea software
BIN=${PROJECT}/bin/gsea/GSEA_Linux_4.2.1
if [ ! -d ${BIN} ]; then
  wget -qO- https://data.broadinstitute.org/gsea-msigdb/gsea/software/desktop/4.2/GSEA_Linux_4.2.1.zip > ${PROJECT}/bin/gsea/GSEA_Linux_4.2.1.zip
  unzip ${PROJECT}/bin/gsea/GSEA_Linux_4.2.1.zip -d ${PROJECT}/bin/gsea/ && rm ${PROJECT}/bin/gsea/GSEA_Linux_4.2.1.zip
fi

# download msigdb relevant datasets
[ -e ${PROJECT}/ref/gsea/h_all_v7.5_symbols.gmt ] || wget -qO- https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5/h.all.v7.5.symbols.gmt > ${PROJECT}/ref/gsea/h_all_v7.5_symbols.gmt
[ -e ${PROJECT}/ref/gsea/c2_cp_kegg_v7.5_symbols.gmt ] || wget -qO- https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5/c2.cp.kegg.v7.5.symbols.gmt > ${PROJECT}/ref/gsea/c2_cp_kegg_v7.5_symbols.gmt
[ -e ${PROJECT}/ref/gsea/c2_cp_reactome_v7.5_symbols.gmt ] || wget -qO- https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5/c2.cp.reactome.v7.5.symbols.gmt > ${PROJECT}/ref/gsea/c2_cp_reactome_v7.5_symbols.gmt
[ -e ${PROJECT}/ref/gsea/c5_go_bp_v7.5_symbols.gmt ] || wget -qO- https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5/c5.go.bp.v7.5.symbols.gmt > ${PROJECT}/ref/gsea/c5_go_bp_v7.5_symbols.gmt

echo ${PROJECT}/ref/gsea/h_all_v7.5_symbols.gmt > ${PROJECT}/tmp/gene_sets.txt
echo ${PROJECT}/ref/gsea/c2_cp_kegg_v7.5_symbols.gmt >> ${PROJECT}/tmp/gene_sets.txt
echo ${PROJECT}/ref/gsea/c2_cp_reactome_v7.5_symbols.gmt >> ${PROJECT}/tmp/gene_sets.txt
echo ${PROJECT}/ref/gsea/c5_go_bp_v7.5_symbols.gmt >> ${PROJECT}/tmp/gene_sets.txt

# run gsea on different contrasts
process() {
  contrast_name=$1
  in=${PROJECT}/results/DE_DESeq2/${contrast_name}.csv
  out=${PROJECT}/results/gsea/${contrast_name}

  # read DE results from that contrast, and extract gene_symbol+log2FC
  tail -n +2 $in | cut -f2,8 -d, | grep -v 'NA$' | tr , "\t" | sed 's/"//g' > ${PROJECT}/tmp/${contrast_name}.rnk
  # run Broad's GSEA software
  ${BIN}/gsea-cli.sh GSEAPreranked -rnk ${PROJECT}/tmp/${contrast_name}.rnk -gmx_list ${PROJECT}/tmp/gene_sets.txt -out $out -rpt_label ${contrast_name}
  # tidy up results folder
  mv ${out}/*/* $out
  mkdir -p ${out}/img
  mv ${out}/*.png ${out}/img
  mkdir -p ${out}/tables
  mv ${out}/*.tsv ${out}/tables
  sed -i "s/src='/src='img\//g" ${out}/*.html
}
export -f process

tail -n +2 ${PROJECT}/contrasts.txt | while read contrast_name contrast mmatrix; do
  process $contrast_name &
done

wait $(jobs -p)
```

## Problematic and blacklisted regions from UCSC

Annotate problematic regions using info from UCSC tables (instead of using ENCODE's blacklisted regions from [Boyle-Lab](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)).
The advantage of UCSC is that they have sub-tracks with other sources of problems, namely:

* The UCSC Unusual Regions subtrack contains annotations collected at UCSC, put together from other tracks, our experiences and support email list requests over the years. For example, it contains the most well-known gene clusters (IGH, IGL, PAR1/2, TCRA, TCRB, etc) and annotations for the GRC fixed sequences, alternate haplotypes, unplaced contigs, pseudo-autosomal regions, and mitochondria. These loci can yield alignments with low-quality mapping scores and discordant read pairs, especially for short-read sequencing data. This data set was manually curated, based on the Genome Browser's assembly description, the FAQs about assembly, and the NCBI RefSeq "other" annotations track data.
* The ENCODE Blacklist subtrack contains a comprehensive set of regions which are troublesome for high-throughput Next-Generation Sequencing (NGS) aligners. These regions tend to have a very high ratio of multi-mapping to unique mapping reads and high variance in mappability due to repetitive elements such as satellite, centromeric and telomeric repeats.
*The GRC Exclusions subtrack contains a set of regions that have been flagged by the GRC to contain false duplications or contamination sequences. The GRC has now removed these sequences from the files that it uses to generate the reference assembly, however, removing the sequences from the GRCh38/hg38 assembly would trigger the next major release of the human assembly. In order to help users recognize these regions and avoid them in their analyses, the GRC have produced a masking file to be used as a companion to GRCh38, and the BED file is available from the GenBank FTP site.

Additionally one can add other genomic features, like centromeres/telomeres from the `gap` table.

Here is the R code used to annotate a GRanges object with such regions:

```R
library(rtracklayer)

# retrieve coordinates of problematic regions
session <- browserSession("UCSC")
#track.names <- trackNames(ucscTableQuery(session))
#table.names <- ucscTables("hg38", track.names[grepl("problematic", track.names)])
#table.names <- ucscTables("hg38", track.names[grepl("gap", track.names)])
ucsc_unusual_regions <- keepStandardChromosomes(track(ucscTableQuery(session, table="comments")), pruning.mode="coarse")
encBlacklist         <- keepStandardChromosomes(track(ucscTableQuery(session, table="encBlacklist")), pruning.mode="coarse")
grcExclusions        <- keepStandardChromosomes(track(ucscTableQuery(session, table="grcExclusions")), pruning.mode="coarse")
gap                  <- keepStandardChromosomes(track(ucscTableQuery(session, table="gap")), pruning.mode="coarse")

# annotate with the following preference: ucsc_unusual_regions, encBlacklist, grcExclusions, gap (centro/telo)
guide_xls <- "GUIDEseq_Lazzarotto2020NBT.xlsx"
guide <- makeGRangesFromDataFrame(as.data.frame(read_xlsx(guide_xls)), keep=TRUE)
guide <- unstrand(guide)

guide$region_annotation <- ""
i <- findOverlaps(guide, ucsc_unusual_regions)
guide$region_annotation[queryHits(i)] <- ucsc_unusual_regions$name[subjectHits(i)]
i <- findOverlaps(guide, encBlacklist)
guide$region_annotation[queryHits(i)] <- encBlacklist$name[subjectHits(i)]
i <- findOverlaps(guide, grcExclusions)
guide$region_annotation[queryHits(i)] <- grcExclusions$name[subjectHits(i)]
i <- findOverlaps(guide, gap)
guide$region_annotation[queryHits(i)] <- gap$type[subjectHits(i)]
```

## Signal summary around regions of interest
The usual metagene plots with histone marks:
```R
metaplot <- function(regions, signal, ext=5000, bins=201, title="") {
  # define regions
  regions <- keepSeqlevels(regions, names(signal), pruning.mode="coarse")
  middle <- mid(regions)
  start(regions) <- middle - ext
  end(regions)   <- middle + ext
  regions.binned <- unlist(tile(regions, n=bins))
  regions.binned$region <- rep(seq_along(regions), each=bins)

  # get binned signal average
  x <- binnedAverage(regions.binned, signal, "score")
  x <- matrix(x$score, ncol=bins, byrow=TRUE)
  x <- lapply(split(seq_len(nrow(x)), regions$scission_profile), function(i) {
    colSums(x[i, ], na.rm=TRUE) / length(i)   # aggregate signal from multiple regions, normalized to the number of regions
  })

  # add some random regions to show the bg
  regions.random <- shift(regions, shift=runif(length(regions), -100000, 100000)) # randomly shift region +/- 1M
  regions.random.binned <- unlist(tile(regions.random, n=bins))                   # tile into bins
  regions.random.binned$region <- rep(seq_along(regions.random), each=bins)
  x.random <- binnedAverage(regions.random.binned, signal, "score") # get binned signal average
  x.random <- matrix(x.random$score, ncol=bins, byrow=TRUE)
  x$random <- colSums(x.random, na.rm=T) / length(regions.random)   # aggregate signal + normalize to the number of regions

  # plot something
  df    <- reshape2::melt(x)
  df$x  <- rep(1:bins, length(x))
  df$L1 <- factor(df$L1, levels=c("blunt", "middle", "staggered", "random"))
  pal   <- c(palette()[1:3], "grey30")

  ggplot(df, aes(x=x, y=value, color=L1)) +
    geom_smooth(span=1/20) +
#    geom_line() +
    geom_vline(xintercept=(1+bins)/2, lty=2) +
    labs(title=title, y="arbitrary units", x="") +
    scale_color_manual("", values=pal) +
    scale_x_continuous(breaks=c(0, (1+bins)/2, bins), labels=c(-ext, "DSB", ext)) +
    theme(legend.position="bottom")
}
```
This is just a simplified version with 1 single coordinate extended up/down stream. A more elaborated version with a gene body and its up/down stream regions requires uneven splitting of the 3 regions (since the gene body will be of different size and the up/down streams of fixed size).

### Other methods to aggregate signal in a bin
A binnedAverage() like function, more general, in a way that it could compute more than just the mean:

```R
binnedView <- function(bins, numvar, varname, na.rm=FALSE, fun=IRanges::viewMeans) {
    if (!is(bins, "GRanges"))
        stop("'x' must be a GRanges object")
    if (!is(numvar, "RleList"))
        stop("'numvar' must be an RleList object")
    if (!identical(seqlevels(bins), names(numvar)))
        stop("'seqlevels(bin)' and 'names(numvar)' must be identical")
    fun2 <- function(v, na.rm = FALSE) {
        if (!isTRUEorFALSE(na.rm))
            stop("'na.rm' must be TRUE or FALSE")
        result <- fun(v, na.rm = na.rm)
        w0 <- width(v)
        v1 <- trim(v)
        w1 <- width(v1)
        if (na.rm) {
            na_count <- sum(is.na(v1))
            w0 <- w0 - na_count
            w1 <- w1 - na_count
        }
        result <- result * w1/w0
        result[w0 != 0L & w1 == 0L] <- 0
        result
    }
    bins_per_chrom <- split(ranges(bins), seqnames(bins))
    result_list <- lapply(names(numvar), function(seqname) {
        v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
        fun2(v, na.rm = na.rm)
    })
    new_mcol <- unsplit(result_list, as.factor(seqnames(bins)))
    mcols(bins)[[varname]] <- new_mcol
    bins
}
```

See the [issue](https://github.com/Bioconductor/GenomicRanges/issues/77) filled in GenomicRange's Github repo.

