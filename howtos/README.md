## Table of Contents

* [Convert BAM to BigWig](#convert-bam-to-bigwig)
* [Deduplicate UMIs](#deduplicate-umis)
* [Downsample](#downsample)
   * [Fastq](#fastq)
   * [BAM](#bam)
* [Extend reads for ChIP/MBD](#extend-reads-for-chipmbd)
* [Generate a STAR index](#generate-a-star-index)
* [Merge BAM files](#merge-bam-files)
* [Merge BigWig tracks](#merge-bigwig-tracks)
* [MISO for isoform quantification](#miso-for-isoform-quantification)
* [Motif discovery](#motif-discovery)
   * [MEME](#meme)
   * [rGadem](#rgadem)
* [Quick and dirty, call peaks with MACS2 and SICER](#quick-and-dirty-call-peaks-with-macs2-and-sicer)
   * [MACS2](#macs2)
   * [SICER](#sicer)
* [Retrieve chromosome sizes](#retrieve-chromosome-sizes)
* [Submit jobs to random machines in the cluster](#submit-jobs-to-random-machines-in-the-cluster)

## Convert BAM to BigWig
This script will loop over the BAM files in a directori, and convert them to BigWig files or visualization in Genome Browsers.

To queue it into LSF, the for loop has to be moved into a master script that will queue with bsub this script passing the BAM file to be processed as a parameter.

**Requirements**

* bedtools
* Picard tools
* bedGraphToBigWig from the UCSC Genome Browser tools
* The cromosome sizes, which can be retrieved from the UCSC Genome Browser tools with the fetchChromSizes. Alternatively, samtools idxstats sample.bam | cut -f1-2 > ./chr.sizes

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

## Extend reads for ChIP/MBD

This script will loop over the BAM files in a directori, and extend 3'end of the reads for $1 bp in order to match the average library's insert size and improve the peak calling.

**Requirements**

* bedtools
* samtools
* The cromosome sizes, which can be retrieved from the UCSC Genome Browser tools with the fetchChromSizes

**Source**
```bash
#!/bin/bash
SIZE=$1
BASEDIR=./mapped
EXEC=/opt
REF=./mm9.chrom.sizes

for f in ${BASEDIR}/*.bam
do
  F=$(basename $f)

  # 1: bam2bed
  # 2: extend 51bp reads to average fragment size 180bp (slopBED)
  # 3: bed2bam
  # 4: samtools sort bam
  # 5: index bam
  echo "${EXEC}/BEDTools/latest/bin/bedtools bamtobed -split -i ${BASEDIR}/${F} | \
        ${EXEC}/BEDTools/latest/bin/bedtools slop -g ${REF} -l 0 -r ${SIZE} -s | \
        ${EXEC}/BEDTools/latest/bin/bedtools bedtobam -ubam -g ${REF} | \
        ${EXEC}/samtools/latest/samtools sort - ${F%.bam}_ext && \
        ${EXEC}/samtools/latest/samtools index ${F%.bam}_ext.bam" | \
  bsub -J ${F%.bam} -o ${F%.bam}_ext.log -cwd $(pwd) -W 2:00 -n 1 -q "testing"
done
```

## Generate a STAR index

This is the STAR command to generate an index file from a fasta reference:

```bash
    STAR --runMode genomeGenerate  \
         --genomeDir ./  \
         --genomeFastaFiles /igenomes_reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
         --runThreadN 12 \
         --sjdbGTFfile /igenomes_reference/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf \
         --sjdbOverhang 50
```

Newer versions of STAR (since 2.4.x) can produce indexes without the read length parameter while generating the index. Though it has to be supplied when running STAR alignment.

```bash
    STAR --runMode genomeGenerate     \
         --genomeDir ${INDEXFOLDER}   \
         --genomeFastaFiles ${GENOME} \
         --runThreadN ${THREADS}      \
         --outTmpDir ${TMPDIR}/indexing/
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
* The cromosome sizes, which can be retrieved from the UCSC Genome Browser tools with the fetchChromSizes. Alternatively, samtools idxstats sample.bam | cut -f1-2 > ./chr.sizes

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

## MISO for isoform quantification

[http://miso.readthedocs.org/en/fastmiso/ MISO] can be used to estimate isoform quantification. It's reported to do a bad job, use it at your own risk!

**Requirements**

* A GFF3 file with the gene annotation. There is extensive documentation in the MISO webpage on how to create your own from UCSC or Ensembl. I usually download directly the GFF version from [http://www.ensembl.org/info/data/ftp/index.html Ensembl] . Older assemblies from Ensembl do not provide GFF gene models; in this case, I download the GTF file and convert it to GFF with the [http://www.sequenceontology.org/software/GAL.html GAL] library.
* Sorted+indexed BAM files

**index the GFF file**

From the MISO suite, run the index_gff tool on the GFF3 file:

```bash
ssayolsp@annotation$ index_gff --index Drosophila_melanogaster.BDGP6.83.gff3 indexed/
```

**Run MISO**
```bash
#!/bin/bash
PROJECT=/project/xxx
ANNOTATION=${PROJECT}/annotation/indexed.BDGP6.83.gff3
CORES=8
READLEN=36

for f in ${PROJECT}/mapped/*.sorted.bam; do
    F=$(basename $f)
    F=${F%_accepted_hits.sorted.bam}
    OUT=${PROJECT}/results/isoform-quantification/${F}

    L1="source /opt/miso/latest/env.sh" 
    L2="miso --run $ANNOTATION $f --output-dir $OUT --read-len $READLEN -p $CORES"
    L3="miso_pack --pack $OUT"
    L4="summarize_miso --summarize-samples $OUT $OUT"
    echo "$L1 && $L2 && $L3 && $L4" | bsub -q short -W 5:00 -R "span[ptile=${CORES}]" -n ${CORES} -app Reserve1G -J miso_${F} -o miso_${F}.out -e miso_${F}.err
done
```

**Parse results**
```perl
#!/bin/perl
# output a 3 column file with gene-transcript-counts columns
# use: $ cat file.miso_summary | perl this_script.pl > file.parsed.miso
while(<>) {

    chomp;

    # split and get gene, isoforms and counts per isoform
    my @l = split /\t/;
    my $g = $l[0];
    my $t = $l[4];
    my $q = $l[6];

    # split isoforms and counts, and get a unique name for the isoform
    my @t = split /,/,$t;
    my @q = split /,/,$q;

    # print registers
    for (my $i=0; $i < scalar @t; $i++) {
        my @ct = $t[$i] =~ m/(FBtr\d+)/;
        my @cq = $q[$i] =~ m/.+:(\d+)/;
        print $g, "\t", $ct[0], "\t", $cq[0], "\n";
    }
}
```

**Compress MISO output folders**

After summarization, save space by compressing the output folders:

```bash
ssayolsp@results$ miso_zip --compress mydata.misozip miso_output/
```

The process can be reverted:

```bash
ssayolsp@results$ miso_zip --uncompress mydata.misozip uncompressed/
```

## Motif discovery

### MEME

```bash
#!/bin/bash
# run meme-chip with the same arguments as calling it from the http://meme-suite.org/
PROJECT=/projects/xxx
OUTDIR=${PROJECT}/GSE51142/results/meme
INDIR=${PROJECT}/GSE51142/results/macs2
BEDTOOLS=/opt/BEDTools/latest/bin/bedtools
MEME=/opt/meme/latest/bin/meme-chip
REF=/igenomes/homo_sapiens/ensembl/grch38/canonical/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
MEMEHUMANDB=/opt/meme/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme
MEMEJASPARDB=/opt/meme/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme

for f in ${INDIR}/*.narrowPeak; do
    if [ ! -e $OUTDIR ]; then
        mkdir $OUTDIR
    fi
    if [ ! -e ${f}.fasta ]; then
        echo "Generating ${f}.fasta..."
        $BEDTOOLS getfasta -fi $REF -bed $f -fo ${f}.fasta
    fi
    echo "$MEME -oc $OUTDIR -time 300 -order 1 -db $MEMEJASPARDB -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 3 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 ${f}.fasta" | bsub -n1 -W5:00 -app Reserve10G -J meme -o meme.out -e meme.err
done
```

### rGadem

```R
library(rGADEM)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

PROJECT <- "/projects/xxx"
TOP <- 50   # take only top peaks (values >50 usually make GADEM cracsh with segfault)

# read in the peaks and call the motif discovery tool
bed <- read.table(paste0(PROJECT, "/GSE51142/results/macs2/zbtb10.vs.igg_peaks.xls"), comment="#", sep="\t", head=T)
sel <- rev(order(bed$X.log10.qvalue.))[1:TOP]
bed$chr <- paste0("chr", bed$chr)
peaks <- with(bed[sel,], RangedData(IRanges(start, end), space=chr))
gadem <- GADEM(peaks, verbose=1, genome=Hsapiens)
```

## Quick and dirty, call peaks with MACS2 and SICER

### MACS2

**targets file**

Create a tab-separated targets file with the IP-input pairs:

```bash
IP                             IP2                            IPname          INPUT                   INPUTname     group
Pol2_IP_rrp6_1.bam       Pol2_IP_rrp6_2.bam       Pol2_rrp6       Input_DNA.bam     Input_DNA     Pol2_rrp6
Pol2_IP_rtt109_1.bam     Pol2_IP_rtt109_2.bam     Pol2_rtt109     Input_DNA.bam     Input_DNA     Pol2_rtt109
Pol2_IP_WT_1.bam         Pol2_IP_WT_2.bam         Pol2_WT         Input_DNA.bam     Input_DNA     Pol2_WT
Pol2_IP_rrp6_1-ss.bam    Pol2_IP_rrp6_2-ss.bam    Pol2_rrp6-ss    Input_DNA-ss.bam  Input_DNA-ss  Pol2_rrp6-ss
Pol2_IP_rtt109_1-ss.bam  Pol2_IP_rtt109_2-ss.bam  Pol2_rtt109-ss  Input_DNA-ss.bam  Input_DNA-ss  Pol2_rtt109-ss
Pol2_IP_WT_1-ss.bam      Pol2_IP_WT_2-ss.bam      Pol2_WT-ss      Input_DNA-ss.bam  Input_DNA-ss  Pol2_WT-ss
```

**MACS2**

Remember to change the thresholds (-m), genome size (-g), cutoff (-q) and add --broad when calling broad peaks:

```bash
#!/bin/bash
MACS2=/opt/macs2/latest
PROJECT=/projects/xxx
TARGETS=${PROJECT}/scripts/targets.txt
tail -n +2 $TARGETS | while read -r TARGET; do
    IP=$(       echo $TARGET | cut -f1 -d" ")
    IP2=$(      echo $TARGET | cut -f2 -d" ")
    IPname=$(   echo $TARGET | cut -f3 -d" ")
    INPUT=$(    echo $TARGET | cut -f4 -d" ")
    INPUTname=$(echo $TARGET | cut -f5 -d" ")
    if [ ! -e "${PROJECT}/results/macs2" ]; then
        mkdir -p "${PROJECT}/results/macs2"
    fi
    L1="source ${MACS2}/env.sh"
    L2="${MACS2}/bin/macs2 callpeak -t ${PROJECT}/mapped/$IP ${PROJECT}/mapped/$IP2 -c ${PROJECT}/mapped/$INPUT -n $IPname -m 5 50 -g mm -q 0.05"
    L3="mv $IPname* ${PROJECT}/results/macs2"
    echo "$L1 && $L2 && $L3" | bsub -n1 -W1:00 -app Reserve2G -J $IPname -o ${IPname}.out -e ${IPname}.err
done
```

### SICER
**targets file**

Create a tab-separated targets file with the IP-input pairs:

```bash
IP      INPUT   name
AE_1_k27m3.bam Input_AE_1_k27m3.bam   AE_1_k27m3
AE_2_k27m3.bam Input_AE_1_k27m3.bam   AE_2_k27m3
E125_1_k27m3.bam       Input_E125_1_k27m3.bam E125_1_k27m3
E125_2_k27m3.bam       Input_E125_1_k27m3.bam E125_2_k27m3
E145_1_k27m3.bam       Input_E145_1_k27m3.bam E145_1_k27m3
E145_2_k27m3.bam       Input_E145_1_k27m3.bam E145_2_k27m3
ISCS_1_k27m3.bam       Input_ISCS_1_k27m3.bam ISCS_1_k27m3
ISCS_2_k27m3.bam       Input_ISCS_1_k27m3.bam ISCS_2_k27m3
```

**SICER**

SICER is potentially good to call broad peaks, like H3K27me3 histone marks. Remember to change SICER's parms:

```bash
#!/bin/bash
BEDTOOLS=/opt/BEDTools/2.25.0/bin
SICER=/opt/sicer/1.1
PROJECT=/projects/xxx
TARGETS=${PROJECT}/peakcalling_SICER/targets.txt

# SICER parms
Species=mm9
redundancy_threshold=1
window_size=200
fragment_size=150
effective_genome_fraction=.74
gap_size=1000
FDR=.01

tail -n +2 $TARGETS | while read -r TARGET; do
        IP=$(   echo $TARGET | cut -f1 -d" ")
        INPUT=$(echo $TARGET | cut -f2 -d" ")
        name=$( echo $TARGET | cut -f3 -d" ")
        WD=${PROJECT}/peakcalling_SICER/${name}
        if [ ! -e "${WD}" ]; then
                mkdir "${WD}"
        fi
        L1="${BEDTOOLS}/bedtools bamtobed -split -i ${PROJECT}/mapping/mapping_results/${IP} > ${WD}/${IP%.bam}.bed"
        L2="${BEDTOOLS}/bedtools bamtobed -split -i ${PROJECT}/mapping/mapping_results/${INPUT} > ${WD}/${INPUT%.bam}.bed"
        L3="source ${SICER}/env.sh && ${SICER}/SICER.sh ${WD} ${IP%.bam}.bed ${INPUT%.bam}.bed ${WD} $Species $redundancy_threshold $window_size $fragment_size $effective_genome_fraction $gap_size $FDR"
        echo "$L1 && $L2 && $L3" | bsub -n1 -W5:00 -app Reserve2G -J $name -o ${WD}/${name}.out -e ${WD}/${name}.err
done
```

## Retrieve chromosome sizes

Use the UCSC toolkit from the deps tree.

The chromosome sizes are often used to generate tracks from a bam file.

**Requirements**

* ucsc tools
* or, alternatively, a mysql client

**Source**

```bash
$ fetchChromSizes hg18 > hg18.chrom.sizes
```

Or:

```bash
$ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  > hg19.genome
```

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

