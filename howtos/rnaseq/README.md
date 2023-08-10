# A collection of R snippets

## Table of Contents

* [RNA-seq workflows](#rna-seq-workflows)
   * [QC](#qc)
      * [dupRadar](#dupradar)
      * [Saturation curves](#saturation-curves)
      * [Dispersion vs. abundance](#dispersion-vs.-abundance)
      * [DGE chromosome overrepresentation](#dge-chromosome-overrepresentation)
   * [Calling DE genes with edgeR](#calling-de-genes-with-edger)
   * [Calling DE genes with DESeq2](#calling-de-genes-with-deseq2)
   * [Calling DE genes with limma voom](#calling-de-genes-with-limma-voom)
   * [Differential Exon Usage](#differential-exon-usage)
      * [DEXSeq](#dexseq)
   * [Differential splicing analysis](#differential-splicing-analysis)
      * [rMATS](#rmats)
   * [MISO for isoform quantification](#miso-for-isoform-quantification)
   * [Plots](#plots-1)
      * [Volcano](#volcano)
   * [Other](#other)
      * [RPK, RPKM, TPM](#rpk-rpkm-tpm)
      * [Generate a STAR index](#generate-a-star-index)
      * [DGE with limma if raw counts are not available](#dge-with-limma-if-raw-counts-are-not-available)

## RNA-seq workflows
### QC
To be done.
#### dupRadar
#### Saturation curves
#### Dispersion vs. abundance
#### DGE chromosome overrepresentation

### Calling DE genes with edgeR
This is just part of a scRNA-seq experiment. ToDo: create the edgeR dataset.

```R
f <- function(group, conts) {
  sce$group <- group

  # prepare edgeR model matrix and contrasts
  de.design <- model.matrix(~ 0 + sce$group + sce$bp, colData(sce))
  colnames(de.design) <- gsub("^sce\\$(group|bp)", "", colnames(de.design))
  de.conts  <- makeContrasts(contrasts=conts, levels=de.design)

  # do the DGE analysis
  y <- convertTo(sce, type="edgeR")
  y <- estimateDisp(y, de.design)

  # QL glm fit
  y.qlfit <- glmQLFit(y, de.design)
  qlf     <- glmQLFTest(y.qlfit, contrast=de.conts)
  res.qlf <- topTags(qlf, n=nrow(qlf), p.value=0.01)
  res.qlf$table
}

res.qlf <- mcMap(
  f,
  group=list(lgr5.1_vs_lgr5.0        =ifelse(sce$Lgr5, "Lgr5.1", "Lgr5.0"),
             lgr5.1_vs_lgr5.0_krt15.0=ifelse(sce$Lgr5, "Lgr5.1",
                                      ifelse(sce$Krt15, "Lgr5.0_Krt15.1", "Lgr5.0_Krt15.0"))),
  conts=list(lgr5.1_vs_lgr5.0        ="Lgr5.1-Lgr5.0",
             lgr5.1_vs_lgr5.0_krt15.0="Lgr5.1-Lgr5.0_Krt15.0"),
  mc.cores=CORES
)
```

### Calling DE genes with DESeq2
To be done.

### Calling DE genes with limma voom
```R
library(limma)

dge <- DGEList(counts(sce))  # sce is a single cell experiment
dge <- calcNormFactors(dge)

design <- model.matrix(~ (cell_type == "Unknown"), colData(sce))

v <- voom(dge, design, plot=TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)
res <- decideTests(fit)

df_limma <- data_frame(log2foldchange = fit$coefficients[,2], 
                       pval = fit$p.value[,2],
                       qval = p.adjust(fit$p.value[,2], method = 'BH'),
                       ensembl_id = rownames(sce),
                       gene_symbol = rowData(sce)$symbol)
```

### Differential exon usage

#### DEXSeq

Do differential exon usage with ''DEXSeq'' and ''Bioconductor''

```R
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
library(DEXSeq)
library(BiocParallel)

PROJECT <- "/project/
GTF <- "/annotation/Homo_sapiens.GRCh38.84.gtf"
FC  <- log2(1.5)    # expect 50% more expression
FDR <- .01

##
## count reads on exons
##
exonicParts   <- exonicParts(makeTxDbFromGFF(GTF))
#intronicParts <- intronicParts(makeTxDbFromGFF(GTF))   # alternatively, one could look at intron retention
bams <- BamFileList(list.files(paste0(PROJECT, "/mapped"), pattern="_read\\.bam$", full=TRUE),
                    index=character(),              # the BAM index file path
                    asMates=TRUE,                   # records should be paired as mates
                    obeyQname=TRUE)                 # BAM file is sorted by ‘qname’

counts <- summarizeOverlaps(exonicParts,
                            bams,
                            mode="Union",           # default htseq union mode
                            singleEnd=FALSE,        # data is paired end
                            inter.feature=FALSE,    # don't discard reads spannings multiple exons
                            fragments=TRUE,         # count also singletons
                            ignore.strand=TRUE,     # it's strand specific, but still ignore the strand
                            BPPARAM=MulticoreParam(workers=6))

##
## DEXSeq
##
colData(counts)$condition <- c(rep("Id2", 3), rep("GFP", 3))
dds <- DEXSeqDataSetFromSE(counts, design= ~ sample + exon + condition:exon )

# normalize, estimate dispersion, test for differential expression and estimate fold changes
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, BPPARAM=MulticoreParam(workers=16))
dds <- testForDEU(dds, BPPARAM=MulticoreParam(workers=16))
dds <- estimateExonFoldChanges(dds, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=16))

# extract significant genes
res <- DEXSeqResults(dds)

res$log2fold_Id2_GFP <- ifelse(is.na(res$log2fold_Id2_GFP), 0, res$log2fold_Id2_GFP)
res$padj <- ifelse(is.na(res$padj), 1, res$padj)
geneids <- unique(res$groupID[abs(res$log2fold_Id2_GFP) > FC & res$padj < FDR])

res2 <- res[res$groupID %in% geneids, ]
x <- do.call(rbind, by(res2, res2$groupID, function(x) data.frame(padj=min(x$padj), log2fold_Id2_GFP=max(x$log2fold_Id2_GFP))))

DEXSeqHTML(res2, path=paste0(PROJECT, "/results/DEXSeq"), FDR=FDR, BPPARAM=MulticoreParam(workers=16), extraCols=x,
           mart=useMart("ensembl",dataset="hsapiens_gene_ensembl"), filter="ensembl_gene_id", attributes="external_gene_name")

save.image(file=paste0(PROJECT, "/results/DEXSeq.RData"))
```

'''Notes:'''
Running on PE inversely stranded protocols, look at [https://support.bioconductor.org/p/65844/ this bioconductor thread].

Some parts of the code need to be adjusted. First, we'll need to define our own invertStrand function:

```R
invertStrand <- function(galp)
{
    ## Using a non-exported helper function and direct slot access is
    ## bad practice and is strongly discouraged. Sorry for doing this here!
    invertRleStrand <- GenomicAlignments:::invertRleStrand
    galp@first <- invertRleStrand(galp@first)
    galp@last <- invertRleStrand(galp@last)
    galp
}
```

Which will then be called from the ''summarizeOverlaps()'' call in ''DEXSeq'', as a preprocessing step:

```R
counts <- summarizeOverlaps(exonicParts,
                            bams,
                            mode="Union",           # default htseq union mode
                            singleEnd=FALSE,        # data is paired end
                            inter.feature=FALSE,    # don't discard reads spannings multiple exons
                            fragments=TRUE,         # count also singletons
                            ignore.strand=FALSE,
                            preprocess.reads=invertStrand, # as suggested in: https://support.bioconductor.org/p/65844/
                            BPPARAM=MulticoreParam(workers=6))
```
### Differential splicing analysis

#### rMATS

[Quick example](https://github.com/Xinglab/rmats-turbo/#examples), starting from FastQ files (includes mapping with STAR with whatever parameters rMATS likes):

```bash
#!/bin/bash
PROJECT="/fsimb/groups/imb-bioinfocf/projects/roukos/imb_roukos_2020_13_sant_rnaseq_aquarius"
GROUP1="${PROJECT}/scripts/rmats_group_aux.txt"
GROUP2="${PROJECT}/scripts/rmats_group_dmso.txt"
STAR_REF="/fsimb/common/genomes/homo_sapiens/gencode/release-26_GRCh38.p10/full/index/star/2.7.3a"
GENESGTF="/fsimb/common/genomes/homo_sapiens/gencode/release-26_GRCh38.p10/full/annotation/gencode.v26.primary_assembly.annotation.gtf"
READLENGTH=42
THREADS=4
TMP="${PROJECT}/tmp"
[[ -n $SLURM_JOB_ID ]] && TMP="/jobdir/$SLURM_JOB_ID"
mkdir -p $TMP

ls ${PROJECT}/rawdata/*.fastq.gz | parallel -j $THREADS "x=\$(basename {}); zcat {} > ${TMP}/\${x%.gz}"
cat $GROUP1 | sed "s~\%path\%~${TMP}~g" > ${TMP}/$(basename $GROUP1)
cat $GROUP2 | sed "s~\%path\%~${TMP}~g" > ${TMP}/$(basename $GROUP2)

ml rmats/4.1.0
python $(which rmats.py) \
  --s1 ${TMP}/$(basename $GROUP1) \
  --s2 ${TMP}/$(basename $GROUP2) \
  --gtf $GENESGTF \
  --bi $STAR_REF \
  -t paired \
  --readLength $READLENGTH \
  --nthread $THREADS \
  --tstat $THREADS \
  --od ${PROJECT}/results/rmats \
  --tmp $TMP
```

The targets files containing the fastq samples to be processed should be included in 2 files:

```bash
$ cat rmats_group_aux.txt 
%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_02_AUX_1.R1.fastq:%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_02_AUX_1.R2.fastq,%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_04_AUX_2.R1.fastq:%path%/imb_roukos_2020_13_sant_r naseq_aquarius_04_AUX_2.R2.fastq,%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_06_AUX_3.R1.fastq:%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_06_AUX_3.R2.fastq,%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_08_AUX_4.R1.fastq :%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_08_AUX_4.R2.fastq,%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_10_AUX_5.R1.fastq:%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_10_AUX_5.R2.fastq
$ cat rmats_group_dmso.txt 
%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_01_DMSO_1.R1.fastq:%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_01_DMSO_1.R2.fastq,%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_03_DMSO_2.R1.fastq:%path%/imb_roukos_2020_13_san t_rnaseq_aquarius_03_DMSO_2.R2.fastq,%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_05_DMSO_3.R1.fastq:%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_05_DMSO_3.R2.fastq,%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_07_DMSO_4.R 1.fastq:%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_07_DMSO_4.R2.fastq,%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_09_DMSO_5.R1.fastq:%path%/imb_roukos_2020_13_sant_rnaseq_aquarius_09_DMSO_5.R2.fastq
```

Notice that we use the pattern %path%, that we will later substitute within the script with the $TMP folder (depending whether the job run interactively or in the cluster)

### MISO for isoform quantification

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

### Plots
#### Volcano
```R
df_limma$is_significant <- df_limma$qval < 0.01

ggplot(df_limma, aes(x = log2foldchange, y = -log10(qval), color = is_significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "darkred"),
                     name = "Significantly differentially expressed") +
  geom_text_repel(data = subset(df_limma, qval < 0.00001), aes(label = gene_symbol), color = 'black', size = 2) +
  labs(x = expression(log[2]~"(fold change)"),
       y = expression(-log[10]~"(q-value)"),
       subtitle = "Differential expression of Unknown vs. Known cell types") +
  theme(legend.position = "bottom",
        legend.box.background = element_rect(linetype = 1, size = 1)) +
  annotate("label", x = 2, y = 0, label = "Upregulated in Unknown") +
  annotate("label", x = -1.2, y = 0, label = "Upregulated in epithelial (sorted) cells")
```

### Other
#### RPK, RPKM, TPM

Zhao, Ye and Stanton (RNA 2020 Aug;26(8):903-909. doi: 10.1261/rna.074922.120) define RPKM and TPM as:

RPKM was initially introduced to facilitate transparent comparison of transcript levels both within and between samples, as it re-scales gene counts to correct for differences in both library sizes and gene length (Mortazavi et al. 2008)

The intended meaning of RPKM is a measure of relative RNA molar concentration (rmc) of a transcript in a sample. If a measure of RNA abundance is proportional to rmc, then their average over genes within a sample should be a constant, namely the inverse of the number of transcripts mapped. Unfortunately, RPKM does not respect this invariance property and thus cannot be an accurate measure of rmc (Wagner et al. 2012).

TPM is unit-less, and it additionally fulfils the invariant average criterion. For a given RNA sample, if you were to sequence one million full length transcripts, a TPM value represents the number of transcripts you would have seen for a given gene or isoform. The average TPM is equal to 10 ^6 (1 million) divided by the number of annotated transcripts in a given annotation, and thus is a constant. TPM is a better unit for RNA abundance since it respects the invariance property and is proportional to the average rmc, and thus adopted by the latest computational algorithms for transcript quantification. 

```R
rpk  <- reads * 1e3 / width(gene)
rpkm <- reads * 1e3 / width(gene) * 1e6 / sum(reads)
tpm  <- 1e6 * (reads / width(gene)) / sum(reads / width(gene))
```

The relationship between TPM and RPKM can be expressed as:

```R
tpm  <- 1e6 * rpkm / sum(rpkm)
```

#### Generate a STAR index

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

#### DGE with limma if raw counts are not available

If raw counts are not available,  proceed as suggested in this post (and in limma's vignette)
  * https://support.bioconductor.org/p/92303/
```
If you really were stuck with nothing but CPM values, then the best approach would be to transform to log2 values:
y <- log2(CPM + 0.1)
and then analyse in the limma package as if it was microarray data, using a limma-trend type analysis.
```

* https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
```
If the sequencing depth is reasonably consistent across the RNA samples, then the simplest and most
robust approach to differential exis to use limma-trend. This approach will usually work well if the
ratio of the largest library size to the smallest is not more than about 3-fold.
In the limma-trend approach, the counts are converted to logCPM values using edgeR’s cpm function
```

