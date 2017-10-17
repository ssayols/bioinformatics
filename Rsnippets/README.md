# A collection of R snippets

## Table of Contents

* [Array analysis](#array-analysis)
   * [Infinium HumanMethylation450 BeadChip](#infinium-humanmethylation450-beadchip)
      * [Read data from raw IDAT files](#read-data-from-raw-idat-files)
      * [Normalization](#normalization)
      * [Quality control](#quality-control)
      * [Differential methylation using limma](#differential-methylation-using-limma)
   * [Affymetrix® Genome-Wide Human SNP Array 6.0 (and probably others)](#affymetrix-genome-wide-human-snp-array-60-and-probably-others)
   * [Affymetrix expression arrays](#affymetrix-expression-arrays)
      * [Read data from raw CEL files](#read-data-from-raw-cel-files)
      * [Normalization](#normalization-1)
      * [Quality control](#quality-control-1)
      * [Differential expression using limma](#differential-expression-using-limma)
* [ChIP-seq workflows](#chip-seq-workflows)
   * [QC](#qc)
      * [FRIP](#frip)
   * [Get blacklisted regions from UCSC](#get-blacklisted-regions-from-ucsc)
   * [Differential binding analysis](#differential-binding-analysis)
      * [Prepare the targets file](#prepare-the-targets-file)
      * [Prepare the contrasts file](#prepare-the-contrasts-file)
      * [DiffBind (Bioconductor)](#diffbind-bioconductor)
   * [Annotate peaks](#annotate-peaks)
   * [Motif discovery](#motif-discovery)
      * [rGADEM](#rgadem)
      * [MEME](#meme)
      * [BCRANK](#bcrank)
   * [Plots](#plots)
      * [Plot signal around features](#plot-signal-around-features)
      * [Distribution of peaks along the genome](#distribution-of-peaks-along-the-genome)
* [Plots](#plots-1)
   * [Color palettes](#color-palettes)
      * [A colorblind-friendly palette](#a-colorblind-friendly-palette)
      * [Color palettes with RColorBrewer](#color-palettes-with-rcolorbrewer)
      * [A very large (255) high contrast color palette](#a-very-large-255-high-contrast-color-palette)
   * [Spaghetti plot: graphs showing regression uncertainty](#spaghetti-plot-graphs-showing-regression-uncertainty)
   * [MA plot](#ma-plot)
   * [Circle](#circle)
   * [Ellipse around the CI95](#ellipse-around-the-ci-95)
   * [Boxplot with average and standard deviation](#boxplot-with-average-and-standard-deviation)
   * [Histogram and density function in the same plot using ''ggplot''](#histogram-and-density-function-in-the-same-plot-using-ggplot)
   * [Placing multiple plots in the same device](#placing-multiple-plots-in-the-same-device)
   * [Placing multiple plots in the same device using ''ggplot'' the ''grid'' package](#placing-multiple-plots-in-the-same-device-using-ggplot-the-grid-package)
   * [Venn Diagrams](#venn-diagrams)
   * [Calling gnuplot](#calling-gnuplot)
   * [Density 2d plots](#density-2d-plots)
      * [Base system](#base-system)
      * [smoothScatter](#smoothscatter)
      * [density2d (contour plot)](#density2d-contour-plot)
      * [ggplot2 density2d](#ggplot2-density2d)
      * [3d plot](#3d-plot)
   * [PCA BiPlot with ggplot2 and Interactive plots in a shiny app](#pca-biplot-with-ggplot2-and-interactive-plots-in-a-shiny-app)
   * [Raster kegg PNG pathways within a coordinate system](#raster-kegg-png-pathways-within-a-coordinate-system)
   * [Plot chromosome ideograms](#plot-chromosome-ideograms)
   * [Plot chromosome ideograms with additional tracks](#plot-chromosome-ideograms-with-additional-tracks)
   * [Gviz plots](#gviz-plots)
   * [Graphical representation of contingency tables](#graphical-representation-of-contingency-tables)
   * [ChIPseq plots](#chipseq-plots)
* [Miscellaneous bioinformatic related stuff](#miscellaneous-bioinformatic-related-stuff)
   * [Barcode design](#barcode-design)
   * [Getting things from Biomart](#getting-things-from-biomart)
   * [Getting things from the Ensembl databases](#getting-things-from-the-ensembl-databases)
   * [Getting things from the UCSC tables and ENCODE tracks](#getting-things-from-the-ucsc-tables-and-encode-tracks)
   * [Drug response curves using the ''drc'' package](#drug-response-curves-using-the-drc-package)
   * [Permutation test to identify if 2 curves are significantly different](#permutation-test-to-identify-if-2-curves-are-significantly-different)
   * [GenomicRanges](#genomicranges)
   * [Converting GTF to GFF3](#converting-gtf-to-gff3)
   * [Flatten a GFF/GTF file by gene_id (and get transcript lengths)](#flatten-a-gffgtf-file-by-gene_id-and-get-transcript-lengths)
   * [Download a dataset from GEO using GEOquery](#download-a-dataset-from-geo-using-geoquery)
   * [Get the genomic sequence of a region and plot its nucleotide content](#get-the-genomic-sequence-of-a-region-and-plot-its-nucleotide-content)
   * [Gene set enrichment analysis](#gene-set-enrichment-analysis)
      * [GSEA of Biological Processes](#gsea-of-biological-processes)
      * [GSEA of Biological Processes (a parallel version)](#gsea-of-biological-processes-a-parallel-version)
      * [GSEA of KEGG or custom db](#gsea-of-kegg-or-custom-db)
      * [GSEA of a custom GO (slim) db](#gsea-of-a-custom-go-slim-db)
      * [GO functional analysis with clusterProfiler](#go-functional-analysis-with-clusterprofiler)
      * [Reducing GO DAGs with Semantic Similarity](#reducing-go-dags-with-semantic-similarity)
      * [Summarizing by scoring frequencies of words](#summarizing-by-scoring-frequencies-of-words)
      * [Summarizing by clustering similarity scores between terms](#summarizing-by-clustering-similarity-scores-between-terms)
         * [Treemap representation of the hierarchy of terms](#treemap-representation-of-the-hierarchy-of-terms)
         * [2D multidimensional scaling of the semantic similarity scores](#2d-multidimensional-scaling-of-the-semantic-similarity-scores)
   * [Precomputing tables of semantic similarities between GO terms](#precomputing-tables-of-semantic-similarities-between-go-terms)
   * [Basic matrix normalization](#basic-matrix-normalization)
      * [Quantile normalization](#quantile-normalization)
      * [Smoothing curves normalization using loess](#smoothing-curves-normalization-using-loess)
   * [Calculate the reverse-complimentari of a sequence](#calculate-the-reverse-complimentari-of-a-sequence)
   * [Read a fasta file](#read-a-fasta-file)
   * [Get transcripts from a Bioconductor's AnnotationDb](#get-transcripts-from-a-bioconductors-annotationdb)
   * [Get information of genomic features](#get-information-of-genomic-features)
* [Statistical analysis](#statistical-analysis)
   * [Principal Component Analysis (PCA)](#principal-component-analysis-pca)
   * [Multidimensional Scaling (MDS)](#multidimensional-scaling-mds)
   * [Affinity Propagation Clustering](#affinity-propagation-clustering)
   * [Other methods for clustering](#other-methods-for-clustering)
   * [Validate the optimal number of clusters](#validate-the-optimal-number-of-clusters)
   * [Evaluate cluster strength](#evaluate-cluster-strength)
   * [Calculate ROC curves from a predictor using the ''ROCR'' package](#calculate-roc-curves-from-a-predictor-using-the-rocr-package)
   * [Impute missing values](#impute-missing-values)
   * [Power and sample size calculation for survival analysis](#power-and-sample-size-calculation-for-survival-analysis)
   * [Model normally distributed data](#model-normally-distributed-data)
   * [Calculate the center of a 2d-distribution](#calculate-the-center-of-a-2d-distribution)
   * [Test the significance of the overlap between 2 lists](#test-the-significance-of-the-overlap-between-2-lists)
   * [Network analysis](#network-analysis)
* [Other stuff that doesn't fit into any other category](#other-stuff-that-doesnt-fit-into-any-other-category)
   * [Supercomputing](#supercomputing)
   * [Parse arguments](#parse-arguments)
   * [Write data into an Excel shit](#write-data-into-an-excel-shit)
   * [Read data from an Excel shit](#read-data-from-an-excel-shit)
   * [SQLite database access](#sqlite-database-access)
   * [Show a progress bar](#show-a-progress-bar)
   * [Parallelize code that generates PDF plots](#parallelize-code-that-generates-pdf-plots)
   * [Parallel by](#parallel-by)

## Array analysis

### Infinium HumanMethylation450 BeadChip

This pipeline uses ''lumi'' and ''IlluminaHumanMethylation450k.db'' packages from ''bioconductor''.

#### Read data from raw IDAT files

```R
 #options(mc.cores=12)
 library("lumi")		#analisi d'arrays illumina
 library("IlluminaHumanMethylation450k.db")

 samplesheet <- read.csv("Curelung_SampleSheet.csv")
 barcodes <- paste(samplesheet$Sentrix_ID,samplesheet$Sentrix_Position,sep="_")
 #idats <- methylumIDAT(barcodes=barcodes,parallel=T,idatPath=paste(getwd(),"idats",sep="/"))
 idats <- methylumIDAT(barcodes=barcodes,idatPath=paste(getwd(),"idats",sep="/"))
 dades <- as(idats,"MethyLumiM")

 # Annotate the samples with some metainfo
 pData(dades)$sampleID <- samplesheet$Sample_Name
 pData(dades)$group    <- factor(samplesheet$Type)
```

#### Normalization
It can be done in a three steps procedure using the ''lumi'' package from ''bioconductor'':

* Color balance adjustment
* Background substraction
* Quantile normalization

```R
 # color balance adjustment: normalization between two color channels
 dades.c.adj <- lumiMethyC(dades)
 # background substraction
 dades.b.adj <- lumiMethyB(dades.c.adj,method="bgAdjust2C")
 # Perform quantile normalization based on color balance adjusted data
 dades.c.quantile <- lumiMethyN(dades.b.adj,method='quantile')
```

#### Quality control

```R
 plotColorBias1D(dades,channel='sum')
 plotColorBias1D(dades.c.quantile,channel='sum',main="Compare density distribution of two color channels - Norm")
 boxplotColorBias(dades,channel='sum')
 boxplotColorBias(dades.c.quantile,channel='sum',main="Boxplots of Red and Green color channels - Norm")
 plotColorBias2D(dades,selSample=1, cex=2)
 plotColorBias2D(dades.c.quantile,selSample=1,cex=2,main="Norm")
```

#### Differential methylation using limma

```R
 ## filter non variant probes by IQR
 eset <- exprs(dades.c.quantile)
 eset.IQR <- varFilter(eset,var.func=IQR,var.cutoff=0.5,filterByQuantile=TRUE)

 # Define design and contrast matrices
 design <- model.matrix(~0 + pData(dades.c.quantile)$grup)
 rownames(design) <-pData(x)$sampleID
 colnames(design) <- c("G1","G2")
 cont.matrix <- makeContrasts(Grup1vsGrup2=G1-G2,levels=design)

 # Get the differentially methylated probes by linear model + bayes stats
 # Model lineal i eBayes
 fit1 <- lmFit(eset.IQR,design)
 fit2 <- contrasts.fit(fit1,cont.matrix)
 fit3 <- eBayes(fit2)

 # Rank the probes
 CG <- topTable(fit3,adjust="fdr",number=nrow(fit3))	# adjust FDR
```

### Affymetrix® Genome-Wide Human SNP Array 6.0 (and probably others)

This is a parallelized version using ''snow'' and 16 cores to analyze SNP arrays with the ''oligo'' package. Can work with the ''MC'' package as well

```R
 library(ff)
 library(foreach)
 library(snow)
 library(doSNOW)
 registerDoSNOW(makeCluster(16,type="SOCK"))
 #library(MC)
 #library(doMC)
 #registerDoMC(16)

 # read arrays
 library(oligo)
 library(pd.genomewidesnp.6)
 #library(human650v3aCrlmm)	# Illumina HumanHap650Y
 #library(human550v3bCrlmm)	# Illumina HumanHap550Y

 fullFilenames <- list.celfiles("arrays",full.names=TRUE)
 outputDir <- file.path(getwd(),"crlmmResults")

 crlmm(fullFilenames,outputDir,pkgname="pd.genomewidesnp.6")
 #crlmm(fullFilenames,outputDir,pkgname="human650v3aCrlmm")
 #crlmm(fullFilenames,outputDir,pkgname="human550v3bCrlmm")
 crlmmOut <- getCrlmmSummaries(outputDir)

 write.csv(calls(crlmmOut),file="SNP.csv")
```

### Affymetrix expression arrays

Tested with Affymetrix Human Genome U133 Plus 2.0 Array, but might work with other from the U133 series and the newer Exon ST and Gene ST. Also might work with other non-Human organisms.

This protocol uses the ''affy'' package from ''bioconductor''

#### Read data from raw CEL files

```R
 ## Load libraries
 library("affy")       #Import and normalize data
 library("limma")      #Differential expression
 library("genefilter") #Gene filtering

 ## Import sample description from targets.txt file
 targets <- readTargets("arrays/targets.txt", row.names="FileName")

 ## Import .CEL files
 data <- ReadAffy(filenames=targets$FileName)
```

#### Normalization

Standard RMA and quantiles normalization using the ''expresso'' function:

```R
 eset <- expresso(
         data,
         bg.correct = TRUE,
         bgcorrect.method="rma",
         normalize = TRUE,
         normalize.method="quantiles",
         pmcorrect.method="pmonly",
         summary.method="medianpolish",
         verbose = TRUE,
         )
```

#### Quality control

Almost useless quality control information:

```R
 par(mfrow = c(1, 2))    #partir la ventana de gráficos en 2 columnas
 # boxplot de datos brutos
 boxplot(
         data,
         main="Boxplot Before Normalization",
         col = "lightgrey")

 # boxplot de datos normalizados         
 exprseset <- as.data.frame(exprs(eset))         
 boxplot(
         data.frame(exprseset),
         main="Boxplot After Normalization (log scale)",
         col = "lightgrey")
```

#### Differential expression using limma

```R
 ## filter non variant probes by IQR
 esetIQR <- varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)

 # Define design and contrast matrices
 design <- model.matrix(~0 + as.factor(targets$Classes))
 colnames(design) <- c("AA","CA","HC")
 rownames(design) <- targets$FileName
 cont.matrix<-makeContrasts(AAvsCA=AA-CA,AAvsHC=AA-HC,CAvsHC=CA-HC,levels=design)

 # Get the differentially methylated probes by linear model + bayes stats
 # Model lineal and eBayes
 fit1 <- lmFit(esetIQR,design)
 fit2 <- contrasts.fit(fit1, cont.matrix)
 fit3 <- eBayes(fit2)

 # List of differentially expressed genes
 toptableIQR <- topTable(fit3, number=nrow(fit3), adjust.method="BH", sort.by="F")
```

## ChIP-seq workflows

### QC
#### FRIP
```R
library(rtracklayer)
library(GenomicRanges)
library(parallel)
CORES <- 4
IN <- "/fsimb/groups/imb-bioinfocf/projects/soshnikova/imb_soshnikova_meta_01_kazakevych_integration/results/chipseq/"

#' @examples
#' library(rtracklayer)
#' frip("chip.bam", import.bed("peaks.bed"))
frip <- function(reads, peaks, singleEnd=T) {
    require(GenomicAlignments)

    if (is.character(reads)) {
        require(Rsamtools)
        reads <- BamFile(reads)
    }

    # find reads in peaks
    overlaps <- summarizeOverlaps(
            peaks,
            reads,
            mode="IntersectionNotEmpty",
            ignore.strand=T,
            singleEnd=singleEnd,
            count.mapped.reads=T
    )

    # sum up all reads in peaks and divide by all mapped reads to get FRiP
    readsInPeaks <- sum(assay(overlaps))
    if (class(reads) %in% c("BamViews", "BamFile")) {
        allReads <- colData(overlaps)$mapped
    } else {
        allReads <- colData(overlaps)$records
    }
    result <- readsInPeaks/allReads

    return(result)
}

f <- list.files(IN,pattern="*.tsv")
x <- mclapply(f,function(f) {
    x <- readLines(paste0(IN,f),n=10)
    x <- unlist(strsplit(x[grepl("# Command line:",x)]," "))
    f.bam <- x[grep("-t",x) + 1]
    peaks <- read.delim(paste0(IN,f),comment.char="#")
    peaks <- with(peaks,GRanges(chr,IRanges(start,end)))
    frip(f.bam,peaks)
},mc.cores=CORES)

write.csv(x,file="./qc/frip.csv",row.names=F)
```

### Get blacklisted regions from UCSC

There are two main tables in UCSC which contain info about abnormal high signal in ChIP-like experiments. They're both derived from the Encode project:

```R
library(rtracklayer)
library(GenomicRanges)

mySession <- browserSession("UCSC")
genome(mySession) <- "hg19"
blcklst <- rbind(getTable(ucscTableQuery(mySession, table="wgEncodeDacMapabilityConsensusExcludable")),
                 getTable(ucscTableQuery(mySession, table="wgEncodeDukeMapabilityRegionsExcludable")))
blcklst <- with(blcklst, GRanges(chrom, IRanges(chromStart, chromEnd)))
```

### Differential binding analysis

#### Prepare the targets file
<pre>
SampleID     Condition  Replicate  bamReads                       ControlID     bamControl                 Peaks                      PeakCaller
AE_1_h2az    AE         1          ./mapped/AE_1_h2az.bam         AE_1_input    ./mapped/AE_1_input.bam    ./results/AE_1_h2az.xls    macs
AE_2_h2az    AE         2          ./mapped/AE_2_h2az.bam         AE_2_input    ./mapped/AE_2_input.bam    ./results/AE_2_h2az.xls    macs
E125_1_h2az  E125       1          ./mapped/E125_1_h2az.bam       E125_1_input  ./mapped/E125_1_input.bam  ./results/E125_1_h2az.xls  macs
E125_2_h2az  E125       2          ./mapped/E125_2_h2az.bam       E125_2_input  ./mapped/E125_2_input.bam  ./results/E125_2_h2az.xls  macs
E145_1_h2az  E145       1          ./mapped/h2az/E145_1_h2az.bam  E145_1_input  ./mapped/E145_1_input.bam  ./results/E145_1_h2az.xls  macs
E145_2_h2az  E145       2          ./mapped/E145_2_h2az.bam       E145_2_input  ./mapped/E145_2_input.bam  ./results/E145_2_h2az.xls  macs
ISC_1_h2az   ISC        1          ./mapped/ISC_1_h2az.bam        ISC_1_input   ./mapped/ISC_1_input.bam   ./results/ISC_1_h2az.xls   macs
ISC_2_h2az   ISC        2          ./mapped/ISC_2_h2az.bam        ISC_2_input   ./mapped/ISC_2_input.bam   ./results/ISC_2_h2az.xls   macs
</pre>

#### Prepare the contrasts file
<pre>
E145vsE125=(E145-E125)
ISCvsE145=(ISC-E145)
AEvsISC=(AE-ISC)
</pre>

#### DiffBind (Bioconductor)
```R
library(DiffBind)

CWD        <- "/fsimb/groups/imb-bioinfocf/projects/soshnikova/imb_soshnikova_meta_01_kazakevych_integration"
BAMS       <- paste0(CWD, "/mapped/chipseq/h2az")
FTARGETS   <- "./scripts/timecourse_h2az_targets.txt"
FCONTRASTS <- "./scripts/timecourse_h2az_contrasts.txt"
pdf("./results/timecourse_h2az.pdf")

setwd(CWD)

# load targets and make analysis
conts   <- read.delim(FCONTRASTS, head=F, comment.char="#")
targets <- read.delim(FTARGETS  , head=T, colClasses="character", comment.char="#", )
h2az <- dba(sampleSheet=targets, config=data.frame(fragmentSize=200, bCorPlot=F))
h2az <- dba.count(h2az)
dba.plotPCA(h2az, DBA_CONDITION, label=DBA_CONDITION)

# parse the formula in cont and do the analysis
h2az.DB <- lapply(conts[, 1], function(cont) {
    cat(cont, fill=T)
    cont.name <- gsub("(.+)=(.+)", "\\1", cont)
    cont.form <- gsub("(.+)=(.+)", "\\2", cont)
    factors   <- unlist(strsplit(cont.form, "\\W"))
    factors   <- factors[factors != ""]

    c1 <- dba.mask(h2az, DBA_CONDITION, factors[1])
    c2 <- dba.mask(h2az, DBA_CONDITION, factors[2])
    db <- dba.contrast(h2az, group1=c1, group2=c2,  name1=factors[1], name2=factors[2], categories=DBA_CONDITION)
    db <- dba.analyze(db)
    dba.plotMA(db)
    dba.plotBox(db)
    dba.plotVenn(db, c1 | c2)

    # retrieve DB sites with p-value <  0.05 and Fold > 2
    dba.report(db, bCalled=T, th=.05, bUsePval=TRUE, fold=2)
})

dev.off()
```

And continue now with the next step: Annotate peaks.

### Annotate peaks
Supose h2az.DB is a list with peaks read from many MACS output files (.xls):

```R
library(WriteXLS)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

h2az.DB <- lapply(list.files(pattern=".xls$"), read.delim, comment.char="#")

txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
ann <- lapply(h2az.DB, annotatePeak, TxDb=txdb, annoDb="org.Mm.eg.db", tssRegion=c(-3000,3000), verbose=T)
lapply(ann, plotAnnoBar)
lapply(ann, plotDistToTSS)
ann <- lapply(ann, as.data.frame)

WriteXLS("ann", ExcelFileName="./results/timecourse_h2az.xls",
         SheetNames=gsub("(.+)=\\((.+)\\)", "\\2", conts[,1]))
```

### Motif discovery

#### rGADEM

```R
library(rGADEM)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

PROJECT <- "/fsimb/groups/imb-bioinfocf/projects/butter/imb_butter_2016_02_bluhm_ChIPseq/"
TOP <- 50   # take only top peaks (values >50 usually make GADEM cracsh with segfault)

# read in the peaks and call the motif discovery tool
bed <- read.table(paste0(PROJECT, "/GSE51142/results/macs2/zbtb10.vs.igg_peaks.xls"), comment="#", sep="\t", head=T)
sel <- rev(order(bed$X.log10.qvalue.))[1:TOP]
bed$chr <- paste0("chr", bed$chr)
peaks <- with(bed[sel,], RangedData(IRanges(start, end), space=chr))
gadem <- GADEM(peaks, verbose=1, genome=Hsapiens)
```

#### MEME

I had little success with rGADEM, being more successful with MEME. Try [[Motif discovery]], or call it from R:

```R
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(ShortRead)

CWD <- "/fsimb/groups/imb-bioinfocf/projects/jgu/jgu_berger_2016_01_kaiser_RNA-Seq/with_UMIs"
MEME_EXE <- "/fsimb/groups/imb-bioinfocf/common-tools/dependencies/meme/4.11.2/bin/meme-chip"
MEME_RES <- paste0(CWD, "/results/meme")
MEME_CORES <- 8
MEME_DB <- c("/fsimb/groups/imb-bioinfocf/common-tools/dependencies/meme/motif_databases/FLY/dmmpmm2009.meme",
             "/fsimb/groups/imb-bioinfocf/common-tools/dependencies/meme/motif_databases/FLY/fly_factor_survey.meme",
             "/fsimb/groups/imb-bioinfocf/common-tools/dependencies/meme/motif_databases/FLY/flyreg.v2.meme",
             "/fsimb/groups/imb-bioinfocf/common-tools/dependencies/meme/motif_databases/FLY/idmmpmm2009.meme",
             "/fsimb/groups/imb-bioinfocf/common-tools/dependencies/meme/motif_databases/FLY/OnTheFly_2014_Drosophila.meme",
             "/fsimb/groups/imb-bioinfocf/common-tools/dependencies/meme/motif_databases/JASPAR/JASPAR_CORE_2016_insects.meme")

#...
# supose we read the peaks as a dataframe
#...

pos <- with(peaks, GRanges(seqnames, ranges=IRanges(start, end)))
s <- getSeq(Dmelanogaster, pos)
names(s) <- sapply(pos, paste, sep="_")
writeFasta(s, outfa <- paste0("./results/starr_rnaseq.motif.", contrast_name, ".fa"))
system(paste(MEME_EXE, "-oc", MEME_RES, paste("-db", MEME_DB, collapse=" "), paste("-meme-p", MEME_CORES),
             "-time 300 -order 1 -meme-mod zoops -meme-minw 5 -meme-maxw 30",
             "-meme-nmotifs 3 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0", outfa))
```

#### BCRANK

Yet another algorithm, this time suggested in a Bioconductor [http://biocluster.ucr.edu/~rkaundal/workshops/R_feb2016/ChIPseq/ChIPseq.html#count-reads-overlapping-the-peak-regions course]:

```R
library(BCRANK)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare sequences
pos <- with(peaks.ann, GRanges(seqnames, ranges=IRanges(start, end)))
s <- getSeq(Hsapiens, pos)
writeFasta(s, "./tmp/motif.fa")

# get motif
BCRANKout <- bcrank("./tmp/motif.fa", restarts=25, use.P1=TRUE, use.P2=TRUE)
toptable(BCRANKout)
topMotif <- toptable(BCRANKout, 1)
weightMatrix <- pwm(topMotif, normalize = FALSE)
weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
seqLogo(weightMatrixNormalized)
```

### Plots
#### Plot signal around features

The ChIPseeker package suggests taking TSS of known genes. Actually, it can be done with any kind of GRanges object.

In this example, we compare colocalization of 2 different peak marks. We plot how it looks like for a ChIP around anoter ChIP.

```R
chip1     <- GRanges(chip1$chr, IRanges(start=chip1$start, end=chip1$end))
chip2     <- GRanges(chip2$chr, IRanges(start=chip2$start, end=chip2$end))
chip2.3kb <- GRanges(seqnames(chip2), IRanges(start=chip2$summit - 1500, end=chip$summit + 1500))
tagMatrix <- getTagMatrix(chip1, windows=chip2.3kb)

image(x=-1500:1500, y=1:nrow(tagMatrix), z=t(tagMatrix), useRaster=TRUE, col=c("black", "yellow"),
      yaxt="n", xaxt="n", ylab="", xlab="")
title("mark 1 on mark 2", cex.main=.5)
Axis(side=1, at=c(-1500, 0, 1500), labels=c("-1500", "mark 2", "+1500"), cex.axis=.5)
```

#### Distribution of peaks along the genome

Maybe see also [Plot chromosome ideograms with additional tracks](#plot-chromosome-ideograms-with-additional-tracks)

'''Read peaks:'''

```R
peaks <- read.xls("./results/diffbind_cpg.xls", sheet=1, header=TRUE)
peaks$seqnames <- paste0("chr", peaks$seqnames)
validchr <- paste0("chr", rev(c(as.character(1:22), "X", "Y")))
```

'''Get centromeres positions from UCSC cytoband table:'''

```R
library(rtracklayer)
mySession <- browserSession("UCSC")
genome(mySession) <- "hg38"
chrlengths  <- getTable(ucscTableQuery(mySession, table="chromInfo"))
chrlengths  <- chrlengths[chrlengths$chrom %in% validchr, ]
centromeres <- getTable(ucscTableQuery(mySession, track="cytoBand"))
centromeres <- centromeres[centromeres$gieStain == "acen" & centromeres$chrom %in% validchr, ]
```

'''Plot the chromosomes and the centromeres:'''

```R
plot(0, type="n", xlim=c(0, max(chrlengths$size)), ylim=c(0, length(validchr)), bty="n", xlab="", ylab="", yaxt="n")
axis(side=2, at=1:length(validchr), labels=validchr, las=2)
x <- sapply(validchr, function(chr) {
    # chromosome sizes
    lines(x=c(0, chrlengths$size[chrlengths$chrom == chr]),
          y=rep(which(validchr == chr), 2))

    # centromeres
    rect(xleft  =min(centromeres$chromStart[centromeres$chrom == chr]),
         xright =max(centromeres$chromEnd  [centromeres$chrom == chr]),
         ybottom=which(validchr == chr) - .2,
         ytop   =which(validchr == chr) + .2,
         border=NA, col="#000000AA")
})
```

'''Draw the peaks on them:'''

```R
library(RColorBrewer)
peaks$x <- (peaks$start + peaks$end) / 2
peaks$y <- sapply(peaks$seqnames, function(chr) which(validchr == chr))
x <- sapply(validchr, function(chr) {
    peaks.chr <- peaks[peaks$seqnames == chr, ]
    if(nrow(peaks.chr) > 0) {
        y <- which(validchr == chr)
        dcols <- densCols(x=peaks.chr$x, colramp=colorRampPalette(brewer.pal(9, "Reds")))
        dcols <- paste0(dcols, "50")
        points(peaks.chr$x, rep(y, length(peaks.chr$x)), col=dcols, pch=16)
    }
})
```

'''Density of peaks along the chromosome:'''

```R
library(ggplot2)
df.peaks <- data.frame(chromosome=factor(peaks$seqnames, levels=rev(validchr)), pos=peaks$x)
df.centromeres <- do.call(rbind, by(centromeres, centromeres$chrom, function(x) data.frame(chromosome=x$chrom[1], start=min(x$chromStart), end=max(x$chromEnd))))
df.centromeres$chromosome <- factor(df.centromeres$chromosome, levels=rev(validchr))
ggplot(df.peaks, aes(pos)) +
   geom_density(fill="grey", alpha=.5, adjust=1/5) +
   geom_segment(data=df.centromeres, mapping=aes(x=start, xend=end, y=0, yend=0), col="red", size=2) +
   facet_wrap(~ chromosome, ncol=2, scales="free", strip.position="left") +
   theme_bw() +
   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
   theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
```

## Plots

### Color palettes

#### A colorblind-friendly palette

To use with ggplot2, it is possible to store the palette in a variable, then use it later.

```R
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
scale_colour_manual(values=cbPalette)
```

#### Color palettes with RColorBrewer

```R
require("RColorBrewer")

# show all available palettes
display.brewer.all()

# show the palette we are planning
display.brewer.pal(9,"Oranges")

# use it as a gradient palette for heatmaps
hmcol   <- colorRampPalette(brewer.pal(9,"Oranges"))(100)
heatmap.2(mat,trace="none",col=rev(hmcol),margin=c(13,13))

# use it to display categorical data
plot(x=x$x[,1],y=x$x[,2],pch=16,cex=.5,
     col=brewer.pal(length(levels(condition)),"Set1")[condition])
```

#### A very large (255) high contrast color palette

```R
pal <-
c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
  "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
  "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
  "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
  "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
  "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
  "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
  "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
  "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
  "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
  "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
  "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
  "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
  "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
  "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
  "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
  "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
  "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
  "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
  "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
  "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
  "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
  "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",
  "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
  "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
  "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
  "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
  "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
  "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
  "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
  "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
  "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
  "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B")
```
### Spaghetti plot: graphs showing regression uncertainty

The idea is taken from [http://andrewgelman.com/2012/08/26/graphs-showing-regression-uncertainty-the-code/ here] and is only partially implemented (only the spaghetti code, not the the shaded CI.

```R
# plot data points
x <- rnorm(100)
y <- 2 * x + rnorm(100, mean=1)^3
plot(x, y, pch=16)#, type="n")

# spaghetti
ci <- replicate(1000, {
    i <- sample(length(x), replace=TRUE)
    newx <- x[i]
    newy <- y[i]
    fit <- loess(newy ~ newx)
    lines(newx[order(newx)], fit$fitted[order(newx)], col="#0000FF10")
    data.frame(x=newx[order(newx)], y=fit$fitted[order(newx)])
}, simplify=FALSE)
ci <- do.call(rbind, ci)

# consensus fit
#fit <- loess(ci$y ~ ci$x)  # takes forever to compute
#lines(ci$x[order(ci$x)], fit$fitted[order(ci$x)], col="red", lwd=2)
lines(lowess(ci$y ~ ci$x, f=1/3), col="red", lwd=2)
```

### MA plot

Visual representation of log ratios (M) and mean average (A) intensities. Useful to compare bias between two samples or in a two-channels array.

```R
 maplot <- function(x,y,...){
 	M=log2(y/x)
 	A=(log2(y)+log2(x))/2
 	plot(A,M,...)
 	fit=loess(M~A,span=1/3)
 	o=order(A)
 	lines(A[o],fit$fitted[o],col=2)
 }

 maplot(y[wh,2],y[wh,3],ylim=c(-5,5))
```

### Circle

There is no primitive function to draw a circle in R.

```R
circle <- function(x,y,r,...) {	# col=fill_colour,border=outline_colour
    polygon(x=x + r*cos( seq(0,2*pi, length.out=360) ),
            y=y + r*sin( seq(0,2*pi, length.out=360) ),
            ...)
}
```

### Ellipse around the CI95

Draw an ellipse describing the CI 95% of some (x,y) points. Useful to describe groups when reducing to only a couple of components (through a PCA or MDS) large datasets with several independent variables describing a sample.

```R
ell95 <- function(df) {
	require(ellipse)
	require(ggplot2)

	ell <- data.frame()
	for(g in levels(df$group)) {
		e <- with(df[df$group==g,],ellipse(cor(x,y),scale=c(sd(x),sd(y)),centre=c(mean(x),mean(y))))
		ell <- rbind(ell,cbind(as.data.frame(e),group=g))
	}
	ggplot(data=df,aes(x=x,y=y,colour=group)) + geom_point() + geom_path(data=ell,aes(x=x,y=y,colour=group))
}

print(ell95(data.frame(x=x,y=y,group=A)))
```

### Boxplot with average and standard deviation

```R
b   <- boxplot(x$gini ~ x$stage, ylab="coef inequality", col=pal, boxwex=.25)
avg  <- tapply(x$gini, x$stage, mean)
sdev <- tapply(x$gini, x$stage, sd)
xi <- 0.1 + seq(b$n)
points(xi, avg, col="orange", pch=18)
arrows(xi, avg - sdev, xi, avg + sdev, code=3, col="orange", angle=75, length=.1)
```
### Histogram and density function in the same plot using ''ggplot''

```R
require(ggplot2)

df <- data.frame(cond=c(rep("before",length(before)),rep("after",length(after))),
                 vals=c(before,after))
p <- ggplot(df,aes(x=vals,fill=cond)) +
	geom_histogram(aes(y=..density..,fill=cond),alpha=.5,position="identity") +
	geom_density(alpha=.2) +
	ggtitle(rownames(ic50)[i])

print(p)
```

### Placing multiple plots in the same device
This solution makes use of the basic R plotting capabilities, through the ''graphics'' package:

```R
# define grid (layout matrix)
layoutRatioWidth <- c(0.75,0.25); layoutRatioHeight <- c(0.25, 0.75)
layout(matrix(c(2,1,0,3),nrow=2),widths=layoutRatioWidth,heights=layoutRatioHeight,respect=F)
# bottom left plot (#1 in layout matrix)
par(mar=c(4,4,1,1))
smoothScatter(x,y,xlab=xlab,ylab=ylab,main=f,xlim=c(0,20),ylim=c(0,20))
# top left plot (#2 in layout matrix)
par(mar=c(1,4,1,1))
z <- density(x,na.rm=T)
plot(z$x,z$y,xlab="",ylab="Density",xaxt="n",yaxt="n",type="l",main="")
# bottom right plot (#3 in layout matrix)
par(mar=c(4,1,1,1))
z <- density(y,na.rm=T)
plot(z$y,z$x,ylab="",xlab="Density",xaxt="n",yaxt="n",type="l",main="")
```

### Placing multiple plots in the same device using ''ggplot'' the ''grid'' package

First option, using the ''grid'' package and creating ''viewports'':

```R
 library("grid")

 vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

 plot1 <- qplot(mtcars,x=wt,y=mpg,geom="point",main="Scatterplot of wt vs. mpg")
 plot2 <- qplot(mtcars,x=wt,y=disp,geom="point",main="Scatterplot of wt vs disp")
 plot3 <- qplot(wt,data=mtcars)
 plot4 <- qplot(wt,mpg,data=mtcars,geom="boxplot")

 # 4 figures arranged in 2 rows and 2 columns
 grid.newpage()
 pushViewport(viewport(layout = grid.layout(2, 2)))
 print(plot1, vp = vplayout(1, 1))
 print(plot2, vp = vplayout(1, 2))
 print(plot3, vp = vplayout(2, 1))
 print(plot4, vp = vplayout(2, 2))
```

Second option, use the ''gridExtra'' package:

```R
 library("gridExtra")

 plot1 <- qplot(...)
 plot2 <- qplot(...)
 sidebysideplot <- grid.arrange(plot1, plot2, ncol=2)
```

### Venn Diagrams
Several Bioconductor packages are available for Venn Diagrams, "VennDiagram" is able to produce high quality proportional diagrams.

```R
require(VennDiagram)
## main argument a list of vectors
## function will calculate overlap between the vectors and produce proportional Venn diagrams.

venn.diagram(list(C = 1700:2500, B = 1:1800, A = 1571:2020),
  fill=RColorBrewer::brewer.pal(3,"Set1"),alpha=.3,lwd=0,cex=1.5,cat.cex=1.5,
  filename="~/Desktop/test2.tif")

# display it on the screen
img <- tiff::readTIFF("~/Desktop/test2.tif",native=T)
plot(1:2,type="n",bty="n",axes=F,xlab="",ylab="")
rasterImage(img,1,1,2,2)
```

To display it on the screen, alternatively one can set to 'filename' to NULL and use the grid package to draw the result:

```R
img <- venn.diagram(list(C = 1700:2500, B = 1:1800, A = 1571:2020),
  fill=RColorBrewer::brewer.pal(3,"Set1"),alpha=.3,lwd=0,cex=1.5,cat.cex=1.5,filename=NULL)
grid.draw(img)
```

Another option is manually calculating the counts table and call the Venn functions from the limma package:

```R
require(limma)

Venn3 <- function(set1, set2, set3, names)
{
  stopifnot( length(names) == 3)

  # Form universe as union of all three sets
  universe <- sort( unique( c(set1, set2, set3) ) )

  Counts <- matrix(0, nrow=length(universe), ncol=3)
  colnames(Counts) <- names

  for (i in 1:length(universe))
  {
    Counts[i,1] <- universe[i] %in% set1
    Counts[i,2] <- universe[i] %in% set2
    Counts[i,3] <- universe[i] %in% set3
  }

  vennDiagram( vennCounts(Counts) )
}

set1 <- letters[1:3]
set2 <- letters[3:5]
set3 <- letters[5:7]

Venn3(set1, set2, set3, c("set1", "set2", "set3"))
```

Yet another option for proportional Venn diagrams with the venneuler package:

```R
plot(venneuler::venneuler(reshape2::melt(list(C=1700:2500, B=1:1800, A=1571:2020))))
```

Colors, labels and other properties can be changed from the VennDiagram object generated by the venneuler::venneuler() call.

Other examples about doing it manually can be found in [https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html this] page.

### Calling gnuplot

From R, one can pipe commands into gnuplot after opening a new session:

```R
library(Rgnuplot)
h1<-gp.init()
gp.cmd(h1,"set terminal pngcairo  transparent enhanced font 'arial,10' fontscale 1.0 size 500, 350")
gp.cmd(h1,"set output 'out.png'")	# set output file
gp.cmd(h1,"set key left box")		# include a boxed legend
gp.cmd(h1,"plot [-10:10] sin(x),atan(x) with points,cos(atan(x)) with impulses")
gp.close(h1)
```

More info in [[Plotting with gnuplot from R]]

### Density 2d plots

#### Base system

```R
dpal <- c("cyan", "blue", "green", "yellow", "red")
dcols <- densCols(x=x, y=y, colramp=colorRampPalette(dpal), nbin=500)
plot(x=x, y=y, col=dcols, pch=20, cex=.25)
```

And now add a legend with the density

```R
dd <- grDevices:::.smoothScatterCalcDensity(data.frame(x=x, y=y), nbin=500)
dens <- as.numeric(dd$fhat)
dens <- dens[dens>0]
colLegend <- data.frame(density=seq(min(dens), max(dens), len=10), color=I(colorRampPalette(dpal)(10)))
par(mar=c(5,0,5,0))
plot(NA, xlim=c(0,10), ylim=c(0,11), type="n", ann=FALSE, axes=FALSE)
rect(0, 1:10, 1, 2:11, border=NA, col=colLegend$col)
text(2, (1:10)+0.5, signif(colLegend$density, 2), adj=0)
```

#### smoothScatter
```R
smoothScatter(log10(DupMat[, "RPK"]), 100 * DupMat[, "dupRate"],
xlab = "expression level (reads/kbp)", ylab = "duplication level (% duplicate reads)",
axes = FALSE, ...)
```

#### density2d (contour plot)
The kernel estimator is borrowed from the MASS package (Modern Applied Statistics with S), usually part of any R installation.

```R
require(MASS)
pal <- colorRampPalette(c('dark blue','blue','light blue','yellow','orange','red','dark red'))
filled.contour(kde2d(log10(DupMat[, "RPK"]), 100 * DupMat[, "dupRate"]),
xlab = "expression level (reads/kbp)", ylab = "duplication level (% duplicate reads)",
axes = FALSE, color.palette=pal, ...)
```

The palette can be modified with RColorBrewer to visually improve the result:
```R
require(RColorBrewer)
pal <- function(n) { colorRampPalette(brewer.pal(9,"Oranges"))(n) }
```

Something easier, leaving the estimation to smoothScatter():
```R
cols <- colorRampPalette(c("black","blue","green","yellow","red"))
smoothScatter(dm,nrpoints=10,nbin=500,colramp=cols)
```

#### ggplot2 density2d
```R
df <- data.frame(x=log10(DupMat[, "RPK"]), y=100 * DupMat[, "dupRate"])
p <- ggplot(df,aes(x,y)) +
stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
scale_x_continuous(breaks=s,labels=as.character(10^s)) +
xlab("expression level (reads/kbp)") + ylab("duplication level (% duplicate reads)") +
theme_bw()
print(p)
```

#### 3d plot
```R
require(MASS)
pal <- colorRampPalette(c('dark blue','blue','light blue','yellow','orange','red','dark red'))
persp(kde2d(log10(DupMat[, "RPK"]), 100 * DupMat[, "dupRate"]),
theta = 30, phi = 30, expand = 0.5, ltheta = 120, shade = 0.75, ticktype = "detailed",
xlab = "expression level (reads/kbp)",
ylab = "duplication level (% duplicate reads)",
zlab = "density") -> res
```

### PCA BiPlot with ggplot2 and Interactive plots in a shiny app

Apart from the inherently ugly builtin biplot that comes with vanilla R, one can use the ggplot2 framework to improve the visual appearance.

The code here can be extended with other [https://github.com/kassambara/factoextra/blob/master/R/fviz_pca.R
http://www.sthda.com/english/wiki/principal-component-analysis-in-r-prcomp-vs-princomp-r-software-and-data-mining PCA visualization ideas].

```R
library(plotly)
library(shiny)

PCbiplot <- function(PC, x="PC1", y="PC2") {
    # PC being a prcomp object
    data <- data.frame(obsnames=row.names(PC$x), PC$x)
    plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
    plot <- plot + geom_hline(aes(yintercept=0), size=.2) + geom_vline(aes(xintercept=0), size=.2)
    datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
    mult <- min(
        (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
        (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
        )
    datapc <- transform(datapc,
            v1 = .7 * mult * (get(x)),
            v2 = .7 * mult * (get(y))
            )
    plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
    plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
    plot + theme_bw()
}

fit <- prcomp(USArrests, scale=T)
shinyApp(ui=bootstrapPage(plotlyOutput("plot")), server=function(input, output) { output$plot <- renderPlotly({ PCbiplot(fit) }) })
```

### Raster kegg PNG pathways within a coordinate system

'''basic plotting system''':

```R
library(png)
library(XML)

pathway <- xmlParse("http://rest.kegg.jp/get/hsa05219/kgml")
pathway <- xmlToList(pathway)
download.file("http://rest.kegg.jp/get/hsa05219/image", "~/Downloads/hsa05219.png")

img <- readPNG("~/Downloads/hsa05219.png")

# draw kegg pict
par(mar=c(0, 0, 0, 0))
plot(0, xlim=c(0, ncol(img)), ylim=c(0, nrow(img)), type="n", bty="n", axes=F, xlab="", ylab="")
rasterImage(img, 0, 0, ncol(img), nrow(img))                                                                                        

# draw something on p53
x <- 796                # center coordinates of the box
y <- nrow(img) - 336    # center coordinates of the box
w <- 46                 # width
h <- 17                 # height
x.ini <- x - w / 2
x.end <- x + w / 2
y.ini <- y - h / 2
y.end <- y + h / 2

# overwrite rectangles and add labels
rect(x.ini, y.ini, x.end, y.end, col="white")     # blank the box
rect(x.ini, y.ini, x.end, y.end, col="#FF000040") # add color depending on the FC
text(x, y, "BRCA1 (2.5, 4.4e-4)\nBRCA2 (1.1, 1.1e-2)", adj=c(.5, .5), cex=.6, bg="black")
```

'''ggplot2''':

```R
library(png)
library(XML)
library(grid)
library(ggplot2)

pathway <- xmlParse("http://rest.kegg.jp/get/hsa05219/kgml")
pathway <- xmlToList(pathway)
download.file("http://rest.kegg.jp/get/hsa05219/image", "~/Downloads/hsa05219.png")

img <- readPNG("~/Downloads/hsa05219.png")

g <- rasterGrob(img, 0, 0, ncol(img), nrow(img), just=c("left","bottom"))
ggplot(data.frame(x=c(0, ncol(img)), y=c(0, nrow(img))), aes(x=x, y=y)) +
    annotation_custom(g, 0, 1, 0, 1) +
    geom_point(data=data.frame(x=796, y=(nrow(img) - 336)), mapping=aes(x=x, y=y), col="red") +
    theme_void()
```

### Plot chromosome ideograms

Download the chromosome bands info from UCSC Table Browser and plot in R with the grid package:

```R
ideogramTab <- read.csv("~/idibell/varis/ideogramTab.csv")

# funcio que torna un objecte ggplot amb l'ideograma del cromosoma
dibuixaIdeograma <- function(chr) {

        ideo <- ideogramTab[ideogramTab$chr == chr,]
        center <- ideo[ideo$braç=="cen","end"]
        chrLength <- max(ideo$end)
        # el viewport comença un 10% abans de la coordenada 0, i acaba un 20% mes tard
        pushViewport(dataViewport(xData=c(0-chrLength * .18,chrLength * 1.32),yData=c(0,6),extension=0,layout.pos.col=1,layout.pos.row=1))      

        # pintem rectangles per a cada banda. El color be definit per:
        #     -si V8=="gpos" --> un gradient de la paleta de grisos [1:100] segons V9
        #     -altrament --> blanc. En aquest cas V9 == NA
        ideo[is.na(ideo[,9]),9] <- 1
        pal <- colorRampPalette(c("white","black"))(100)
    for(i in seq(along = ideo[,1])) {
                grid.rect(x=ideo[i,6],
                y=2,
                width=ideo[i,7]-ideo[i,6],
                height=2,
                gp=gpar(col=pal[ideo[i,9]],
                fill=pal[ideo[i,9]]),
                default.units="native",
                just=c("left","bottom"))
        }

        # center == centromer (on creuarem els dos braços)
        # requadre al braç p
        grid.lines(c(0,center-500000),c(4,4),default.units = "native")
        grid.lines(c(0,center-500000),c(2,2),default.units = "native")
        grid.lines(c(0,0),c(2,4),default.units = "native")
        # requadre al braç q
        grid.lines(c(center+500000,chrLength),c(4,4),default.units = "native")
        grid.lines(c(center+500000,chrLength),c(2,2),default.units = "native")
        grid.lines(c(chrLength,chrLength),c(2,4),default.units = "native")
        # creuament al centromer
        grid.lines(c(center-500000,center+500000),c(4,2),default.units = "native")
        grid.lines(c(center-500000,center+500000),c(2,4),default.units = "native")
        popViewport()
}
```

### Plot chromosome ideograms with additional tracks

Ideograms aka karyograms may contain additional information like p-Values or coverage information.   
https://stackoverflow.com/questions/44003072/annotate-karyogram-with-granges-track/44043471#44043471 points to an option how to annotate karyograms with additional data.


```R
library("ggbio")
library("GenomicRanges")
# needs biovizBase as well

# test data look like this
    CHROM    POS  fisher_het
 1:    10 134775 0.299587633
 2:    10 135237 1.000000000
 3:    10 135277 0.483198279
 4:    10 135331 0.224587437
 5:    10 135334 0.068035761
 6:    10 135656 0.468998144
 7:    10 135708 0.746611845
 8:    10 135801 0.242257762
 9:    10 135853 0.001234701
10:    10 137186 0.774670848

# load banding data
data(hg19IdeogramCyto, package = "biovizBase")
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))

# create a test GRanges object
# from the test data given above
test.granges <- GRanges(seqnames = paste0("chr", df.test.data$CHROM),
                        ranges=IRanges(start = df.test.data$POS,
                                       end = df.test.data$POS),
                        strand = "*",
                        fisher_het = df.test.data$fisher_het)

# attach chromosome lengths
data(hg19Ideogram, package = "biovizBase")
seqlengths(test.granges) <- seqlengths(hg19Ideogram)[names(seqlengths(test.granges))]

```

- the x-coordinate from the GRanges object to be used would be `start`
- `ylim` defines the size of the point-subplots
- `size` needs to be reduced a little as default point size would produce a nearly unreadable plot
- `geom` understands all ggplot2 geom_... , e.g. `line` and `point`

```R
# plot karyogram with p-Values as point ploy on top of chromosomes

ggplot(hg19) +
     layout_karyogram(cytoband = TRUE) +
     layout_karyogram(data = test.granges,
                      geom = "point",
                      aes(x=start, y=fisher_het),
                      ylim = c(10,50),
                      color = "black",
                      size = 0.4
                      )

```



### Gviz plots

The nice [http://bioconductor.org/packages/release/bioc/html/Gviz.html GViz] package for visualization of genomic data. More examples [https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/ here] and [http://www.sthda.com/english/wiki/gviz-visualize-genomic-data here]:

```R
#################################
##
## Gviz visualization of KMT2A coordinates: chr11:118,307,205-118,397,539 (hg19)
##
#################################
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb   <- TxDb.Hsapiens.UCSC.hg19.knownGene
CHR    <- "chr11"
START  <- 118307205
END    <- 118397539
GENOME <- "hg19"
WD     <- "/fsimb/groups/imb-bioinfocf/projects/roukos/imb_roukos_2016_03_roukos_KMT2A"
setwd(WD)

# read tracks in
bw <- list.files(path="./tracks", pattern="\\.bw$", full.names=TRUE)
bw.tracks <- lapply(bw, import.bw, which=GRanges(CHR, IRanges(START, END)))
names(bw.tracks) <- gsub("\\.bw$", "", basename(bw))

# genome browser
tracks <- c(list(IdeogramTrack(genome=GENOME, chromosome=CHR),
                 GenomeAxisTrack(),
                 GeneRegionTrack(txdb, chromosome=CHR, start=START, end=END, showId=TRUE, name="Gene Annotation")),
            mapply(function(track, name) DataTrack(range=track, name=name, chromosome=CHR, start=START, end=END, genome=GENOME),
                   bw.tracks, names(bw.tracks), SIMPLIFY=FALSE))

plotTracks(tracks, type="histogram")
```

### Graphical representation of contingency tables

```R
# fisher test of the DB which have CpG affinity in WT vs. ZBTB48 KD vs. all promoters
m <- matrix(c(length(unique(i_WT[, 1])), length(pos_WT) - length(unique(i_WT[, 1])),
              length(unique(i_KD[, 1])), length(pos_KD) - length(unique(i_KD[, 1]))),
            nrow=2, ncol=2, dimnames=list(c("cpg", "no_cpg"), c("ZBTB48 targets", "ZBTB48 KO")))

x <- chisq.test(m)
main <- if(x$p.value < .01) "p<0.01" else paste0("p=", round(x$p.value))

mosaicplot(t(m), color=TRUE, main="contingency table", sub=main)
assocplot(m, main="deviation from independence of\ncpg islands and ZBTB48 targets")
```

### ChIPseq plots

See the section [Plots](#plots) for specific ChIP-seq plots.

## Miscellaneous bioinformatic related stuff

### Barcode design
Design custom and robust DNA barcode sets capable of correcting substitution errors or insertion, deletion, and substitution errors. Use the [http://bioconductor.org/packages/release/bioc/html/DNABarcodes.html DNABarcodes] package from Bioconductor. More [http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0036852 here].

```R
library("DNABarcodes")
mySet <- create.dnabarcodes(5)
## 1) Creating pool ...  of size 592
## 2) Conway closing...  done
show(mySet)
##  [1] "GAGAA" "AGCAA" "CCTAA" "CAAGA" "ACGGA" "GTCGA" "TGTGA" "GGACA"
##  [9] "CTGCA" "TACCA" "CGAAG" "TCGAG" "GTTAG" "ATAGG" "AAGCG" "GAATG"
## [17] "TGCTG" "ACTTG" "ACAAC" "CACAC" "TAGGC" "CTTGC" "TTACC" "GATCC"
## [25] "AGGTC" "GCCAT" "TCAGT" "AACGT" "TGGCT" "CAGTT"
analyse.barcodes(mySet)
##                   Description  hamming   seqlev levenshtein
## 1               Mean Distance 5.242908 3.656915    4.560284
## 2             Median Distance 5.000000 4.000000    5.000000
## 3            Minimum Distance 3.000000 1.000000    2.000000
## 4            Maximum Distance 7.000000 7.000000    7.000000
## 5 Guaranteed Error Correction 1.000000 0.000000    0.000000
## 6  Guaranteed Error Detection 2.000000 0.000000    1.000000
```

### Getting things from Biomart

Biomart is the data access interface for Ensembl databases. The ''biomaRt'' package from ''bioconductor'' provides an API to access the datamart.

```R
#######################################
##
## Biomart access
##
## listDatasets(useMart="ensembl") --> returns the available organisms
## listFilters(mart) --> returns the filters that can be applied on the dataset
## listAttributes(mart) --> returns a list of available fields on the dataset
##
#######################################
library(biomaRt)

## normal connection (ensembl genes, hg19)
gene.ann <- getBM(attributes=c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position"),
                  filter="biotype", values="protein_coding",
                  mart=useMart("ensembl",dataset="hsapiens_gene_ensembl"))

## connect to HG18
#mart  <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",
#	host="may2009.archive.ensembl.org",
#	path="/biomart/martservice",archive=FALSE)	# Ensembl54 (hg18/NCBI36))

## connect to Variation (SNPs)
# mart <- biomaRt::useMart("snp",dataset="hsapiens_snp")
# snp  <- biomaRt::getBM(attributes=c("refsnp_id","chr_name","chrom_start"),
#	filters="snp_filter",values=snp,mart=mart)
```

### Getting things from the Ensembl databases
Alternatively, one can get the same stuff directly from the Ensembl databases. Also an alternitive from their messy perl API...

```R
library(RMySQL)
con   <- dbConnect(MySQL(),user="anonymous",host="ensembldb.ensembl.org",dbname="homo_sapiens_core_79_38")
query <- paste("SELECT g.stable_id, t.stable_id, p.stable_id, x.display_label",
               "  FROM gene AS g, transcript AS t, translation AS p, xref AS x",
               " WHERE t.transcript_id = g.canonical_transcript_id",
               "   AND t.transcript_id = p.transcript_id",
               "   AND x.xref_id = g.display_xref_id")
x     <- dbGetQuery(conn=con,statement=query)
```

### Getting things from the UCSC tables and ENCODE tracks

Using the '''rtracklayer''' package. Examples adapted from there:

```R
library(rtracklayer)
mySession <- browserSession("UCSC")
genome(mySession) <- "hg19"
```

Discovering which tracks and tables are available from UCSC:

```R
track.names <- trackNames(ucscTableQuery(mySession))
tableNames(ucscTableQuery(mySession,track=track.names[1]))
```

Identify repeat-masked regions in and around the transcription
start site (TSS) of the human E2F3 gene, in hg19:

```R
e2f3.tss.grange <- GRanges("chr6", IRanges(20400587, 20403336))
tbl.rmsk <- getTable(ucscTableQuery(mySession,track="rmsk",range=e2f3.tss.grange,table="rmsk"))
```

Get DNaseI hypersensitivity regions in the K562 Cell Line from the ENCODE project:

```R
track.name <- "wgEncodeUwDgf"
table.name <- "wgEncodeUwDgfK562Hotspots"
e2f3.grange <- GRanges("chr6", IRanges(20400587, 20403336))
tbl.k562.dgf.e2f3 <- getTable(ucscTableQuery(mySession,track=track.name,range=e2f3.grange,table=table.name))
tbl.k562.dgf.hg19 <- getTable(ucscTableQuery(mySession,track=track.name,table=table.name))
```

### Drug response curves using the ''drc'' package

```R
 library(drc)

 # llegir dades
 x <- read.csv("x.csv")
 x$dose <- 2^x$concen

 # fit into a sigmoidal model
 mock <- drm(mock ~ dose,data=x,
 			fct=LL.4(names=c("Slope","Lower Limit","Upper Limit", "ED50")))

 oe   <- drm(oe ~ dose,data=x,
 			fct=LL.4(names=c("Slope","Lower Limit","Upper Limit", "ED50")))

 # report the IC50 with its confidence intervals
 ED(mock,c(5,50,95),interval="delta")
 ED(oe  ,c(5,50,95),interval="delta")

 # plot
 plot(mock,broken=T,type="all",col="red",ylab="viability")
 plot(oe  ,broken=T,type="all",col="blue",add=T)
 legend("center",fill=c("red","blue"),legend=c("mock","oe"))
```

### Permutation test to identify if 2 curves are significantly different

First method is based on anova multiple sample test described in Elso et al, 2004 and implemented in the ''statmod'' package:

```R
 ##
 ## Input data has the folowing structure:
 ##
 ## +---------------+---------+-----+---------+-------+
 ## |  concentracio | x.clon1 | ... | x.clonN | Group |
 ## +---------------+---------+-----+---------+-------+
 ## |           1uM |   0.66  |     |   0.69  |  MOCK |
 ## |           2uM |   0.68  |     |   0.71  |  MOCK |
 ## |  ...                                            |
 ## |           1uM |   0.56  |     |   0.59  | SCRAM |
 ## |           2uM |   0.58  |     |   0.51  | SCRAM |
 ## |  ...                                            |
 ## |           1uM |   0.54  |     |   0.69  | SKMEL |
 ## |           2uM |   0.55  |     |   0.61  | SKMEL |
 ## |  ...                                            |
 ## +---------------+---------+-----+---------+-------+
 ##
 library(statmod)
 library(reshape2)
 library(parallel)

 # permutation test
 tperm <- function(f,nperm=1000) {

 	# read data an normalize every clon with its first concentration value
 	df <- read.csv(f)
 	vals <- grep("^x",colnames(df))
 	df.norm <- as.list(by(df[,vals],df$Group,function(df) { apply(df,2,function(x) x / x[1]) }))
 	A  <- rep(names(df.norm),each=length(vals))
 	dades <- Reduce(cbind,df.norm)

 	# permutation test
 	s <- compareGrowthCurves(A,t(as.matrix(dades)),nsim=nperm)
 	write.csv(s,file=paste0("stats.",f))

 	# plot the dosage curves (be careful, the regression is polinomic (loess) and not logistic)
 	pdf(paste0("stats.",f,".pdf"))
 	df.norm <- lapply(df.norm,as.data.frame)
 	for(i in 1:length(df.norm)) {
 		df.norm[[i]]$uM <- log2(unique(df$uM))
 		df.norm[[i]]$Group <- names(df.norm)[i]
 	}
 	dades <- Reduce(rbind,df.norm)
 	df2 <- melt(dades,id.vars=c("uM","Group"))
 	p <- ggplot(df2,aes(x=uM,y=value,colour=Group)) + geom_smooth() + geom_point()
 	print(p)
 	dev.off()
 }

 mclapply(list.files(pattern="*.csv"),tperm,nperm=1000,mc.cores=4)
```

Second method is based on the F-score from an anova test, were we compare the original score to those obtained by permutation, using an empirical cumulative distribution to get the pvalue:

```R
 tperm <- function(f,nperm) {
 	dades <- read.csv(paste(f,".csv",sep=""))
 	dades <- reshape2::melt(dades,id.vars=c("Concentration","Group"))
 	dades$Concentration <- as.factor(dades$Concentration)
 	dades$Group <- as.factor(dades$Group)
 	d <- anova(lm(value ~ Group*Concentration,data=dades))["Group:Concentration","F value"]

 	# get nperm scores (permutation test)
 	require(parallel)

 	d.null <- mclapply(1:nperm,function(x) {
 		# permute (only permute the values inside every concentration)
 		perm <- by(dades,dades[,"Concentration"],function(x) {
 						x$value <- x$value[sample(1:nrow(x))];
 						return(x) })
 		perm <- do.call(rbind,perm)
 		d <- anova(lm(value ~ Group*Concentration,data=perm))["Group:Concentration","F value"]
 	},mc.cores=4)

 	# build an empirical distribution, and output the pvalue to get more extreme D
 	p.null <- ecdf(unlist(d.null))
 	print(paste(f,"pval:",1 - p.null(d)))
 	return(p.null)
 }

 s <- sapply(c("akt","sft31","shikonin"),tperm,10000)
```

### GenomicRanges

Example using the ''GenomicRanges'' package from ''bioconductor'' to find overlapping features:

```R
## more things can be done (GAP, union, intersect, setdiff...)
library(GenomicRanges)

# read the two tables containing positions
dmr  <- read.csv("~/Desktop/DMR.csv")
gens <- read.csv("~/Desktop/gens.csv")

# build the GRanges objects (the columns, in this order: CHR,INI,FI,STRAND,FEATURES...)
rd1 <- with(dmr, GRanges(CHR,IRanges(INI,FI),STRAND,nom=NOM))
rd2 <- with(gens,GRanges(CHR,IRanges(INI,FI),STRAND,nom=NOM,pval=PVAL,diff=DIFF))

# do the matching; maxgap distance accepted
mm <- findOverlaps(rd1,rd2,maxgap=2000)

# add the columns with the feature names from the query
solucio <- data.frame(as.data.frame(rd1[as.matrix(mm)[,1],]),as.data.frame(rd2[as.matrix(mm)[,2],]))
write.csv(solucio,file="solucio.csv")
```

### Converting GTF to GFF3

Use the ''rtracklayer'' package from ''Bioconductor'' to convert between GTF <--> GFF formats:

```R
library(rtracklayer)
GTF <- "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/rnaseq/ref-chr1/genes+biotypes-chr1.gtf"
GFF <- "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/rnaseq/ref-chr1/genes+biotypes-chr1.gff"
export.gff(import.gff(GTF,format="gtf"),GFF,format="gff3")
```

### Flatten a GFF/GTF file by gene_id (and get transcript lengths)

Taking a standard GFF/GTF and splitting exons by gene, and then flattening them in order to merge overlapping exons and get rid of duplicated coordinates:

```R
# grab the genes from the igenomes annotation
library(GenomicRanges)
library(rtracklayer)
gtf <- import.gff("./test/genes.gtf", format="gtf", feature.type="exon")
gtf.flat <- unlist(reduce(split(gtf, elementMetadata(gtf)$gene_id)))
gene.lengths <- tapply(width(gtf.flat), names(gtf.flat), sum)
```

Tag a data.frame with the just produced transcript lengths. (watch out to use the names of the length matrix as the last parameter in match!)

```R
DE.data$transcript_length <- gene.lengths[match(DE.data$gene_id, names(gene.lengths))]

```

### Download a dataset from GEO using GEOquery

Get a dataset from GEO directly quering from R and the ''GEOquery'' package from ''bioconductor'':

```R
 library(Biobase)
 library(GEOquery)

 ##
 ## Retina i Retina Detachment
 ##
 # expression arrays
 gset <- getGEO("GSE28133",GSEMatrix=T,AnnotGPL=T)	# be careful if the dataset uses more than 1 platform or is divided in several parts
 ex   <- exprs(gset[[1]])
 colnames(ex) <- c(rep("Ret",19),rep("RD",19))

 # platform annotation
 gpl    <- annotation(gset[[1]])	# GPL platform name
 platf  <- getGEO(gpl,AnnotGPL=TRUE)
 ncbifd <- data.frame(attr(dataTable(platf),"table"))

 # average expression of the several probes that interrogate a gene
 ex <- ex[rownames(ex) %in% ncbifd[ncbifd$Gene.symbol != "","ID"],]
 ex <- apply(ex,2,function(x,gens) { tapply(x,gens,mean,na.rm=T) },
 			sapply(rownames(ex),function(x) ncbifd[ncbifd$ID == x,"Gene.symbol"]))

 # save tables
 write.csv(ex[,grep("Ret",colnames(ex))],file="Ret.csv")
 write.csv(ex[,grep("RD",colnames(ex))],file="RD.csv")
```

### Get the genomic sequence of a region and plot its nucleotide content

This code uses the ''BSgenome'' package from ''bioconductor'' and calculates de nucleotide content in a ''w'' bp window with an ''w/5'' stepsize:

```R
 # get the sequence
 library(BSgenome)
 library(BSgenome.Hsapiens.UCSC.hg19)
 library(ggplot2)

 ini <- 66217200
 fi  <- 66221000
 w   <- 50
 st  <- w/5	# step window
 gen <- "HMGA2"
 s <- getSeq(Hsapiens,"chr12",ini,fi,as.character=TRUE)

 # plot the A/T and C/G content (as the average of w bp)
 pdf(paste(gen,".pdf",sep=""),width=12)

 f <- function(i,x,s,w) { length(gregexpr(x,substr(s,i-w,i+w))[[1]]) / (w * 2) }
 A=sapply(seq(from=w,to=nchar(s) - w,by=st),f,"A",s,w)
 T=sapply(seq(from=w,to=nchar(s) - w,by=st),f,"T",s,w)
 C=sapply(seq(from=w,to=nchar(s) - w,by=st),f,"C",s,w)
 G=sapply(seq(from=w,to=nchar(s) - w,by=st),f,"G",s,w)

 # A/T content
 df <- data.frame(x=rep(1:length(A),2),y=c(A,T),class=c(rep("A",length(A)),rep("T",length(T))))
 p <- ggplot(df,aes(x=x,y=y,colour=class)) + geom_line(size=1.05) + ggtitle(gen) +
 #	 scale_x_discrete(breaks=df$x,labels=as.character(ini + df$x * w)) +
 	 theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))
 print(p)

 # C/G content
 df <- data.frame(x=rep(1:length(C),2),y=c(C,G),class=c(rep("C",length(C)),rep("G",length(G))))
 p <- ggplot(df,aes(x=x,y=y,colour=class)) + geom_line(size=1.05) + ggtitle(gen) +
 #	 scale_x_discrete(breaks=df$x,labels=as.character(ini + df$x * w)) +
 	 theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))
 print(p)

 dev.off()
```

### Gene set enrichment analysis

#### GSEA of Biological Processes

This code calculates for enrichment in biological processes from Gene Ontology, using the ''org.Hs.eg.db'' and ''GOstats'' packages from ''bioconductor'':

```R
library(org.Hs.eg.db)
library(GOstats)

# list of official gene names
geneids <- read.table("list.csv",head=F,row.names=NULL,stringsAsFactors=F)[,1]
geneids <- geneids[!is.na(geneids)]
genes   <- unlist(mget(geneids,ifnotfound=NA,revmap(org.Hs.egSYMBOL)))
#be careful here if one Symbol which is returned is associated to more than one entrez id
# there will be entrez ids added to the list of genes and the names of the genes vector are changed
# by adding 1,2 to the name! Check for your symbol and correct. An easy solution for ensembl ids (with length 11) is
#genes <-genes[!duplicated(substr(names(genes), 1, 11))]
#names(genes) <- substr(names(genes), 1, 11)
univ    <- Lkeys(org.Hs.egGO)
param2  <- new("GOHyperGParams",geneIds=genes,universeGeneIds=univ,annotation="org.Hs.eg.db",ontology="BP")
hyp     <- hyperGTest(param2)

# take the categories with +2 genes
result<-summary(hyp,categorySize=2)
result<-data.frame(result,FDR=p.adjust(result$Pvalue,"fdr"))

# complement the GO term with the genes from the list it contains
result$genes <- sapply(as.list(geneIdUniverse(hyp)[result$GOBPID]),function(x) {
   paste(names(genes)[genes %in% x],collapse="; ")
})
```

Also, gene SYMBOLS can be retrieved from org.Hs.eg.db with the entrez gene Ids:
```R
result$genes <- sapply(as.list(geneIdUniverse(hyp)[result$GOBPID]),function(x) {
   # get the entrezid of the genes in the GO category which are also in our list
   gene_symbols <- select(org.Hs.eg.db, columns="SYMBOL", keys=x[x %in% entrezid], keytype="ENTREZID")$SYMBOL
   paste(sort(unique(gene_symbols)), collapse="; ") # gene symbol
})
```

Alternatively genes names can be obtained from the SQL database (in case the previous step doesn't work):
```R
library(DBI)
library(RSQLite)

drv <- dbDriver("SQLite")
con <- dbConnect(drv,dbname="/usr/local/lib/R/site-library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite")

 genes <- apply(result,1,function(x) {
		res <- dbGetQuery(con,paste("select distinct a.symbol",
 					    "  from gene_info a,go_bp_all b",
 					    " where a._id = b._id",
 					    "   and go_id = ?"),
 				  data.frame(go_id=x["GOBPID"]))
 		res <- intersect(unique(geneids),res[,1])
 		return(paste(res,collapse=";"))
 })

 # save the results
 write.csv(data.frame(result,genes=genes),file="llista_verda.GO.csv")
```

#### GSEA of Biological Processes (a parallel version)

A parallel version based on "parallel" (not "multicore", "snow", "foreach") can be prepared to speed up when multiple GO tests have to be performed. Notice the use of the package "parallel", which is the only one allowing parallel access to SQLite databases.

```R
# functions for parallel GO testing
makeGoParms <- function(geneids) {
    genes  <- unlist(mget(geneids,ifnotfound=NA,revmap(org.Mm.egSYMBOL)))
    #be careful here if one Symbol which is returned is associated to more than one entrez id
    # there will be entrez ids added to the list of genes and the names of the genes vector are changed
    # by adding 1,2 to the name! Check for your symbol and correct. An easy solution for ensembl ids (with length 11) is
    #genes <-genes[!duplicated(substr(names(genes), 1, 11))]
    #names(genes) <- substr(names(genes), 1, 11)
    univ   <- Lkeys(org.Mm.egGO)
    param2 <- new("GOHyperGParams",geneIds=genes,universeGeneIds=univ,annotation="org.Mm.eg.db",ontology="BP")
    list(genes=genes, param2=param2)
}

GoTest <- function(parms) {
    hyp <- hyperGTest(parms$param2)
    result <- summary(hyp,categorySize=2)
    result <- data.frame(result,FDR=p.adjust(result$Pvalue,"fdr"))
    result$genes <- sapply(as.list(geneIdUniverse(hyp)[result$GOBPID]),function(x) {
        paste(names(parms$genes)[parms$genes %in% x],collapse="; ")
    })
    result
}

# calculate overlaps as in venn
ABC <- Reduce(intersect, x)
AB  <- setdiff(intersect(x[[1]], x[[2]]), ABC)
AC  <- setdiff(intersect(x[[1]], x[[3]]), ABC)
BC  <- setdiff(intersect(x[[2]], x[[3]]), ABC)
A   <- Reduce(setdiff, list(x[[1]], AB, AC, ABC))
B   <- Reduce(setdiff, list(x[[2]], AB, BC, ABC))
C   <- Reduce(setdiff, list(x[[3]], AC, BC, ABC))

# make cluster
cl <- makeCluster(CORES)
x  <- clusterEvalQ(cl, library(GOstats))
x  <- clusterEvalQ(cl, library(org.Mm.eg.db))

# calculate GO
res <- parLapply(cl, lapply(list(A, B, C, AB, AC, BC, ABC), makeGoParms), GoTest)
stopCluster(cl)
names(res) <- c("E145 only", "ISC only", "AE only", "E145+ISC", "E145+AE", "ISC+AE", "E145+ISC+AE")
WriteXLS("res", ExcelFileName="./results/timecourse_rna.GO.xls")
```

#### GSEA of KEGG or custom db

```R
library(org.Dm.eg.db)
library(GOstats)
library(GSEABase)

# create gene set collection from KEGG.db
frame <- toTable(org.Dm.egPATH)
keggframeData <- data.frame(frame$path_id, frame$gene_id)
keggFrame <- KEGGFrame(keggframeData, organism="Drosophila melanogaster")
gsc  <- GeneSetCollection(keggFrame, setType=KEGGCollection())

# generate parms and call to GOstats

g  <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id"), filters="flybase_transcript_id",
            values=unique(unlist(strsplit(x, ";"))), mart=mart)
xx <- unique(unlist(mget(g$ensembl_gene_id, ifnotfound=NA, revmap(org.Dm.egENSEMBL))))
#now we have to get rid of the genes which have multiple entrez ids and just take one
# since mget is adding 0,1,2 to the name we can filter by the ensembl id length (11)
xx <-xx[!duplicated(substr(names(xx), 1, 11))]
names(xx) <- substr(names(xx), 1, 11)
univ <- Lkeys(org.Dm.egGO)
param <- GSEAKEGGHyperGParams(name="KEGG", geneSetCollection=gsc, geneIds=xx,
                              universeGeneIds=univ, pvalueCutoff=0.05, testDirection="over")
hyp    <- hyperGTest(param$param)
result <- summary(hyp, categorySize=catsize)
result <- data.frame(result, FDR=p.adjust(result$Pvalue, "fdr"))
```

#### GSEA of a custom GO (slim) db

```R
library(GOstats)
library(GSEABase)
goslim <- read.delim("./analysis/db/go_slim_mapping.tab", comment="#")
goslim$evidence <- "EXP"    # fake experimental evidence code
goframe <- GOAllFrame(GOFrame(goslim[, c("GOid", "evidence", "gene")], organism="Saccharomyces cerevisiae"))
gsc  <- GeneSetCollection(goframe, setType=GOCollection())

# and generate parms and call to GOstats
# ...
```

#### GO functional analysis with clusterProfiler

Use visualization functions provided by ''clusterProfiler'' to plot GO enrichment analysis calculated through DAVID. ''clusterProfiler'' also allows calculating the GO analysis through up-to-date R packages instead of using the old and unmantained DAVID.

```R
require(DOSE)
require(clusterProfiler)
require(biomaRt)

##
## GO ENRICHMENT + PLOTS
##
mart <- useDataset("mmusculus_gene_ensembl",useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org"))
GOenrich <- function(g) {
     g <- getBM(attributes="ensembl_gene_id",filters="external_gene_name",values=g,mart=mart)
     x <- enrichDAVID(g$ensembl_gene_id,idType="ENSEMBL_GENE_ID",listType="Gene",
                      annotation="GOTERM_BP_ALL",pvalueCutoff=0.01)
     x
}
david1 <- GOenrich(genes1)	# official gene symbols
david2 <- GOenrich(genes2)	# official gene symbols

# some plots
barplot(david1)
cnetplot(david1)

##
## Compare the GO enrichment of the 2 gene lists
##
GOcompare <- function(g) {
    g <- lapply(g,function(g) {
        cat("connecting to biomart...\n")
        x <- getBM(attributes="entrezgene",filters="external_gene_name",values=g,mart=mart)
        x$entrezgene
    })
    cat("connecting to david...\n")
#   x <- compareCluster(g,fun="enrichDAVID",idType="ENTREZ_GENE_ID",listType="Gene",
#                       annotation="GOTERM_BP_ALL",pvalueCutoff=0.01,minGSSize=5)
    x <- compareCluster(g,fun="enrichGO",organism="mouse",ont="BP")
    x
}

plot(GOcompare(list(genes1,genes2)),showCategory=25)
```

#### Reducing GO DAGs with Semantic Similarity

* '''Problem''': Many redundant terms in GSEA
* '''Solution''':
  * Use slimmed down ontologies --> not available for most organisms
  * Several methods proposed based on Information Content criteria
  * Or graph based methods for quantitative comparison between nodes of a DAG
  * A Bioconductor package implements all these methods described before [http://bioconductor.org/packages/release/bioc/html/GOSemSim.html GOSemSim]

How large would you like the resulting list to be?
* Large (allowed similarity=0.9)
* Medium (0.7)
* Small (0.5)
* Tiny (0.4)

```R
# goterms is the results of gostats as a dataframe or something similar
# (it will work as long as the Pvalue and GOBPID columns are present)
goterms <- read.csv("goterms.csv")

# load the parallel backend
library(parallel)
library(Rcpp)
library(GOSemSim)

MAXTERMS <- 64  # reduce the top MAXTERMS significant terms
CORES <- 16
cl    <- makeCluster(CORES)
clusterEvalQ(cl, library(GOSemSim))

# calculate the semantic similarity with previous terms (from bottom (less signif) to top)
goterms <- goterms[order(goterms$Pvalue),]
goterms <- goterms[1:min(c(MAXTERMS, nrow(goterms))),]
drop <- parLapply(cl, nrow(goterms):2, function(i, goterms) {
    semsim <- lapply((i-1):1, function(j) {
        goSim(goterms$GOBPID[i], goterms$GOBPID[j], semData=godata("org.Hs.eg.db", ont="BP"), measure="Rel")
    })
    # drop if there was a similar term with lower FDR
    any(semsim > .7, na.rm=T)
}, goterms)

drop <- c(FALSE, rev(unlist(drop)))
goterms.reduced <- goterms[!drop,]  # take just the columns you want...

stopCluster(cl)
```

Since it's an extremely CPU intensive and parallelizable problem, it's advisable to use parLapply instead of the first lapply loop and set a maximum of the number of terms to be processed (by min(c(MAXTERMS,nrow(cl):2)) in the lapply loop).

#### Summarizing by scoring frequencies of words

An original way to visually summarize GO terms by displaying variable sized words depending on their frequency:

```R
library(wordcloud)
x <- read.csv("categorical_onlyE_apclusters_GO.csv")
x[1:5,]
#      GOBPID       Pvalue OddsRatio   ExpCount Count Size                                                                 Term
#1 GO:0001112 5.290072e-06  143.6133 0.03804965     3    7                 DNA-templated transcriptional open complex formation
#2 GO:0001113 5.290072e-06  143.6133 0.03804965     3    7 transcriptional open complex formation at RNA polymerase II promoter
#3 GO:0001120 5.290072e-06  143.6133 0.03804965     3    7                                       protein-DNA complex remodeling
#4 GO:0034367 5.290072e-06  143.6133 0.03804965     3    7                                    macromolecular complex remodeling
#5 GO:0001109 8.431172e-06  114.8812 0.04348532     3    8                promoter clearance during DNA-templated transcription
wordcloud(x$Term, colors=brewer.pal(8, "Dark2"))
```

Which is just a shortcut for the text mining analysis:

```R
crude <- Corpus(VectorSource(x$Term))
crude <- tm_map(crude, PlainTextDocument)
crude <- tm_map(crude, stripWhitespace)
crude <- tm_map(crude, removePunctuation)
crude <- tm_map(crude, function(x)removeWords(x, stopwords()))
tdm <- TermDocumentMatrix(crude)
m <- as.matrix(tdm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
wordcloud(d$word, d$freq, min.freq=2, colors=brewer.pal(8, "Dark2"))
```

#### Summarizing by clustering similarity scores between terms

First one needs to calculate the similarity scores between terms, and cluster them based on the distance.

```R
# recalculate the semantic similarity scores of the reduced terms for the core proteome
library(parallel)  # create your cluster of workers
library(GOSemSim)

scores <- {
    x <- goterms$GOBPID[order(goterms$Pvalue)][1:min(c(MAXTERMS, length(goterms$GOBPID)))]
    parSapply(cl, 1:length(x), function(i, x) {
        c(rep(NA, i-1),
            sapply(i:length(x), function(j, i) {
                goSim(x[i], x[j], ont="BP", measure="Wang", organism="fly")
            }, i)
        )
    }, x)
}
rownames(scores) <- colnames(scores) <- goterms$GOBPID

# cluster based on distance between terms, and calculate the cluster representative as the broadest term
h <- cutree(hclust(as.dist(1-scores)), h=.7)            # calculate distances, cluster and cut into subclusters
goterms$cluster <- h[match(goterms$GOBPID, names(h))]
clrep <- c(by(goterms, goterms$cluster, function(x) {   # calculate the representative
    x$Term[which.max(x$Size)]
}))
goterms$clrep <- clrep[goterms$cluster]
```

##### Treemap representation of the hierarchy of terms

Split a rectangle into multiple rectangles, which area is proportional to a value (odds ratio, pvalue, genes in the term...), and colour according to a hierarchy of the data.

```R
library(treemap)

# plot treemap, the color i based on the term and the size on the oddsRatio
pal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pal2<- colorRampPalette(pal)(length(unique(h)))
treemap(goterms, index=c("clrep", "Term"), vSize="Size", type="index", title="", fontsize.labels=c(14, 11),
        overlap.labels=1, fontcolor.labels=c("#FFFFFFDD", "#00000080"), bg.labels=0, border.col="#00000080",
        palette=pal2)
```

##### 2D multidimensional scaling of the semantic similarity scores

Using the scores previously calculated, do a multidimensional scaling of the distance matrix and plot the hierarchy of clusters as a scatter plot.

```R
library(ggplot2)

fit <- cmdscale(as.matrix(as.dist(1-scores)), eig=TRUE, k=2)

df  <- as.data.frame(fit$points)
df$term  <- goterms$Term[match(rownames(df), goterms$GOBPID)]
df$size  <- goterms$Size[match(rownames(df), goterms$GOBPID)]
df$clrep <- goterms$clrep[match(rownames(df), goterms$GOBPID)]
df2 <- df[df$clrep == df$term,]

p <- ggplot(df, aes(x=V1, y=V2, label=clrep)) +
        geom_point(aes(color=clrep, size=size), alpha=.5) +
        geom_label_repel(data=df2, aes(x=V1, y=V2, label=clrep, color=clrep), box.padding=unit(2,"lines")) +
        scale_color_manual("cluster", values=pal2, guide=F) +
        scale_size_continuous(guide=F, range=c(0,25)) +
        scale_x_continuous(name="") +
        scale_y_continuous(name="") +
        theme_minimal() +
        theme(axis.text.x=element_blank(), axis.text.y=element_blank())
print(p)
```

### Precomputing tables of semantic similarities between GO terms

Pre-render the table of all GO terms by semantic similarity. This is a one-time and extremely time consuming task that is worth having it pre-calculated before reducing any list of GO terms.

```R
library(parallel)
library(org.Dm.eg.db)
set.seed(666)
setwd("/fsimb/groups/imb-buttergr/=EVOREG/drosophila_Development/tmp")
CORES=48
LOG="semsim.Dm.log"

# make cluster
cl <- makeCluster(CORES)
x  <- clusterEvalQ(cl, library(GOSemSim))

# get all GO terms from the annotation and calculate the Wang similarity scores
goterms <- Rkeys(org.Dm.egGO)
cat(length(goterms), file=LOG, append=F, fill=T)
x <- clusterExport(cl, "goterms")
reduced <- sapply(1:length(goterms), function(i) {
    cat(".", file=LOG, append=T)
    c(rep(NA, i-1),
        parSapply(cl, i:length(goterms), function(j, i) {
            goSim(goterms[i], goterms[j], ont="BP", measure="Wang", organism="fly")
        }, i)
    )
})

# save result
colnames(reduced) <- rownames(reduced) <- goterms
write.csv(reduced, file="semsim.Dm.csv")

stopCluster(cl)
```

### Basic matrix normalization

#### Quantile normalization

* sort the values per sample (column)
* calculate the median per row
* rank the values
* linear interpolation using a linear model y ~ x

```R
 f.qn <- function(x) {
 	xm <- apply(x,2,sort)
 	xm <- apply(xm,1,median,na.rm=T)
 	xr <- c(apply(x,2,rank))
 	return(array(approx(1:nrow(x),xm,xr)$y,dim(x),dimnames(x)))
 }
```

#### Smoothing curves normalization using loess

* calculate the median per row
* for every sample, fit a polinomic model (loess) with the real values (y) and the medians (x)
* substract from the real value the distance between the fitted value and the median

```R
library(parallel)
cl <- makeCluster(CORES)
f.smooth <- function(x) {
    xm <- apply(x,1,median,na.rm=T)
    parApply(cl,x,2,function(x,xm) {
        smooth.fit <- fitted.values(loess(x ~ xm))
        dev <- smooth.fit - xm
        x - dev
    },xm)
}
stopCluster(cl)
```

### Calculate the reverse-complimentari of a sequence

```R
RC <- function(s) {
	chartr("ACTG","TGAC",paste(rev(unlist(strsplit(s,""))),collapse=""))
}
```

### Read a fasta file

Can be done with the ShortRead package:

```R
library("ShortRead")
fasta <- readFasta("file.fa")
for(i in 1:length(fasta)) {
    se <- as.character(sread(fasta[i]))	# sequence
    id <- as.character(id(fasta[i]))	# id
}
```

Or with the Biostrings package:

```R
library("Biostrings")
fasta <- readAAStringSet("file.fa")     # for peptides. Try also readDNA or readRNA
for(i in 1:length(fasta)) {
    se <- as.character(fasta[i])	# sequence
    id <- names(as.character(fasta[i]))	# id
}
```

### Get transcripts from a Bioconductor's AnnotationDb

```R
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
columns(txdb)
x <- select(txdb, columns=c("TXNAME","CDSCHROM","TXSTART","TXEND","TXSTRAND"),
            keys=keys(txdb,"TXID"), keytype=c("TXID"))
```

### Get information of genomic features

Taken with a grain of salt, these numbers mean the total number of non-redundant bases the features have. For introns its specially sensitive, since some isoforms may have exons on them (which are not subtracted). Thus, at least for introns, may be better to use transcripts minus exons. '''It's still all a rough estimate'''

```R
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

# generic function to flatten and report total length
flatten <- function(x) {
    x.flat <- reduce(unlist(x))
    sum(width(x.flat))
}

# sum bases per feature
df <- data.frame(bp=c(flatten(promoters(txdb, upstream=1000, downstream=1000)),
                      flatten(transcriptsBy(txdb, by="gene")),
                      flatten(exonsBy(txdb, by="tx")),
                      flatten(cdsBy(txdb, by="tx")),
                      flatten(intronsByTranscript(txdb)),
                      flatten(threeUTRsByTranscript(txdb)),
                      flatten(fiveUTRsByTranscript(txdb)),
                      sum(as.numeric(seqlengths(txdb)))))
df$perc <- round(100 * df$bp / df$bp[8], 2) # [8] is the total genome size
rownames(df) <- c("promoters", "transcripts", "exons", "cds", "introns", "threeUTRs", "fiveUTRs", "genome")
```

## Statistical analysis

### Principal Component Analysis (PCA)

```R
# Verify variance is uniform
plot(apply(x,1,var))

# and scale otherwise
x <- data.frame(scale(x))

#calculate the PCA
pc <- princomp(x)
plot(pc,type='l') # use elbow rule to select how many components

# See which components dominate. What are the loadings?
summary(pc)  # components summary
loadings(pc) # loadings of original variables onto components

# Get principal component vectors using prcomp instead of princomp
# samples as rows & measurements/variables in columns
pc <- prcomp(x)

# Plot first 2 PC
plot(data.frame(pc$x[,1:2]),pch=16,col=rgb(0,0,0,0.5))
```

### Multidimensional Scaling (MDS)

```R
 # Classical MDS
 # N rows (objects) x p columns (variables)
 # each row identified by a unique row name

 d   <- dist(mydata) # euclidean distances between the rows
 fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
 fit # view results

 # plot solution
 x <- fit$points[,1]
 y <- fit$points[,2]
 plot(x,y,xlab="Coordinate 1",ylab="Coordinate 2",main="Metric MDS",type="n")
 text(x,y,labels=row.names(mydata),cex=.7)
```

### Affinity Propagation Clustering

Autoclustering method implemented as described in the paper by Brendan J. Frey and Delbert Dueck (2007). Clustering by passing messages between data points. Science 315:972-977. DOI: 10.1126/science.1136800. The method automatically determines the best number of subgroups in a group, and optimizes the distances of its members from the centroid.

```R
require(apcluster)

 ##
 ## k trees defined by affinity propagation, similarites as negative squared Euclidean distances
 ##
 f <- function(s,type,x) {
         h <- apcluster(s)       # determine number of clusters
         a <- aggExCluster(s,h)  # merge clusters

         plot(h,s,main=paste("similarities",type,"matrix"))      # plot similarities heatmap labeling the found clusters
         plot(a,main=paste("similarities",type,"matrix"))        # plot the cluster
         # plot the first 2 PCA components and the cluster grouping
         for(i in length(h@clusters):1) {
                 plot(a,x,k=i,main=paste("similarities",type,",",i,"clusters")) }

         # write csv with cluster samples
         r <- data.frame()
         for(i in 1:length(h@clusters)) {
                 r <- rbind(r,data.frame(sample=names(h@clusters[[i]]),cluster=i)) }
         write.csv(r,file=paste("lung.AP.",type,".csv",sep=""))
 }

 # for graphical representation, we calculate a PCA to get x,y coordinate based on the first 2 components
 pca <- prcomp(t(betas))
 x   <- t(rbind(pca$x[,1],pca$x[,2]))

 # distnaces can be calculated as correlations, euclidean or manhattan
 f((as.matrix(as.dist(cor(betas)))*100)^2,"correlation",x)
 #f(negDistMat(t(betas),r=2),"euclidean",x)
 #f(negDistMat(t(betas),r=2,method="manhattan"),"manhattan",x)
```

### Other methods for clustering

'''k-means clustering'''
```R
NCLUST=9
k <- lapply(1:NCLUST, function(i) {                                                                                       
    kmeans(models$lfq, centers=i, iter.max=100, nstart=25)                                                                
})
```

'''Mfuzz clustering'''
```R
NCLUST=9
exprs <- new("ExpressionSet", exprs=models$lfq)                                
k <- lapply(2:NCLUST, function(i) {                                          
   mfuzz(exprs, c=i, m=1.25)                                                  
})
```

### Validate the optimal number of clusters
```R
dissE <- daisy(models$lfq)                                                                                                
lapply(k, function(k) {                                                                                                   
    sk <- silhouette(k$cl, dissE)                                                                                         
    plot(sk)                                                                                                              
    abline(v=summary(sk)$avg.width, col="red", lty=2)                                                                     
    cat(i, "-->", sum(summary(sk)$clus.avg.widths > summary(sk)$avg.width), fill=T)                                       

    x <- prcomp(models$lfq)                                                                                               
    plot(x$x[,1], x$x[,2], type="n")                                                                                      
    points(x$x[,1], x$x[,2], col=paste0(brewer.pal(12, "Set1")[as.factor(k$cl)], "80"), pch=20, cex=.5)                   
}
```

### Evaluate cluster strength

Evaluate cluster strength by calculating p-values for hierarchical clustering via multiscale bootstrap resampling:

```R
 require("pvclust")
 require("snow")

 load("data.RData")

 cl <- makeCluster(12,type="SOCK")

 result <- parPvclust(cl,data,nboot=10000) #method.dist="euclidean",method.hclust="complete",nboot=1000,init.rand=F)
 stopCluster(cl)
```

### Calculate ROC curves from a predictor using the ''ROCR'' package

```R
 library(ROCR)

 x <- read.csv("weka.out.csv")

 # R> head(x)
 #   inst. actual predicted
 # 1     1      1         1
 # 2     2      1         1
 # 3     3      0         0
 # 4     4      1         0
 # 5     5      1         1
 # 6     6      0         0

 pred <- prediction(x$predicted,x$actual)
 perf <- performance(pred,"tpr","fpr")
 plot(perf,colorize=T,main="ROC plot")
```

### Impute missing values

Impute missing values in a table using the method of the K-Neares Neighbours and the ''imputation'' package

```R
 library(imputation)

 k <- cv.kNNImpute(data)$k       #cross validate by artificially erasing data and get the number of neighbours
 impdata <- kNNImpute(data,k)$x  #run with the previously obtained number of neighbours
```

### Power and sample size calculation for survival analysis
The calculations follow the method described in Rosner B. (2006), which was proposed on Freedman, L.S. (1982), implemented in the R package powerSurvEpi:

```R
library(powerSurvEpi)
library(survival)

## nomes del training (29 mostres HR=.03, mirar la KM)
dat.clin <- read.csv("../RSF.nosurgery/analisi.rsf.v3.MI.csv")
dat.clin$days  <- round(dat.clin$days / 30)
dat.clin$days <- ifelse(dat.clin$days > 60,60,dat.clin$days)
dat.clin$grup <- ifelse(dat.clin$grup == "baix","C","E")
ssizeCT(formula = Surv(days,status) ~ grup, dat = dat.clin,
        power = 0.8,k = 0.35, RR = 0.03, alpha = 0.05)
powerCT(formula = Surv(days,status) ~ grup, dat = dat.clin,
        nE = 7, nC = 22, RR = 0.03, alpha = 0.05)
```

### Model normally distributed data

Most of the thoughts, text and formulas are extracted from [http://www.statmethods.net/stats/power.html statmethods].

The following four quantities have an intimate relationship:
* sample size
* effect size = total part of the observed between group variance which can be explained as an effect of the group
* significance level = P(Type I error) = probability of finding an effect that is not there
* power = 1 - P(Type II error) = probability of finding an effect that is there
Given any three, we can determine the fourth.

Specifying an [http://en.wikipedia.org/wiki/Effect_size effect size] can be a daunting task. ES formulas and Cohen's suggestions (based on social science research) are provided below. Cohen's suggestions should only be seen as very rough guidelines. Your own subject matter experience should be brought to bear.

```R
models <- data.frame(pval=apply(lfq, 1, function(x) anova(lm(x ~ group))["group", "Pr(>F)"]))
models$fdr  <- p.adjust(models$pval, method="fdr")
models$es   <- apply(lfq, 1, function(x) {
    sqrt(sum(tapply(x, group, function(xi) (length(xi) / length(x)) * (mean(xi) - mean(x))^2)) / sd(x)^2)
})
```

For t-tests, the effect size is assessed as d = |u_1 - u_2| / sd, where u_1=mean group 1, u_2 mean group 2, sd square root of the variance.

For ANOVA, the effect size is defined as f = sqrt( sum(p_i * (u_i - u)^2) / sd^2), where p_i is the percentage of samples in group i.

Cohen suggests that d values of 0.2, 0.5, and 0.8 represent small, medium, and large effect sizes respectively.

### Calculate the center of a 2d-distribution

```R
library(KernSmooth)

optimalBandwidth <- function(x) { # return a more or less useful bandwith (from densCols)
    bandwidth <- diff(apply(x, 2, quantile, probs=c(0.05,0.95), na.rm=T, names=F))/25
    bandwidth[bandwidth == 0] <- 1
    bandwidth
}

xy <- xy.coords(x,y)
select <- is.finite(xy$x) & is.finite(xy$y)
xy <- cbind(xy$x, xy$y)[select,]
dens <- bkde2D(xy,bandwidth=optimalBandwidth(xy))
maxpos <- which(dens$fhat == max(dens$fhat), arr.ind=TRUE)
points(dens$x1[maxpos[1]], dens$x2[maxpos[2]],cex=3,pch="+")
```

### Test the significance of the overlap between 2 lists

E.g. testing for overrepresentation of differentially expressed genes (RNA-seq) in a set of targets from a ChIP-seq experiment.

Use the hypergeometric test.

The background is usually the number of genes which are expressed, or all genes. We're calculating the probability to get a bigger overlap between the two lists (ChIP targets and DE genes) with as many trials as DE genes we have:

```R
bg <- nrow(rnaseq)
A  <- sum(rnaseq$fdr < 0.01)
B  <- sum(peaks$fdr < 0.01)
AB <- length(intersect(rnaseq$gene[rnaseq$fdr < 0.01], peaks$gene[peaks$fdr < 0.01]))

pval <- sum(dhyper(AB:A, B, bg - B, A))
```

### Network analysis

Calculate the similarity between PPI profiles of a set of proteins, and build a network. The assumtion underneath is that proteins which interact with a similar set of proteins are likely to have a similar function. Thus, in the network, whese proteins (nodes) will be connected by an edge.

The edge can be directed (an arrow, showing directionality of the interaction) and weighted (showing the strength of the interaction (correlation, in this case). In order to keep the example simple, the graph in not "directed" nor the edges "weighted".

```R
# calculate similarities between interaction profiles in our RBP
pcc <- cor(t(cellmap[rownames(cellmap) %in% genes, ]), use="pairwise.complete.obs")

# build the adjancency matrix
pcc.binary <- ifelse(abs(pcc) > CORRELATION_THRESHOLD, 1, 0)  # only link nodes with an score above the threshold
out <- apply(pcc.binary, 1, function(x) all(x == 0))
pcc.binary <- pcc.binary[!out, !out]
colnames(pcc.binary) <- rownames(pcc.binary) <- annotate(colnames(pcc.binary))

# build network and detect communities
net <- graph_from_adjacency_matrix(pcc.binary, mode="undirected")
clp <- cluster_label_prop(net)
V(net)$color <- "#FF000050" # could also color nodes depending on a factor

plot(clp, net, vertex.frame.color="white", col=V(net)$color, mark.col="#00000010", mark.border="#000000AA",
     vertex.size=10, vertex.label.family="Helvetica", vertex.label.cex=3/4, vertex.label.color="black",
     layout=layout_nicely)#with_kk) # Kamada-Kawai network layout (aka the spring-embedded network layout)

# calculate GO on the communities detected
# ...
```

## Other stuff that doesn't fit into any other category

### Supercomputing
Make use of several machines cooperating together. Data may reside in the main node which is then shared via sockets to the slaves:

```R
primary <- '192.168.1.235'
machineAddresses <- list(
  list(host=primary,user='johnmount',ncore=4),
  list(host='192.168.1.70',user='johnmount',ncore=4)
)

spec <- lapply(machineAddresses,function(machine) {
  rep(list(list(host=machine$host,user=machine$user)),machine$ncore)
})
spec <- unlist(spec,recursive=FALSE)

parallelCluster <- parallel::makeCluster(type='PSOCK',master=primary,spec=spec)
print(parallelCluster)
## socket cluster with 8 nodes on hosts
##                   ‘192.168.1.235’, ‘192.168.1.70’
```

### Parse arguments
This is just a suggestion on how to deal with input parms:

```R
parseArgs <- function(args,string,default=NULL,convert="as.character") {

    if(length(i <- grep(string,args,fixed=T)) == 1)
        return(do.call(convert,list(gsub(string,"",args[i]))))

    if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(T)
ftargets     <- parseArgs(args,"targets=","targets.txt")     # file describing the targets
fcontrasts   <- parseArgs(args,"contrasts=","contrasts.txt") # file describing the contrasts
mmatrix      <- parseArgs(args,"mmatrix=","~0+group")        # model matrix (or ~0+group for multiple comparisons)
filter.genes <- parseArgs(args,"filter=",TRUE,convert="as.logical") # filter invariant genes?
pre          <- parseArgs(args,"prefix=","")    # prefix to remove from the sample name
suf          <- parseArgs(args,"suffix=","_readcounts.tsv")    # suffix to remove from the sample name
cwd          <- parseArgs(args,"cwd=","./")     # current working directory
out          <- parseArgs(args,"out=","DE.edgeR") # output filename

runstr <- "Rscript DE.edgeR.R [targets=targets.txt] [contrasts=contrasts.txt] [mmatrix=~0+group] [filter=TRUE] [prefix=RE] [suffix=RE] [cwd=.] [out=DE.edgeR]"                                                       
if(!file.exists(ftargets))   stop("File",ftargets,"does NOT exist. Run with:\n",runstr)
if(!file.exists(fcontrasts)) stop("File",fcontrasts,"does NOT exist. Run with:\n",runstr)
if(!file.exists(cwd))        stop("Dir",cwd,"does NOT exist. Run with:\n",runstr)
if(is.na(filter.genes))      stop("Filter (filter invariant genes) has to be either TRUE or FALSE. Run with:\n",runstr)
```

### Write data into an Excel shit

Exporting dataframes into Excel shits can be done through the ''WriteXLS'' package. Suports old and new Excel format, which is decided by the function based on the extension of the putput file (.xls or .xlsx)

```R
require("WriteXLS")

df1 <- as.data.frame(matrix(rnorm(100  ),ncol=10))
df2 <- as.data.frame(matrix(rnorm(1000 ),ncol=10))
df3 <- as.data.frame(matrix(rnorm(10000),ncol=10))

out <- list(df1,df2,df3)
WriteXLS("out",ExcelFileName="R.xls",SheetNames=c("dataframe1","dataframe2","dataframe3"))
```

Keep in mind the limitations of both formats (xls=64K rows, 256 columns, xlsx= 1M rows, 16K columns).

### Read data from an Excel shit

Data can be the taken back into R with the ''readxl'' package:

```R
library(readxl)
all <- read_excel("./analysis/SUMMARY.xlsx", sheet="ALL")
```

Yet another option, using the gdata package:

```R
library(gdata)
df <- read.xls("myfile.xlsx"), sheet=1, header=TRUE)
```

### SQLite database access

Access a local SQLite database using the ''DBI'' and ''SQLite'' packages:

```R
 library("DBI")
 library("RSQLite")

 ## database connection
 drv <- dbDriver("SQLite")
 con <- dbConnect(drv, dbname="knowngenes.db")

 ## database access
 obtenirRegio <- function(chr,pos,strand) {

     chr    <- paste("chr",chr,sep="")
     strand <- ifelse(strand == "F","+","-")

     query <- paste("SELECT name,txStart,txEnd from genes",
                    " WHERE chrom=?",
                    "   AND ? BETWEEN txStart-1500 AND txEnd+1500",
                    "   AND strand != ?")

     regio <- dbGetQuery(con,query,data.frame(chr,pos,strand))

     return(c(paste(regio$name,collapse=";"),paste(regio$txStart,collapse=";"),paste(regio$txEnd,collapse=";")))
 }


 regio <- apply(fData450k[,c("chr","pos","strand")],1,obtenirRegio)
```

### Show a progress bar
The trick here is to print the "\r" special character to get to the first column of the current row of the screen, thus overwritting what has just been displayed:

```R

# progress bar
progress.old <<- 0
progress <- function(current) {
	if(current > progress.old) {	# just print every 1% increase, not for every step
		cat(paste(c("\r[",rep("=",current),rep(" ",100-current),"]",current,"%"),collapse=""))
		flush.console()
		progress.old <<- current
	}
}

s <- rnorm(1000000)

for(i in 1:length(s)) {

	# print the progress bar
	progress(round(100 * i / length(s)))

	# (your code here)
}
```

### Parallelize code that generates PDF plots

PDF files contain instructions for the drawing engine, thus forcing to be writen sequentially. To speed up parallel code which generates graphics, a trick is to generate high resolution TIFF images using the TIFF driver, and eventually merge them in a PDF (sequentially):

```R
library(parallel)
library(tiff)
CORES <- 32

# calculations part
img <- mclapply(1:1000000,function(i) {
    # some very heavy computation
    tmp   <- paste0("./tmp/",as.character(as.hexmode(round(rnorm(1) * 2^28))),".tiff")
    tiff(tmp)
    # plot whatever
    dev.off()
    tmp
},mc.cores=CORES)

# plotting part
pdf("plot.pdf")
lapply(img,function(tmp) {
    img <- readTIFF(tmp,native=T)
    plot(1:2,type="n",bty="n",axes=F,xlab="",ylab="")
    rasterImage(img,1,1,2,2)
    file.remove(tmp)
})
dev.off()
```

### Parallel by

''by()'' is often used as a replacement for ''tapply()'' to split-map-reduce like treatment of data frames.

```R
res <- do.call(rbind, by(x, x$key, function(x) {
    x[which.min(x$distance), ]
}))
```

The code is slow for relatively large datasets, both the ''do.call(rbind, ...'' and the single thread execution of ''by()''. The solutions comes to use multiple threads from the ''parallel'' package, use the ''data.table::rbindlist'' to reduce the final list:

```R
library(parallel)

res <- data.table::rbindlist(mclapply(split(x, x$key), function(x) {
    x[which.min(x$distance), ]
}, mc.cores=getOption("mc.cores", 4L)))
```
