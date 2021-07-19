# A collection of R snippets

## Table of Contents

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

f <- list.files(IN, pattern="*.tsv")
x <- mclapply(f, function(f) {
    x <- readLines(paste0(IN, f), n=10)
    x <- unlist(strsplit(x[grepl("# Command line:", x)], " "))
    f.bam <- x[grep("-t", x) + 1]
    peaks <- read.delim(paste0(IN, f), comment.char="#")
    peaks <- with(peaks, GRanges(chr, IRanges(start, end)))
    frip(f.bam, peaks)
}, mc.cores=CORES)

write.csv(x, file="./qc/frip.csv", row.names=F)
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
ann <- lapply(h2az.DB, annotatePeak, TxDb=txdb, annoDb="org.Mm.eg.db", tssRegion=c(-3000, 3000), verbose=T)
lapply(ann, plotAnnoBar)
lapply(ann, plotDistToTSS)
ann <- lapply(ann, as.data.frame)

WriteXLS("ann", ExcelFileName="./results/timecourse_h2az.xls",
         SheetNames=gsub("(.+)=\\((.+)\\)", "\\2", conts[, 1]))
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
peaks <- with(bed[sel, ], RangedData(IRanges(start, end), space=chr))
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

The ChIPseeker package suggests taking TSS of known genes. Actually, it can be done with any kind of GRanges object. ChIPseeker provides a not very flexible `plotHeatmap()` function which does the trick, but here I suggest a more customizable representation of the tags matrix.

**As a heatmap:**

In this example, we compare colocalization of 2 different peak marks. We plot how it looks like for a ChIP around anoter ChIP.

```R
chip1     <- GRanges(chip1$chr, IRanges(start=chip1$start, end=chip1$end))
chip2     <- GRanges(chip2$chr, IRanges(start=chip2$start, end=chip2$end))
chip2.3kb <- GRanges(seqnames(chip2), IRanges(start=chip2$summit - 1500, end=chip$summit + 1500))
tagMatrix <- getTagMatrix(chip1, windows=chip2.3kb)
tagMatrix <- t(apply(tagMatrix, 1, function(x) x/max(x))) # normalize, as in ChIPseeker:::peakHeatmap.internal
tagMatrix <- tagMatrix[order(rowSums(tagMatrix)), ]       # sort by signal

image(x=-1500:1500, y=1:nrow(tagMatrix), z=t(tagMatrix), useRaster=TRUE, col=c("black", "yellow"),
      yaxt="n", xaxt="n", ylab="", xlab="")
title("mark 1 on mark 2", cex.main=.5)
Axis(side=1, at=c(-1500, 0, 1500), labels=c("-1500", "mark 2", "+1500"), cex.axis=.5)
```

**As density lines:**

```R
getBamSignal <- function(f, gr, bins=101, region=500, filetype="bam") {
  # extend regions and read bam file. Center around the middle of the window
  gr.center <- start(gr) + floor(width(gr) / 2)
  gr2 <- GRanges(seqnames(gr), IRanges(gr.center - region, gr.center + region), strand=strand(gr))
  cvg <- switch(filetype,
                bam =coverage(readGAlignments(f, param=ScanBamParam(which=gr2))),
                bw  =coverage(rtracklayer::import.bw(f)),
                bed =coverage(rtracklayer::import.bed(f))
  )
  gr2 <- gr2[gr2 %within% reduce(GRanges(names(unlist(ranges(cvg))), unlist(ranges(cvg))))]

  # calcualte the coverage
  x <- do.call(rbind, mclapply(gr2, function(x) {
    # get the coverage around the motif
    if(as.character(strand(x)) == "+")
      x <- unlist(cvg[x], use.names=FALSE)
    else
      x <- rev(unlist(cvg[x], use.names=FALSE))

    # uncompress the coverage and split the region in 101 windows, and report the avg coverage per bin
    x <- tapply(decode(x), cut(1:length(x), bins), mean, na.rm=TRUE)     # aggregate reads per bin
  }))

  # aggregate the signal per bin and express it relative to the max
  avg <- apply(x, 2, sum)
  avg / max(avg)
}

# call the function for all bams in a folder
f <- list.files(pattern="\\.bam$")
motif_signal <- mcMap(function(f, ft) getBamSignal(f, motif[strand(motif) == "+"], filetype=ft),
                      f, gsub(".+(bed|bw|bam)$", "\\1", f), mc.cores=4)

# and plot
df <- melt(lapply(motif_signal, as.matrix))
ggplot(df) +
  geom_line(aes(x=as.numeric(Var1), y=value, color=L1, group=L1)) + #, lty=L1)) +
  scale_color_brewer("", palette="Set1") +
  scale_x_continuous(breaks=c(1, 26, 51, 76, 101), labels=c("-500", "-250", "CTCF", "+250", "+500")) +
  xlab("") + ylab("relative signal") + ggtitle("relative signal of CTCF positive BLISS hotspots") +
  theme_bw()

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
