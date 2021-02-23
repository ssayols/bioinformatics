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

## Array analysis

### Infinium HumanMethylation450 BeadChip

This pipeline uses ''lumi'' and ''IlluminaHumanMethylation450k.db'' packages from ''bioconductor''.

#### Read data from raw IDAT files

```R
 #options(mc.cores=12)
 library("lumi")    #analisi d'arrays illumina
 library("IlluminaHumanMethylation450k.db")

 samplesheet <- read.csv("Curelung_SampleSheet.csv")
 barcodes <- paste(samplesheet$Sentrix_ID, samplesheet$Sentrix_Position, sep="_")
 #idats <- methylumIDAT(barcodes=barcodes, parallel=T, idatPath=paste(getwd(), "idats", sep="/"))
 idats <- methylumIDAT(barcodes=barcodes, idatPath=paste(getwd(), "idats", sep="/"))
 dades <- as(idats, "MethyLumiM")

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
 dades.b.adj <- lumiMethyB(dades.c.adj, method="bgAdjust2C")
 # Perform quantile normalization based on color balance adjusted data
 dades.c.quantile <- lumiMethyN(dades.b.adj, method='quantile')
```

#### Quality control

```R
 plotColorBias1D(dades, channel='sum')
 plotColorBias1D(dades.c.quantile, channel='sum', main="Compare density distribution of two color channels - Norm")
 boxplotColorBias(dades, channel='sum')
 boxplotColorBias(dades.c.quantile, channel='sum', main="Boxplots of Red and Green color channels - Norm")
 plotColorBias2D(dades, selSample=1, cex=2)
 plotColorBias2D(dades.c.quantile, selSample=1, cex=2, main="Norm")
```

#### Differential methylation using limma

```R
 ## filter non variant probes by IQR
 eset <- exprs(dades.c.quantile)
 eset.IQR <- varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)

 # Define design and contrast matrices
 design <- model.matrix(~0 + pData(dades.c.quantile)$grup)
 rownames(design) <-pData(x)$sampleID
 colnames(design) <- c("G1", "G2")
 cont.matrix <- makeContrasts(Grup1vsGrup2=G1-G2, levels=design)

 # Get the differentially methylated probes by linear model + bayes stats
 # Model lineal i eBayes
 fit1 <- lmFit(eset.IQR, design)
 fit2 <- contrasts.fit(fit1, cont.matrix)
 fit3 <- eBayes(fit2)

 # Rank the probes
 CG <- topTable(fit3, adjust="fdr", number=nrow(fit3))  # adjust FDR
```

### Affymetrix® Genome-Wide Human SNP Array 6.0 (and probably others)

This is a parallelized version using ''snow'' and 16 cores to analyze SNP arrays with the ''oligo'' package. Can work with the ''MC'' package as well

```R
 library(ff)
 library(foreach)
 library(snow)
 library(doSNOW)
 registerDoSNOW(makeCluster(16, type="SOCK"))
 #library(MC)
 #library(doMC)
 #registerDoMC(16)

 # read arrays
 library(oligo)
 library(pd.genomewidesnp.6)
 #library(human650v3aCrlmm)  # Illumina HumanHap650Y
 #library(human550v3bCrlmm)  # Illumina HumanHap550Y

 fullFilenames <- list.celfiles("arrays", full.names=TRUE)
 outputDir <- file.path(getwd(), "crlmmResults")

 crlmm(fullFilenames, outputDir, pkgname="pd.genomewidesnp.6")
 #crlmm(fullFilenames, outputDir, pkgname="human650v3aCrlmm")
 #crlmm(fullFilenames, outputDir, pkgname="human550v3bCrlmm")
 crlmmOut <- getCrlmmSummaries(outputDir)

 write.csv(calls(crlmmOut), file="SNP.csv")
```

### Affymetrix expression arrays

Tested with Affymetrix Human Genome U133 Plus 2.0 Array, but might work with other from the U133 series and the newer Exon ST and Gene ST. Also might work with other non-Human organisms.

This protocol uses the ''affy'' package from ''bioconductor''

#### Read data from raw CEL files

```R
 # Load libraries
 library("affy")       #Import and normalize data
 library("limma")      #Differential expression
 library("genefilter") #Gene filtering

 # Import sample description from targets.txt file
 targets <- readTargets("arrays/targets.txt", row.names="FileName")

 # Import .CEL files
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
 # filter non variant probes by IQR
 esetIQR <- varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)

 # Define design and contrast matrices
 design <- model.matrix(~0 + as.factor(targets$Classes))
 colnames(design) <- c("AA", "CA", "HC")
 rownames(design) <- targets$FileName
 cont.matrix<-makeContrasts(AAvsCA=AA-CA, AAvsHC=AA-HC, CAvsHC=CA-HC, levels=design)

 # Get the differentially methylated probes by linear model + bayes stats
 # Model lineal and eBayes
 fit1 <- lmFit(esetIQR, design)
 fit2 <- contrasts.fit(fit1, cont.matrix)
 fit3 <- eBayes(fit2)

 # List of differentially expressed genes
 toptableIQR <- topTable(fit3, number=nrow(fit3), adjust.method="BH", sort.by="F")
```

