# A collection of R snippets

## Table of Contents

* [Plots](#plots-2)
   * [Color palettes](#color-palettes)
      * [Display a palette](#display-a-palette)
      * [A colorblind-friendly palette](#a-colorblind-friendly-palette)
      * [Color palettes with RColorBrewer](#color-palettes-with-rcolorbrewer)
      * [Color palettes with colorspace](#color-palettes-with-colorspace)
      * [ggplo2 color palette](#ggplo2-color-palette)
   * [Spaghetti plot: graphs showing regression uncertainty](#spaghetti-plot-graphs-showing-regression-uncertainty)
   * [Smooth spline ggplot2](#smooth-spline-ggplot2)
   * [MA plot](#ma-plot)
   * [Circle](#circle)
   * [Ellipse around the CI95](#ellipse-around-the-ci-95)
   * [Boxplot with average and standard deviation](#boxplot-with-average-and-standard-deviation)
   * [Histogram and density function in the same plot using ''ggplot''](#histogram-and-density-function-in-the-same-plot-using-ggplot)
   * [Pairs plot (nice)](#pairs-plot-nice)
   * [Placing multiple plots in the same device](#placing-multiple-plots-in-the-same-device)
   * [Placing multiple plots in the same device using ''ggplot'' the ''grid'' package](#placing-multiple-plots-in-the-same-device-using-ggplot-the-grid-package)
   * [Heatmap + Upset](#heatmap-+-upset)
   * [Venn Diagrams](#venn-diagrams)
   * [Venn Diagrams of GRanges objects](#venn-diagrams-of-granges-objects)
   * [Calling gnuplot](#calling-gnuplot)
   * [Density 2d plots](#density-2d-plots)
      * [Base system](#base-system)
      * [smoothScatter](#smoothscatter)
      * [density2d (contour plot)](#density2d-contour-plot)
      * [ggplot2 density2d](#ggplot2-density2d)
      * [ggplot2 color points by density](#ggplot2-color-points-by-density)
      * [3d plot](#3d-plot)
   * [PCA BiPlot with ggplot2 and Interactive plots in a shiny app](#pca-biplot-with-ggplot2-and-interactive-plots-in-a-shiny-app)
   * [Raster kegg PNG pathways within a coordinate system](#raster-kegg-png-pathways-within-a-coordinate-system)
   * [Plot chromosome ideograms](#plot-chromosome-ideograms)
   * [Plot chromosome ideograms with additional tracks](#plot-chromosome-ideograms-with-additional-tracks)
   * [Gviz plots](#gviz-plots)
   * [Circos plots](#circos-plots)
   * [Graphical representation of contingency tables](#graphical-representation-of-contingency-tables)
   * [Sequence motif](#sequence-motif)
* [Miscellaneous bioinformatic related stuff](#miscellaneous-bioinformatic-related-stuff)
   * [Barcode design](#barcode-design)
   * [Getting things from Biomart](#getting-things-from-biomart)
   * [Getting things from the Ensembl databases](#getting-things-from-the-ensembl-databases)
   * [Getting things from the UCSC tables and ENCODE tracks](#getting-things-from-the-ucsc-tables-and-encode-tracks)
   * [Drug response curves using the ''drc'' package](#drug-response-curves-using-the-drc-package)
   * [Permutation test to identify if 2 curves are significantly different](#permutation-test-to-identify-if-2-curves-are-significantly-different)
   * [GenomicRanges](#genomicranges)
   * [Converting seqlevel styles](#converting-seqlevel-styles)
   * [LiftOver coordinates between different assemblies](#liftover-coordinates-between-different-assemblies)
   * [Converting GTF to GFF3](#converting-gtf-to-gff3)
   * [Flatten a GFF/GTF file by gene_id (and get transcript lengths)](#flatten-a-gffgtf-file-by-gene_id-and-get-transcript-lengths)
   * [Get genome wide distribution of features](#get-genome-wide-distribution-of-features)
   * [Download a dataset from GEO using GEOquery](#download-a-dataset-from-geo-using-geoquery)
   * [Download a dataset from Gene Expression Atlas](#download-a-dataset-from-gene-expression-atlas)
   * [Get the genomic sequence of a region and plot its nucleotide content](#get-the-genomic-sequence-of-a-region-and-plot-its-nucleotide-content)
   * [Gene set enrichment analysis](#gene-set-enrichment-analysis)
      * [GSEA of Biological Processes](#gsea-of-biological-processes)
      * [GSEA of Biological Processes (a parallel version)](#gsea-of-biological-processes-a-parallel-version)
      * [GSEA of KEGG or custom db](#gsea-of-kegg-or-custom-db)
      * [GSEA of a custom GO (slim) db](#gsea-of-a-custom-go-slim-db)
      * [GO functional analysis with clusterProfiler](#go-functional-analysis-with-clusterprofiler)
      * [GO functional analysis with clusterProfiler and DAVID](#go-functional-analysis-with-clusterprofiler-and-david)
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
   * [Convert BAM to BigWig](#convert-bam-to-bigwig)
   * [Shannon index for the nt diversity of a sequence](#shannon-index-for-the-nt-diversity-of-a-sequence)
* [Statistical analysis](#statistical-analysis)
   * [Principal Component Analysis (PCA)](#principal-component-analysis-pca)
   * [High Level PCA](#hihg-level-pca)
   * [Project a new vector onto PCA space](#project-a-new-vector-onto-pca-space)
   * [Multidimensional Scaling (MDS)](#multidimensional-scaling-mds)
   * [Affinity Propagation Clustering](#affinity-propagation-clustering)
   * [Other methods for clustering](#other-methods-for-clustering)
   * [Validate the optimal number of clusters](#validate-the-optimal-number-of-clusters)
   * [Evaluate cluster strength](#evaluate-cluster-strength)
   * [Evaluate cluster similarity](#evaluate-cluster-similarity)
   * [Calculate ROC curves from a predictor using the ''ROCR'' package](#calculate-roc-curves-from-a-predictor-using-the-rocr-package)
   * [Impute missing values](#impute-missing-values)
   * [Geometric mean](#geometric-mean)
   * [Power analysis](#power-analysis)
   * [Power and sample size calculation for survival analysis](#power-and-sample-size-calculation-for-survival-analysis)
   * [Calculate the center of a 2d-distribution](#calculate-the-center-of-a-2d-distribution)
   * [Test the significance of the overlap between 2 lists](#test-the-significance-of-the-overlap-between-2-lists)
   * [Network analysis](#network-analysis)
   * [Split data for CV](#split-data-for-cv)
   * [Compute PI-value based on FC and p-value ](#compute-pi-value-based-on-fc-and-p-value)
* [Other stuff that doesn't fit into any other category](#other-stuff-that-doesnt-fit-into-any-other-category)
   * [Supercomputing](#supercomputing)
   * [Parse arguments](#parse-arguments)
   * [Write data into an Excel shit](#write-data-into-an-excel-shit)
   * [Read data from an Excel shit](#read-data-from-an-excel-shit)
   * [SQLite database access](#sqlite-database-access)
   * [Other databases](#other-databases)
   * [In memory database](#in-memory-database)
   * [Show a progress bar](#show-a-progress-bar)
   * [Parallelize code that generates PDF plots](#parallelize-code-that-generates-pdf-plots)
   * [Parallel by](#parallel-by)
   * [Split and Join PDF files](#split-and-join-pdf-files)

## Plots

### Color palettes

#### Display a palette

Very manual approach:

```R
show_palette <- function(colors) {
  image(1:length(colors), 1, as.matrix(1:length(colors)), col=colors, 
    xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
}

show_palette(rainbow(6))
```

Alternatively, using the scales package:

```R
palette("Classic Tableau")
scales::show_col(palette())
```

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
display.brewer.pal(9, "Oranges")

# use it as a gradient palette for heatmaps
hmcol   <- colorRampPalette(brewer.pal(9, "Oranges"))(100)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))

# use it to display categorical data
plot(x=x$x[, 1], y=x$x[, 2], pch=16, cex=.5,
     col=brewer.pal(length(levels(condition)), "Set1")[condition])
```

#### Color palettes with colorspace

```R
require("colorspace")

# show all available palettes
hcl_palettes(plot = TRUE)

# Usage with base graphics
q4 <- qualitative_hcl(4, "Dark 3")
plot(log(EuStockMarkets), plot.type = "single", col = q4, lwd = 2)
legend("topleft", colnames(EuStockMarkets), col = q4, lwd = 3, bty = "n")

# Usage with ggplot2
library("ggplot2")
ggplot(iris, aes(x = Sepal.Length, fill = Species)) + geom_density(alpha = 0.6) +
  scale_fill_discrete_qualitative(palette = "Dark 3")
  
# Palette visualization and assessment
demoplot(q4, "bar")
hclplot(q4)
specplot(q4, type = "o")

s9 <- sequential_hcl(9, "Purples 3")
demoplot(s9, "heatmap")
hclplot(s9)
specplot(s9, type = "o")
```

#### ggplo2 color palette

```R
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
```

### Spaghetti plot: graphs showing regression uncertainty

The idea is taken from [http://andrewgelman.com/2012/08/26/graphs-showing-regression-uncertainty-the-code/ here] and is only partially implemented (only the spaghetti code, not the the shaded CI).

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

### Smooth spline ggplot2

Use a different fit function from the predefined ones in `geom_smooth()`.

```R
# From Hadley, https://groups.google.com/forum/#!topic/ggplot2/FJ36CJH-ODo

smooth.spline2 <- function(formula, data, ...) {
  mat <- model.frame(formula, data)
  smooth.spline(mat[, 2], mat[, 1])
}

predictdf.smooth.spline <- function(model, xseq, se, level) {
  pred <- predict(model, xseq)
  data.frame(x = xseq, y = pred$y)
}

qplot(mpg, wt, data = mtcars) + geom_smooth(method="smooth.spline2", se=F)
```

### MA plot

Visual representation of log ratios (M) and mean average (A) intensities. Useful to compare bias between two samples or in a two-channels array.

```R
 maplot <- function(x, y, ...){
   M=log2(y/x)
   A=(log2(y)+log2(x))/2
   plot(A, M, ...)
   fit=loess(M~A, span=1/3)
   o=order(A)
   lines(A[o], fit$fitted[o], col=2)
 }

 maplot(y[wh, 2], y[wh, 3], ylim=c(-5, 5))
```

### Circle

There is no primitive function to draw a circle in R.

```R
circle <- function(x, y, r, ...) {  # col=fill_colour, border=outline_colour
    polygon(x=x + r*cos( seq(0, 2*pi, length.out=360) ),
            y=y + r*sin( seq(0, 2*pi, length.out=360) ),
            ...)
}
```

### Ellipse around the CI95

Draw an ellipse describing the CI 95% of some (x, y) points. Useful to describe groups when reducing to only a couple of components (through a PCA or MDS) large datasets with several independent variables describing a sample.

```R
ell95 <- function(df) {
  require(ellipse)
  require(ggplot2)

  ell <- data.frame()
  for(g in levels(df$group)) {
    e <- with(df[df$group==g, ], ellipse(cor(x, y), scale=c(sd(x), sd(y)), centre=c(mean(x), mean(y))))
    ell <- rbind(ell, cbind(as.data.frame(e), group=g))
  }
  ggplot(data=df, aes(x=x, y=y, colour=group)) + geom_point() + geom_path(data=ell, aes(x=x, y=y, colour=group))
}

print(ell95(data.frame(x=x, y=y, group=A)))
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

df <- data.frame(cond=c(rep("before", length(before)), rep("after", length(after))),
                 vals=c(before, after))
p <- ggplot(df, aes(x=vals, fill=cond)) +
  geom_histogram(aes(y=..density.., fill=cond), alpha=.5, position="identity") +
  geom_density(alpha=.2) +
  ggtitle(rownames(ic50)[i])

print(p)
```

### Pairs plot (nice)

A version with color density:

```R
#' panel.smooth for a pairs() plot
#' Modified version of the base function panel.smooth, in order to accomodate a
#' density color palette.
#' @param x vector of numbers corresponding to the x coordinates
#' @param y vector of numbers corresponding to the y coordinates
#' @param ... rest of the parameters that go to points() and lines()
#'
#' @examples
#'   pairs(x, panel=panel.smooth.dens, diag.panel=panel.hist, pch=16, col="#00000050")
panel.smooth.dens <- function (x, y, bg=NA, pch=par("pch"), cex=1,
                               col.smooth="red", span=2/3, iter=3, ...)
{
  ok <- is.finite(x) & is.finite(y)
  if(any(ok)) {
    dcols <- densCols(x=x[ok], y=y[ok], colramp=viridis, nbin=100)
    points(x[ok], y[ok], pch=pch, col=dcols, bg=bg, cex=cex)
    lines(stats::lowess(x[ok], y[ok], f=span, iter=iter), col=col.smooth, ...)
  }
}

#' panel.hist for a pairs() plot
#' Taken literally from examples(pairs)
#' @param x vector of numbers
panel.hist <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, ...)
}
```

Yet another version, with a regular scatter and the correlation in the lower diagonal, taken from [here](http://www.sthda.com/english/wiki/scatter-plot-matrices-r-base-graphs):

```R
my_cols <- c("#00AFBB", "#E7B800", "#FC4E07")

# Correlation panel
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits=2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)  # text is proportional to the correlations
}

# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = my_cols[iris$Species])
}

# Create the plots
pairs(iris[,1:4], 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
```

The correlation panel could be integrated into the upper panel, and display the `panel.hist()` instead.

### Placing multiple plots in the same device

This solution makes use of the basic R plotting capabilities, through the `graphics` package:

```R
# define grid (layout matrix)
layoutRatioWidth <- c(0.75, 0.25); layoutRatioHeight <- c(0.25, 0.75)
layout(matrix(c(2, 1, 0, 3), nrow=2), widths=layoutRatioWidth, heights=layoutRatioHeight, respect=F)
# bottom left plot (#1 in layout matrix)
par(mar=c(4, 4, 1, 1))
smoothScatter(x, y, xlab=xlab, ylab=ylab, main=f, xlim=c(0, 20), ylim=c(0, 20))
# top left plot (#2 in layout matrix)
par(mar=c(1, 4, 1, 1))
z <- density(x, na.rm=T)
plot(z$x, z$y, xlab="", ylab="Density", xaxt="n", yaxt="n", type="l", main="")
# bottom right plot (#3 in layout matrix)
par(mar=c(4, 1, 1, 1))
z <- density(y, na.rm=T)
plot(z$y, z$x, ylab="", xlab="Density", xaxt="n", yaxt="n", type="l", main="")
```

### Placing multiple plots in the same device using ''ggplot'' the ''grid'' package

First option, using the ''grid'' package and creating ''viewports'':

```R
 library("grid")

 vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

 plot1 <- qplot(mtcars, x=wt, y=mpg, geom="point", main="Scatterplot of wt vs. mpg")
 plot2 <- qplot(mtcars, x=wt, y=disp, geom="point", main="Scatterplot of wt vs disp")
 plot3 <- qplot(wt, data=mtcars)
 plot4 <- qplot(wt, mpg, data=mtcars, geom="boxplot")

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

### Heatmap + Upset

Plot a heatmap annotated with an upset plot.

```R
#' Heatmap + Upset function
#'
#' Plot a heatmap annotated with an upset plot
#'
#' @param counts Counts matrix
#' @param gs List of gene sets
#' @param main Plot title
#'
#' @example
#' # Create some demo data
#' counts <- matrix(rpois(600, 10), ncol=6) # counts matrix with 100 genes, 2 conditions and 3 replicates
#' colnames(counts) <- c("A1", "A2", "A3", "B1", "B2", "B3")
#' rownames(counts) <- paste("gene", 1:100)
#'
#' genesets <- list(set1=c("gene 1" , "gene 2" , "gene 34", "gene 11", "gene 66"),
#'                  set2=c("gene 2" , "gene 3" , "gene 35"),
#'                  set3=c("gene 34", "gene 66", "gene 67", "gene 68"))
#'
#' heatmap_upset(counts, genesets, "test")
heatmap_upset <- function(counts, gs, main="") {

    require(ggplot2)
    require(grid)
    require(reshape2)

    # get the gene expression per gene set
    x <- lapply(gs, function(genes) {
        counts[intersect(rownames(counts), genes), ]
    })

    # remove empty gene sets and melt
    x     <- x[sapply(x, nrow) > 0]
    df    <- melt(x)
    colnames(df) <- c("gene", "sample", "counts", "geneset")

    # calculate how the rows (genes) cluster
    x <- do.call(rbind, x)
    x <- unique(x, MARGIN=1)
    h <- hclust(dist(x))
    df$gene <- factor(df$gene, levels=h$labels[h$order])

    # plot heatmap
    p1 <- ggplot(df, aes(sample, gene)) +
            geom_tile(aes(fill=counts), colour="white") +
            scale_fill_distiller(palette="RdBu", direction=1, limits=c(0, max(df$counts))) +
            ggtitle(main) +
            xlab("") +
            theme_minimal() +
            theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
            theme(legend.position="left")

    # plot upset
    p2 <- ggplot(df, aes(geneset, gene)) +
            geom_point() +
            ggtitle("gene belongs to set") +
            xlab("") +
            theme_minimal() +
            theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
            theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

    grid.newpage()
    grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size="last"))
}
```

### Venn Diagrams

Several Bioconductor packages are available for Venn Diagrams, "VennDiagram" is able to produce high quality proportional diagrams.

```R
require(VennDiagram)
# main argument a list of vectors
# function will calculate overlap between the vectors and produce proportional Venn diagrams.

venn.diagram(list(C = 1700:2500, B = 1:1800, A = 1571:2020),
  fill=RColorBrewer::brewer.pal(3, "Set1"), alpha=.3, lwd=0, cex=1.5, cat.cex=1.5,
  filename="~/Desktop/test2.tif")

# display it on the screen
img <- tiff::readTIFF("~/Desktop/test2.tif", native=T)
plot(1:2, type="n", bty="n", axes=F, xlab="", ylab="")
rasterImage(img, 1, 1, 2, 2)
```

To display it on the screen, alternatively one can set to 'filename' to NULL and use the grid package to draw the result:

```R
img <- venn.diagram(list(C = 1700:2500, B = 1:1800, A = 1571:2020),
  fill=RColorBrewer::brewer.pal(3, "Set1"), alpha=.3, lwd=0, cex=1.5, cat.cex=1.5, filename=NULL)
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
    Counts[i, 1] <- universe[i] %in% set1
    Counts[i, 2] <- universe[i] %in% set2
    Counts[i, 3] <- universe[i] %in% set3
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

### Venn diagrams of GRanges objects

The easiest is to use `makeVennDiagram()` from the `ChIPpeakAnno` package:

```R
library(ChIPpeakAnno)
peaks1 <- GRanges("chr1", IRanges(seq(1, 100, 5), width=2), "+")
peaks2 <- GRanges("chr1", IRanges(seq(2, 20, 3), width=2), "+")
peaks3 <- GRanges("chr1", IRanges(seq(10, 50, 4), width=2), "+")
res <- makeVennDiagram(Peaks=list(peaks1, peaks2, peaks3),
                       NameOfPeaks=c("TF1", "TF2", "TF3"))
```

Or:

```R
peaks1$type<-"TF1"
peaks2$type<-"TF2"
peaks3$type<-"TF3"
gr <- c(peaks1, peaks2, peaks3) # like your data

grl <- splitAsList(gr, gr$type)
res <- makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl))
```

Otherwise, one can calculate overlaps with `GenomicRanges` and manually plot with `venneuler`:

```R
library(venneuler)

x <- reduce(c(peaks1, peaks2))   # common+unique peakset
x$peaks1 <- x %over% peaks1
x$peaks2 <- x %over% peaks2

AB <- sum( x$peaks1 &  x$peaks2)
A  <- sum( x$peaks1 & !x$peaks2)
B  <- sum(!x$peaks1 &  x$peaks2)

v  <- venneuler(x[, c("peaks1", "peaks2")])
names(v$labels) <- v$labels
v$labels["peaks1"] <- paste(A , "peaks in cond1")
v$labels["peaks2"] <- paste(B , "peaks in cond2")
# `v` can be manipulated to include additional labels (eg. in AB). See `v` for x,y coordinates of the circles
plot(v, sub=paste("peaks overlap=", AB))
```

### Calling gnuplot

From R, one can pipe commands into gnuplot after opening a new session:

```R
library(Rgnuplot)
h1<-gp.init()
gp.cmd(h1, "set terminal pngcairo  transparent enhanced font 'arial, 10' fontscale 1.0 size 500, 350")
gp.cmd(h1, "set output 'out.png'")  # set output file
gp.cmd(h1, "set key left box")    # include a boxed legend
gp.cmd(h1, "plot [-10:10] sin(x), atan(x) with points, cos(atan(x)) with impulses")
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
par(mar=c(5, 0, 5, 0))
plot(NA, xlim=c(0, 10), ylim=c(0, 11), type="n", ann=FALSE, axes=FALSE)
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
pal <- colorRampPalette(c('dark blue', 'blue', 'light blue', 'yellow', 'orange', 'red', 'dark red'))
filled.contour(kde2d(log10(DupMat[, "RPK"]), 100 * DupMat[, "dupRate"]),
xlab = "expression level (reads/kbp)", ylab = "duplication level (% duplicate reads)",
axes = FALSE, color.palette=pal, ...)
```

The palette can be modified with RColorBrewer to visually improve the result:
```R
require(RColorBrewer)
pal <- function(n) { colorRampPalette(brewer.pal(9, "Oranges"))(n) }
```

Something easier, leaving the estimation to smoothScatter():
```R
cols <- colorRampPalette(c("black", "blue", "green", "yellow", "red"))
smoothScatter(dm, nrpoints=10, nbin=500, colramp=cols)
```

#### ggplot2 density2d

```R
df <- data.frame(x=log10(DupMat[, "RPK"]), y=100 * DupMat[, "dupRate"])
p <- ggplot(df, aes(x, y)) +
stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
scale_x_continuous(breaks=s, labels=as.character(10^s)) +
xlab("expression level (reads/kbp)") + ylab("duplication level (% duplicate reads)") +
theme_bw()
print(p)
```

#### ggplot2 color points by density

Slightly different approach, using the native `geom_point` colored by density. Equivalent to using `densCols()` from base.

```R
# from Kamil Slowikowski (http://slowkow.com/notes/ggplot2-color-by-density/)
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
# @examples
#   dat$density <- get_density(dat$x, dat$y)
#   ggplot(dat) + geom_point(aes(x, y, color = density)) + scale_color_viridis()
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
```

#### 3d plot

```R
require(MASS)
pal <- colorRampPalette(c('dark blue', 'blue', 'light blue', 'yellow', 'orange', 'red', 'dark red'))
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
        (max(data[, y]) - min(data[, y])/(max(datapc[, y])-min(datapc[, y]))),
        (max(data[, x]) - min(data[, x])/(max(datapc[, x])-min(datapc[, x])))
        )
    datapc <- transform(datapc,
            v1 = .7 * mult * (get(x)),
            v2 = .7 * mult * (get(y))
            )
    plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
    plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2, "cm")), alpha=0.75, color="red")
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

g <- rasterGrob(img, 0, 0, ncol(img), nrow(img), just=c("left", "bottom"))
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

        ideo <- ideogramTab[ideogramTab$chr == chr, ]
        center <- ideo[ideo$braç=="cen", "end"]
        chrLength <- max(ideo$end)
        # el viewport comença un 10% abans de la coordenada 0, i acaba un 20% mes tard
        pushViewport(dataViewport(xData=c(0-chrLength * .18, chrLength * 1.32), yData=c(0, 6), extension=0, layout.pos.col=1, layout.pos.row=1))

        # pintem rectangles per a cada banda. El color be definit per:
        #     -si V8=="gpos" --> un gradient de la paleta de grisos [1:100] segons V9
        #     -altrament --> blanc. En aquest cas V9 == NA
        ideo[is.na(ideo[, 9]), 9] <- 1
        pal <- colorRampPalette(c("white", "black"))(100)
    for(i in seq(along = ideo[, 1])) {
                grid.rect(x=ideo[i, 6],
                y=2,
                width=ideo[i, 7]-ideo[i, 6],
                height=2,
                gp=gpar(col=pal[ideo[i, 9]],
                fill=pal[ideo[i, 9]]),
                default.units="native",
                just=c("left", "bottom"))
        }

        # center == centromer (on creuarem els dos braços)
        # requadre al braç p
        grid.lines(c(0, center-500000), c(4, 4), default.units = "native")
        grid.lines(c(0, center-500000), c(2, 2), default.units = "native")
        grid.lines(c(0, 0), c(2, 4), default.units = "native")
        # requadre al braç q
        grid.lines(c(center+500000, chrLength), c(4, 4), default.units = "native")
        grid.lines(c(center+500000, chrLength), c(2, 2), default.units = "native")
        grid.lines(c(chrLength, chrLength), c(2, 4), default.units = "native")
        # creuament al centromer
        grid.lines(c(center-500000, center+500000), c(4, 2), default.units = "native")
        grid.lines(c(center-500000, center+500000), c(2, 4), default.units = "native")
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
                      ylim = c(10, 50),
                      color = "black",
                      size = 0.4
                      )

```



### Gviz plots

The nice [http://bioconductor.org/packages/release/bioc/html/Gviz.html GViz] package for visualization of genomic data. More examples [https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/ here] and [http://www.sthda.com/english/wiki/gviz-visualize-genomic-data here]:

```R
#################################
##
## Gviz visualization of KMT2A coordinates: chr11:118, 307, 205-118, 397, 539 (hg19)
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

### Circos plots

The nice [https://jokergoo.github.io/circlize_book/book/ Circlize] package for visualization of genomic data. More examples in his book.

Very basic usage, with a barplot showing number of breaks at some genomic loci:

```R
library(circlize)

# load counts of each sample to known AsiSI cut sites
x <- read.delim("./results/asisi_sites.txt")
x$seqnames <- paste0("chr", x$seqnames)

#  seqnames   start     end width strand asisi_ctrl_1 asisi_ctrl_2 asisi_4oth_1
#1        1 1110919 1110926     8      *            0            0            4
#2        1 2063312 2063319     8      *            0            0            0
#3        1 2152681 2152688     8      *            0            0            0
#4        1 2322833 2322840     8      *            0            0            0
#5        1 2985156 2985163     8      *            0            0            1
#6        1 3103043 3103050     8      *            0            0            1

for(i in grep("^asisi", colnames(x))) {  #  samples start with asisi
    df <- x[x[[i]] > 0, c(1:3, i)]       # columns 1:5 are chr, start, end, width, strand
     
    circos.initializeWithIdeogram(species="hg19")
        circos.genomicTrack(df, panel.fun=function(region, value, ...) {
        circos.genomicLines(region, value, type="h", col=c("#FF000080"))
          })
    text(0, 0, colnames(x)[i], cex = 0.6)
}
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

### Sequence motif
Conventional motifs can be computed using the [motifStack](https://bioconductor.org/packages/release/bioc/html/motifStack.html) package. Arbitrary motifs (that use arbitrary numbers instead of PCM, PWM, ICM, etc.) can be plotted using [ggseqlogo](https://omarwagih.github.io/ggseqlogo/):

```R
# Create a custom matrix 
set.seed(123)
custom_mat = matrix( rnorm(20), nrow=4, dimnames=list(c('A', 'T', 'G', 'C')))

# Generate sequence logo
ggseqlogo(custom_mat, method='custom', seq_type='dna') + ylab('my custom height')
```

## Miscellaneous bioinformatic related stuff

### Barcode design

Design custom and robust DNA barcode sets capable of correcting substitution errors or insertion, deletion, and substitution errors. Use the [http://bioconductor.org/packages/release/bioc/html/DNABarcodes.html DNABarcodes] package from Bioconductor. More [http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0036852 here].

```R
library("DNABarcodes")
mySet <- create.dnabarcodes(5)
# 1) Creating pool ...  of size 592
# 2) Conway closing...  done
how(mySet)
#  [1] "GAGAA" "AGCAA" "CCTAA" "CAAGA" "ACGGA" "GTCGA" "TGTGA" "GGACA"
#  [9] "CTGCA" "TACCA" "CGAAG" "TCGAG" "GTTAG" "ATAGG" "AAGCG" "GAATG"
# [17] "TGCTG" "ACTTG" "ACAAC" "CACAC" "TAGGC" "CTTGC" "TTACC" "GATCC"
# [25] "AGGTC" "GCCAT" "TCAGT" "AACGT" "TGGCT" "CAGTT"
nalyse.barcodes(mySet)
#                   Description  hamming   seqlev levenshtein
# 1               Mean Distance 5.242908 3.656915    4.560284
# 2             Median Distance 5.000000 4.000000    5.000000
# 3            Minimum Distance 3.000000 1.000000    2.000000
# 4            Maximum Distance 7.000000 7.000000    7.000000
# 5 Guaranteed Error Correction 1.000000 0.000000    0.000000
# 6  Guaranteed Error Detection 2.000000 0.000000    1.000000
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
gene.ann <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"),
                  filter="biotype", values="protein_coding",
                  mart=useMart("ensembl", dataset="hsapiens_gene_ensembl"))

## connect to HG18
#mart  <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
#  host="may2009.archive.ensembl.org",
#  path="/biomart/martservice", archive=FALSE)  # Ensembl54 (hg18/NCBI36))

## connect to Variation (SNPs)
# mart <- biomaRt::useMart("snp", dataset="hsapiens_snp")
# snp  <- biomaRt::getBM(attributes=c("refsnp_id", "chr_name", "chrom_start"),
#  filters="snp_filter", values=snp, mart=mart)
```

### Getting things from the Ensembl databases

Alternatively, one can get the same stuff directly from the Ensembl databases. Also an alternitive from their messy perl API...

```R
library(RMySQL)
con   <- dbConnect(MySQL(), user="anonymous", host="ensembldb.ensembl.org", dbname="homo_sapiens_core_79_38")
query <- paste("SELECT g.stable_id, t.stable_id, p.stable_id, x.display_label",
               "  FROM gene AS g, transcript AS t, translation AS p, xref AS x",
               " WHERE t.transcript_id = g.canonical_transcript_id",
               "   AND t.transcript_id = p.transcript_id",
               "   AND x.xref_id = g.display_xref_id")
x     <- dbGetQuery(conn=con, statement=query)
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
tableNames(ucscTableQuery(mySession, track=track.names[1]))
```

Identify repeat-masked regions in and around the transcription
start site (TSS) of the human E2F3 gene, in hg19:

```R
e2f3.tss.grange <- GRanges("chr6", IRanges(20400587, 20403336))
tbl.rmsk <- getTable(ucscTableQuery(mySession, track="rmsk", range=e2f3.tss.grange, table="rmsk"))
```

Get DNaseI hypersensitivity regions in the K562 Cell Line from the ENCODE project:

```R
track.name <- "wgEncodeUwDgf"
table.name <- "wgEncodeUwDgfK562Hotspots"
e2f3.grange <- GRanges("chr6", IRanges(20400587, 20403336))
tbl.k562.dgf.e2f3 <- getTable(ucscTableQuery(mySession, track=track.name, range=e2f3.grange, table=table.name))
tbl.k562.dgf.hg19 <- getTable(ucscTableQuery(mySession, track=track.name, table=table.name))
```

### Drug response curves using the ''drc'' package

```R
 library(drc)

 # llegir dades
 x <- read.csv("x.csv")
 x$dose <- 2^x$concen

 # fit into a sigmoidal model
 mock <- drm(mock ~ dose, data=x,
       fct=LL.4(names=c("Slope", "Lower Limit", "Upper Limit", "ED50")))

 oe   <- drm(oe ~ dose, data=x,
       fct=LL.4(names=c("Slope", "Lower Limit", "Upper Limit", "ED50")))

 # report the IC50 with its confidence intervals
 ED(mock, c(5, 50, 95), interval="delta")
 ED(oe  , c(5, 50, 95), interval="delta")

 # plot
 plot(mock, broken=T, type="all", col="red", ylab="viability")
 plot(oe  , broken=T, type="all", col="blue", add=T)
 legend("center", fill=c("red", "blue"), legend=c("mock", "oe"))
```

### Permutation test to identify if 2 curves are significantly different

First method is based on anova multiple sample test described in Elso et al, 2004 and implemented in the ''statmod'' package:

```R
 #
 # Input data has the folowing structure:
 #
 # +---------------+---------+-----+---------+-------+
 # |  concentracio | x.clon1 | ... | x.clonN | Group |
 # +---------------+---------+-----+---------+-------+
 # |           1uM |   0.66  |     |   0.69  |  MOCK |
 # |           2uM |   0.68  |     |   0.71  |  MOCK |
 # |  ...                                            |
 # |           1uM |   0.56  |     |   0.59  | SCRAM |
 # |           2uM |   0.58  |     |   0.51  | SCRAM |
 # |  ...                                            |
 # |           1uM |   0.54  |     |   0.69  | SKMEL |
 # |           2uM |   0.55  |     |   0.61  | SKMEL |
 # |  ...                                            |
 # +---------------+---------+-----+---------+-------+
 #
 library(statmod)
 library(reshape2)
 library(parallel)

 # permutation test
 tperm <- function(f, nperm=1000) {

   # read data an normalize every clon with its first concentration value
   df <- read.csv(f)
   vals <- grep("^x", colnames(df))
   df.norm <- as.list(by(df[, vals], df$Group, function(df) { apply(df, 2, function(x) x / x[1]) }))
   A  <- rep(names(df.norm), each=length(vals))
   dades <- Reduce(cbind, df.norm)

   # permutation test
   s <- compareGrowthCurves(A, t(as.matrix(dades)), nsim=nperm)
   write.csv(s, file=paste0("stats.", f))

   # plot the dosage curves (be careful, the regression is polinomic (loess) and not logistic)
   pdf(paste0("stats.", f, ".pdf"))
   df.norm <- lapply(df.norm, as.data.frame)
   for(i in 1:length(df.norm)) {
     df.norm[[i]]$uM <- log2(unique(df$uM))
     df.norm[[i]]$Group <- names(df.norm)[i]
   }
   dades <- Reduce(rbind, df.norm)
   df2 <- melt(dades, id.vars=c("uM", "Group"))
   p <- ggplot(df2, aes(x=uM, y=value, colour=Group)) + geom_smooth() + geom_point()
   print(p)
   dev.off()
 }

 mclapply(list.files(pattern="*.csv"), tperm, nperm=1000, mc.cores=4)
```

Second method is based on the F-score from an anova test, were we compare the original score to those obtained by permutation, using an empirical cumulative distribution to get the pvalue:

```R
 tperm <- function(f, nperm) {
   dades <- read.csv(paste(f, ".csv", sep=""))
   dades <- reshape2::melt(dades, id.vars=c("Concentration", "Group"))
   dades$Concentration <- as.factor(dades$Concentration)
   dades$Group <- as.factor(dades$Group)
   d <- anova(lm(value ~ Group*Concentration, data=dades))["Group:Concentration", "F value"]

   # get nperm scores (permutation test)
   require(parallel)

   d.null <- mclapply(1:nperm, function(x) {
     # permute (only permute the values inside every concentration)
     perm <- by(dades, dades[, "Concentration"], function(x) {
             x$value <- x$value[sample(1:nrow(x))];
             return(x) })
     perm <- do.call(rbind, perm)
     d <- anova(lm(value ~ Group*Concentration, data=perm))["Group:Concentration", "F value"]
   }, mc.cores=4)

   # build an empirical distribution, and output the pvalue to get more extreme D
   p.null <- ecdf(unlist(d.null))
   print(paste(f, "pval:", 1 - p.null(d)))
   return(p.null)
 }

 s <- sapply(c("akt", "sft31", "shikonin"), tperm, 10000)
```

### GenomicRanges

Example using the ''GenomicRanges'' package from ''bioconductor'' to find overlapping features:

```R
# more things can be done (GAP, union, intersect, setdiff...)
library(GenomicRanges)

# read the two tables containing positions
dmr  <- read.csv("~/Desktop/DMR.csv")
gens <- read.csv("~/Desktop/gens.csv")

# build the GRanges objects (the columns, in this order: CHR, INI, FI, STRAND, FEATURES...)
rd1 <- with(dmr, GRanges(CHR, IRanges(INI, FI), STRAND, nom=NOM))
rd2 <- with(gens, GRanges(CHR, IRanges(INI, FI), STRAND, nom=NOM, pval=PVAL, diff=DIFF))

# do the matching; maxgap distance accepted
mm <- findOverlaps(rd1, rd2, maxgap=2000)

# add the columns with the feature names from the query
solucio <- data.frame(as.data.frame(rd1[as.matrix(mm)[, 1], ]), as.data.frame(rd2[as.matrix(mm)[, 2], ]))
write.csv(solucio, file="solucio.csv")
```

### Converting seqlevel styles

eg:UCSC to Ensembl, using GenomeInfoDb:

```R
library(GenomeInfoDb)

txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
seqlevels(txdb)
##  [1] "chr2L"     "chr2R"     "chr3L"     "chr3R"     "chr4"      "chrX"
##  [7] "chrU"      "chrM"      "chr2LHet"  "chr2RHet"  "chr3LHet"  "chr3RHet"
## [13] "chrXHet"   "chrYHet"   "chrUextra"

genomeStyles("Drosophila melanogaster")
##    circular   sex  auto  NCBI      UCSC                   Ensembl
## 1     FALSE FALSE  TRUE    2L     chr2L                        2L
## 2     FALSE FALSE  TRUE    2R     chr2R                        2R
## 3     FALSE FALSE  TRUE    3L     chr3L                        3L
## ...

x <- mapSeqlevels(seqlevels(txdb), "Ensembl")
##     chr2L     chr2R     chr3L     chr3R      chr4      chrX      chrU
##      "2L"      "2R"      "3L"      "3R"       "4"       "X"      "Un"
##      chrM  chr2LHet  chr2RHet  chr3LHet  chr3RHet   chrXHet   chrYHet
##      "MT"   "2LHet"   "2Rhet"   "3LHet"   "3RHet"    "Xhet"    "Yhet"
## chrUextra
##        NA

seqlevels(txdb) <- x[seqlevels(txdb)]
seqlevels(txdb)
##  [1] "2L"     "2R"     "3L"     "3R"     "4"      "X"
##  [7] "U"      "M"      "2LHet"  "2RHet"  "3LHet"  "3RHet"
## [13] "XHet"   "YHet"   "Uextra"
```

Or even faster, with the same library but all steps implicit:

```R
seqlevelsStyle(gr) <- "Ensembl"
```

### LiftOver coordinates between different assemblies

Use the ''rtracklayer'' package from ''Bioconductor''. Careful about overlapping regions in the converted file: `export.bw` will fail to export the object, and need to be cleared out first.

```R
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)

# read
x.hg19 <- import.bw("sample1.hg19.bw")
genome(x.hg19) <- genome(BSgenome.Hsapiens.UCSC.hg19)

# convert
x.hg38 <- unlist(liftOver(x.hg19, import.chain("hg19ToHg38.over.chain")))
genome(x.hg38)     <- genome(BSgenome.Hsapiens.UCSC.hg38)
seqlevels(x.hg38)  <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
seqlengths(x.hg38) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)

# drop the overlapping ranges
hits   <- findOverlaps(x.hg38, drop.self=TRUE)
x.hg38 <- x.hg38[-queryHits(hits)]

# and export
export.bw(x.hg38, "sample1.hg38.bw")
```

### Converting GTF to GFF3

Use the ''rtracklayer'' package from ''Bioconductor'' to convert between GTF <--> GFF formats:

```R
library(rtracklayer)
GTF <- "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/rnaseq/ref-chr1/genes+biotypes-chr1.gtf"
GFF <- "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/rnaseq/ref-chr1/genes+biotypes-chr1.gff"
export.gff(import.gff(GTF, format="gtf"), GFF, format="gff3")
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

### Get genome wide distribution of features

Based on the info annotated in TxDb.Hsapiens.UCSC.hg19.knownGene.

Approximation. Projected into one single strand.

```R
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

pr <- local({
  x <- promoters(txdb)
  strand(x) <- "*"
  reduce(x)
})

ex <- local({
  x <- exons(txdb)
  strand(x) <- "*"
  reduce(setdiff(reduce(x), pr))
})

it <- local({
  x <- unlist(intronsByTranscript(txdb))
  strand(x) <- "*"
  reduce(setdiff(reduce(setdiff(reduce(x), ex)), pr))
})

ig <- gaps(reduce(c(ex, it, pr)))
ig <- ig[strand(ig) == "*"]  # take only one strand

# add location info and merge regions
pr$type <- "promoter"
ex$type <- "exon"
it$type <- "intron"
ig$type <- "intergenic"
regions <- c(pr, ex, it, ig)
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
 gset <- getGEO("GSE28133", GSEMatrix=T, AnnotGPL=T)  # be careful if the dataset uses more than 1 platform or is divided in several parts
 ex   <- exprs(gset[[1]])
 colnames(ex) <- c(rep("Ret", 19), rep("RD", 19))

 # platform annotation
 gpl    <- annotation(gset[[1]])  # GPL platform name
 platf  <- getGEO(gpl, AnnotGPL=TRUE)
 ncbifd <- data.frame(attr(dataTable(platf), "table"))

 # average expression of the several probes that interrogate a gene
 ex <- ex[rownames(ex) %in% ncbifd[ncbifd$Gene.symbol != "", "ID"], ]
 ex <- apply(ex, 2, function(x, gens) { tapply(x, gens, mean, na.rm=T) },
       sapply(rownames(ex), function(x) ncbifd[ncbifd$ID == x, "Gene.symbol"]))

 # save tables
 write.csv(ex[, grep("Ret", colnames(ex))], file="Ret.csv")
 write.csv(ex[, grep("RD", colnames(ex))], file="RD.csv")
```

### Download a dataset from Gene Expression Atlas

Interesting, as it includes severals of the major experiments (ENCODE, GTEx, NIH Epigenome Roadmap, etc.).
Web sources: [Gene Expression Atlas (bulk)](https://www.ebi.ac.uk/gxa/home) and [single cell](https://www.ebi.ac.uk/gxa/sc/home).

```R
study_table <- rbind(c("E-MTAB-3871", 64   , "Homo Sapiens", "NIH Roadmap Epigenomics Mapping Consortium"),
                     c("E-MTAB-513" , 16   , "Homo Sapiens", "Illumina Body Map"),
                     c("E-MTAB-5214", 18736, "Homo Sapiens", "GTEx"),
                     c("E-MTAB-4344", 25   , "Homo Sapiens", "Michael Snyder's lab (ENCODE)"),
                     c("E-MTAB-2836", 200  , "Homo Sapiens", "Uhlen's lab 32 tissue samples of 122 individuals"))
colnames(study_table) <- c("ENA id", "Assays", "Species", "Experiment")
study_table <- as.data.frame(study_table)

library(ExpressionAtlas)
library(SummarizedExperiment)
library(httr)
library(XML)

# get all experiments matching our query
response <- httr::GET('https://www.ebi.ac.uk/arrayexpress/xml/v2/experiments?gxa=TRUE&species=Homo%20sapiens')
stopifnot(status_code(response) == 200)
parsedXML   <- xmlParse(content(response))
allExpsNode <- xmlRoot(parsedXML)
allExperiments <- xmlElementsByTagName(allExpsNode, "experiment")
atlasExperimentsList <- lapply(allExperiments, function(experimentNode) {
  c(accession=xmlValue(xmlElementsByTagName(experimentNode, "accession")$accession),
    title    =xmlValue(xmlElementsByTagName(experimentNode, "name")$name),
    species  =xmlValue(xmlElementsByTagName(experimentNode, "organism")$organism),
    expType  =xmlValue(xmlElementsByTagName(experimentNode, "experimenttype")$experimenttype))
})
atlasExperiments <- as.data.frame(do.call(rbind, atlasExperimentsList), row.names=1)

# get data from those experiments (RNA-seq only)
atlasData <- getAtlasData(study_table$`ENA id`)
atlasData <- lapply(atlasData, function(x) x[[1]])
stopifnot(all(sapply(atlasData, function(x) rownames(x) == rownames(atlasData[[1]]))))

# extract samples from relevant tissues only
atlasData <- lapply(atlasData, function(x) {
  x[, grepl("(kidney|renal|liver|pancreas)", colData(x)$"organism_part")]
})
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
 st  <- w/5  # step window
 gen <- "HMGA2"
 s <- getSeq(Hsapiens, "chr12", ini, fi, as.character=TRUE)

 # plot the A/T and C/G content (as the average of w bp)
 pdf(paste(gen, ".pdf", sep=""), width=12)

 f <- function(i, x, s, w) { length(gregexpr(x, substr(s, i-w, i+w))[[1]]) / (w * 2) }
 A=sapply(seq(from=w, to=nchar(s) - w, by=st), f, "A", s, w)
 T=sapply(seq(from=w, to=nchar(s) - w, by=st), f, "T", s, w)
 C=sapply(seq(from=w, to=nchar(s) - w, by=st), f, "C", s, w)
 G=sapply(seq(from=w, to=nchar(s) - w, by=st), f, "G", s, w)

 # A/T content
 df <- data.frame(x=rep(1:length(A), 2), y=c(A, T), class=c(rep("A", length(A)), rep("T", length(T))))
 p <- ggplot(df, aes(x=x, y=y, colour=class)) + geom_line(size=1.05) + ggtitle(gen) +
 #   scale_x_discrete(breaks=df$x, labels=as.character(ini + df$x * w)) +
    theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))
 print(p)

 # C/G content
 df <- data.frame(x=rep(1:length(C), 2), y=c(C, G), class=c(rep("C", length(C)), rep("G", length(G))))
 p <- ggplot(df, aes(x=x, y=y, colour=class)) + geom_line(size=1.05) + ggtitle(gen) +
 #   scale_x_discrete(breaks=df$x, labels=as.character(ini + df$x * w)) +
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
geneids <- read.table("list.csv", head=F, row.names=NULL, stringsAsFactors=F)[, 1]
geneids <- geneids[!is.na(geneids)]
genes   <- unlist(mget(geneids, ifnotfound=NA, revmap(org.Hs.egSYMBOL)))
#be careful here if one Symbol which is returned is associated to more than one entrez id
#there will be entrez ids added to the list of genes and the names of the genes vector are changed
#by adding 1, 2 to the name! Check for your symbol and correct. An easy solution for ensembl ids (with length 11) is
#genes <-genes[!duplicated(substr(names(genes), 1, 11))]
#names(genes) <- substr(names(genes), 1, 11)
univ    <- Lkeys(org.Hs.egGO)
param2  <- new("GOHyperGParams", geneIds=genes, universeGeneIds=univ, annotation="org.Hs.eg.db", ontology="BP")
hyp     <- hyperGTest(param2)

# take the categories with +2 genes
result<-summary(hyp, categorySize=2)
result<-data.frame(result, FDR=p.adjust(result$Pvalue, "fdr"))

# complement the GO term with the genes from the list it contains
result$genes <- sapply(as.list(geneIdUniverse(hyp)[result$GOBPID]), function(x) {
   paste(names(genes)[genes %in% x], collapse="; ")
})
```

Also, gene SYMBOLS can be retrieved from org.Hs.eg.db with the entrez gene Ids:

```R
result$genes <- sapply(as.list(geneIdUniverse(hyp)[result$GOBPID]), function(x) {
   # get the entrezid of the genes in the GO category which are also in our list
   gene_symbols <- select(org.Hs.eg.db, columns="SYMBOL", keys=x[x %in% entrezid], keytype="ENTREZID")$SYMBOL
   paste(sort(unique(gene_symbols)), collapse="; ") # gene symbol
})
```

Alternatively genes names can be obtained from the SQL database (in case the previous step doesn't work):

```R
library(DBI)
library(RSQLite)

drv <- dbDriver("kQLite")
con <- dbConnect(drv, dbname="/usr/local/lib/R/site-library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite")

genes <- apply(result, 1, function(x) {
   res <- dbGetQuery(con, paste("select distinct a.symbol",
              "  from gene_info a, go_bp_all b",
              " where a._id = b._id",
              "   and go_id = ?"),
          data.frame(go_id=x["GOBPID"]))
  res <- intersect(unique(geneids), res[, 1])
  return(paste(res, collapse=";"))
})

# save the results
write.csv(data.frame(result, genes=genes), file="llista_verda.GO.csv")
```

#### GSEA of Biological Processes (a parallel version)

A parallel version based on "parallel" (not "multicore", "snow", "foreach") can be prepared to speed up when multiple GO tests have to be performed. Notice the use of the package "parallel", which is the only one allowing parallel access to SQLite databases.

```R
# functions for parallel GO testing
makeGoParms <- function(geneids) {
    genes  <- unlist(mget(geneids, ifnotfound=NA, revmap(org.Mm.egSYMBOL)))
    #be careful here if one Symbol which is returned is associated to more than one entrez id
    # there will be entrez ids added to the list of genes and the names of the genes vector are changed
    # by adding 1, 2 to the name! Check for your symbol and correct. An easy solution for ensembl ids (with length 11) is
    #genes <-genes[!duplicated(substr(names(genes), 1, 11))]
    #names(genes) <- substr(names(genes), 1, 11)
    univ   <- Lkeys(org.Mm.egGO)
    param2 <- new("GOHyperGParams", geneIds=genes, universeGeneIds=univ, annotation="org.Mm.eg.db", ontology="BP")
    list(genes=genes, param2=param2)
}

GoTest <- function(parms) {
    hyp <- hyperGTest(parms$param2)
    result <- summary(hyp, categorySize=2)
    result <- data.frame(result, FDR=p.adjust(result$Pvalue, "fdr"))
    result$genes <- sapply(as.list(geneIdUniverse(hyp)[result$GOBPID]), function(x) {
        paste(names(parms$genes)[parms$genes %in% x], collapse="; ")
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
# since mget is adding 0, 1, 2 to the name we can filter by the ensembl id length (11)
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

First, convert gene ids from entrez to ensembl and gene symbols:

```R
library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Hs.eg.db)
head(gene.df)
##   ENTREZID         ENSEMBL SYMBOL
## 1     4312 ENSG00000196611   MMP1
## 2     8318 ENSG00000093009  CDC45
## 3    10874 ENSG00000109255    NMU
## 4    55143 ENSG00000134690  CDCA8
## 5    55388 ENSG00000065328  MCM10
## 6      991 ENSG00000117399  CDC20
```

Then do GO over-representation test:

```R
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
head(ego)
```

and/or GSEA analysis:

```R
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
```

How to perform other analysis with `clusterProfiler` using different databases, including KEGG and MsigDB, can be found in the [clusterProfiler book](https://yulab-smu.github.io/clusterProfiler-book/index.html).

#### GO functional analysis with clusterProfiler and DAVID

Use visualization functions provided by ''clusterProfiler'' to plot GO enrichment analysis calculated through DAVID. ''clusterProfiler'' also allows calculating the GO analysis through up-to-date R packages instead of using the old and unmantained DAVID.

```R
require(DOSE)
require(clusterProfiler)
require(biomaRt)

##
## GO ENRICHMENT + PLOTS
##
mart <- useDataset("mmusculus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))
GOenrich <- function(g) {
     g <- getBM(attributes="ensembl_gene_id", filters="external_gene_name", values=g, mart=mart)
     x <- enrichDAVID(g$ensembl_gene_id, idType="ENSEMBL_GENE_ID", listType="Gene",
                      annotation="GOTERM_BP_ALL", pvalueCutoff=0.01)
     x
}
david1 <- GOenrich(genes1)  # official gene symbols
david2 <- GOenrich(genes2)  # official gene symbols

# some plots
barplot(david1)
cnetplot(david1)

##
## Compare the GO enrichment of the 2 gene lists
##
GOcompare <- function(g) {
    g <- lapply(g, function(g) {
        cat("connecting to biomart...\n")
        x <- getBM(attributes="entrezgene", filters="external_gene_name", values=g, mart=mart)
        x$entrezgene
    })
    cat("connecting to david...\n")
#   x <- compareCluster(g, fun="enrichDAVID", idType="ENTREZ_GENE_ID", listType="Gene",
#                       annotation="GOTERM_BP_ALL", pvalueCutoff=0.01, minGSSize=5)
    x <- compareCluster(g, fun="enrichGO", organism="mouse", ont="BP")
    x
}

plot(GOcompare(list(genes1, genes2)), showCategory=25)
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

See [gist](https://gist.github.com/ssayols/296a1de802440f2db6c065f968aecd24).

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
goterms <- goterms[order(goterms$Pvalue), ]
goterms <- goterms[1:min(c(MAXTERMS, nrow(goterms))), ]
semdata <- godata("org.Hs.eg.db", ont="BP")
suppressMessages(clusterExport(cl, c("goterms", "semdata")))
drop <- parLapply(cl, nrow(goterms):2, function(i, goterms) {
    semsim <- lapply((i-1):1, function(j) {
        goSim(goterms$GOBPID[i], goterms$GOBPID[j], semData=semdata, measure="Rel")
    })
    # drop if there was a similar term with lower FDR
    any(semsim > .7, na.rm=T)
})

drop <- c(FALSE, rev(unlist(drop)))
goterms.reduced <- goterms[!drop, ]  # take just the columns you want...

stopCluster(cl)
```

Since it's an extremely CPU intensive and parallelizable problem, it's advisable to use parLapply instead of the first lapply loop and set a maximum of the number of terms to be processed (by min(c(MAXTERMS, nrow(cl):2)) in the lapply loop).

#### Summarizing by scoring frequencies of words

An original way to visually summarize GO terms by displaying variable sized words depending on their frequency:

```R
library(wordcloud)
x <- read.csv("categorical_onlyE_apclusters_GO.csv")
x[1:5, ]
#      GOBPID       Pvalue OddsRatio   ExpCount Count Size                                                                 Term
#1 GO:0001112 5.290072e-06  143.6133 0.03804965     3    7 DNA-templated open complex formation
#2 GO:0001113 5.290072e-06  143.6133 0.03804965     3    7 transcriptional open complex formation
#3 GO:0001120 5.290072e-06  143.6133 0.03804965     3    7 protein-DNA complex remodeling
#4 GO:0034367 5.290072e-06  143.6133 0.03804965     3    7 macromolecular complex remodeling
#5 GO:0001109 8.431172e-06  114.8812 0.04348532     3    8 promoter clearance during transcription
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
v <- sort(rowSums(m), decreasing=TRUE)
d <- data.frame(word = names(v), freq=v)
wordcloud(d$word, d$freq, min.freq=2, colors=brewer.pal(8, "Dark2"))
```

#### Summarizing by clustering similarity scores between terms

First one needs to calculate the similarity scores between terms, and cluster them based on the distance.

```R
# recalculate the semantic similarity scores of the reduced terms for the core proteome
library(parallel)  # create your cluster of workers
library(GOSemSim)

semdata <- godata("org.Hs.eg.db", ont="BP")
scores <- {
    x <- goterms$GOBPID[order(goterms$Pvalue)][1:min(c(MAXTERMS, length(goterms$GOBPID)))]
    parSapply(cl, 1:length(x), function(i, x, semdata) {
        c(rep(NA, i-1),
            sapply(i:length(x), function(j, i) {
                goSim(x[i], x[j], semData=semdata, measure="Wang")
            }, i)
        )
    }, x, semdata)
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
df2 <- df[df$clrep == df$term, ]

p <- ggplot(df, aes(x=V1, y=V2, label=clrep)) +
        geom_point(aes(color=clrep, size=size), alpha=.5) +
        geom_label_repel(data=df2, aes(x=V1, y=V2, label=clrep, color=clrep), box.padding=unit(2, "lines")) +
        scale_color_manual("cluster", values=pal2, guide=F) +
        scale_size_continuous(guide=F, range=c(0, 25)) +
        scale_x_continuous(name="") +
        scale_y_continuous(name="") +
        theme_minimal() +
        theme(axis.text.x=element_blank(), axis.text.y=element_blank())
print(p)
```

### Precomputing tables of semantic similarities between GO terms

Pre-render the table of all GO terms by semantic similarity. This is a one-time and extremely time consuming task that is worth having it pre-calculated before reducing any list of GO terms.

See [gist](https://gist.github.com/ssayols/296a1de802440f2db6c065f968aecd24).

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

# get all GO terms from the annotation 
goterms <- Rkeys(org.Dm.egGO)
semdata <- godata("org.Hs.eg.db", ont="BP")
cat(length(goterms), file=LOG, append=F, fill=T)
suppressMessages(clusterExport(cl, c("goterms", "semdata")))

# calculate the Wang similarity scores
# for some reason, parallelizing the inner loop is ~30% faster
reduced <- sapply(1:length(goterms), function(i) {
    cat(".", file=LOG, append=T)
    c(rep(NA, i-1),
        parSapply(cl, i:length(goterms), function(j, i) {
            goSim(goterms[i], goterms[j], semData=semdata, measure="Wang")
        }, i)
    )
})

# save result
colnames(reduced) <- rownames(reduced) <- goterms
write.csv(reduced, file="semsim.Dm.csv")

stopCluster(cl)
```

***Note***: the Rcpp code is fast and efficient, and tables can be computed on the way. See gist, or basically this code:

```R
library(GOSemSim)

ORGDBs <- c("org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db", "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db", "org.Sc.sgd.db", "org.Ss.eg.db", "org.Xl.eg.db")
ONTs <- c("BP", "MF", "CC")
METHODs <- c("Resnik", "Lin", "Rel", "Jiang", "Wang")

tables <- lapply(ORGDBs, function(ORGDB) {

  if(!require(ORGDB, char=TRUE)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(ORGDB)
    require(ORGDB, char=TRUE)
  }
  orgDb <- eval(parse(text=ORGDB))

  lapply(ONTs, function(ONT) {
    lapply(METHODs, function(METHOD) {
      goAnno  <- suppressMessages(
        select(orgDb,
               keys=keys(orgDb, keytype = "ENTREZID"),
               keytype="ENTREZID",
               columns=c("GO", "ONTOLOGY"))
      )
      goAnno  <- goAnno[!is.na(goAnno$GO), ]
      goterms <- unique(goAnno$GO[goAnno$ONTOLOGY == ONT])
      semdata <- godata(ORGDB, ont=ONT)
      goSim(goterms, goterms, semData=semdata, measure=METHOD)
    })
  })
})
```

### Basic matrix normalization

#### Quantile normalization

* sort the values per sample (column)
* calculate the median per row
* rank the values
* linear interpolation using a linear model y ~ x

```R
 f.qn <- function(x) {
   xm <- apply(x, 2, sort)
   xm <- apply(xm, 1, median, na.rm=T)
   xr <- c(apply(x, 2, rank))
   return(array(approx(1:nrow(x), xm, xr)$y, dim(x), dimnames(x)))
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
    xm <- apply(x, 1, median, na.rm=T)
    parApply(cl, x, 2, function(x, xm) {
        smooth.fit <- fitted.values(loess(x ~ xm))
        dev <- smooth.fit - xm
        x - dev
    }, xm)
}
stopCluster(cl)
```

Alternatively, call `affy::normalize.loess`. It probably does the real loess normalization right.

**NOTE:** `GEOquery::getGEO()` most likely come already normalized (under responsability of the submitter):

```
Value measurements for each Sample within a GEO DataSet (GDS) are assumed to be 
calculated in an equivalent manner, that is, considerations such as background 
processing and normalization are consistent across the dataset. Information 
reflecting experimental design is provided through GDS subsets.
```

Also, from GEO submitter's spreadsheet info, we get that:

* **Matrix table:** The matrix table is a spreadsheet containing the final, normalized values that are comparable across rows and Samples, and preferably processed as described in any accompanying manuscript.
* **Raw data files:** In addition to the normalized data provided in the Matrix table, submitters are required to provide raw data, usually in the form of supplementary raw data files. This facilitates the unambiguous interpretation of the data and potential verification of the conclusions as described in the MIAME and MINSEQE standards.

The **Matrix table** is downloaded with `GEOquery::getGEO()`, while the **Raw data files:** are downloaded from the web interface, with the `Supplementary file` link at the bottom (usually a tar of rawdata files).

### Calculate the reverse-complimentari of a sequence

```R
RC <- function(s) {
  chartr("ACTG", "TGAC", paste(rev(unlist(strsplit(s, ""))), collapse=""))
}
```

### Read a fasta file

Can be done with the ShortRead package:

```R
library("ShortRead")
fasta <- readFasta("file.fa")
for(i in 1:length(fasta)) {
    se <- as.character(sread(fasta[i]))  # sequence
    id <- as.character(id(fasta[i]))  # id
}
```

Or with the Biostrings package:

```R
library("Biostrings")
fasta <- readAAStringSet("file.fa")     # for peptides. Try also readDNA or readRNA
for(i in 1:length(fasta)) {
    se <- as.character(fasta[i])  # sequence
    id <- names(as.character(fasta[i]))  # id
}
```

### Get transcripts from a Bioconductor's AnnotationDb

```R
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
columns(txdb)
x <- select(txdb, columns=c("TXNAME", "CDSCHROM", "TXSTART", "TXEND", "TXSTRAND"),
            keys=keys(txdb, "TXID"), keytype=c("TXID"))
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
### Convert BAM to BigWig

This method doesn't include normalization to total coverage.

```R
rtracklayer::export.bw(
  GenomicAlignments::coverage(                        # convert to coverages
    GenomicAlignments::readGAlignments("sample.bam")  # read in BAM file (use readGAlignmentPairs for paired-end files)
  ), con="sample.bw"                                  # export to bigwig
)
```

### Shannon index for the nt diversity of a sequence
In the Shannon index, p is the proportion (n/N) of individuals of one 
particular species found (n) divided by the total number of individuals found 
(N), and s is the number of species.

It is calculated as:

$$ Shannon Index (H) = - \sum_{i=1}^{s} p_{i} \ln p_{i} $$

It can be easilly adjusted to sequence analysis, by assuming the alphabet $s=\{A, C, T, G\}$.

```R
shannon <- function(x, s=c("A", "C", "G", "T")) {
  x <- unlist(strsplit(x, ""))
  -sum(sapply(s, function(s) { p <- sum(x == s) / length(x); p * log(p) }), na.rm=TRUE)
}
```

## Statistical analysis

### Principal Component Analysis (PCA)

```R
library(scatterplot3d)

# The teapot data is from https://www.khanacademy.org/computer-programming/3d-teapot/971436783
teapot <- read.csv("data/teapot.csv", header = FALSE)
dat <- matrix(unlist(teapot), ncol = 3, byrow = TRUE)
head(dat)
##      [,1] [,2] [,3]
## [1,] 40.6 28.3 -1.1
## [2,] 40.1 30.4 -1.1
## [3,] 40.7 31.1 -1.1
## [4,] 42.0 30.4 -1.1
## [5,] 43.5 28.3 -1.1
## [6,] 37.5 28.3 14.5
 
# teapot in 3D
scatterplot3d(dat[, 1], dat[, 2], dat[, 3], highlight.3d=TRUE, angle=75,
              pch=19, lwd=15, xlab="", ylab="", zlab="", main="teapot in 3D")

# PCA
eigenvec <- eigen(cov(dat))$vectors  # equivalent to prcomp(dat)$rotation
PCA_2 <- dat %*% eigenvec[ , 1:2]  # take only the two eigenvectors with biggest eigenvalues
plot(PCA_2, ylim=c(50, -55), pch=19, lwd=35, col="blue", xlab="", ylab="", main="teapot in 2D")

# Variance explained by each component
eigenval <- eigen(cov(dat))$values   # equivalent to prcomp(dat)$sdev
round(cumsum(eigenval)/sum(eigenval) * 100, 2) # cumulative percentage of retrieved information
## [1]  61.18  83.49 100.00
```

### High Level PCA

```R
# Verify variance is uniform
plot(apply(x, 1, var))

# and scale otherwise
x <- data.frame(scale(x))

#calculate the PCA
pc <- princomp(x)
plot(pc, type='l') # use elbow rule to select how many components

# See which components dominate. What are the loadings?
summary(pc)  # components summary
loadings(pc) # loadings of original variables onto components

# Get principal component vectors using prcomp instead of princomp
# samples as rows & measurements/variables in columns
pc <- prcomp(x)

# Plot first 2 PC
plot(data.frame(pc$x[, 1:2]), pch=16, col=rgb(0, 0, 0, 0.5))
```

### Project a new vector onto PCA space

As described [here](https://stats.stackexchange.com/questions/2592/how-to-project-a-new-vector-onto-pca-space), Either use

```R
predict.prcomp
getS3method("predict", "prcomp")
```

Which is basically:

```R
# perform principal components analysis
pca <- prcomp(data) 

# project new data onto the PCA space (newdata is a matrix with 1 row and 'n' cols
scale(newdata, pca$center, pca$scale) %*% pca$rotation
```

### Multidimensional Scaling (MDS)

```R
# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

d   <- dist(mydata) # euclidean distances between the rows
fit <- cmdscale(d, eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
x <- fit$points[, 1]
y <- fit$points[, 2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="n")
text(x, y, labels=row.names(mydata), cex=.7)
```

### Affinity Propagation Clustering

Autoclustering method implemented as described in the paper by Brendan J. Frey and Delbert Dueck (2007). Clustering by passing messages between data points. Science 315:972-977. DOI: 10.1126/science.1136800. The method automatically determines the best number of subgroups in a group, and optimizes the distances of its members from the centroid.

```R
require(apcluster)

##
## k trees defined by affinity propagation, similarites as negative squared Euclidean distances
##
f <- function(s, type, x) {
        h <- apcluster(s)       # determine number of clusters
        a <- aggExCluster(s, h)  # merge clusters

        plot(h, s, main=paste("similarities", type, "matrix"))      # plot similarities heatmap labeling the found clusters
        plot(a, main=paste("similarities", type, "matrix"))        # plot the cluster
        # plot the first 2 PCA components and the cluster grouping
        for(i in length(h@clusters):1) {
                plot(a, x, k=i, main=paste("similarities", type, ", ", i, "clusters")) }

        # write csv with cluster samples
        r <- data.frame()
        for(i in 1:length(h@clusters)) {
                r <- rbind(r, data.frame(sample=names(h@clusters[[i]]), cluster=i)) }
        write.csv(r, file=paste("lung.AP.", type, ".csv", sep=""))
}

# for graphical representation, we calculate a PCA to get x, y coordinate based on the first 2 components
pca <- prcomp(t(betas))
x   <- t(rbind(pca$x[, 1], pca$x[, 2]))

# distnaces can be calculated as correlations, euclidean or manhattan
f((as.matrix(as.dist(cor(betas)))*100)^2, "correlation", x)
#f(negDistMat(t(betas), r=2), "euclidean", x)
#f(negDistMat(t(betas), r=2, method="manhattan"), "manhattan", x)
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
    plot(x$x[, 1], x$x[, 2], type="n")
    points(x$x[, 1], x$x[, 2], col=paste0(brewer.pal(12, "Set1")[as.factor(k$cl)], "80"), pch=20, cex=.5)
}
```

### Evaluate cluster strength

Evaluate cluster strength by calculating p-values for hierarchical clustering via multiscale bootstrap resampling:

```R
 require("pvclust")
 require("snow")

 load("data.RData")

 cl <- makeCluster(12, type="SOCK")

 result <- parPvclust(cl, data, nboot=10000) #method.dist="euclidean", method.hclust="complete", nboot=1000, init.rand=F)
 stopCluster(cl)
```

### Evaluate cluster similarity

Evaluate how similar are two different clusterings of the same data.
This can be done using the [Rand index](https://en.wikipedia.org/wiki/Rand_index), which is defined as the frequency of occurrence of agreements over the total pairs: 

Ri = (TP + TN) / (TP + TN + FP + FN)

It is described in this [blog post](https://davetang.org/muse/2017/09/21/the-rand-index/) and implemented in the [fossil package](https://cran.r-project.org/web/packages/fossil/index.html) as:

```R
rand.index <- function (group1, group2) {
    # For each clustering, create a matrix of all possibles pairs.
    # A 0 means both are in the same cluster, >0 otherwise.
    x <- abs(sapply(group1, function(x) x - group1))
    x[x > 1] <- 1
    y <- abs(sapply(group2, function(x) x - group2))
    y[y > 1] <- 1
    
    # sum to get all the disagreements (divide by two because the pairs are counted twice)
    sg <- sum(abs(x - y))/2
    
    # get the total number of pairs
    bc <- choose(dim(x)[1], 2)
    
    # calculate the rand index
    ri <- 1 - sg/bc
    return(ri)
}
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

 pred <- prediction(x$predicted, x$actual)
 perf <- performance(pred, "tpr", "fpr")
 plot(perf, colorize=T, main="ROC plot")
```

### Impute missing values

Impute missing values in a table using the method of the K-Nearest Neighbours and the ''imputation'' package

```R
 library(imputation)

 k <- cv.kNNImpute(data)$k       #cross validate by artificially erasing data and get the number of neighbours
 impdata <- kNNImpute(data, k)$x  #run with the previously obtained number of neighbours
```

Also, other methods like Random Forests have a fancy way to impute missing values.

### Geometric mean

```R
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
```

### Power analysis

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

### Power and sample size calculation for survival analysis

The calculations follow the method described in Rosner B. (2006), which was proposed on Freedman, L.S. (1982), implemented in the R package powerSurvEpi:

```R
library(powerSurvEpi)
library(survival)

## nomes del training (29 mostres HR=.03, mirar la KM)
dat.clin <- read.csv("../RSF.nosurgery/analisi.rsf.v3.MI.csv")
dat.clin$days  <- round(dat.clin$days / 30)
dat.clin$days <- ifelse(dat.clin$days > 60, 60, dat.clin$days)
dat.clin$grup <- ifelse(dat.clin$grup == "baix", "C", "E")
ssizeCT(formula = Surv(days, status) ~ grup, dat = dat.clin,
        power = 0.8, k = 0.35, RR = 0.03, alpha = 0.05)
powerCT(formula = Surv(days, status) ~ grup, dat = dat.clin,
        nE = 7, nC = 22, RR = 0.03, alpha = 0.05)
```

### Calculate the center of a 2d-distribution

```R
library(KernSmooth)

optimalBandwidth <- function(x) { # return a more or less useful bandwith (from densCols)
    bandwidth <- diff(apply(x, 2, quantile, probs=c(0.05, 0.95), na.rm=T, names=F))/25
    bandwidth[bandwidth == 0] <- 1
    bandwidth
}

xy <- xy.coords(x, y)
select <- is.finite(xy$x) & is.finite(xy$y)
xy <- cbind(xy$x, xy$y)[select, ]
dens <- bkde2D(xy, bandwidth=optimalBandwidth(xy))
maxpos <- which(dens$fhat == max(dens$fhat), arr.ind=TRUE)
points(dens$x1[maxpos[1]], dens$x2[maxpos[2]], cex=3, pch="+")
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

### Split data for CV

Taken from [here](https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation).

```R
#Randomly shuffle the data
yourData<-yourData[sample(nrow(yourData)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)

#Perform 10 fold cross validation
for(i in 1:10){
    #Segement your data by fold using the which() function 
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- yourData[testIndexes, ]
    trainData <- yourData[-testIndexes, ]
    #Use the test and train data partitions however you desire...
}
```

### Compute PI-value based on FC and p-value 

PI-value is a score suggested in REF, which combines log2FC and FDR to assess significance.
It basically transforms the pvalue (-log10) by penalizing or enhancing it based on the log2FC.
Check the publication and supplemetary materials for a more detailed explanation.
REF: https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btr671

```R
# @param pval a vector of p-values
# @param log2fc a vector of log2 tranformed fc
PIscore <- function(pval, log2fc) {
  list(score=log2fc * -log10(pval),   # the score
       pval =pval ^ log2fc)           # the transformed pval derived from the new score
}
```

## Other stuff that doesn't fit into any other category

### Supercomputing

Make use of several machines cooperating together. Data may reside in the main node which is then shared via sockets to the slaves:

```R
primary <- '192.168.1.235'
machineAddresses <- list(
  list(host=primary, user='johnmount', ncore=4),
  list(host='192.168.1.70', user='johnmount', ncore=4)
)

spec <- lapply(machineAddresses, function(machine) {
  rep(list(list(host=machine$host, user=machine$user)), machine$ncore)
})
spec <- unlist(spec, recursive=FALSE)

parallelCluster <- parallel::makeCluster(type='PSOCK', master=primary, spec=spec)
print(parallelCluster)
## socket cluster with 8 nodes on hosts
##                   ‘192.168.1.235’, ‘192.168.1.70’
```

### Parse arguments

This is just a suggestion on how to deal with input parms:

```R
parseArgs <- function(args, string, default=NULL, as.is=TRUE) {
  if(length(i <- grep(string, args, fixed=T)) == 1) {
    x <- gsub(string, "", args[i])
    return(if(as.is) x else eval(parse(text=x)))
  }
  if(!is.null(default)) default else NA
}

args   <- commandArgs(T)
input  <- parseArgs(args, "--input=")
output <- parseArgs(args, "--output=")
cutoff <- parseArgs(args, "--cutoff=")
color  <- parseArgs(args, "--color=" , "blue")
lwd    <- parseArgs(args, "--lwd=" , 1, as.is=FALSE)
smooth_method <- parseArgs(args, "--smooth=", "none")
smoothness    <- parseArgs(args, "--smoothness=" , 2/3, as.is=FALSE)

runstr <- paste("Call with: Rscript plotProfile.R <arguments>",
                "Mandatory arguments:",
                "  --input=<matrix_A.txt>",
                "  --output=<matrix_A.pdf>",
                "  --cutoff=<10>",
                "Optional arguments:",
                "  --color=[blue]     : line color (hex code or color name)",
                "  --lwd=[1]          : line width. Range [0..Inf], defaults to 1",
                "  --smooth=[method]  : draw smoothed line instead, using [method]. Supported methods: none, lowess, tukey",
                "                       see ?lowess and ?smooth in R for details on this methods. Defaults to none"
                "  --smoothness=[2/3] : smoothness parameter. Defaults to 2/3 for method lowess",
                sep="\n")
if(any(grepl("^-h|^--help", args))) stop(runstr)
if(is.na(input )) stop(runstr)
if(is.na(output)) stop(runstr)
if(is.na(cutoff)) stop(runstr)
```

### Write data into an Excel shit

Exporting dataframes into Excel shits can be done through the ''WriteXLS'' package. Suports old and new Excel format, which is decided by the function based on the extension of the putput file (.xls or .xlsx)

```R
require("WriteXLS")

df1 <- as.data.frame(matrix(rnorm(100  ), ncol=10))
df2 <- as.data.frame(matrix(rnorm(1000 ), ncol=10))
df3 <- as.data.frame(matrix(rnorm(10000), ncol=10))

out <- list(df1, df2, df3)
WriteXLS("out", ExcelFileName="R.xls", SheetNames=c("dataframe1", "dataframe2", "dataframe3"))
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

### [[SQLite]] database access

Access a local SQLite database using the ''DBI'' and ''SQLite'' packages:

```R
 library("DBI")
 library("RSQLite")

 ## database connection
 drv <- dbDriver("SQLite")
 con <- dbConnect(drv, dbname="knowngenes.db")

 ## database access
 obtenirRegio <- function(chr, pos, strand) {

     chr    <- paste("chr", chr, sep="")
     strand <- ifelse(strand == "F", "+", "-")

     query <- paste("SELECT name, txStart, txEnd from genes",
                    " WHERE chrom=?",
                    "   AND ? BETWEEN txStart-1500 AND txEnd+1500",
                    "   AND strand != ?")

     regio <- dbGetQuery(con, query, data.frame(chr, pos, strand))

     return(c(paste(regio$name, collapse=";"), paste(regio$txStart, collapse=";"), paste(regio$txEnd, collapse=";")))
 }


 regio <- apply(fData450k[, c("chr", "pos", "strand")], 1, obtenirRegio)
```

### Other databases

MySQL, for intsance. All drivers use a similar (same?) API:

```R
library(RMySQL)
mydb <- dbConnect(MySQL(), 
                  user='reader',
                  password='reader',
                  dbname='homo_sapiens_core_70_37',
                  host='vieciaepg.eu.boehringer.com')
dbListTables(mydb)
x <- dbSendQuery(mydb, 'select * from gene;')
x <- fetch(x, n=5)
# or
x <- dbGetQuery(mydb, 'select * from gene;')
```

### In memory database

Put a data.frame in a in-memory SQLite database, and write queries like a pro!

```R
library(DBI)
db <- dbConnect(RSQLite::SQLite(), ":memory:")
dbWriteTable(db, "mtcars", mtcars)
head(dbReadTable(db, "mtcars"))
#    mpg cyl disp  hp drat    wt  qsec vs am gear carb
# 1 21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
# 2 21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
# 3 22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
# 4 21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
# 5 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
# 6 18.1   6  225 105 2.76 3.460 20.22  1  0    3    1
dbGetQuery(db, "SELECT COUNT(*) FROM mtcars WHERE mpg < 20")
#   COUNT(*)
# 1       18
dbDisconnect(db)
```

NOTE: for very big data.frames or intensive queries, no index are create with `dbWriteTable()`.

```R
dbGetQuery(db, "EXPLAIN QUERY PLAN SELECT COUNT(*) FROM mtcars WHERE mpg < 20")
#   id parent notused            detail
# 1  3      0       0 SCAN TABLE mtcars   <-- it's going to crawl through whole the table
```

An option is to pre-create the table + indexes with `dbExecute()` prior to dumping the data.frame on the table (preventing the function to create the table for us).

```R
dbCreateTable(db, "mtcars2", mtcars)
dbExecute(db, "CREATE INDEX idx_mpg ON mtcars2(mpg)")
# [1] 0
dbWriteTable(db, "mtcars2", mtcars, append=TRUE)
dbGetQuery(db, "EXPLAIN QUERY PLAN SELECT COUNT(*) FROM mtcars2 WHERE mpg < 20")
#   id parent notused                                                    detail
# 1  3      0       0 SEARCH TABLE mtcars2 USING COVERING INDEX idx_mpg (mpg<?)
```

### Show a progress bar

The trick here is to print the "\r" special character to get to the first column of the current row of the screen, thus overwritting what has just been displayed:

```R

# progress bar
progress.old <<- 0
progress <- function(current) {
  if(current > progress.old) {  # just print every 1% increase, not for every step
    cat(paste(c("\r[", rep("=", current), rep(" ", 100-current), "]", current, "%"), collapse=""))
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
img <- mclapply(1:1000000, function(i) {
    # some very heavy computation
    tmp   <- paste0("./tmp/", as.character(as.hexmode(round(rnorm(1) * 2^28))), ".tiff")
    tiff(tmp)
    # plot whatever
    dev.off()
    tmp
}, mc.cores=CORES)

# plotting part
pdf("plot.pdf")
lapply(img, function(tmp) {
    img <- readTIFF(tmp, native=T)
    plot(1:2, type="n", bty="n", axes=F, xlab="", ylab="")
    rasterImage(img, 1, 1, 2, 2)
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

### Split and Join PDF files

```R
# Load pdftools
library(pdftools)

# extract some pages
pdf_subset('https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf',
  pages = 1:3, output = "subset.pdf")

# Should say 3
pdf_length("subset.pdf")
Similarly pdf_combine() is used to join several pdf files into one.

# Generate another pdf
pdf("test.pdf")
plot(mtcars)
dev.off()

# Combine them with the other one
pdf_combine(c("test.pdf", "subset.pdf"), output = "joined.pdf")

# Should say 4
pdf_length("joined.pdf")
```
