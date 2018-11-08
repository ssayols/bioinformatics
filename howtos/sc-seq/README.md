# sc-seq

Here we forge the tools to analyze single cell sequencing experiments.

## scRNA-seq: our current pipeline

As the libraries generated so far were always produced with the [Smart-seq2 kit](http://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2639.html), we usually end up generating small datasets (<500 cells) that can (and ideally should be) processed with our standard RNAseq pipeline for bulk data.

This will usually include:
- FastQC
- Mapping to the genome using STAR
- Quantification with featureCounts (Subread)

If other methods based on CEL/Drop/MARS which are more high throughput in terms of cells sequenced, consider doing alignment free quantification such as Salmon or Kallisto.

From here on, take the Rmd file in ./rnaseq/sc_analysis.Rmd to analyze the count data. This step could of course be automatized and run from the pipeline, but so far some things that are not clear yet need to be parametrized. Some of these things include thresholds for discarding cells, parameters to use in the regression model, and marker genes.

The Rmd file use, among others, the following tools and methods:
- QC: the [scater](http://bioconductor.org/packages/release/bioc/html/scater.html) package.
- Normalization: the [scran](http://bioconductor.org/packages/release/bioc/html/scran.html) package.
- Differential expression analysis: the [scde](http://bioconductor.org/packages/release/bioc/html/scde.html) package.
- Trajectory analysis (pseudotime): the [monocle](https://bioconductor.org/packages/release/bioc/html/monocle.html) package.

## Resources
- ./rnaseq/sc_analysis.Rmd is mainly based on the workflow from Lun ATL, McCarthy DJ and Marioni JC. A step-by-step workflow for low-level analysis of single-cell RNA-seq data [version 1; referees: 5 approved with reservations]. F1000Research 2016, 5:2122 (doi: 10.12688/f1000research.9501.1
- QC: the [scater](http://bioconductor.org/packages/release/bioc/html/scater.html) package.
- Normalization: the [scran](http://bioconductor.org/packages/release/bioc/html/scran.html) package.
- Differential expression analysis: the [scde](http://bioconductor.org/packages/release/bioc/html/scde.html) package.
- Trajectory analysis (pseudotime): the [monocle](https://bioconductor.org/packages/release/bioc/html/monocle.html) package.
- A tutorial: https://hemberg-lab.github.io/scRNA.seq.course/index.html
- [Notes of interest](https://hemberg-lab.github.io/scRNA.seq.course/ideal-scrnaseq-pipeline-as-of-oct-2017.html) from the last guys

