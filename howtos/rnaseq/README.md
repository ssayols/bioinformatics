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
   * [Plots](#plots-1)
     * [Volcano](#volcano)
   * [Other](#other)
     * [RPK, RPKM, TPM](#rpk-rpkm-tpm)

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

