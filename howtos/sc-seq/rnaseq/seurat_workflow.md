## Seurat v5 workflow

### Preamble

```R
library(Seurat)
library(parallel)
library(SingleR)
library(celldex)
#library(SoupX)
#library(DoubletFinder)   # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force=TRUE)
library(ggplot2)

CORES <- 10
PROJECT <- "~/prj/aa_ober_reynolds/RD_DataScience.initiative.oberreynoldsaadataset"
options(mc.cores=CORES)
setwd(PROJECT)
```

### Load MTX file + Seurat workflow + cell type annotation
Process through Suerat and annotate cell types

```R
seu <- mcMap(
  sort(list.files("data", pattern=".mtx.gz$", full.names=TRUE)),
  sort(list.files("data", pattern=".features.tsv.gz$", full.names=TRUE)),
  sort(list.files("data", pattern=".barcodes.tsv.gz$", full.names=TRUE)),
  sub(".+_(.+)_.+", "\\1", sort(list.files("data", pattern=".mtx.gz$", full.names=TRUE))),
  f=function(m, f, b, s) {
    cells <- read.delim(b, head=FALSE)[[1]]
    seu <- ReadMtx(mtx     =m,
                   features=f,
                   cells   =b) |>
           CreateSeuratObject(project=paste0("Ober_Reynolds_", s),
                              min.cells=3,
                              min.features=200,
                              meta.data=data.frame(sample   =rep(s, length(cells)),
                                                   group    =rep(sub("^(..)(\\d)", "\\1", s), length(cells)),
                                                   replicate=rep(sub("^(..)(\\d)", "\\2", s), length(cells)),
                                                   row.names=cells))
    
   # Seurat workflow
   seu <- NormalizeData(seu) 
   seu <- FindVariableFeatures(seu) 
   seu <- ScaleData(seu)
   seu <- RunPCA(seu)
   seu <- RunUMAP(seu, dims=1:10)
   seu <- FindNeighbors(seu, dims=1:10)
   seu <- FindClusters(seu, resolution=0.5)
   
   ## Remove RNA ambient contaminants with SoupX
   #sc <- SoupChannel(as.matrix(seu@assays$RNA$counts), as.matrix(seu@assays$RNA$scale.data))
   #sc <- setClusters(sc, seu$seurat_clusters)
   #sc <- autoEstCont(sc)
   #seu <- SetAssayData(seu, slot="counts", new.data=adjustCounts(sc))
   #
   ### Remove doublets with DoubletFinder
   #nExp_poi <- round(0.5 * nrow(seu@meta.data))
   #seu <- doubletFinder(seu, PCs=1:10, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN=FALSE, sct=FALSE)
   #seu <- subset(seu, DF.classifications == "Singlet")
   
   # annotate cell types
   ref <- celldex::HumanPrimaryCellAtlasData()
   seu$HPCA_label      <- SingleR(test=GetAssayData(seu, slot="data"), ref=ref, labels=ref$label.main)$labels
   seu$HPCA_label_fine <- SingleR(test=GetAssayData(seu, slot="data"), ref=ref, labels=ref$label.fine)$labels
   
   ref <- celldex::BlueprintEncodeData()
   seu$BlueprintEncode_label      <- SingleR(test=GetAssayData(seu, slot="data"), ref=ref, labels=ref$label.main)$labels
   seu$BlueprintEncode_label_fine <- SingleR(test=GetAssayData(seu, slot="data"), ref=ref, labels=ref$label.fine)$labels
   
   seu
  }
)
names(seu) <- sub(".+_(.+)_.+", "\\1", names(seu))
```

### Plot number of cells / cell types detected in each sample

```R
df <- lapply(seu, function(x) table(x$BlueprintEncode_label)) |> reshape2::melt()
df$Var1 <- factor(df$Var1, levels=sort(levels(df$Var1)))
df$L2 <- factor(sub("\\d$", "", df$L1),
                levels=c("AA", "PB", "SD"),
                labels=c("Alopecia Areata", "Healthy (punch biopsy)", "Healthy (derma surgery)"))

ggplot(df, aes(x=Var1, y=value, fill=L2)) +
  stat_summary(fun=sum, geom="bar", position="dodge") +
  scale_fill_manual("", values=c("#F28E2B", "#2EAAAA", "#2EAAFF")) +
  labs(title="Number of cells / cell type", x="", y="Number of cells") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="bottom")
```

### Calculate pseudobluk counts from sc data

```R
seu_pseudo_bulks <- mclapply(seu, function(x) {
  #AggregateExpression(x, group.by="BlueprintEncode_label", scale.factor=1e6)$RNA |>
  PseudobulkExpression(
    object  =x,
    assay   ="RNA",
    normalization.method="LogNormalize",
    scale.factor = 1e6,
    group.by="BlueprintEncode_label",
    method  ="average")$RNA |>
  as.matrix()
})
```

### Plot CD70 average expression in every sample

```R
# prepare data on CD70
df <- lapply(seu_pseudo_bulks, function(x) {
  x <- reshape2::melt(x["CD70", ])
  x$cell_type <- rownames(x)
  x
}) |> reshape2::melt()
df$L2 <- factor(sub("\\d$", "", df$L1),
                levels=c("AA", "PB", "SD"),
                labels=c("Alopecia Areata", "Healthy (punch biopsy)", "Healthy (derma surgery)"))

ggplot(df, aes(x=cell_type, y=value, fill=L2)) +
  stat_summary(fun=mean, geom="bar", position="dodge") +
  scale_fill_manual("", values=c("#F28E2B", "#2EAAAA", "#2EAAFF")) +
  labs(title="CD70 expression", x="", y="average expression (arbitrary units)") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
```
