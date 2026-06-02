# Trichoderma_atroviride_genomics

In this github we share code for transcriptomic and correlation networks analysis of Trichoderma atroviride RNA-seq data

This is a WGCNA tutorial for gene expression data. We will use the R package WGCNA to construct a gene co-expression network and identify modules of co-expressed genes. 
The tutorial assumes that you have already performed gene quantification and have a count matrix ready for analysis.

## We load the required libraries and read the gene quantification table in R 

```{r}
library("limma")
library("edgeR")
library("statmod")
library("DESeq2")
library("WGCNA")

red = read.table("counts.txt")
options(stringsAsFactors = FALSE)

```

# Normalize (VST with DESeq2)

```{r}
dds <- DESeqDataSetFromMatrix(
  countData = round(red),
  colData = data.frame(row.names = colnames(red)),
  design = ~1
)
vsd <- vst(dds, blind = FALSE)
expr <- assay(vsd)
vsd <- vst(dds, blind = FALSE)
expr <- assay(vsd)
```

#  Format for WGCNA
```{r}
datExpr <- t(expr)
```


# Quality control

```{r}
gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
#7. Detection of outliers
sampleTree <- hclust(dist(datExpr), method = "average")

plot(sampleTree, main = "Sample clustering")

# court
clust <- cutreeStatic(sampleTree, cutHeight = 320, minSize = 1)
```


# Remove outliers (only if applicable)
```{r}
keepSamples <- (clust != 3)  # ajustar según resultado
datExpr <- datExpr[keepSamples, ]
```

# Recalculate tree (IMPORTANT)
```{r}
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Clean sample clustering")

#Selecction of soft-threshold

powers <- 1:20

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

par(mfrow = c(1,2))
```

# Scale-free topology
```{r}
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab="Power", ylab="Scale-free R^2",
     type="n")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers, col="red")

abline(h=0.7, col="red")
```

# Mean connectivity
```{r}
plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab="Power", ylab="Mean connectivity",
     type="n")
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels=powers, col="red")
```  



# Network construction

```{r}
net <- blockwiseModules(
  datExpr,
  corType = "pearson",              
  power = 9,                     
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "WGCNA_TOM",
  verbose = 3
)


moduleColors <- labels2colors(net$colors)
```

# Relate modules to treatment

```{r}
MEs0 <- net$MEs
MEs0 <- orderMEs(MEs0)
module_order <- names(MEs0) %>% gsub("ME", "", .)
MEs0$treatment <- rownames(MEs0)
library(tidyr)
library(dplyr)
library(ggplot2)

mME <- MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )
```

# Heatmap
```{r}
ggplot(mME, aes(x = treatment, y = name, fill = value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limits = c(-1,1)
  ) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(
    title = "Module Eigengene Expression",
    y = "Modules",
    fill = "Corr"
  )
```

# Calculate hub genes by connectivity---structural hub

```{r}
colors <- labels2colors(net$colors)
topHubs <- chooseTopHubInEachModule(
  datExpr,
  colors,
  omitColors = "grey",
  type = "signed"
)
print(topHubs)
```


# Finally we export the network to cytoscape
```{r}
TOM = TOMsimilarityFromExpr(datExpr, power = 9, corType = "pearson", networkType = "signed");
probes = colnames(datExpr)


cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("Network_Tatroviride_0.33_edges", ".txt", sep=""),
                               nodeFile = paste("Network_Tatroviride_0.33_nodes", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.33,
                               nodeNames = probes,
                               nodeAttr = moduleColors);

```


