# Trichoderma_atroviride_genomics

In this github we share code for transcriptomic and correlation networks analysis of Trichoderma atroviride RNA-seq data

## Install toolkit from https://www.ncbi.nlm.nih.gov/home/tools/ and download data from NCBI


```{bash}
## download in .sra format

prefetch --option-file SRR_Acc_List.txt

#change to fastq format
fastq-dump *.sra
```

## Eliminate low quality reads and trim adapters with Trimmomatic

Loop for single-end reads

```{bash}
for filename in /path/to/your/fastq/files/*.fastq.gz
do
    # Extract the base name of the file (without the path)
    basename=$(basename "$filename" .fastq.gz)

    echo "Trimming $filename"
    
    # Run Trimmomatic, this code will name the output file as "trimmed" for all fastq output
    java -jar /path/to/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 "$filename" "trimmed-${basename}.fastq.gz" ILLUMINACLIP:/path/to/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    echo "Finished $filename"
done
```

Loop for paired-end reads

```{bash}
for forward_file in /path/to/your/fastq/files/*_1.fq.gz
do
    # Extract the base name (without _1.fastq.gz)
    basename=$(basename "$forward_file" _1.fq.gz)
    
    # Construct the reverse file name
    reverse_file="/path/to/your/fastq/files/${basename}_2.fq.gz"
    
    echo "Trimming paired-end files: $forward_file and $reverse_file"
    
    # Run Trimmomatic for paired-end reads
    java -jar /path/to/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
        "$forward_file" "$reverse_file" \
        "trimmed-${basename}_1_paired.fastq.gz" "trimmed-${basename}_1_unpaired.fastq.gz" \
        "trimmed-${basename}_2_paired.fastq.gz" "trimmed-${basename}_2_unpaired.fastq.gz" \
        ILLUMINACLIP:/path/to/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    echo "Finished $basename"
done
```

# Mapping reads to the reference genome and quantifying gene expression

Create the index using Hisat2 for the reference genome 

```{bash}
hisat2-build reference_genome.fa reference_index
```

Code for mapping and converting to BAM format for single-end reads
```{bash}
#!/bin/bash

# Path to the directory containing your FASTQ files
FASTQ_DIR="/path/to/your/fastq/files"

# Path to the HISAT2 index (without the .ht2 extension)
HISAT2_INDEX="/path/to/your/hisat2/index"

# Path to the output directory for BAM files
OUTPUT_DIR="/path/to/your/output/directory"

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through all fastq files in the directory
for FASTQ_FILE in $FASTQ_DIR/*.fastq.gz; do
   
    # Extract the base name of the file (without extension)
    BASENAME=$(basename "$FASTQ_FILE" .fastq.gz)

    # Define the output files
    SAM_FILE="$OUTPUT_DIR/${BASENAME}.sam"
    BAM_FILE="$OUTPUT_DIR/${BASENAME}.bam"
    SORTED_BAM="$OUTPUT_DIR/${BASENAME}_sorted.bam"

    # Run HISAT2 alignment for single-end reads
    echo "Aligning $BASENAME to genome (single-end mode)..."
    hisat2 -x $HISAT2_INDEX \
           -U "$FASTQ_FILE" \
           -S "$SAM_FILE"

    # Check if HISAT2 ran successfully
    if [ $? -ne 0 ]; then
        echo "Error: HISAT2 failed for $BASENAME"
        continue
    fi

    # Convert SAM to BAM and sort
    echo "Converting $SAM_FILE to BAM and sorting..."
    samtools view -Sb "$SAM_FILE" | samtools sort -o "$SORTED_BAM"

    # Index the BAM file
    echo "Indexing $SORTED_BAM..."
    samtools index "$SORTED_BAM"

    # Remove the intermediate SAM file to save space
    rm "$SAM_FILE"

    echo "Finished processing $BASENAME. Output saved to $SORTED_BAM"
done

echo "All single-end files processed!"
```

Code for mapping and converting to BAM format for paired-end reads

```{bash}
#!/bin/bash

# Path to the directory containing your FASTQ files
FASTQ_DIR="/path/to/your/fastq/files"

# Path to the HISAT2 index (without the .ht2 extension)
HISAT2_INDEX="/path/to/your/hisat2/index"

# Path to the output directory for BAM files
OUTPUT_DIR="/path/to/your/output/directory"

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Get list of unique sample names 
SAMPLES=$(ls $FASTQ_DIR/trimmed-*_1_paired.fastq.gz | sed 's/_1_paired\.fastq\.gz//' | sed 's/.*\///')

# Loop through each sample
for SAMPLE in $SAMPLES; do
    # Define input files
    FASTQ_R1="${FASTQ_DIR}/${SAMPLE}_1_paired.fastq.gz"
    FASTQ_R2="${FASTQ_DIR}/${SAMPLE}_2_paired.fastq.gz"
    
    # Define output files
    SAM_FILE="$OUTPUT_DIR/${SAMPLE}.sam"
    BAM_FILE="$OUTPUT_DIR/${SAMPLE}.bam"

    echo "Aligning $SAMPLE to genome (paired-end mode)..."
    echo "Input R1: $FASTQ_R1"
    echo "Input R2: $FASTQ_R2"
    
    # Run HISAT2 alignment for paired-end reads
    hisat2 -x $HISAT2_INDEX \
           -1 $FASTQ_R1 \
           -2 $FASTQ_R2 \
           -S $SAM_FILE

    # Check if HISAT2 ran successfully
    if [ $? -ne 0 ]; then
        echo "Error: HISAT2 failed for $SAMPLE"
        continue
    fi

    # Convert SAM to BAM and sort
    echo "Converting $SAM_FILE to BAM and sorting..."
    samtools view -Sb $SAM_FILE | samtools sort -o $BAM_FILE

    # Check if conversion was successful
    if [ $? -ne 0 ]; then
        echo "Error: SAM to BAM conversion failed for $SAMPLE"
        continue
    fi

    # Index the BAM file
    echo "Indexing $BAM_FILE..."
    samtools index $BAM_FILE

    # Remove the intermediate SAM file to save space
    rm $SAM_FILE

    echo "Finished processing $SAMPLE. Output saved to $BAM_FILE"
done

echo "All paired-end files processed!"
```

## Quantify gene expression with featureCounts using the R package Rsubread

The function featureCounts uses the bam or sam files and the GFF3 annotation file

```{r}

library (Rsubread)

#this is an example of how to use the function
# first put the path to the bam/sam file
# then you write the path to the annotation file (gff or gtf or gff3 file)
#    the function featureCounts contains several parameters, for example it can count     #    multimapping reads if we enable "countMultiMappingReads = TRUE",
#    "GTF.featureType" here you specify the feature: gene, mRNA, lncRNA, 5 prime/3 prime  #    UTR, intron, exon, etc. 
#    you can also specify the strand you want to quantify using the option                #    "strandSpecific", 1 = stranded, 2 = reversely stranded or 0 unstranded (default)
#    "fraction" means that multimapping reads are quantified in fractional counts
#    take note that Multimapping counting is a valuable feature for this function
#    "isPairedEnd": TRUE for paired end reads or FALSE for single end reads
#    this process can take several minutes depending on the number of reads, in this case #    it would be fast but if you want to accelerate the process you can specify the number #    of threads 


counts<-featureCounts("/path/to/your/bam/file.bam",
         annot.ext = "/path/to/your/annotation/file.gff",
         isGTFAnnotationFile = TRUE,
         GTF.featureType = "mRNA", # specify the feature type you want to quantify
         GTF.attrType = "ID", # here you specify the gene id 
         countMultiMappingReads = TRUE,
         strandSpecific="0",
         fraction = TRUE,
         nthreads=4,
         allowMultiOverlap = TRUE,
         isPairedEnd = TRUE #specify if your data is paired end or single end
         ) 
```


## WGCNA code for Trichoderma atroviride RNA-seq data

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


