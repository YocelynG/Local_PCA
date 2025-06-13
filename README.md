# 🚀 Local PCA

Scripts to run local PCA using
[lostruct](https://doi.org/10.1534/genetics.118.301747) library

genotype call files (bcf) are separated by chromosome (for example:
chromosome_1.bcf, chromosome_2.bcf, chromosome_3.bcf . . .)

## 🧭 Overview

This pipeline includes the following steps:

- **Step 1:** Library requirements  
- **Step 2:** Genomic window sizes and PCA
- **Step 3:** localPCA

------------------------------------------------------------------------

## 🔧 Library requirements

``` r
install.packages("data.table")
devtools::install_github("petrelharp/local_pca/lostruct")

library(tidyverse)
library(lostruct)
library(SNPRelate)
library(ggplot2)
library(viridis)
library(gridExtra)
```

------------------------------------------------------------------------

## VCF files to window size

For each chromosome, genomic windows of fixed size (window_size = 100
SNPs) are defined using vcf_windower() Later, a genetic PCA on each
window using eigen_windows() is generated, computing pairwise distances
between windows using pc_dist() based on the PCA result. Run MDS on the
non-missing distances

Pop = name of the population file. For example: 1_Manaus.bcf

``` r
process_chromosome <- function(chrom_number, window_size = 100, k_kept = 2) {
  bcf.file <- paste0(chrom_number, "_Pop.bcf")
  sites <- vcf_positions(bcf.file)
  win.fn.snp <- vcf_windower(bcf.file, size = window_size, type = "snp", sites = sites)

  cat("Processing chromosome", chrom_number, "\n")

  system.time(snp.pca <- eigen_windows(win.fn.snp, k = 2))
  system.time(pcdist <- pc_dist(snp.pca))

  na.inds <- is.na(pcdist[, 1])
  if (sum(na.inds) == length(na.inds)) {
    na.inds <- is.na(pcdist[, 2])
  }

  mds <- cmdscale(pcdist[!na.inds, !na.inds], eig = TRUE, k = k_kept)
  mds.coords <- mds$points
  colnames(mds.coords) <- paste0("MDS coordinate ", 1:ncol(mds.coords))

  win.regions <- region(win.fn.snp)()
  win.regions$n <- 1:nrow(win.regions)
  win.regions <- win.regions[!na.inds, ]
  win.regions <- win.regions %>%
    mutate(mid = (start + end) / 2)

  for (i in 1:k_kept) {
    win.regions[[paste0("mds", str_pad(i, 2, pad = "0"))]] <- mds.coords[, i]
  }

  output_file <- paste0("Pop_Chrm", chrom_number, ".txt")
  write.table(win.regions, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Loop across the 19 chromosomes
for (chrom in 1:19) {
  process_chromosome(chrom)
}
```

📊 Tables generated in previous step

``` r
file_list <- list.files(pattern = "*.txt", full.names = T)

plot_list <- list()
# Loop 
for (file in file_list) {
  
  win.regions <- read.table(file, header=TRUE)
  
  p <- print(
    win.regions %>%
        gather(., mds, value, colnames(win.regions)[c(6)]) %>% 
        ggplot(.,aes(x=mid,y=value)) + geom_point(shape=21, colour="black", fill="lightseagreen",size=2) + facet_grid(mds~.,scales = "free") + theme(plot.title = element_text(hjust = 0.5, size=12))
  )
  
  plot_list[[length(plot_list) + 1]] <- p
}

grid.arrange(grobs = plot_list, nrow = 3, ncol = 6)
```
