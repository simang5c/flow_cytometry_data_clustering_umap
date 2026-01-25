# Flow Cytometry Data Processing Pipeline (R)

This repository contains an R script for automated processing, gating, clustering, and visualization of flow cytometry data. The workflow includes quality control, variance stabilization, gating, FlowSOM clustering, UMAP/tSNE dimensionality reduction, and interactive visualizations.

## **Requirements**

- R version ≥ 4.3
- Packages:
  - Bioconductor: `flowCore`, `flowWorkspace`, `openCyto`, `flowAI`, `flowVS`, `flowClust`
  - CRAN: `ggcyto`, `ggplot2`, `dplyr`, `tidyr`, `reshape2`, `FlowSOM`, `ConsensusClusterPlus`, `umap`, `Rtsne`, `plotly`, `htmlwidgets`, `pheatmap`, `sessioninfo`

## **Input**

- `sampleplan.csv`: A table with at least two columns:
  - `Samples`: Sample names
  - `Filenames`: Corresponding FCS file paths

- FCS files for all samples.

## **Processing Steps**

1. **Data Import**  
   FCS files are imported into R as a `flowSet` using `flowCore`.

2. **Quality Control**  
   Automatic quality control is performed with `flowAI` to remove low-quality events.  

3. **Variance Stabilization**  
   Fluorescence intensities are variance-stabilized using `flowVS` with manually defined cofactors, applying an arcsinh transformation.

4. **Gating**  
   Sequential gating is applied:
   - Non-debris cells (FSC vs SSC)
   - Singlets (FSC-H vs FSC-A)
   - Viable cells (APC.Cy7 vs SSC)
   - CD45+ leukocytes (FlowClust per-sample clustering)

5. **FlowSOM Clustering**  
   - Marker expression is used for clustering with FlowSOM.
   - Meta-clustering is performed using consensus clustering.

6. **Dimensionality Reduction**  
   - UMAP and tSNE coordinates are calculated for visualization.
   - Interactive plots are saved as HTML files.

7. **Visualization**  
   - Density plots for CD45+ clusters
   - FlowSOM heatmaps
   - UMAP and tSNE plots (static and interactive)

8. **Output Files**  
   - `fs_qc.RDS`: Quality-controlled flowSet  
   - `expr_all_modified.txt`: Final expression matrix  
   - `ALLMARKERS_without_filter.csv`: Raw marker expression  
   - `FlowSOM_Metacluster_Heatmap.pdf`  
   - `UMAP_clusters.html` / `UMAP_samplewise.html`  
   - `TSNE_clusters.html` / `TSNE_samplewise.html`  
   - `cell_counts.txt`: Summary of cells after each gating step  
   - `FlowSOM_formatted_annotations.csv`: Median marker expression per meta-cluster  

## **Usage**

```r
# Clone the repository
git clone https://github.com/simang5c/flowcytometry-pipeline.git
setwd("flowcytometry-pipeline")

# Run the script
source("flow_cytometry_pipeline.R")
