# Note:
# This R script performs the initial round of data processing and clustering.
# Subsequent analysis includes multiple iterative rounds of reclustering and backgating
# to further refine population structure and systematically remove problematic or low-quality populations.
# Backgating, as well as manual cluster identification and confirmation, were performed
# in collaboration with an experienced flow cytometry expert.
######################################################################################################

# Load all libraries 
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(flowAI)
library(flowVS)
library(flowClust)

library(ggcyto)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

library(FlowSOM)
library(ConsensusClusterPlus)

library(umap)
library(Rtsne)

library(plotly)
library(htmlwidgets)
library(pheatmap)


# Load the sampleplan 
sampleplan <- read.csv("sampleplan.csv", sep = ",", header = TRUE)

# Load the flowcore data
if (all(file.exists(sampleplan$Filenames))) {
  fs <- read.flowSet(
    files = sampleplan$Filenames,
    path = ".",
    truncate_max_range = FALSE,
    alter.names = TRUE
  )
} else {
  stop("One or more FCS files listed in sampleplan do not exist.")
}


# Update filename 
for (i in seq_along(fs)) {
  # Update filename inside flowFrame
  fs[[i]]@description[["FILENAME"]] <- sampleplan$Samples[i]
}

# Update the flowSet names
sampleNames(fs) <- sampleplan$Samples
# quality check after compensation
fs_qc<-flow_auto_qc(fs,
                    remove_from = "all",
                    mini_report = TRUE,
                    html_report = TRUE)


saveRDS(fs_qc, "fs_qc.RDS")

# markers that will be used for clustering
marker_channels <- c("FITC.A","BV421.A","PE.Cy7.A","APC.A",
                     "BV605.A","PE.A","Alexa.700.A","PerCP.Cy5.5.A")
# check epressions
fcs_1<-exprs(fs_qc[[1]])
fcs_1[, marker_channels]

stats_table <- data.frame(
  Min    = apply(fcs_1, 2, min),
  Mean   = apply(fcs_1, 2, mean),
  Median = apply(fcs_1, 2, median),
  Max    = apply(fcs_1, 2, max)
)

# see the table
stats_table

# The ararcsinh transformation formula is:
#       y=arcsinh(x/c)
  
# x = raw fluorescence intensity
# c = cofactor
# y = transformed intensity

# The goal of the cofactor is to:
# Keep the negative and low-intensity background near zero roughly linear
# Avoid compressing the positive population excessively
# Stabilize variance across events


# see column information
column_informations<-pData(parameters(fs_qc[[1]]))

# manual cofactors for scaling
manualfactors <- c(
  "FITC.A"         = 211,    # CD15
  "BV421.A"        = 840,   # CD11b
  "PE.Cy7.A"       = 167,    # CD8
  "APC.A"          = 100,    # CD3
  "BV605.A"        = 350,   # CD4
  "PE.A"           = 140,    # CD56
  "Alexa.700.A"    = 77,    # CD19
  "PerCP.Cy5.5.A"  = 809    # CD45
)

# variance transformation
fcs_trans <- transFlowVS(fs_qc, channels = names(manualfactors), cofactors = manualfactors)

sampleNames(fcs_trans) <- sampleplan$Samples

# Then create GatingSet
gs <- GatingSet(fcs_trans)

#####################################################################
# GATING
## None Debris (FSC.A vs SSC.A)
#gs_pop_remove(gs, "noneDebris_gate")
fsc_min <- 20000   # forward scatter minimum
ssc_max <- 280000     # side scatter minimum

nonDebris_gate <- rectangleGate(
  filterId = "noneDebris_gate",
  .gate = matrix(c(fsc_min, Inf, 0, ssc_max), ncol = 2,
                 dimnames = list(NULL, c("FSC.A", "SSC.A")))
)

# Apply the gate to all samples in the GatingSet
fs_data <- gs_pop_get_data(gs)

for (i in seq_along(fs_data)) {
  gh <- gs[[i]]
  gs_pop_add(gh, nonDebris_gate, parent = "root", name = "noneDebris_gate")
}

# Recompute populations
recompute(gs)

#  Visualize
p<-autoplot(gs, x = "FSC.A", y = "SSC.A", "noneDebris_gate") +
  geom_gate() +
  ggtitle("Non-debris Cells") +
  theme_bw()

ggsave(
  filename = "Non_debris_cells.png",
  plot = p,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)

#gs_pop_remove(gs, "singlets")
slope <- 0.82         # slope of the line (FSC-H / FSC-A)
width <- 20000        # total perpendicular width in FSC units
intensity_min <- 5000
fr_nd <- gs_pop_get_data(gs, "noneDebris_gate")[[1]]
FSC.A <- exprs(fr_nd)[, "FSC.A"]

# Horizontal range
fsc_min <- intensity_min
fsc_max <- max(FSC.A)

# Unit vector along slope
dx <- 1 / sqrt(1 + slope^2)
dy <- slope / sqrt(1 + slope^2)

# Perpendicular unit vector
px <- -dy
py <- dx

# Perpendicular offset
ox <- px * width / 2
oy <- py * width / 2

# Coordinates of rotated rectangle (slope in middle)
polygon_coords <- matrix(
  c(fsc_min - ox, fsc_min * slope - oy,   # bottom-left
    fsc_max - ox, fsc_max * slope - oy,   # bottom-right
    fsc_max + ox, fsc_max * slope + oy,   # top-right
    fsc_min + ox, fsc_min * slope + oy),  # top-left
  ncol = 2, byrow = TRUE,
  dimnames = list(NULL, c("FSC.A", "FSC.H"))
)

# Create gate
singlet_gate <- polygonGate(filterId = "singlets", .gate = polygon_coords)

# Add to gatingSet
#gs_pop_remove(gs, "singlets")
gs_pop_add(gs, singlet_gate, parent = "noneDebris_gate", name = "singlets")
recompute(gs)

# Plot to confirm
p <- autoplot(gs, x = 'FSC.A', y = 'FSC.H', "singlets", bins = 256)


ggsave(
  filename = "Non_debris_singlets_cells.png",
  plot = p,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)

#gs_pop_remove(gs, "Viable")
fs_data <- gs_pop_get_data(gs, "singlets")
g_manual <- rectangleGate(
  filterId = "Viable",
  .gate = matrix(c(-200, 9000, 0, Inf), ncol = 2, byrow = FALSE,
                 dimnames = list(NULL, c("APC.Cy7.A", "SSC.A"))),
  dimension = c("APC.Cy7.A", "SSC.A")
)
gs_pop_add(gs, g_manual, parent = "singlets", name="Viable")
recompute(gs)
chnl <- "APC.Cy7.A"

# Plot x vs y and overlay the rectangle gate
p<-autoplot(gs, x = "APC.Cy7.A", y = "SSC.A", "Viable") +
  geom_gate(g_manual)

ggsave(
  filename = "Non_debris_singlets_Viable_cells.png",
  plot = p,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)


## CD45+ Cells Gating (PerCP.Cy5.5.A vs SSC.A)

#gs_pop_remove(gs, "CD45pos")
set.seed(123)

# Extract Viable cells from GatingSet
fs_viable <- gs_pop_get_data(gs, "Viable")

# FlowClust 2D clustering per sample
res_list <- lapply(fs_viable, function(fr){
  tryCatch(
    flowClust(fr, varNames = c("PerCP.Cy5.5.A", "SSC.A"), K = 2),
    error = function(e) {
      message("flowClust failed: ", e$message)
      return(NULL)
    }
  )
})

# Process all samples to identify CD45+ and prepare per-sample gates
sample_names <- sampleNames(fs_viable)
num_samples <- length(fs_viable)

# Data frame for summary stats (optional, for validation)
summary_df <- data.frame(
  sample = sample_names,
  mean_pos = NA,
  mean_neg = NA,
  thresh = NA,  # We'll store the midpoint threshold here
  stringsAsFactors = FALSE
)

# Create a list of gates (one per sample)
gate_list <- list()

for (i in 1:num_samples) {
  fit <- res_list[[i]]
  if (is.null(fit)) {
    message("Skipping sample ", sample_names[i], ": flowClust failed.")
    next
  }
  
  fr <- fs_viable[[i]]
  expr_data <- exprs(fr)
  clusters <- factor(fit@label)
  
  df <- data.frame(
    CD45 = expr_data[, "PerCP.Cy5.5.A"],
    SSC = expr_data[, "SSC.A"],
    cluster = clusters
  )
  
  # Identify CD45+ based on higher mean
  cluster_means <- tapply(df$CD45, df$cluster, mean)
  cd45_pos_cluster <- names(which.max(cluster_means))
  
  # Calculate midpoint threshold for manual gate
  mean_pos <- cluster_means[cd45_pos_cluster]
  mean_neg <- cluster_means[setdiff(names(cluster_means), cd45_pos_cluster)]
  thresh <- (mean_neg + mean_pos) / 2
  
  # Create manual rectangle gate for this sample
  g_manual <- rectangleGate("PerCP.Cy5.5.A" = c(thresh, Inf), filterId = "Cd45Positive")
  gate_list[[sample_names[i]]] <- g_manual
  
  # Update summary
  summary_df$mean_pos[i] <- mean_pos
  summary_df$mean_neg[i] <- mean_neg
  summary_df$thresh[i] <- thresh
  
  # Optional: Print for this sample
  cat("\nSample:", sample_names[i], "\n")
  cat("CD45+ mean:", mean_pos, "| CD45- mean:", mean_neg, "| Threshold:", thresh, "\n")
}

# View summary (optional)
cat("\nSummary across samples:\n")
print(summary_df)



for (i in seq_along(gate_list)) {
  g_manual <- gate_list[[i]]
  gh <- gs[[names(gate_list)[i]]]
  
  # Remove existing gate if it exists
  if ("CD45Pos" %in% getNodes(gh)) gs_pop_remove(gh, "CD45Pos")
  
  # Add the manual gate under Viable
  gs_pop_add(gh, g_manual, parent = "Viable", name = "CD45Pos")
}

# Recompute the GatingSet
recompute(gs)

# Checkin Gating using CD45 (PerCP.Cy5.5.A) Density Plots
for (i in 1:num_samples) {
  fit <- res_list[[i]]
  if (is.null(fit)) {
    message("Skipping sample ", sample_names[i], ": flowClust failed.")
    next
  }
  
  fr <- fs_viable[[i]]
  expr_data <- exprs(fr)
  clusters <- factor(fit@label)
  
  df <- data.frame(
    CD45 = expr_data[, "PerCP.Cy5.5.A"],
    SSC = expr_data[, "SSC.A"],
    cluster = clusters
  )
  
  # Density plot
  p_dens <- ggplot(df, aes(x = CD45, color = factor(cluster))) +
    geom_density() +
    labs(title = paste("Density of CD45 by Cluster -", sample_names[i]),
         x = "CD45 (PerCP.Cy5.5.A)", y = "Density") +
    theme_minimal()
  
  # Display the plot 
  ggsave(
  filename = paste0("P", i, ".png"),
  plot = p_dens,
  width = 6,
  height = 4,
  dpi = 300,
  bg = "white"
)


}

# plot the hierarchy
#gs_pop_remove(gs, "nonDebris")
plot(gs)
##########################################################################################################
# Plot the hierarchy
png("hierarchy.png", width = 1800, height = 1200, res = 300)
plot(gs)
dev.off()
##########################################################################################################
# Cell counts after each step
data_counts<-gs_pop_get_stats(gs)
counts<-reshape2::dcast(data_counts, sample~pop)
colnames(counts)[1:dim(counts)[2]]<-c("Samples","After_debris_removal","Singlets","Viable_cells","CD45pos","After_QC")
counts<-counts[,c("Samples","After_QC",colnames(counts)[2:5])]

write.table(counts, "cell_counts.txt",sep = "\t", row.names = F)
####################################################################################################
# Cleaning and extracting
clean_frames <- lapply(seq_along(gs), function(i) {
  gh <- gs[[i]]
  samp <- sampleNames(gs)[i]

  # Extract CD45+ under Viable
  fr_keep <- tryCatch({
    getData(gh, "/noneDebris_gate/singlets/Viable/CD45Pos")
  }, error = function(e) NULL)

  if (is.null(fr_keep) || nrow(fr_keep) == 0) return(NULL)

  # Store sample ID in keyword (metadata)
  keyword(fr_keep)[["SAMPLE_ID"]] <- samp

  fr_keep
})

# Remove empty samples
clean_frames <- clean_frames[!sapply(clean_frames, is.null)]

# Assign flowFrame names using SAMPLE_ID keyword
names(clean_frames) <- sapply(clean_frames, function(fr) keyword(fr, "SAMPLE_ID"))

fs_clean <- flowSet(clean_frames)
gs_clean <- GatingSet(fs_clean)

# creating expr_all with selected markers
expr_list <- lapply(sampleNames(fs_clean), function(sn) {
  fr <- fs_clean[[sn]]
  marker_use <- marker_channels[marker_channels %in% colnames(exprs(fr))]
  out <- exprs(fr)[, marker_use, drop = FALSE]
  rownames(out) <- paste0(sn, "_", seq_len(nrow(out)))
  out
})

expr_all <- do.call(rbind, expr_list)

# creating all markers data
expr_list_all <- lapply(sampleNames(fs_clean), function(sn) {
  fr <- fs_clean[[sn]]

  out <- exprs(fr)                 # take ALL markers
  rownames(out) <- paste0(sn, "_", seq_len(nrow(out)))

  out
})

# All markers table
expr_all_markers_all <- do.call(rbind, expr_list_all)
write.csv(expr_all_markers_all, "ALLMARKERS_without_filter.csv", row.names = TRUE, sep="\t", quote=F)


temp1 <- data.frame(interest_cols = colnames(expr_all))
temp1 <- left_join(temp1, column_informations, by = c("interest_cols" = "name"))
temp1[is.na(temp1)] <- "SSC.A"
expr_all_modified <- expr_all
expr_all_modified<-as.data.frame(expr_all_modified)
colnames(expr_all_modified) <- temp1$desc

# Keep sample ID in a separate column
expr_all_modified<-as.data.frame(expr_all_modified)
expr_all_modified["SampleID"] <- gsub("\\.", "_", rownames(expr_all_modified))

# Convert to matrix for FlowSOM

fsom_mat <- as.matrix(expr_all_modified[,temp1$desc])


#save the RDS data
#run flowsompipeline
flowsom <- FlowSOM(
  input = fsom_mat,
  transform = FALSE,    # Already done via transFlowVS
  scale = FALSE,        # Already scaled via manual cofactors
  colsToUse = 1:ncol(fsom_mat),
  xdim = 10,
  ydim = 10,
  seed = 100
)

# Build MST 
flowsom <- BuildMST(flowsom)

plot_flowsom_elbow <- function(flowsom, k_min = 5, k_max = 15, reps = 50, seed = 100,
                               out_dir = "consensus_cluster") {


# Example usage
result <- plot_flowsom_elbow(flowsom, k_min = 3, k_max = 30, reps = 50, seed = 99)

# Meta clustering
flowsom$metaclustering <- metaClustering_consensus(flowsom$map$codes, k = 15, seed = 100)
# Get metaclusters
clusters_flowsom <- GetMetaclusters(flowsom)

# marker expression values
expr_all_modified$SOM_cluster <- flowsom$map$mapping[,1]

expr_all_modified$MetaCluster <- flowsom$metaclustering[expr_all_modified$SOM_cluster]
expr_all_modified$MetaCluster <- as.factor(expr_all_modified$MetaCluster)


expr_all_modified <- expr_all_modified %>%
  separate_wider_delim(
    SampleID,
    delim = "_",
    names = c("Group", "CellID"),
    too_many = "merge"
  )


expr_all_modified<-as.data.frame(expr_all_modified)

# plot cluster distribution
count_FlowSOM <- table(expr_all_modified[,c("MetaCluster", "Group")])
count_FlowSOM <- as.data.frame(count_FlowSOM)

count_FlowSOM$SampleID <- as.character(count_FlowSOM$Group)


# UMAP

# Downsample
#max_cells_per_cluster <- 10000
 
#downsampled_data <- expr_all_modified %>%
#   group_by(Group, MetaCluster) %>% 
#   slice_sample(n=max_cells_per_cluster) %>%
#   ungroup()

#pca <- prcomp(
#  downsampled_data[, c(1:8)],
#  center = TRUE,
#  scale. = FALSE
#)

# UMAP
#UMAP_results <- umap(
#  pca$x[, 1:5],
#  n_neighbors = 30,
#  min_dist = 0.3,
#  metric = "euclidean",
#  nn_method = "annoy",
#  random_state = 42
#)

# Add back metadata
#downsampled_data$UMAP1 <- UMAP_results$layout[, 1]
#downsampled_data$UMAP2 <- UMAP_results$layout[, 2]

 UMAP_results <- umap(expr_all_modified[,c(1:8)],
                     n_neighbors = 30,
                     min_dist = 0.3,
                     metric = "euclidean",
                     random_state= 42
 )

saveRDS(UMAP_results, "UMAP_result.RDS")

#UMAP_results<-readRDS("UMAP_result.RDS")

umap_df <- data.frame(UMAP1 = UMAP_results$layout[,1],
                      UMAP2 = UMAP_results$layout[,2])

expr_all_modified <- cbind(expr_all_modified, umap_df)

# Assuming downsampled_data has columns: UMAP1, UMAP2, MetaCluster, SampleID

# Static ggplots
p1 <- ggplot(expr_all_modified, aes(x = UMAP1, y = UMAP2, color = MetaCluster)) +
  geom_point(size = 0.5) +
  theme_minimal()

p2 <- ggplot(expr_all_modified, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 0.5) +
  theme_minimal()
library(plotly)
# Convert to interactive Plotly plots
p1_plotly <- ggplotly(p1)
p2_plotly <- ggplotly(p2)

# Display
#p1_plotly
#p2_plotly

library(htmlwidgets)

# Save the interactive plot as an HTML file
saveWidget(as_widget(p1_plotly), "UMAP_clusters.html")
saveWidget(as_widget(p2_plotly), "UMAP_samplewise.html")

# tSNE
 tSNE_results <- Rtsne(expr_all_modified[,c(1:8)], perplexity = 30,
                       theta = 0.5, 
                       dims = 2, 
                       normalize = TRUE 
 )
saveRDS(tSNE_results, "tSNE_result.RDS")

#tSNE_results<-readRDS("tSNE_result.RDS")
tSNE_coords <- as.data.frame(tSNE_results$Y)
colnames(tSNE_coords) <- c("tSNE1", "tSNE2")
expr_all_modified <- cbind(expr_all_modified, tSNE_coords)

p1<-ggplot(expr_all_modified, aes(x=tSNE1, y = tSNE2, color= factor(MetaCluster))) + 
  geom_point(alpha = 0.6, size = 1) + # add transparency for dense plots
  scale_color_viridis_d() + # better color scheme for many clusters
  theme_minimal() 

p2<-ggplot(expr_all_modified, aes(x=tSNE1, y = tSNE2, color= factor(Group))) + 
  geom_point(alpha = 0.6, size = 1) + # add transparency for dense plots
  scale_color_viridis_d() + # better color scheme for many clusters
  theme_minimal() 

library(plotly)
p1_plotly <- ggplotly(p1)
p2_plotly <- ggplotly(p2)
#p1_plotly
#p2_plotly
library(htmlwidgets)

# Save the interactive plot as an HTML file
saveWidget(as_widget(p1_plotly), "TSNE_clusters.html")
saveWidget(as_widget(p2_plotly), "TSNE_samplewise.html")
write.table(expr_all_modified, "expr_all_modified.txt",sep="\t", row.names=T)

# Annotation
marker_cols <- colnames(expr_all_modified)[1:8]

cluster_markers <- expr_all_modified %>%
  group_by(MetaCluster) %>%
  summarize(across(all_of(marker_cols),
                   median,
                   .names = "median_{.col}"))


plot_data<-as.data.frame(cluster_markers)
plot_data$MetaCluster<-paste0("Cluster_",rownames(plot_data))
plot_data<-tibble::column_to_rownames(plot_data, "MetaCluster")

formatted_numbers <- matrix(
  sprintf("%.2f", t(plot_data)),
  nrow = nrow(t(plot_data)),
  ncol = ncol(t(plot_data))
)

# Call pheatmap with the formatted numbers
pheatmap::pheatmap(
  t(plot_data),
  color = colorRampPalette(c("#2166ac", "white", "#d7191c"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  angle_col = 45,
  main = "FlowSOM Metaclusters: Marker Expression",
  fontsize_row = 14,
  fontsize_col = 14,
  border_color = "grey30",
  cellwidth = 40,
  cellheight = 35,
  gaps_col = c(5),
  # Set display_numbers to the formatted matrix
  display_numbers = formatted_numbers, 
  legend = TRUE,
  filename = "FlowSOM_Metacluster_Heatmap.pdf"
)


formatted_numbers_matrix <- matrix(
  sprintf("%.2f", t(plot_data)),
  nrow = nrow(t(plot_data)),
  ncol = ncol(t(plot_data))
)
# Convert the character matrix back to a data frame for easier CSV saving (optional, but robust)
formatted_df <- as.data.frame(formatted_numbers_matrix)
colnames(formatted_df) <- colnames(t(plot_data)) # Restore column names if necessary
rownames(formatted_df) <- rownames(t(plot_data)) # Restore row names if necessary

# Save the formatted data to a file
write.csv(formatted_df, "FlowSOM_formatted_annotations.csv", row.names = TRUE)



# Session Info

sessioninfo::session_info()
