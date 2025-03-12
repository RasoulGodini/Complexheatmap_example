library(ComplexHeatmap)
library(tidyverse)
library(truncnorm)
library(circlize)
library(factoextra)

#-------- Generate random data ---------
set.seed(123)  # For reproducibility

n_participants <- 200
n_genes <- 100

sample_id <- sample(paste0("Sample_", sprintf("%03d", 1:n_participants)), n_participants, replace = FALSE)
sex <- sample(c("male", "female"), n_participants, replace = TRUE) # 0=male, 1=female
Age <- sample(c(30:70), n_participants, replace = TRUE)
Weight <- sample(c(40:90), n_participants, replace = TRUE)
Ethnicity <- sample(c("A", "B", "C", "D", "E"), n_participants, replace = TRUE)
blood_sugar <- sample(c(5:12), n_participants, replace = TRUE)


# Initialize matrix for gene expression
gene_data <- matrix(nrow = n_participants, ncol = n_genes)

# Generate gene expression values based on Sex
for (i in 1:n_participants) {
  if (sex[i] == "male") {
    gene_data[i, ] <- rtruncnorm(n_genes, a = 0, b = 300, mean = 70, sd = 60)
  } else {
    gene_data[i, ] <- rtruncnorm(n_genes, a = 0, b = 200, mean = 50, sd = 30)
  }
}

# Convert to dataframe and add Sex column
colnames(gene_data) <- paste0("Gene_", 1:n_genes)  # Name genes
df <- data.frame(Sampe_ID = sample_id, Sex = sex, Age = Age, 
                 Weight = Weight, Ethnicity = Ethnicity, 
                 Blood_sugar = blood_sugar, gene_data)  # Make the datframe

genes_log2 <- log2(df[,7:ncol(df)]+1)
df_log2 <- cbind(df[1:6], genes_log2)
rownames(df_log2) <- df_log2$Sampe_ID
#df_log2 <- df_log2[,-1]


# Prepare gene annotation data
gene_names <- colnames(gene_data)
gene_annotation <- data.frame(Gene_ID = gene_names, 
                              Gene_family = sample(c("Family_A", "Family_B", 
                                                     "Family_C", "Family_D"), 
                                                   n_genes, replace = TRUE))


#Prepare files for the heatmap
# This part is required for row annotation later
df_t <- t(df_log2[,7:ncol(df_log2)])

df_log_scale <- scale(t(df_t))
df_log_scale_clinical <- merge(df_log2[,c(1:5)], df_log_scale, 
                               by.x = "row.names", 
                               by.y = "row.names")


#Arrange participants according to sex (if columns are not clustured)
df_log_scale_clinical_ordered <- df_log_scale_clinical %>% arrange(desc(Sex))
rownames(df_log_scale_clinical_ordered) <- df_log_scale_clinical_ordered$Sampe_ID
df_log_scale_clinical_ordered <- df_log_scale_clinical_ordered[,-1]
df_log_scale_clinical_ordered$Sex <- factor(df_log_scale_clinical_ordered$Sex, 
                                                   levels = c("male", "female"))

# To prepare gene annotation (make sure the order of genes in the heatmap data and annotation data are the same)
df_heatmap <- as.matrix(t(df_log_scale_clinical_ordered[,6:ncol(df_log_scale_clinical_ordered)]))
df_heatmap_gene_family <- merge(gene_annotation, df_heatmap, 
                                by.x = "Gene_ID", by.y = "row.names", sort = FALSE)

max(df_heatmap)
min(df_heatmap)


# -------------- Prepare the annotations -------------- 
Sex_stat <- data.frame(Sex_stat = as.factor(df_log_scale_clinical_ordered$Sex))
col_sex <- structure(c("royalblue4", "darkorange2"), names = c("male", "female"))
.Sex_stat <- anno_simple(Sex_stat, simple_anno_size = unit(1.2, "mm"),
                         which = "col", col = col_sex)

Ethnicity <- data.frame(Ethnicity = as.factor(df_log_scale_clinical_ordered$Ethnicity))
col_Ethnicity <- structure(c("red", "green", "yellow", "purple", "blue"),
                           names = c("A", "B", "C", "D", "E"))

.Ethnicity <- anno_simple(Ethnicity, simple_anno_size = unit(1.2, "mm"),
                          which = "col", col = col_Ethnicity)

.age = anno_lines(df_log_scale_clinical_ordered$Age, height = unit(1, "cm"), 
                  pt_gp = gpar(col = "skyblue2"), pch = 1,
                  axis_param = list(
                    side = "left",
                    gp = gpar(cex = 0.5)
                  ),
                  size = unit(0.6, "mm"), axis = TRUE, smooth = TRUE)


.weight = anno_lines(df_log_scale_clinical_ordered$Weight, height = unit(1, "cm"),
                  pt_gp = gpar(col = "hotpink"), pch = 1,
                  axis_param = list(
                    side = "left",
                    gp = gpar(cex = 0.5)
                  ),
                  size = unit(0.6, "mm"), axis = TRUE, smooth = TRUE)


bottom_ann = HeatmapAnnotation("Sex status" = .Sex_stat,
                               "Ethnicity " = .Ethnicity,
                               "Age" = .age,
                               "Weight" = .weight,
                               annotation_name_gp = gpar(fontsize = 5),
                               show_annotation_name = TRUE ,show_legend = TRUE,
                               annotation_name_side = "left")


#Left annotation of the genes, multiple can be added like above
#Note: you need to define the columns above
col1 <- structure(c("darkred", "blue", 
                    "green", "orange1"),
                  names = c("Family_A", "Family_B", 
                            "Family_C", "Family_D"))

anno_df <- df_heatmap_gene_family$Gene_family #Get the annotation for each gene
.Gene_family <- anno_simple(x = anno_df, col = col1, which = "row")


left_ann = HeatmapAnnotation("Gene family" = .Gene_family, which = "row",
                             annotation_name_gp = gpar(fontsize = 10),
                             annotation_width = unit(2, "mm"),
                             show_annotation_name = FALSE,
                             show_legend = TRUE)

#Color of the heatmap
colors_heatmap <- hcl.colors(11, palette = "Purple-Green")
col_fun_heatmap = colorRamp2(c(-6, -3, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3), colors_heatmap)


# To add coordinates to splite columns
col_split <- factor(c("male" = rep("male", 106), "female" = rep("female", 94)), 
                    levels = c("male", "female"))


# -------------- Cluster and visualize the heatmap -------------- 
set.seed(123)#To set seeds for km repeats
p_heaatmap <- Heatmap(df_heatmap, col = col_fun_heatmap,
                      border = TRUE, border_gp = gpar(col = "black", lwd = 0.5),
                      cluster_columns = TRUE, show_column_dend = FALSE,
                      clustering_distance_columns = "euclidean",
                      clustering_method_columns = "complete",
                      cluster_rows = TRUE, show_row_dend = TRUE, 
                      clustering_distance_rows = "euclidean", 
                      clustering_method_rows = "average",
                      row_title = NULL, column_title = NULL, 
                      show_row_names = FALSE, row_names_gp = gpar(fontsize = 5), row_names_side = "right",
                      show_column_names = FALSE, column_names_gp = gpar(fontsize = 2),
                      row_km = 4,# Number of cluster obtained from fviz_nbclust() bellow
                      row_km_repeats = 1000, #Row clustering repeats
                      column_split = col_split, #To splite columns based on values defined above
                      row_dend_width = unit(0.7, "cm"), gap = unit(0.1, "cm"), 
                      show_parent_dend_line = FALSE,
                      bottom_annotation = bottom_ann,
                      left_annotation = left_ann,
                      width = unit(6, "cm"), 
                      height = unit(4, "cm"),
                      show_heatmap_legend = FALSE)

#Define the cosyume legends for the heatmap
legend_sex <- Legend(title = "Sex status", labels = c("Male", "Female"),
                     legend_gp = gpar(fill = c("royalblue4", "darkorange2")), direction = "horizontal", 
                     labels_gp = gpar(fontsize = 5), grid_height = unit(1, "mm"), title_position = "topcenter",
                     grid_width = unit(2, "mm"), , legend_width = unit(2, "mm"), nrow = 1,
                     title_gp = gpar(fontsize = 5, fontface ="bold"))


legend_log2_expression <- Legend(title = "Z-score", col_fun = col_fun_heatmap, at = c(-6, -3, -1.5, 0, 1, 2, 3), 
                             labels = c("-6", "-3", "-1.5", "0", "1", "2", "3"), direction = "horizontal", 
                             labels_gp = gpar(fontsize = 5), grid_height = unit(2, "mm"), 
                             title_position = "topcenter", 
                             legend_height = unit(10, "mm"), legend_width = unit(20, "mm"),
                             grid_width = unit(1, "mm"), title_gp = gpar(fontsize = 5, fontface ="bold"))


legend_gene_family <- Legend(title = "Gene family", labels = c("Family_A", "Family_B", 
                                                                "Family_C", "Family_D"),
                              legend_gp = gpar(fill = c("darkred", "blue", 
                                                        "green", "orange1")), direction = "horizontal", title_position = "topcenter",
                              labels_gp = gpar(fontsize = 5), grid_height = unit(1, "mm"), nrow = 4,
                              grid_width = unit(2, "mm"), legend_width = unit(2, "mm"),
                              title_gp = gpar(fontsize = 5, fontface ="bold"))


legend_Ethni <- Legend(title = "Ethnicity", labels = c("A", "B", "C", "D", "E"),
                       legend_gp = gpar(fill = c("red", "green", "yellow", "purple", "blue")), 
                       direction = "horizontal", 
                       labels_gp = gpar(fontsize = 5), grid_height = unit(1, "mm"), nrow = 3,
                       grid_width = unit(2, "mm"), title_gp = gpar(fontsize = 5, fontface ="bold"))

# -------------- Visualize the heatmap and legends -------------- 
pdf("Complexheatmap_example.pdf", width = 5, height = 5, bg = "transparent")

draw(p_heaatmap, 
     annotation_legend_list = list(legend_log2_expression, 
                                   legend_gene_family, 
                                   legend_sex, 
                                   legend_Ethni), 
     heatmap_legend_side = "right", annotation_legend_side = "top", #Heatmap and annotation legens coordinates 
     merge_legend = TRUE)


dev.off()

# To get the optimum number of clusters based on the appearance of the plot 
p_elbow <- fviz_nbclust(df_heatmap, kmeans, method = "wss")
plot(p_elbow)
