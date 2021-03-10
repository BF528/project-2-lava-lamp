#biologist analysis
#Divya Sundaresan 

library(ggplot2)
library(dplyr)
library(pheatmap)

#-----------------------------------7.3
#Read Data
fpkm_0 <- read.table('/projectnb/bf528/users/lava_lamp/project_2/P0_1_cufflinks/genes.fpkm_tracking', header = TRUE)
fpkm_0 <- fpkm_0  %>%  select(FPKM, tracking_id, gene_short_name)
colnames(fpkm_0) = gsub("FPKM", "P0_1", colnames(fpkm_0))
fpkm_matrix <- read.csv('/project/bf528/project_2/data/fpkm_matrix.csv', sep =  "\t")
list_de_genes <- read.delim("/projectnb/bf528/users/lava_lamp/project_2/cuffdiff_out/gene_exp.diff",
                            header = TRUE, stringsAsFactors = FALSE, quote = "", sep = "\t")
length()
#combine fpkm p0_1 into rest and rename columns
colnames(fpkm_matrix) = gsub("_FPKM", "", colnames(fpkm_matrix))
fpkm_combined <- merge(fpkm_0, fpkm_matrix, by = 'tracking_id')
fpkm_combined <- fpkm_combined %>% relocate(P0_1, .after = Ad_2)
fpkm_combined

#get significant genes, get top 100 deferentially expressed genes 
list_de_genes  <- list_de_genes %>% arrange(q_value)  %>% slice_head(n=1000)
list_de_genes  <- list_de_genes[list_de_genes$significant =='yes',]
top  <- list_de_genes %>% arrange(q_value)  %>% slice_head(n=100)
deg <- top$gene

#subset to top 100 genes 
fpkm_combined <- fpkm_combined[fpkm_combined$gene_short_name %in% deg, ]

 
#reformat create final_fpkm_combined and make genes row names, remove rows with all 0

final_fpkm_combined <- fpkm_combined[-1:-2]
rownames(final_fpkm_combined) <- make.names(fpkm_combined[,2], unique = TRUE)
row.names(final_fpkm_combined) <- fpkm_combined$gene_short_name
final_fpkm_combined <- final_fpkm_combined[apply(final_fpkm_combined[,-1], 1, function(x) !all(x==0)),]

#convert to matrix for heatmap
a <- as.matrix(final_fpkm_combined[,])

# colors for heatmap 
my_colors = brewer.pal(n = 11, name = "RdBu")
my_colors = colorRampPalette(my_colors)(50)
my_colors = rev(my_colors)

#saveheatmap
my_heatmap <- pheatmap(a, scale = "row", color = my_colors,fontsize_row = 4,border_color = NA, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean")

#function for save file
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#save heatmap to current directory
save_pheatmap_png(my_heatmap, "my_heatmap.png")
