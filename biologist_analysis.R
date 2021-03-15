#biologist analysis
#Divya Sundaresan 

library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(tibble)

#Read Data
fpkm_0 <- read.table('/projectnb/bf528/users/lava_lamp/project_2/P0_1_cufflinks/genes.fpkm_tracking', header = TRUE)
fpkm_0 <- fpkm_0  %>%  select(FPKM, tracking_id, gene_short_name)
colnames(fpkm_0) = gsub("FPKM", "P0_1", colnames(fpkm_0))
fpkm_matrix <- read.csv('/project/bf528/project_2/data/fpkm_matrix.csv', sep =  "\t")
list_de_genes <- read.delim("/projectnb/bf528/users/lava_lamp/project_2/cuffdiff_out/gene_exp.diff",
                            header = TRUE, stringsAsFactors = FALSE, quote = "", sep = "\t")

#combine fpkm p0_1 into rest and rename columns
colnames(fpkm_matrix) = gsub("_FPKM", "", colnames(fpkm_matrix))
fpkm_combined <- merge(fpkm_0, fpkm_matrix, by = 'tracking_id')
fpkm_combined <- fpkm_combined %>% relocate(P0_1, .after = Ad_2)
fpkm_combined
#-----------------------------------7.1
sarc_list <- c('Pdlim5', 'Pygm', 'Myoz2', 'Des', 'Csrp3', 'Tcap', 'Cryab')
mito_list <- c("Mpc1","Prdx3","Acat1","Echs1","Slc25a11","Phyh")
cell_list <- c("Cdc7","E2f8","Cdk7","Cdc26","Cdc6","Cdc27",
            "E2f1","Cdc45","Rad51","Aurkb","Cdc23")

f <- fpkm_combined[fpkm_combined$gene_short_name %in% gene_list, ][-1]
f <- f %>% relocate(Ad_1, .after = P7_2)
f <- f %>% relocate(Ad_2, .after = Ad_1)  %>% select(gene_short_name, P0_1, P4_1, P7_1, Ad_1) #not sure to average or not
e <- f[-1]
row.names(e) <- f$gene_short_name
e <- as.data.frame(t(e))
e$Samples <- row.names(e)
row.names(e) <- NULL
line_plot <- e %>% 
  mutate(Samples = factor(Samples, levels = unique(Samples))) %>%
  gather(Genes, FPKM, -Samples) %>% 
  ggplot(aes(Samples, FPKM)) + 
  geom_line(aes(color = Genes, group = Genes)) +
  geom_point() +
  scale_y_log10() + 
  ggtitle("Sarcomere") +
  scale_color_brewer(palette = "Dark2")  +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))
png("sarc_line.png")
line_plot
dev.off()

#-----------------------------------7.2
upreg_genes <- read.csv("UpReg_genes.txt",header=F,stringsAsFactors=F)

names(top_up_clus) <- top_up_clus[2,]
top_up_clus <- top_up_clus %>%
  filter(Category %in% c("GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"))


top_up_clus <- read.csv("/projectnb/bf528/users/saxophone/project2/upreg_paper_reference.csv",header=F,stringsAsFactors=F)
colnames(top_up_clus) <- top_up_clus[2,]
#-----------------------------------7.3

#average duplicates
fpkm_combined <- fpkm_combined[-1]
fpkm_combined_1 <- aggregate(.~gene_short_name, data=fpkm_combined, mean)

#get significant genes, get top 100 deferentially expressed genes 
#list_de_genes  <- list_de_genes %>% arrange(q_value)  %>% slice_head(n=1000)
sig_de_genes  <- list_de_genes[list_de_genes$significant =='yes',]
top  <- sig_de_genes %>% arrange(q_value)  %>% slice_head(n=120)
deg <- top$gene
length(unique(deg))

#subset to top 100 genes 
fpkm_combined_sub <- fpkm_combined_1[fpkm_combined_1$gene_short_name %in% deg, ]

 
#reformat create final_fpkm_combined and make genes row names, remove rows with all 0

rownames( fpkm_combined_sub ) <- NULL
fpkm_combined_sub <- data.frame(column_to_rownames(fpkm_combined_sub, var = "gene_short_name"))

#convert to matrix for heatmap
a <- as.matrix(fpkm_combined_sub[,])

# colors for heatmap 
my_colors = brewer.pal(n = 11, name = "RdBu")
my_colors = colorRampPalette(my_colors)(50)
my_colors = rev(my_colors)

#saveheatmap
my_heatmap <- pheatmap(a, scale = "row", color = my_colors,fontsize_row = 4,border_color = NA, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean")

#function for save file
save_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#save heatmap to current directory
save_png(my_heatmap, "my_heatmap.png")
