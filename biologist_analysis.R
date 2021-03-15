#biologist analysis
#Divya Sundaresan 

library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(tibble)
require(gridExtra)
library(cowplot)
require(grid)
library(patchwork)



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

#-----------------------------------7.1
sarc_list <- c('Pdlim5', 'Pygm', 'Myoz2', 'Des', 'Csrp3', 'Tcap', 'Cryab')
mito_list <- c("Mpc1","Prdx3","Acat1","Echs1","Slc25a11","Phyh")
cell_list <- c("Cdc7","E2f8","Cdk7","Cdc26","Cdc6","Cdc27",
            "E2f1","Cdc45","Rad51","Aurkb","Cdc23")

#list of lists for for loop
lists <- list(sarc_list,mito_list,cell_list)

#function to reformat datframe for graphing
reform_df <- function(list, samp_set, fpkm_m) {
  f <- fpkm_combined[fpkm_combined$gene_short_name %in% list, ]
  f <- f %>% relocate(Ad_1, .after = P7_2)
  if (samp_set == 1) {
    sample <- f %>% relocate(Ad_2, .after = Ad_1)  %>% select(gene_short_name, P0_1, P4_1, P7_1, Ad_1)
  }
  else if (samp_set == 2) {
    sample <- f %>% relocate(Ad_2, .after = Ad_1)  %>% select(gene_short_name, P0_2, P4_2, P7_2, Ad_2)
  }
  rownames(sample) <- NULL
  sample <- as.data.frame(column_to_rownames(sample, var = "gene_short_name"))
  sample <- as.data.frame(t(sample))
  sample$Samples <- row.names(sample)
  return (sample)
}

#plotting function
plot_lines <- function (df) {
  line_plot <- df %>% 
    mutate(Samples = factor(Samples, levels = unique(Samples))) %>%
    gather(Genes, FPKM, -Samples) %>% 
    ggplot(aes(Samples, FPKM)) + 
    geom_line(aes(color = Genes, group = Genes)) +
    geom_point() +
    scale_y_log10() + 
    scale_color_brewer(palette = "Paired")  +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))
  return (line_plot)
}

#for loop for plotting all graphs and saving into png
w <- 3
titles <- c('Sarcomere', 'Mitochodria', 'Cell Cycle')
titles[w]
for (i in lists) {
  file_name <- paste("line_plot", w, ".png", sep="")
  p1 <- plot_lines(reform_df(i, 1, fpkm_combined))
  p2 <- plot_lines(reform_df(i, 2, fpkm_combined))
  png(file_name)
  l <- grid.arrange(p1, p2, nrow = 1)
  w <- w +1
}

# Additional plot code used for formatting using patchwork for layout however 
#has issues running in for loop therefore kept the raw graphs above in loop 
# p1 + p2 + plot_spacer() + plot_spacer() + plot_annotation(title = titles[w],tag_levels = 'i', tag_suffix = ')')
                    

dev.off()

#-----------------------------------7.2
#up regulated ------
#read in reference paper analysis
paper_up_reg <- read.csv("upreg_paperreference.csv", stringsAsFactors=FALSE, header = FALSE)  
colnames(paper_up_reg) <- paper_up_reg[2,] 
paper_up_reg <- paper_up_reg %>%
  filter(Category %in% c("GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"))

#get go terms 
terms <- paper_up_reg$Term

#read in our analysis  
my_up_reg <- read.csv("Up_Reg_Clusters.txt",header=F,stringsAsFactors=F, sep = '\t')
colnames(my_up_reg) <- my_up_reg[2,]
my_up_reg <- my_up_reg %>%
  filter(Category %in% c("GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"))

my_up_reg$Paper_ref_overlap <- my_up_reg$Term %in% terms
write.csv(my_up_reg,"Up_regulated_extended.csv")


#down regulated ---------
#read in reference paper analysis
paper_down_reg <- read.csv("upreg_paperreference.csv", stringsAsFactors=FALSE, header = FALSE)  
colnames(paper_down_reg) <- paper_down_reg[2,] 
paper_down_reg <- paper_down_reg %>%
  filter(Category %in% c("GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"))


#get go terms 
terms <- paper_down_reg$Term

#read our data
my_down_reg <- read.csv("Down_Reg_clusters.txt",header=F,stringsAsFactors=F, sep = '\t')
colnames(my_down_reg) <- my_down_reg[2,]
my_down_reg <- my_down_reg %>%
  filter(Category %in% c("GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"))

my_down_reg$Paper_ref_overlap <- my_down_reg$Term %in% terms
write.csv(my_down_reg,"Down_regulated_extended.csv")
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
