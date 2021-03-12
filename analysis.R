install.packages('ggplot2')
install.packages('ggpubr')
library('ggplot2')
library('ggpubr')

setwd("/projectnb2/bf528/users/lava_lamp/project_2")

# Loading in expression file and selecting 10 most significant genes
gene_diff <- read.delim("cuffdiff_out/gene_exp.diff")
gene_diff <- gene_diff[order(gene_diff$q_value,decreasing = FALSE),]
top_ten_genes <- gene_diff[1:10,c(3,8,9,10,12,13)]
write.table(top_ten_genes,'top_ten_diffexp_genes.txt',sep = '\t')

#histogram of log2 fold change
all_genes_hist <- ggplot(gene_diff, aes(x = log2.fold_change.)) + geom_histogram(binwidth = 0.3, fill = 'dodgerblue3', color = 'black') + 
                    theme_light() + xlab('Log2 Fold Change') + ylab('Count') + scale_y_log10() +
                    ggtitle('All Genes') + theme(plot.title = element_text(size=11))


# new df of only significant genes and histogram
sig_genes <- subset(gene_diff,significant == 'yes')
sig_genes_hist <- ggplot(sig_genes, aes(x = log2.fold_change.)) + geom_histogram(binwidth = 0.3, fill = 'darkorange1', color = 'black') + 
                    theme_light() + xlab('Log2 Fold Change') + ylab('Count') + scale_y_log10() +
                    ggtitle('Significant Genes')+ theme(plot.title = element_text(size=11))

#Combined histograms
pdf('DiffExp_Hist.pdf',height = 6,width = 5)
combined_hists <- ggarrange(all_genes_hist,sig_genes_hist, ncol = 1)
annotate_figure(combined_hists, top=text_grob('Differential Expression of Genes Between \nP0 and Adult Mouse Cardiac Myocytes', size = 14))
dev.off()

#subsetting into upregulated and downregulated
up_reg <- sig_genes[sig_genes$log2.fold_change. > 0,] #1091 up regulated genes
down_reg <- sig_genes[sig_genes$log2.fold_change. < 0,] #1032 down regulated genes

#writing gene names to files
write(up_reg$gene, file = 'UpReg_genes.txt')
write(down_reg$gene, file = 'DownReg_genes.txt')
