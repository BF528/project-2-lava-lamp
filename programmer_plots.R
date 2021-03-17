#!/usr/bin/env/Rscript
#libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library("gplots")
library(tidyverse)


#paths
setwd('/projectnb/bf528/users/lava_lamp/project_2/esaake_pr2')
fpkmpath="/projectnb/bf528/users/lava_lamp/project_2/esaake_pr2/P0_1_cufflinks/genes.fpkm_tracking"

#data load
genedata=read.table(fpkmpath, header=TRUE)

#view data briefly
head(genedata)
nrow(genedata)

#viewing columnames
names(genedata)

#removing genes with zero fpkm values
genefpkm_trimmed <- genedata[genedata$FPKM!=0, ]  # Remove zero-rows

#number of genes after deletion of genes with FPKM=0
nrow(genefpkm_trimmed)

set.seed(25)

#checking FPKM field data
fpkmplotdata<-as.data.frame(genefpkm_trimmed)
fpkmplotdata<-fpkmplotdata %>% select(FPKM,gene_short_name,gene_id)
fpkmplotdata2<-fpkmplotdata

#removing duplicates
fpkmplotdata2<-fpkmplotdata2 %>% distinct(fpkmplotdata2$gene_short_name, .keep_all = TRUE)

#count of genes after removal of duplicates
nrow(fpkmplotdata2)

#viewing rownames and select out only two columns FPKM and geneshortname
#NB: Didn't use this processed data as intended...it was to be for hist gene vs fpkm
row.names(fpkmplotdata2)<-(fpkmplotdata2$gene_short_name)
fpkmplotdata2<-fpkmplotdata2 %>% select(FPKM,gene_short_name)


#Histogram with overlay of density lines
  #if you want to save plots uncomment png and dev() ---currently commented to avoid duplication

#png("./Qseqcplots/genefpkm_density.png")
fpkmplotdata<-as.numeric(genefpkm_trimmed$FPKM, center=TRUE)
max(fpkmplotdata)
fpkm_sd = sd(fpkmplotdata)
fpkm_mean = mean(fpkmplotdata)
fpkm_median = median(fpkmplotdata)
hist(fpkmplotdata,probability=T,xlab="FPKM values",main=paste(c("Mean=",fpkm_mean,";","SD=",fpkm_sd),collapse=""),border="blue")
lines(density(fpkmplotdata,bw=10),col='red')
#dev.off()

#png("./Qseqcplots/genefpkmhist_rnormed.png")
fpkmplotdata_norm<-rnorm(as.numeric(genefpkm_trimmed$FPKM, center=TRUE))
#computing sd, mean, median for main() label
fpkmn_sd = sd(fpkmplotdata_norm)
fpkmn_mean = mean(fpkmplotdata_norm)
fpkmn_median = median(fpkmplotdata_norm)
write(x=c("Name","Mean","Median","sd"), sep="	", file=stdout(),ncolumns=4)
write(c(out_file,fpkmn_mean,fpkmn_median,fpkmn_sd),sep="	", file=stdout(),ncolumns=4)
hist(fpkmplotdata_norm,probability=T,breaks=100,xlab="FPKM values",main=paste(c("Mean=",fpkmn_mean,";","SD=",fpkmn_sd),collapse=""),border="blue")
lines(density(fpkmplotdata_norm,),col='red')
#dev.off()


#Plots from QSEQC RUN
----------------------
  
#1. GeneBodyCoverage
  accepted_hits <- c(0.0,0.09573223137350666,0.162482517593406,0.20489697669329596,0.24827554492780282,0.27730426840290884,0.30360319932004937,0.3330505491299232,0.35501812461744847,0.379637793078077,0.4155762296355731,0.4318178214304969,0.4538868820900871,0.4644904968650596,0.4933503426710263,0.5184473085815228,0.5305343511450382,0.538797464139263,0.5623792881449462,0.5837601524814711,0.6041808719479127,0.620793518903199,0.6420118483938386,0.6529262615082599,0.6717224253687559,0.6867660164215694,0.6966838134321968,0.7229637158795752,0.7268479022697794,0.7270207442034524,0.7173368388003184,0.7303982975862386,0.7384830186764431,0.7564728511307668,0.7727041358379028,0.7807674498996249,0.7861255498434908,0.8169524637111225,0.8117870270235985,0.8185381695246371,0.8056891953177279,0.8131245698772981,0.8244013167701075,0.8341058363488172,0.8506114481616912,0.8521979468281127,0.8553971083068784,0.8806042807714142,0.8968569725070326,0.9079402632905933,0.9040291199015594,0.9033456806959346,0.9051795494692643,0.9207495948521647,0.9170897858345728,0.9171492498025796,0.910760441079929,0.9129669507194347,0.9087148805804952,0.9253735130043734,0.9280430487414253,0.9337262183770618,0.9477303792691165,0.9465482355851413,0.9696115337898051,0.984951651829746,0.9893115499640045,0.988453683118893,0.9737280260816892,0.9607157241759879,0.9635921945217035,0.9532407069710799,0.9674486310601712,0.9652762140956561,0.9555431518123032,0.9641352987628323,0.9722533196751206,0.9753811243922782,0.9893123428169113,1.0,0.9940980029620985,0.9816073982690435,0.9725062397523762,0.9702370947332367,0.9647156670905787,0.9637047796344631,0.9416404760923135,0.9372076354906332,0.9331902498120939,0.9091580853553726,0.8921061978897427,0.8939765378967832,0.8801610759965368,0.863324844521545,0.8242744603050264,0.7802949095671975,0.738035849637032,0.693598822772004,0.608829368540231,0.43582569287415523)
#png("./Qseqcplots/rseqcgB.geneBodyCoverage.curves.png")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,accepted_hits,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
#dev.off()

#2. InnerDistance
out_file = 'rseqcindist'
#png('./Qseqcplots/rseqcindist.inner_distance_plot.png')
fragsize=rep(c(-247.5,-242.5,-237.5,-232.5,-227.5,-222.5,-217.5,-212.5,-207.5,-202.5,-197.5,-192.5,-187.5,-182.5,-177.5,-172.5,-167.5,-162.5,-157.5,-152.5,-147.5,-142.5,-137.5,-132.5,-127.5,-122.5,-117.5,-112.5,-107.5,-102.5,-97.5,-92.5,-87.5,-82.5,-77.5,-72.5,-67.5,-62.5,-57.5,-52.5,-47.5,-42.5,-37.5,-32.5,-27.5,-22.5,-17.5,-12.5,-7.5,-2.5,2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5,82.5,87.5,92.5,97.5,102.5,107.5,112.5,117.5,122.5,127.5,132.5,137.5,142.5,147.5,152.5,157.5,162.5,167.5,172.5,177.5,182.5,187.5,192.5,197.5,202.5,207.5,212.5,217.5,222.5,227.5,232.5,237.5,242.5,247.5),times=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,49,632,549,590,634,693,781,1071,1614,2463,3791,5819,8920,13199,18950,24236,30986,36202,40520,44134,48253,48530,49772,48764,48112,45552,43596,40775,37738,34111,30639,28150,24899,22332,19970,17170,15149,13423,11940,10270,9032,7871,6936,5958,5358,4584,3968,3520,2952,2734,2353,2090,1933,1698,1578,1453,1346,1191,1083))
frag_sd = sd(fragsize)
frag_mean = mean(fragsize)
frag_median = median(fragsize)
write(x=c("Name","Mean","Median","sd"), sep="	", file=stdout(),ncolumns=4)
write(c(out_file,frag_mean,frag_median,frag_sd),sep="	", file=stdout(),ncolumns=4)
hist(fragsize,probability=T,breaks=100,xlab="mRNA insert size (bp)",main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse=""),border="blue")
lines(density(fragsize,bw=10),col='red')
#dev.off()