setwd("C:\\Users\\jasonxu\\Desktop\\Loquat_weight")
library(BioCircos)
library(tidyverse)
library(htmltools)
library(shiny)
library(readxl)
source("calcDensityGFF.r")
mycolors <- rev(c("#9370DB","#483D8B","#B0C4DE",
                  "#4169E1","#4682B4","#5F9EA0","#00CED1","#008080","#3CB371","#BDB76B","#FFD700",
                  "#FFA500","#800000",'#3A97CA','#5CB8EA','#F1E648',"#F14A48"))

eqtl <- read.table("05_pQTL/F1_pQTL_results.txt",header=T)
eqtl <- eqtl %>%
  mutate(CHR=str_extract(SNP,"\\w+:")) %>%
  mutate(CHR=gsub(":","",CHR)) %>%
  mutate(BP=str_extract(SNP,":(\\w+)")) %>%
  mutate(BP=gsub(":","",BP)) %>%
  mutate(BP=as.numeric(BP)) %>%
  dplyr::select(CHR,BP,p.value) %>%
  mutate(value=-log10(p.value)) %>%
  mutate(colors=as.character(factor(CHR,labels=mycolors))) %>%
  arrange(CHR)
mycolors2 <- eqtl$colors
names(mycolors2) <- eqtl$CHR
Fst <- read.table("03_vcf//Fst//large_small.windowed.weir.fst",head=T)
Fst <- Fst %>%
    dplyr::select(CHROM,BIN_START,BIN_END,MEAN_FST)
for(genome in c('SS')){
  karyotype <- read.table(paste0("genome/",genome,".karyotype"))
  myGenome <- lapply(1:nrow(karyotype), function(x){
    karyotype$V3[x]
  })
  tracks <- BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.99, maxRadius = 1.03,
                                              borderColors = "white", 
                                              chrPad=0,
                                              borderSize = 0, fillColors = "white") 
  names(myGenome) <- karyotype$V1
  geneDensity <- calcDensityGFF(paste0("genome/",genome,".gff3.gz"),karyotype_info=karyotype,window=300000)
  tracks <- tracks+BioCircosHeatmapTrack("heatmap1", 
                                         geneDensity$Chr,
                                         geneDensity$Start - 1, 
                                         geneDensity$End,
                                         geneDensity$Count,
                                         color = c('#F0E68C','#3CB371'),
                                         minRadius = 0.69, maxRadius = 0.78)
  tracks <- tracks+BioCircosSNPTrack("snp_track1",
                    chromosomes = eqtl$CHR,
                    positions = eqtl$BP,
                    values = -eqtl$value,
                    maxRadius = 0.79,
                    minRadius = 1,
                    colors = mycolors2,
                    size=1
  )
  homologs <- read_excel("05_pQTL/snp-gene-gene.xlsx")
  homologs <- homologs[,-1]

  chrs <- unique(homologs$gene1Chromosomes)
  for(chr in chrs){
    inx <- which(chrs==chr)
    tmp <- subset(homologs,gene1Chromosomes==chr)
    links <- BioCircosLinkTrack("myLinkTrack", 
                                gene1Chromosomes = tmp$gene1Chromosomes, 
                                gene1Starts = tmp$gene1Starts, 
                                gene1Ends = tmp$gene1Ends,
                                axisPadding = 6, 
                                color = mycolors[inx], 
                                gene2Chromosomes = tmp$gene2Chromosomes,
                                gene2Starts = tmp$gene2Starts, 
                                gene2Ends = tmp$gene2Ends,
                                labels = "",
                                labelSize = "0.5em",
                                maxRadius = 0.28,
                                width = "0.1em")

    tracks <- tracks + links
  }
  pic <- BioCircos(genome = myGenome, tracks,
                   genomeFillColor = mycolors,
                   genomeTicksDisplay = T,genomeLabelOrientation = -90,
                   chrPad = 0.02,genomeBorderSize = 3,genomeBorderColor = "white",
                   genomeLabelTextSize = "8pt",genomeLabelDy = 20,
                   genomeTicksTextSize = "0.3em", 
                   genomeTicksScale = 5e+6)
}
pic
