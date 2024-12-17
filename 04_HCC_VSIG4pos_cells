hcc <- subset(hcc, subset = leiden %in% c("17", "21", "23"), invert = TRUE)
phago3 <- read.gmt("GOBP_PHAGOCYTOSIS_RECOGNITION.v2024.1.Hs.gmt")
phago3<-phago3$gene
phagocytosis<-as.data.frame(phago3)
phagocytosis<-as.list(phagocytosis)
hcc <- AddModuleScore(hcc,features = phagocytosis,name = 'Phagocytosis')
ggviolin(hcc@meta.data, x = "VSIG4_status", y = "Phagocytosis1",
         color = "VSIG4_status",add = 'mean_sd',fill = 'VSIG4_status',
         add.params = list(color = "black"))+
  scale_color_manual(values = anno_col) + 
  scale_fill_manual(values = anno_col) +
  theme(axis.text.x.bottom = element_text(angle = 45,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')+
  stat_pvalue_manual(stat.test, label = "p")


#基因交集
gene_vsig_low <- c("ALAS1", "ANPEP", "PGD", "SLC43A3", "LPL", "RASGRP3", 
                 "TXNRD1", "DUSP4", "TRAF1", "EEPD1", "ANKRD37", 
                 "SEMA4D", "CHST11", "PTGER2")
gene_vsig_hi <- c("VSIG4","NDUFAF4", "SERPING1", "ITGA9", "IFI16", "MAF", "CP")





