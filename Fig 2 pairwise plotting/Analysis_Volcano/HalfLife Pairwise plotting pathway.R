######### Half-life pair-wise comparison 

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggrepel)

df_byPro = read.csv('../output_7/Peptide_fit_statistics.csv')

pathway = read.csv('../../CrossAnalysis_CC/Peroxisome/PeroxisomeDB_Mus_musculus_geneNames.csv')
#pathway = pathway[pathway$Peroxisome == 'True', 'Entry']
pathway = pathway$X

POI = subset(df_byPro, Genes %in% pathway)
#write.csv(POI,'Peptide_fit_statistics_Peroxisome.csv')

## Volcano Plots
### Hypoxia vs Normoxia
FC_Cutoff = 2
Pval_Cutoff = 0.05

#Sig = subset(df_byPro, abs(Oxy8vs21_Log2FC)>log2(FC_Cutoff) & Oxy8vs21_FDR<Pval_Cutoff)

pdf('Oxy8vs21_Volcano.pdf', width = 4, height = 4)
ggplot(data = df_byPro, aes(x=Oxy8vs21_Log2FC, y =-log10(Oxy8vs21_FDR))) +
  geom_point(size = 1, color = "grey", alpha = 0.6) +
  theme_classic() +
  xlab("Log2 Fold-Change") + ylab("-log(FDR)") +
  theme(plot.title = element_text(colour = "black", size = 16, hjust = 0.5),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
        axis.ticks = element_line(color = "black", size = 0.5, linetype = "solid"),
        legend.title = element_blank()) +
  geom_hline(yintercept = -log10(Pval_Cutoff), linetype="dashed", color = "grey", size = 0.2) +
  geom_vline(xintercept = c(-log2(FC_Cutoff),0, log2(FC_Cutoff)), linetype="dashed", color = "grey", size = 0.2) +
  geom_point(POI, mapping = aes(x=Oxy8vs21_Log2FC, y =-log10(Oxy8vs21_FDR)), 
             size = 1, color = "#21618C", alpha = 0.8) +
  geom_label_repel(POI, mapping = aes(x=Oxy8vs21_Log2FC, y =-log10(Oxy8vs21_FDR)), 
                   label = POI$Genes, color = "#21618C",max.overlaps = 40) +
  xlim(-5,5) +
  ylim(0,20)
dev.off()

### Hyperoxia vs Normoxia
FC_Cutoff =2
Pval_Cutoff = 0.05

#Sig = subset(df_byPro, abs(Oxy60vs21_Log2FC)>log2(FC_Cutoff) & Oxy60vs21_FDR<Pval_Cutoff)

pdf('Oxy60vs21_Volcano_Sarcomere.pdf', width = 4, height = 4)
ggplot(data = df_byPro, aes(x=Oxy60vs21_Log2FC, y =-log10(Oxy60vs21_FDR))) +
  geom_point(size = 1, color = "grey") +
  theme_classic() +
  xlab("Log2 Fold-Change") + ylab("-log(pval)") +
  theme(plot.title = element_text(colour = "black", size = 16, hjust = 0.5),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
        axis.ticks = element_line(color = "black", size = 0.5, linetype = "solid"),
        legend.title = element_blank()) +
  geom_hline(yintercept = -log10(Pval_Cutoff), linetype="dashed", color = "grey", size = 0.2) +
  geom_vline(xintercept = c(0, log2(FC_Cutoff),-log2(FC_Cutoff)), linetype="dashed", color = "grey", size = 0.2) +
  geom_point(POI, mapping = aes(x=Oxy60vs21_Log2FC, y =-log10(Oxy60vs21_FDR)), 
             size = 1, color = "#7B241C") +
  geom_label_repel(POI, mapping = aes(x=Oxy60vs21_Log2FC, y =-log10(Oxy60vs21_FDR)), 
                   label = POI$Genes, color = "#7B241C",max.overlaps = 40) +
  xlim(-5,5) +
  ylim(0,20)


dev.off()
