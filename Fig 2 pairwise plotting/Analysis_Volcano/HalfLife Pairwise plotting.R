######### Half-life pair-wise comparison 

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggrepel)

df_byPro = read.csv('../output_7/Peptide_fit_statistics.csv')

## Volcano Plots
### Hypoxia vs Normoxia
FC_Cutoff = 1.3
Pval_Cutoff = 0.05

Sig_fast = subset(df_byPro, (Oxy8vs21_Log2FC>log2(FC_Cutoff)) & (Oxy8vs21_FDR<Pval_Cutoff))
Sig_slow = subset(df_byPro, (Oxy8vs21_Log2FC<(-log2(FC_Cutoff)))& (Oxy8vs21_FDR<Pval_Cutoff))


## Volcano
pdf('Oxy8vs21_Volcano.pdf', width = 3.5, height = 3.5)
ggplot(data = df_byPro, aes(x=Oxy8vs21_Log2FC, y =-log10(Oxy8vs21_FDR))) +
  geom_point(size = 0.5, color = "grey") +
  theme_classic() +
  xlab("Log2(Fold-Change)") + ylab("-log10(FDR)") +
  theme(plot.title = element_text(colour = "black", size = 16, hjust = 0.5),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
        axis.ticks = element_line(color = "black", size = 0.5, linetype = "solid"),
        legend.title = element_blank()) +
  geom_hline(yintercept = -log10(Pval_Cutoff), linetype="dashed", color = "grey", size = 0.2) +
  geom_vline(xintercept = c(-log2(FC_Cutoff),0, log2(FC_Cutoff)), linetype="dashed", color = "grey", size = 0.2) +
  geom_point(Sig_fast, mapping = aes(x=Oxy8vs21_Log2FC, y =-log10(Oxy8vs21_FDR)), 
             size = 0.5, color = "#9d2128") +
  geom_point(Sig_slow, mapping = aes(x=Oxy8vs21_Log2FC, y =-log10(Oxy8vs21_FDR)), 
             size = 0.5, color = "#305679") +
  xlim(-5,5) +
  ylim(0,20)
dev.off()

## Pie chart
f = sum(df_byPro$Protein.Group %in% Sig_fast$Protein.Group)
s = sum(df_byPro$Protein.Group %in% Sig_slow$Protein.Group)
df = data.frame(
  group =c('Fast','Others','Slow'),
  value =c(f,dim(df_byPro)[1]-f-s, s)
)
df$fraction = df$value/dim(df_byPro)[1] *100
df$lab.ypo = cumsum(df$fraction) - 0.5*df$fraction

pdf('Oxy8vs21_Pie.pdf', width = 3, height = 3)
ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity",alpha = 0.9)+
  coord_polar("y", start=0)+
  scale_fill_manual(values = c( "#9d2128",'grey',"#305679")) +
  theme_void()
dev.off()

### Hyperoxia vs Normoxia
FC_Cutoff =1.3
Pval_Cutoff = 0.05

Sig_fast = subset(df_byPro, (Oxy60vs21_Log2FC>log2(FC_Cutoff)) & (Oxy60vs21_FDR<Pval_Cutoff))
Sig_slow = subset(df_byPro, (Oxy60vs21_Log2FC<(-log2(FC_Cutoff)))& (Oxy60vs21_FDR<Pval_Cutoff))


pdf('Oxy60vs21_Volcano.pdf', width = 3.5, height = 3.5)
ggplot(data = df_byPro, aes(x=Oxy60vs21_Log2FC, y =-log10(Oxy60vs21_FDR))) +
  geom_point(size = 0.5, color = "grey") +
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
  geom_point(Sig_fast, mapping = aes(x=Oxy60vs21_Log2FC, y =-log10(Oxy60vs21_FDR)), 
             size = 0.5, color = "#9d2128") +
  geom_point(Sig_slow, mapping = aes(x=Oxy60vs21_Log2FC, y =-log10(Oxy60vs21_FDR)), 
             size = 0.5, color = "#305679") +
  xlim(-5,5) +
  ylim(0,20)
dev.off()


## Pie chart
f = sum(df_byPro$Protein.Group %in% Sig_fast$Protein.Group)
s = sum(df_byPro$Protein.Group %in% Sig_slow$Protein.Group)
df = data.frame(
  group =c('Fast','Others','Slow'),
  value =c(f,dim(df_byPro)[1]-f-s, s)
)
df$fraction = df$value/dim(df_byPro)[1] *100
df$lab.ypo = cumsum(df$fraction) - 0.5*df$fraction

pdf('Oxy60vs21_Pie.pdf', width = 3, height = 3)
ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity",alpha = 0.9)+
  coord_polar("y", start=0)+
  scale_fill_manual(values = c( "#9d2128",'grey',"#305679")) +
  theme_void()
dev.off()

