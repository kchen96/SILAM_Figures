library(dplyr)
library(reshape2)
library(ggplot2)

####### 
# Data loading and preprosessing
dfl = read.csv('../Lung/output_7/Peptide_fit_statistics.csv')
dfh = read.csv('../Heart/output_7/Peptide_fit_statistics.csv')
dfb = read.csv('../Brain/output_7/Peptide_fit_statistics.csv')


dfl_s = melt(data = dfl[,c('Protein.Group','Halflife_8','Halflife_21','Halflife_60')], 
             id.var = c('Protein.Group'), variable.name = 'Oxygen') %>% na.omit()
dfl_s['Tissue'] = 'Lung'

dfh_s = melt(data = dfh[,c('Protein.Group','Halflife_8','Halflife_21','Halflife_60')], 
             id.var = c('Protein.Group'), variable.name = 'Oxygen') %>% na.omit()
dfh_s['Tissue'] = 'Heart'

dfb_s = melt(data = dfb[,c('Protein.Group','Halflife_8','Halflife_21','Halflife_60')], 
             id.var = c('Protein.Group'), variable.name = 'Oxygen') %>% na.omit()
dfb_s['Tissue'] = 'Brain'


df_s = rbind(dfl_s, dfh_s, dfb_s)
df_s['log_value'] = log2(df_s$value)

df_s$Tissue = factor(df_s$Tissue, levels = c('Lung','Heart','Brain'))


####### 
# Half life plotting
pdf('./HalfLife_plots/HalfLife_ByTissue_desnityPlots.pdf',width = 5, height=3)
ggplot(data = df_s, mapping= aes(x = log_value,  fill = Oxygen)) +
  scale_fill_manual(values=c('#305679',"#d4ebe6","#df3d46")) +
  geom_density(aes(y=..density..),alpha=.35) +
  facet_grid(Tissue~.)+
  xlim(0,6) + 
  xlab('Log2 Half-life') +
  theme_classic()
dev.off()

pdf('./HalfLife_plots/HalfLife_ByOxygen_desnityPlots.pdf',width = 5, height=5)
ggplot(data = df_s, mapping= aes(x = log_value, color = Tissue, fill = Tissue)) +
  geom_density(aes(y=..density..),alpha=.2) +
  facet_grid(Oxygen~.)+
  xlim(0,6) + 
  xlab('Log2 Half-life') +
  theme_classic()
dev.off()

pdf('./HalfLife_plots/HalfLife_ByOxygen_density_Plots_Heart.pdf',width = 5, height=5)
ggplot(data = df_s[df_s$Tissue == 'Heart',], mapping= aes(x = value,fill= Oxygen )) +
  geom_density(aes(y=..density..),alpha=.2) +
  facet_grid(Oxygen~.)+
  xlim(0,20) + 
  xlab('Half-life') +
  theme_classic()
dev.off()

pdf('./HalfLife_plots/HalfLife_ByOxygen_density_Plots_21.pdf',width = 7, height=3)
ggplot(data = df_s[df_s$Oxygen == 'Halflife_21',], mapping= aes(x = value,fill= Tissue )) +
  geom_density(aes(y=..density..),alpha=.4) +
  scale_fill_manual(values=c('#305679',"#9d2128","#e8d8a5")) +
  xlim(0,20) + 
  xlab('Half-life (Days)') +
  theme_classic()
dev.off()

pdf('./HalfLife_plots/HalfLife_Boxplot_Heart.pdf',width = 3.2, height=2.4)
ggplot(data = df_s[df_s$Tissue == 'Heart',], mapping= aes(x = Oxygen, y = log_value, fill= Oxygen)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values=c('#305679',"#d4ebe6",'#df3d46')) +
  ylim(0,6) + 
  ylab('Log2 Half-life (Days)') +
  theme(legend.position = "None") +
  theme_classic()
dev.off()
