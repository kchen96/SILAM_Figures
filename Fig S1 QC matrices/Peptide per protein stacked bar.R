library(ggplot2)
library(dplyr)
library(openxlsx)
library(reshape2)

df = read.xlsx('Peptide_per_protein.xlsx')

df_sub = df[df$Oxygen == 'Total',c(1,3:6)]
df_sub$`1-2.peptides` = df_sub$`1-2.peptides`/df_sub$Sum
df_sub$`3-5.peptides` = df_sub$`3-5.peptides`/df_sub$Sum
df_sub$`>5.peptides` = df_sub$`>5.peptides`/df_sub$Sum
df_sub = df_sub[,1:4]

df_s = melt(df_sub, id.var = c('Tissue'),
            variable.name = 'var')
df_s$var = factor(df_s$var, levels = c('>5.peptides',
                                       '3-5.peptides',
                                       '1-2.peptides'))
df_s$Tissue = factor(df_s$Tissue,levels=c('Lung','Heart','Brain'))

pdf('Peptide_per_protein.pdf',width = 4,height = 3.6)
ggplot(data = df_s, mapping = aes(x = Tissue, y = value, fill = var)) +
  geom_bar(position = 'stack',stat='identity') +
  scale_fill_manual(values=c('#305679','#d4ebe6',"#7193a7")) +
  ylab('Proportion')+
  xlab('') + 
  theme_classic()
dev.off()
