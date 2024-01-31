library(ggplot2)
library(dplyr)
library(GGally)

dfl = read.csv('../Lung/output_7/peptide_fit.csv', row.names = 1)
dfh = read.csv('../Heart/Output_7/peptide_fit.csv', row.names = 1)
dfb = read.csv('../Brain/Output_7/peptide_fit.csv', row.names = 1)

intersect(dfl['Protein.Group'], dfh['Protein.Group']) %>% dim()  #2066
intersect(dfl['Protein.Group'], dfb['Protein.Group']) %>% dim()  #2491
intersect(dfh['Protein.Group'], dfb['Protein.Group']) %>% dim()  #1928


data = inner_join(dfl[c('Protein.Group','Half_life','Oxygen')],dfh[c('Protein.Group','Half_life','Oxygen')],
            by = c('Protein.Group','Oxygen'))%>% 
  inner_join(dfb[c('Protein.Group','Half_life','Oxygen')], by = c('Protein.Group','Oxygen'))
colnames(data) = c('Protein.Group','Lung','Oxygen','Heart','Brain')
data$Oxygen = as.character(data$Oxygen)
data$Lung = log2(data$Lung)
data$Heart = log2(data$Heart)
data$Brain = log2(data$Brain)

dfsub = data[data$Oxygen ==21, c('Lung','Heart','Brain')]

pdf('./HalfLife_plots/CorrPlots_AllOxygen.pdf')
ggpairs(dfsub,ggplot2::aes( alpha = 0.6), upper = 6, lower=0)
dev.off()

pdf('./HalfLife_plots/CorrPlots_HeartVSBrain_21.pdf')
ggplot(data = dfsub, aes(x = Heart, y = Brain)) +
  geom_point(color='grey40') +
  geom_smooth(method='lm') +
  xlim(0,6) +
  ylim(0,6) +
  xlab('Log2 Half-life (Days)') +
  ylab('Log2 Half-life (Days)') +
  theme_minimal()
dev.off()

l <- lm(dfsub$Brain~ dfsub$Heart)
summary(l)

ggpairs(data, columns = c(2,4,5), ggplot2::aes(colour=Oxygen, alpha = 0.6)) 
