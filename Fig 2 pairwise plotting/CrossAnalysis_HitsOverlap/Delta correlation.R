library(ggplot2)
library(dplyr)
library(reshape2)


####### 
# Data loading and preprosessing
dfl = read.csv('../Lung/output_7/Peptide_fit_statistics.csv')
dfh = read.csv('../Heart/output_7/Peptide_fit_statistics.csv')
dfb = read.csv('../Brain/output_7/Peptide_fit_statistics.csv')


data = inner_join(dfl[c('Protein.Group','Oxy8vs21_Log2FC','Oxy60vs21_Log2FC')],
                  dfh[c('Protein.Group','Oxy8vs21_Log2FC','Oxy60vs21_Log2FC')],
                  by = c('Protein.Group'))%>% 
  inner_join(dfb[c('Protein.Group','Oxy8vs21_Log2FC','Oxy60vs21_Log2FC')], 
             by = c('Protein.Group'))

colnames(data) = c('Protein.Group','Lung8','Lung60',
                   'Heart8','Heart60',
                   'Brain8','Brain60')

## Single tissue
pdf('./Delta-corr/Delta-corr_Lg.pdf', width = 5, height = 4)
ggplot(data = data, mapping= aes(x = Lung8, y = Lung60)) +
  geom_point(alpha = 0.6, size = 0.3) +
  geom_smooth(method=lm) +
  xlab('Log2FC 8vs21') +
  ylab('Log2FC 60vs21') +
  theme_minimal()
dev.off()

pdf('./Delta-corr/Delta-corr_Ht.pdf', width = 5, height = 4)
ggplot(data = data, mapping= aes(x = Heart8, y = Heart60)) +
  geom_point(alpha = 0.6, size = 0.3) +
  geom_smooth(method=lm) +
  xlab('Log2FC 8vs21') +
  ylab('Log2FC 60vs21') +
  theme_minimal()
dev.off()

pdf('./Delta-corr/Delta-corr_Br.pdf', width = 5, height = 4)
ggplot(data = data, mapping= aes(x = Brain8, y = Brain60)) +
  geom_point(alpha = 0.6, size = 0.3) +
  geom_smooth(method=lm) +
  xlab('Log2FC 8vs21') +
  ylab('Log2FC 60vs21') +
  theme_minimal()
dev.off()

## Across organs
pdf('./Delta-corr/Delta-corr_BrLg_8vs21.pdf', width = 5, height = 4)
ggplot(data = data, mapping= aes(x = Brain8, y = Lung8)) +
  geom_point(alpha = 0.6, size = 0.3) +
  xlab('Brain Log2FC 8vs21') +
  ylab('Lung Log2FC 8vs21') +
  theme_minimal()
dev.off()

pdf('./Delta-corr/Delta-corr_HtLg_8vs21.pdf', width = 5, height = 4)
ggplot(data = data, mapping= aes(x = Heart8, y = Lung8)) +
  geom_point(alpha = 0.6, size = 0.3) +
  xlab('Heart Log2FC 8vs21') +
  ylab('Lung Log2FC 8vs21') +
  theme_minimal()
dev.off()

pdf('./Delta-corr/Delta-corr_BrHt_8vs21.pdf', width = 5, height = 4)
ggplot(data = data, mapping= aes(x = Brain8, y = Heart8)) +
  geom_point(alpha = 0.6, size = 0.3) +
  xlab('Brain Log2FC 8vs21') +
  ylab('Heart Log2FC 8vs21') +
  theme_minimal()
dev.off()

pdf('./Delta-corr/Delta-corr_BrLg_60vs21.pdf', width = 5, height = 4)
ggplot(data = data, mapping= aes(x = Brain60, y = Lung60)) +
  geom_point(alpha = 0.6, size = 0.3) +
  xlab('Brain Log2FC 60vs21') +
  ylab('Lung Log2FC 60vs21') +
  theme_minimal()
dev.off()

pdf('./Delta-corr/Delta-corr_HtLg_60vs21.pdf', width = 5, height = 4)
ggplot(data = data, mapping= aes(x = Heart60, y = Lung60)) +
  geom_point(alpha = 0.6, size = 0.3) +
  xlab('Heart Log2FC 60vs21') +
  ylab('Lung Log2FC 60vs21') +
  theme_minimal()
dev.off()

pdf('./Delta-corr/Delta-corr_BrHt_60vs21.pdf', width = 5, height = 4)
ggplot(data = data, mapping= aes(x = Brain60, y = Heart60)) +
  geom_point(alpha = 0.6, size = 0.3) +
  xlab('Brain Log2FC 60vs21') +
  ylab('Heart Log2FC 60vs21') +
  theme_minimal()
dev.off()
