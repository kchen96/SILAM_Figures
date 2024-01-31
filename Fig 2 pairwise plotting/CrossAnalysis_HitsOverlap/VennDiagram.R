library(ggplot2)
library(dplyr)
library(ggVennDiagram)

lg = read.csv('../Lung/output_7/Peptide_fit_statistics.csv')
ht = read.csv('../Heart/Output_7/Peptide_fit_statistics.csv')
br = read.csv('../Brain/Output_7/Peptide_fit_statistics.csv')

p_cutoff = 0.05
Log2FC_cutoff =  log2(1.3)

CommonProteins = intersect(lg$Protein.Group,ht$Protein.Group)
CommonProteins = intersect(CommonProteins,br$Protein.Group)

Hypo_slow = list(Lung = lg %>% filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy8vs21_Log2FC<(-Log2FC_cutoff)  & Oxy8vs21_FDR<p_cutoff)  %>% 
                   select(Protein.Group) %>% unlist,
                 Brain = br %>% filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy8vs21_Log2FC<(-Log2FC_cutoff)  & Oxy8vs21_FDR<p_cutoff) %>% 
                   select(Protein.Group)  %>% unlist,
                 Heart = ht %>% filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy8vs21_Log2FC<(-Log2FC_cutoff ) &Oxy8vs21_FDR<p_cutoff)  %>% 
                   select(Protein.Group) %>% unlist)
pdf('./Venn/Hypo_slow.pdf', width=4.5, height = 3.6)
ggVennDiagram(Hypo_slow) +
  scale_fill_gradient(low = "white", high = "#305679")+
  scale_color_manual(values = c('#00121a','#e8d8a5','#9d2128'))
dev.off()

Hypo_fast = list(Lung = lg %>% filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy8vs21_Log2FC>Log2FC_cutoff  & Oxy8vs21_FDR<p_cutoff)  %>% 
                   select(Protein.Group) %>% unlist,
                 Brain = br %>% filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy8vs21_Log2FC>Log2FC_cutoff  & Oxy8vs21_FDR<p_cutoff) %>% 
                   select(Protein.Group)  %>% unlist,
                 Heart = ht %>% filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy8vs21_Log2FC>Log2FC_cutoff  &Oxy8vs21_FDR<p_cutoff)  %>% 
                   select(Protein.Group) %>% unlist)
pdf('./Venn/Hypo_fast.pdf', width=4.5, height = 3.6)
ggVennDiagram(Hypo_fast) +
  scale_fill_gradient(low = "white", high = "#a6373d") +
  scale_color_manual(values = c('#00121a','#e8d8a5','#9d2128'))
dev.off()

Hypr_slow = list(Lung = lg %>%filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy60vs21_Log2FC<(-Log2FC_cutoff)  & Oxy60vs21_FDR<p_cutoff)  %>% 
                   select(Protein.Group) %>% unlist,
                 Brain = br %>%filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy60vs21_Log2FC<(-Log2FC_cutoff ) & Oxy60vs21_FDR<p_cutoff) %>% 
                   select(Protein.Group)  %>% unlist,
                 Heart = ht %>%filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy60vs21_Log2FC<(-Log2FC_cutoff ) &Oxy60vs21_FDR<p_cutoff)  %>% 
                   select(Protein.Group) %>% unlist)
pdf('./Venn/Hypr_slow.pdf', width=4.5, height = 3.6)
ggVennDiagram(Hypr_slow) +
  scale_fill_gradient(low = "white", high = "#305679") +
  scale_color_manual(values = c('#00121a','#e8d8a5','#9d2128'))
dev.off()

Hypr_fast = list(Lung = lg %>% filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy60vs21_Log2FC>Log2FC_cutoff  & Oxy60vs21_FDR<p_cutoff)  %>% 
                   select(Protein.Group) %>% unlist,
                 Brain = br %>% filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy60vs21_Log2FC>Log2FC_cutoff  & Oxy60vs21_FDR<p_cutoff) %>% 
                   select(Protein.Group)  %>% unlist,
                 Heart = ht %>% filter(Protein.Group %in% CommonProteins) %>% 
                   filter(Oxy60vs21_Log2FC>Log2FC_cutoff  &Oxy60vs21_FDR<p_cutoff)  %>% 
                   select(Protein.Group) %>% unlist)

pdf('./Venn/Hypr_fast.pdf', width=4.5, height = 3.6)
ggVennDiagram(Hypr_fast) +
  scale_fill_gradient(low = "white", high = "#a6373d") +
  scale_color_manual(values = c('#00121a','#e8d8a5','#9d2128'))
dev.off()
