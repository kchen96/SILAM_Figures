library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(openxlsx)

df = read.xlsx('HalfLife_Percentage.xlsx')

df_scaled = df[df$Tissue == 'Brain',3:8]

## hierarchical clustering
#heatmap(as.matrix(df_log), scale = "row")

pdf('./HalfLife_plots/HalfLife_percentage_ht_Brain.pdf', width = 4, height =4)
ht = Heatmap(as.matrix(t(df_scaled)),
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             show_row_names = FALSE,
             show_column_names = FALSE,
             col = colorRamp2(c(0,0.1,0.4), c("white", "#e8d8a5",'#9d2128')),
             rect_gp = gpar(col = "white", lwd = 1))

ht = draw(ht) 
dev.off()


