print(5)
gplots::heatmap.2(binary_mat,trace='none',Rowv = dend,Colv = FALSE,dendrogram = 'row',
                  RowSideColors = color_df$color)
#browser()