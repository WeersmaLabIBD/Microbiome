library(pheatmap)
library(randomcoloR)


tmp.direction <- tmp.pvalue 
tmp.direction[tmp.direction < 0.05] = "*"
tmp.direction[tmp.direction > 0.05] = ""

# set colors
breaksList = seq(-max(abs(tmp.coef)), +max(abs(tmp.coef)),by=2 * max(abs(tmp.coef))/55)
cols <- colorRampPalette(c("#EECC66","White","#6699CC"))(length(breaksList))
cols[seq(length(cols)/2-1,length(cols)/2+1)] <- "White"

anno_color=data.frame("Metabolite"=colnames(tmp.coef))
anno_color=merge(anno_color,annotate_fecal,by="Metabolite",all=F)
anno_col=data.frame("Class"=(unique(anno_color$Class)),"color"=distinctColorPalette(length(unique(anno_color$Class))))
rownames(anno_color)=anno_color$Metabolite
anno_color$Metabolite=NULL

Var1        <- anno_col$color
names(Var1) <- anno_col$Class
anno_colors <- list(Var1 = Var1)

pheatmap(tmp.coef,legend = T,show_rownames = T, cluster_rows = T,angle_col = 90,display_numbers = tmp.direction,show_colnames = T,cluster_cols = T,
         fontsize_number = 10,border_color = "#EEEEEE",na_col = "white",
         color = cols,legend_labels = "log(p-value)",
         fontsize_col = 5,fontsize_row = 5,breaks = breaksList,
         annotation_col = anno_color,annotation_row = anno_color,annotation_colors = (anno_colors))

