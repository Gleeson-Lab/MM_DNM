library(ggplot2)
library(ggpubr)
library(forcats)
setwd("/Users/jinyha/Documents/NTDnew/Analysis/Network/propagation")



df = read.csv("Submodule_GO_plot_e8.txt", sep = '\t')

p1 <- ggplot(subset(df, df$Cluster == "Rho GTPase")) + 
  geom_point(aes(x=-log10(FDR.value), y=fct_reorder(description, -log10(FDR.value)), color = -log10(FDR.value), size = X..genes), fill = "black") + 
  theme_bw() +
  xlim(c(7,14)) + 
  scale_color_continuous_sequential(palette = "Reds", l1 = 20, c2 = 70, p1 = 1) + 
  theme(legend.position='bottom', 
        aspect.ratio = 3,
        axis.title.y.right = element_text( angle = 90,size = 12),
        axis.text.x = element_text(angle = -45, hjust = 0.01,size = 10, color = "black", face = "bold"),
        axis.text.y = element_text(size = 10, color = "black", face = "bold")) + 
  ggtitle("GTPase mediated actin cytoskeleton")
p1
p2 <- ggplot(subset(df, df$Cluster == "Cillum")) + 
  geom_point(aes(x=-log10(FDR.value), y=fct_reorder(description, -log10(FDR.value)), color = -log10(FDR.value), size = X..genes)) + 
  theme_bw() +
  xlim(c(8,17))+
  scale_color_continuous_sequential(palette = "Reds", l1 = 20, c2 = 70, p1 = 1) + 
  theme(legend.position='bottom', 
        aspect.ratio = 3,
        axis.title.y.right = element_text( angle = 90,size = 12),
        axis.text.x = element_text(angle = -45, hjust = 0.01,size = 10, color = "black", face = "bold"),
        axis.text.y = element_text(size = 10, color = "black", face = "bold")) + 
  ggtitle("Microtubule cytoskeleton")
p2
combined <- ggarrange(p1, p2, common.legend = T, ncol = 1)
combined
ggsave("GO_module.pdf", plot= combined, width = 4, height = 10)
 
