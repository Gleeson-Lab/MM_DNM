setwd("/Users/jinyha/Documents/NTDnew/Analysis/Merfish")
library(ggplot2)
install.packages("reshape")
library(reshape)

df = read.csv("36gene_expression.txt", sep = '\t')
head(df)
df$X
df.melt = melt(df,id = "X")
head(df.melt)

celltype = c("Indeterminate", "Blood","Dorsal Root Ganglia","Mesoderm","Neural Crest", "Pre-EMT-NCP", "Neural Progenitor", "Neuron"  )

df.melt$X= factor(df.melt$X, levels=celltype)
ggplot(df.melt, aes(fill=X, y=value, x=variable)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_x_discrete(limits = rev) + 
  scale_fill_manual(values= c("#A4A1A6", "#253659", "#698EBF", "#B9BF04", "#D97652", "#D9C24E",  "#D0E2F2","#734870" )) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8, vjust = 0.4, color = "black", face = "bold"),
        axis.text.y = element_text(size = 8,color = "black", face = "bold"))+ 
  ggtitle("Expression of 36 damaging DNM genes in MERFISH") + 
  labs(fill='DNM gene', x='', y='Expression (%)') 

ggsave("Merfish_expression.pdf", plot= last_plot(), units = "cm")
