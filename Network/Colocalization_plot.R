library(ggplot2)
library(ggridges)
library(tidyverse)
library(reshape)

setwd("/Users/jinyha/Documents/NTDnew//Analysis/Network/Connectivity")

#Colocalization - Fig.2b
df = read.csv("ntd_coloc_100k_iter_0.8contsize.txt", sep = '\t')
df$cohort = factor(df$cohort, levels = c("Control", "MM") )
head(df)
dim(df)
head(df)

ggplot(df) + 
  geom_histogram(aes(x = connection, y = ..density..,fill = cohort), 
                 alpha = 0.9, binwidth = 6, position = "identity")+ 
  scale_fill_manual(values = c( '#285059', '#A8361D')) + 
  geom_vline(data = subset(df, df$cohort == 'MM'), aes(xintercept=median(connection), color= 'white'),
             linetype="dashed")+
  geom_vline(data = subset(df, df$cohort == 'Control'), aes(xintercept=median(connection), color= 'white'),
             linetype="dashed")+
  theme_classic() + 
  xlim(c(430,860))+
  labs(fill='Cohort', x='Number of edges', y='Density') +
  theme(legend.position='right', 
        aspect.ratio = 0.5,
        axis.title.x = element_text( size = 9,color = "black", face = "bold"),
        axis.title.y = element_text( angle = 90,size = 9,color = "black", face = "bold"),
        axis.text.x = element_text( size = 7, color = "black", face = "bold"),
        axis.text.y = element_text(size = 7,color = "black", face = "bold" ))

ggsave("Network_colocalization.pdf", plot= last_plot(), width = 4, height = 3)



df_mm = subset(df, df$cohort == 'MM')
df_ssc = subset(df, df$cohort == 'Control')


wilcox.test(df_mm$connection,df_ssc$connection)


##Edge per node Fig.2c
df_pernode = read.csv("per_node_ntd_100k.txt", sep = '\t')
df_pernode = subset(df_pernode, df_pernode$edge_per_node < 11)
head(df_pernode)

df_pernode$cohort = factor(df_pernode$cohort, levels = c("Control", "MM") )

df_pernode_w = data.frame()
for (i in 0:max(df_pernode$edge_per_node)){
  df_sub = subset(df_pernode, df_pernode$edge_per_node == i)
  var = var.test(subset(df_sub, df_sub$cohort== 'MM')$edge_cnt , subset(df_sub, df_sub$cohort == 'Control')$edge_cnt)
  wt = wilcox.test(subset(df_sub, df_sub$cohort== 'MM')$edge_cnt , subset(df_sub, df_sub$cohort == 'Control')$edge_cnt, alternative = "greater")
  print (wt)
  #Getting quartiles 
  mm_1st = as.numeric(quantile(subset(df_sub, df_sub$cohort== 'MM')$edge_cnt,0.25))
  mm_3rd = as.numeric(quantile(subset(df_sub, df_sub$cohort== 'MM')$edge_cnt,0.75))
  ssc_1st = as.numeric(quantile(subset(df_sub, df_sub$cohort== 'Control')$edge_cnt,0.25))
  ssc_3rd = as.numeric(quantile(subset(df_sub, df_sub$cohort== 'Control')$edge_cnt,0.75))
  mm_median = as.numeric(quantile(subset(df_sub, df_sub$cohort== 'MM')$edge_cnt,0.5))
  mm_lower = mm_1st
  mm_upper = mm_3rd
  ssc_median =  as.numeric(quantile(subset(df_sub, df_sub$cohort== 'Control')$edge_cnt,0.5))
  ssc_lower = ssc_1st
  ssc_upper = ssc_3rd
  result = c(i,   "MM", mm_median, mm_1st, mm_3rd, wt$p.value)
  df_pernode_w =rbind(df_pernode_w, result)
  result2 = c(i,"Control", ssc_median, ssc_1st, ssc_3rd, wt$p.value)
  df_pernode_w =rbind(df_pernode_w, result2)
}
df_pernode_w
colnames(df_pernode_w) = c("Turn",  "cohort", "median", "Q1", "Q3", "wilcox pval")

df_pernode_w



df_pernode_w$pvalue = as.numeric(df_pernode_w$pvalue)
df_pernode_w$median = as.numeric(df_pernode_w$median)
df_pernode_w$Q1 = as.numeric(df_pernode_w$Q1)
df_pernode_w$Q3 = as.numeric(df_pernode_w$Q3)

df_pernode_w$Turn = factor(df_pernode_w$Turn, levels= c("0","1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

df_pernode_w$median = as.numeric(df_pernode_w$median)


ggplot(data=df_pernode_w, aes(x=Turn, y= as.numeric(median), color = cohort, shape = cohort)) + 
  geom_pointrange(aes(ymin=as.numeric(Q1), ymax=as.numeric(Q3)), position = position_dodge(0.8)) + 
  scale_color_manual(values = c( '#285059', '#A62F14')) + 
 # ylim(c(0,16)) + 
  labs(fill='Cohort', x='Edge per node', y='Occurence') +
  theme_bw() +
  theme(legend.position='right', 
        aspect.ratio = 0.5,
        axis.title.x = element_text( size = 9,color = "black", face = "bold"),
        axis.title.y = element_text( angle = 90,size = 9,color = "black", face = "bold"),
        axis.text.x = element_text( size = 7, color = "black", face = "bold"),
        axis.text.y = element_text(size = 7,color = "black", face = "bold" )) + 
  ggtitle("Edge count per node")
  
ggsave("Colocalization_edge_per_node.pdf", plot= last_plot(), width = 5, height = 5)



