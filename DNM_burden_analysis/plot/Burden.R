library(ggplot2)
library(ggpubr)
library(ggsignif)    


#setwd("/Users/jinyha/Documents/NTDnew/Analysis/DNM/Burden/ForPlot/")
df_all_stat <- read.csv("Stat.txt", header= T, sep ='\t')
df_rr <- read.csv("RR_total.txt", header= T, sep ='\t')
total_size = 59281518 #size of coding region (hg38)
head(df_all_stat)


figure_a = c('LGD','DMis','DMisMC2', 'DamagingS', 'DamagingMC2', 'OtherMissense', 'Synonymous')
label_a = c('LGD','D-Mis (high-conf.)', 'D-Mis (all)','Damaging (high-conf.)','Damaging (all)','Tolerant Mis', 'Synonymous')

df_all_stat$category <- factor(df_all_stat$category, levels = figure_a)
df_all_stat$cohort <- factor(df_all_stat$cohort, levels = c("SSC", "MM"))
df_all_stat$left = as.numeric(df_all_stat$mean) - as.numeric(df_all_stat$lower)
df_all_stat$right = as.numeric(df_all_stat$upper) - as.numeric(df_all_stat$mean)
 
#Constrained genes
p_con <- ggplot(subset(df_all_stat, df_all_stat$category %in% figure_a & Gene =="Constrained genes"), aes(category, as.numeric(mean))) +
  geom_errorbar(aes(x = category, ymin=lower, ymax = upper,group = cohort,color = category),
                position = position_dodge(0.5), width = 0.2,size = .3) + 
  geom_point(data = subset(df_all_stat, df_all_stat$category %in% figure_a & Gene =="Constrained genes"),
             aes(category, as.numeric(mean), group = cohort, color = category, shape = cohort),
             size =2.5,
             position =position_dodge(0.5)) + 
  scale_color_manual(values = c("#642E44", "#A65D7C",  "#E5A3B1", "#D94929", "#A62F14", "#C2C0A6", "#889C9B", "#B2BEBF")) +
  scale_x_discrete(breaks=figure_a,labels=label_a) + 
  scale_y_continuous(sec.axis = sec_axis(~ . / (100000000 / 2) * total_size, name = 'Theoretical Variant Rate (per child)'), limits= c(-0.005,0.14)) +
  labs(fill='Cohort', x='', y='Variant rate (x 10e-8)') +
  theme_classic() +
  theme(legend.position='bottom', 
        aspect.ratio = 0.6,
        axis.title.y.right = element_text( angle = 90,size = 12),
        axis.text.x = element_text(angle = -45, hjust = 0.01,size = 10, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12))+ 
  geom_text(data =subset(df_rr,  category %in% figure_a & Gene == 'Constrained'), aes(x = category, y = 0.13,
     label = paste(paste("RR=",round(as.numeric(RR),2), sep=""), 
     paste("p=",ifelse(pvalue>0.01, round(pvalue,2), formatC(pvalue, format = "e", digits = 1)),sep="" ), sep='\n')), size = 3) + 
  geom_text(data =subset(df_rr,  category %in% figure_a & Gene == 'Constrained'), 
            aes(x = category, y =0.1, label = ifelse(pvalue > 0.05, "NS", ifelse(pvalue>0.01, "*", ifelse(pvalue > 0.001, "**", ifelse(pvalue<= 0.001, "***",""))))), size =3) + 
  ggtitle("Constrained genes (pLI ≥ 0.9)")

p_con

#All genes
p_all <- ggplot() +
  geom_errorbar(data = subset(df_all_stat, df_all_stat$category %in% figure_a & Gene =="All genes"),aes(x = category, ymin=lower, ymax = upper,group = cohort,color = category),
                position = position_dodge(0.5), width = 0.2,size = .3) + 
  geom_point(data = subset(df_all_stat, df_all_stat$category %in% figure_a & Gene =="All genes"),
             aes(category, as.numeric(mean), group = cohort, color = category, shape = cohort),
             size =2.5,
             position =position_dodge(0.5)) + 
  scale_color_manual(values = c("#642E44", "#A65D7C",  "#E5A3B1", "#D94929", "#A62F14", "#C2C0A6", "#889C9B", "#B2BEBF")) +
  scale_x_discrete(breaks=figure_a,labels=label_a) + 
  scale_y_continuous(sec.axis = sec_axis(~ . / (100000000 / 2) * total_size, name = 'Theoretical Variant Rate (per child)'),limits = c(-0.005,0.45)) +
  labs(fill='Cohort', x='', y='Variant rate (x 10e-8)') +
  theme_classic() +
  theme(legend.position='bottom', 
        aspect.ratio = 0.6,
        axis.title.y.right = element_text( angle = 90,size = 12),
        axis.text.x = element_text(angle = -45, hjust = 0.01,size = 10, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12))+ 
  geom_text(data =subset(df_rr,  category %in% figure_a & Gene == 'All'), aes(x = category, y = 0.4,
     label = paste(paste("RR=",round(as.numeric(RR),2), sep=""), 
     paste("p=",ifelse(pvalue>0.01, round(pvalue,2), formatC(pvalue, format = "e", digits = 1)),sep="" ), sep='\n')), size = 3) + 
  geom_text(data =subset(df_rr,  category %in% figure_a & Gene == 'All'), 
            aes(x = category, y =0.35, label = ifelse(pvalue > 0.05, "NS", ifelse(pvalue>0.01, "*", ifelse(pvalue > 0.001, "**", ifelse(pvalue<= 0.001, "***",""))))), size = 3) + 
 
  ggtitle("All genes")
p_all

p_combined <- ggarrange(p_all, p_con, common.legend = TRUE)
p_combined

ggsave("Rate_burden.pdf", plot= p_combined, width = 12, height = 6)



