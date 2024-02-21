setwd("/Users/jinyha/Documents/NTDnew/Analysis/Validation/Raw_data/") 
library(ggplot2)
library(reshape2)
library(plyr)


#TNK2
df_tnk2_IP = read.csv("TNK2_IP_Kinase.txt", sep = '\t')
df_tnk2_IP.format = aggregate(data=df_tnk2_IP, CPM ~ Type, mean)
df_tnk2_IP.format$sd = aggregate(data=df_tnk2_IP, CPM ~ Type, sd)
df_tnk2_IP.format$length = aggregate(data=df_tnk2_IP, CPM ~ Type, length)
df_tnk2_IP.format$sem = df_tnk2_IP.format$sd$CPM / sqrt(df_tnk2_IP.format$length$CPM)
df_tnk2_IP.format
table(df_tnk2_IP$Type)
condition = c("Mock", "A156T TNK2", "WT TNK2", "P168L TNK2")
label=c("Mock", "A156T TNK2","WT TNK2", "P168L TNK2")
df_tnk2_IP.format$Type <- factor(df_tnk2_IP.format$Type, levels = condition)
df_tnk2_IP.format
#Statistical test
df_tnk2_IP_vector= subset(df_tnk2_IP, df_tnk2_IP$Type == 'Mock')
df_tnk2_IP_wt= subset(df_tnk2_IP, df_tnk2_IP$Type == 'WT TNK2')
df_tnk2_IP_mt= subset(df_tnk2_IP, df_tnk2_IP$Type == 'P168L TNK2')
df_tnk2_IP_kd= subset(df_tnk2_IP, df_tnk2_IP$Type == 'A156T TNK2')
#Not possible, data point only two each


p <- ggplot(df_tnk2_IP.format, aes(x = Type, y = CPM, fill = Type)) +
  geom_bar(position=position_dodge(width=0.9) , stat="identity") +
  geom_jitter(data= df_tnk2_IP[df_tnk2_IP$Type=='Mock',], fill="black", size=1, alpha = 0.4) +
  geom_jitter(data= df_tnk2_IP[df_tnk2_IP$Type=='WT TNK2',], colour="black", size=1, alpha = 0.4, ) +
  geom_jitter(data= df_tnk2_IP[df_tnk2_IP$Type=='P168L TNK2',], colour="black", size=1, alpha = 0.4)+
  geom_jitter(data= df_tnk2_IP[df_tnk2_IP$Type=='A156T TNK2',], colour="black", size=1, alpha = 0.4)+
  geom_errorbar(aes(ymax=CPM+sem, ymin=CPM-sem), position=position_dodge(width=0.9) , width=0.6, size = 0.5)+
  scale_fill_manual(values = c("#AFAFB2","#4B9D96", "#4B704D", "#C8972E" )) + 
  scale_color_manual(values = c("#AFAFB2","#4B9D96", "#4B704D", "#C8972E" )) + 
  scale_x_discrete(breaks=condition,labels=label) + 
  theme_classic() +
  #ylim(c(0,3)) + 
  theme(legend.position = "none",
        aspect.ratio = 1.4,
        axis.text.x = element_text(angle = 45, size = 6, vjust = 0.4, color = "black", face = "bold"),
        axis.text.y = element_text(size = 6,color = "black", face = "bold"))+ 
  ggtitle("WASP IP-Kinase reaction") + 
  labs(fill='Cohort', x='', y='cpm-blank') 
p

ggsave("Main_tnk2_IP.pdf", plot= p, width = 6, height = 9, units = "cm")



df <- ddply(df, .(col1, col2, col3), transform, colSEM = colSD/sqrt(length(na.omit(colSD))))

####################
#Tiam1 actin 
####################

####################

df_tiam1_phallo = read.csv("TIAM1_Rac1_actin.txt", sep = '\t')
colnames(df_tiam1_phallo) <- c("Set", "Phallo", "Type")
df_tiam1_phallo.format = aggregate(data=df_tiam1_phallo, Phallo ~ Type, mean)
df_tiam1_phallo.format$sd = aggregate(data=df_tiam1_phallo, Phallo ~ Type, sd)
df_tiam1_phallo.format$length = aggregate(data=df_tiam1_phallo, Phallo ~ Type, length)
df_tiam1_phallo.format$sem = df_tiam1_phallo.format$sd$Phallo / sqrt(df_tiam1_phallo.format$length$Phallo)
df_tiam1_phallo.format
condition = c("Vector", "WT Tiam1", "H1149P Tiam1")
label=c("Vector", "WT Tiam1", "H1149P Tiam1")
df_tiam1_phallo.format$Type <- factor(df_tiam1_phallo.format$Type, levels = condition)

#Statistical test
df_tiam1_phallo_vector= subset(df_tiam1_phallo, df_tiam1_phallo$Type == 'Vector')
df_tiam1_phallo_wt= subset(df_tiam1_phallo, df_tiam1_phallo$Type == 'WT Tiam1')
df_tiam1_phallo_mt= subset(df_tiam1_phallo, df_tiam1_phallo$Type == 'H1149P Tiam1')
shapiro.test(df_tiam1_phallo_vector$Phallo) # p-value =  2.933e-14
shapiro.test(df_tiam1_phallo_wt$Phallo) #p-value = 2.392e-10
shapiro.test(df_tiam1_phallo_mt$Phallo) #p-value = 0.0001924
#Outlier estimation by Hampel filter, consists of considering as outliers the values outside the interval (I) formed by the median, plus or minus 3 median absolute deviations.
upper_vector = median(df_tiam1_phallo_vector$Phallo) + 3 * mad(df_tiam1_phallo_vector$Phallo, constant = 1)
upper_wt = median(df_tiam1_phallo_wt$Phallo) + 3 * mad(df_tiam1_phallo_wt$Phallo, constant = 1)
upper_mt = median(df_tiam1_phallo_mt$Phallo) + 3 * mad(df_tiam1_phallo_mt$Phallo, constant = 1)


kruskal.test(data=df_tiam1_phallo, Phallo ~ Type)
pairwise.wilcox.test(df_tiam1_phallo$Phallo, df_tiam1_phallo$Type, p.adjust.method = "bonferroni")

p <- ggplot(df_tiam1_phallo.format, aes(x = Type, y = Phallo, fill = Type)) 
position_dodge(width=0.9) 
limits <- aes(ymax=value+sd, ymin=value-sd) #Set up the error bars
df_tiam1_phallo.format$mean_se

p <- ggplot(df_tiam1_phallo.format, aes(x = Type, y = Phallo, fill = Type)) +
  geom_bar(position=position_dodge(width=0.9) , stat="identity") +
  geom_jitter(data= subset(df_tiam1_phallo, df_tiam1_phallo$Phallo <= upper_vector)[df_tiam1_phallo$Type=='Vector',], fill="black", size=0.5, alpha = 0.4) +
  geom_jitter(data= subset(df_tiam1_phallo, df_tiam1_phallo$Phallo <= upper_wt)[df_tiam1_phallo$Type=='WT Tiam1',], colour="black", size=0.5, alpha = 0.4, ) +
  geom_jitter(data= subset(df_tiam1_phallo, df_tiam1_phallo$Phallo <= upper_mt)[df_tiam1_phallo$Type=='H1149P Tiam1',], colour="black", size=0.5, alpha = 0.4)+
  geom_errorbar(aes(ymax=Phallo+sem, ymin=Phallo-sem), position=position_dodge(width=0.9) , width=0.6, size = 0.5)+
  scale_fill_manual(values = c("#AFAFB2","#4B704D", "#C8972E")) + 
  scale_color_manual(values = c("#AFAFB2","#4B704D", "#C8972E")) + 
  scale_x_discrete(breaks=condition,labels=label) + 
  theme_classic() +
  #ylim(c(0,3)) + 
  ylim(c(0,3)) + 
  theme(legend.position = "none",
        aspect.ratio = 1.4,
        axis.text.x = element_text(angle = 45, size = 6, vjust = 0.4, color = "black", face = "bold"),
        axis.text.y = element_text(size = 6,color = "black", face = "bold"))+ 
  ggtitle("Filamentous actin") + 
  labs(fill='Cohort', x='', y='Phalloidin Intensity (A.U.)') 
p
ggsave("Tiam1_actin.pdf", plot= p, width = 6, height = 9, units = "cm")


####################
#Tiam1 Rac1 activation (FRET) with C.A. Src
####################
df_tiam1_fret = read.csv("TIAM1_Rac1_FRET.txt", sep = '\t')

colnames(df_tiam1_fret) <- c("Set", "FRET.Donor", "Type")
df_tiam1_fret.format = aggregate(data=df_tiam1_fret, FRET.Donor ~ Type, mean)
df_tiam1_fret.format$sd = aggregate(data=df_tiam1_fret, FRET.Donor ~ Type, sd)
df_tiam1_fret.format$length = aggregate(data=df_tiam1_fret, FRET.Donor ~ Type, length)
df_tiam1_fret.format$sem = df_tiam1_fret.format$sd$FRET.Donor / sqrt(df_tiam1_fret.format$length$FRET.Donor)
df_tiam1_fret.format
condition = c("Vector + Src CA", "WT Tiam1 +Src CA", "H1149P Tiam1 + Src CA")
label=c("Vector", "WT Tiam1", "H1149P Tiam1")
df_tiam1_fret.format$Type <- factor(df_tiam1_fret.format$Type, levels = condition)

#Statistical test
df_tiam1_fret_vector= subset(df_tiam1_fret, df_tiam1_fret$Type == 'Vector + Src CA')
df_tiam1_fret_wt= subset(df_tiam1_fret, df_tiam1_fret$Type == 'WT Tiam1 +Src CA')
df_tiam1_fret_mt= subset(df_tiam1_fret, df_tiam1_fret$Type == 'H1149P Tiam1 + Src CA')
shapiro.test(df_tiam1_fret_vector$FRET.Donor) # p-value = 0.5127
shapiro.test(df_tiam1_fret_wt$FRET.Donor) #p-value = 8.404e-08
shapiro.test(df_tiam1_fret_mt$FRET.Donor) #p-value = 0.6559

kruskal.test(data=df_tiam1_fret, FRET.Donor ~ Type)
#post hoc anlaysis
pairwise.wilcox.test(df_tiam1_fret$FRET.Donor, df_tiam1_fret$Type, p.adjust.method = "bonferroni", correct=FALSE)

upper_vector = median(df_tiam1_fret_vector$FRET.Donor) + 3 * mad(df_tiam1_fret_vector$FRET.Donor, constant = 1)
upper_wt = median(df_tiam1_fret_wt$FRET.Donor) + 3 * mad(df_tiam1_fret_wt$FRET.Donor, constant = 1)
upper_mt = median(df_tiam1_fret_mt$FRET.Donor) + 3 * mad(df_tiam1_fret_mt$FRET.Donor, constant = 1)

p <- ggplot(df_tiam1_fret.format, aes(x = Type, y = FRET.Donor, fill = Type)) +
  geom_bar(position=position_dodge(width=0.9) , stat="identity") +
  geom_jitter(data= subset(df_tiam1_fret, df_tiam1_fret$FRET.Donor <= upper_vector)[df_tiam1_fret$Type=='Vector + Src CA',], fill="black", size=0.5, alpha = 0.4) +
  geom_jitter(data= subset(df_tiam1_fret, df_tiam1_fret$FRET.Donor <= upper_wt)[df_tiam1_fret$Type=='WT Tiam1 +Src CA',], colour="black", size=0.5, alpha = 0.4, ) +
  geom_jitter(data= subset(df_tiam1_fret, df_tiam1_fret$FRET.Donor <= upper_mt)[df_tiam1_fret$Type=='H1149P Tiam1 + Src CA',], colour="black", size=0.5, alpha = 0.4)+
  geom_errorbar(aes(ymax=FRET.Donor+sem, ymin=FRET.Donor-sem), position=position_dodge(width=0.9) , width=0.6, size = 0.5)+
  scale_fill_manual(values = c("#AFAFB2","#4B704D", "#C8972E")) + 
  scale_color_manual(values = c("#AFAFB2","#4B704D", "#C8972E")) + 
  scale_x_discrete(breaks=condition,labels=label) + 
  theme_classic() +
  ylim(c(0,2.5)) + 
  theme(legend.position = "none",
        aspect.ratio = 1.4,
        axis.text.x = element_text(angle = 45, size = 6, vjust = 0.4, color = "black", face = "bold"),
        axis.text.y = element_text(size = 6,color = "black", face = "bold"))+ 
  ggtitle("Rac1 activation with Src C.A.") + 
  labs(fill='Cohort', x='', y='Normalized Rac1 Activity\n(FRET/Donor)') 
p

ggsave("Main_Tiam1_Fret.pdf", plot= p, width = 6, height = 9, units = "cm")




#Xenopus
####################

setwd("/Users/jinyha/Documents/NTDnew/Analysis/Validation/Raw_data/")
library(DescTools)

df_xeno = read.csv("Xenopus_main_dist.txt", sep = '\t')
df_xeno$condition = paste(df_xeno$treatment, " ",df_xeno$Dose,"ng", sep="")

label=c("Control", "Whamm")
dose_label=c("0","5","10")
#df_xeno = subset(df_xeno, df_xeno$Dose == "5"| df_xeno$Dose == "10")
df_xeno$treatment <- factor(df_xeno$treatment, levels = label)
df_xeno$Dose <- factor(df_xeno$Dose, levels = dose_label)

df_xeno.format = aggregate(data=df_xeno, Distance..mm. ~ condition, mean)
df_xeno.format$sd = aggregate(data=df_xeno, Distance..mm. ~ condition, sd)
df_xeno.format$length = aggregate(data=df_xeno, Distance..mm. ~ condition, length)
df_xeno.format$sem = df_xeno.format$sd$Distance..mm. / sqrt(df_xeno.format$length$Distance..mm.)

df_xeno.format %>% separate_wider_delim(df_xeno.format$condition, delim = " ", names = c("treatment", "Dose"))
df_xeno.format
df_xeno.format$treatment = str_split_i(df_xeno.format$condition, " ",1)
df_xeno.format

#Statistical test
df_xeno_vector= subset(df_xeno, df_xeno$condition == 'Control 0ng')
df_xeno_W5= subset(df_xeno, df_xeno$condition == 'Whamm 5ng')
df_xeno_W10= subset(df_xeno, df_xeno$condition == 'Whamm 10ng')
shapiro.test(df_xeno_vector$Distance..mm.) # p-value = 0.1432
shapiro.test(df_xeno_W5$Distance..mm.) #p-value = 0.2436
shapiro.test(df_xeno_W10$Distance..mm.) #p-value = 0.6095
#ANOVA
model <- aov(data = df_xeno, Distance..mm. ~ condition)
summary(model) # 2.37e-12 ***
m = aov(Distance..mm. ~ condition, data = df_xeno)
PostHocTest(m, method='newmankeuls')

condition = c("Control 0ng", "Whamm 5ng",  "Whamm 10ng")
df_xeno.format$condition <- factor(df_xeno.format$condition, levels = rev(condition))
df_xeno$condition <- factor(df_xeno$condition, levels = rev(condition))
p <- ggplot() + 
  geom_bar(position=position_dodge(width=0.9) , stat="identity") +
  geom_jitter(data=df_xeno,aes(x=condition, y = Distance..mm., color = treatment), size=1.5, shape=18) +
  geom_point(data=df_xeno.format,aes(x=condition, y = Distance..mm.), fill = "black") + 
  geom_errorbar(data = df_xeno.format, aes(x = condition, ymax=Distance..mm.+sem, ymin=Distance..mm.-sem), color = "black", position=position_dodge(width=0.9) , width=0.6, size = 0.5) + 
  scale_fill_manual(values = c("#AFAFB2","#F27C38", "#8C4820", "#724673")) + 
  scale_color_manual(values = c("#AFAFB2","#F27C38", "#8C4820", "#724673")) + 
  scale_x_discrete(breaks=condition) + 
  theme_classic() +
  #ylim(c(0,3)) + 
  theme(legend.position = "none",
        aspect.ratio = 0.4,
        axis.text.x = element_text(angle = 45, size = 6, vjust = 0.4, color = "black", face = "bold"),
        axis.text.y = element_text(size = 6,color = "black", face = "bold"))+ 
  labs(fill='Cohort', x='', y='Distance (mm)') + 
  coord_flip()
p
ggsave("Main_Xenopus_quanti.pdf", plot= p, width = 12, height = 9, units = "cm")


