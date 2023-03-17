
########################################################################
##                                                                    ##
##   Gut microbiota across the early stages of alpha-synucleinopathy  ##
##                                                                    ##
########################################################################

library(compositions)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ape)

########################### Microbiota Community Composition ##########################

library(car)
library(rstatix)

metadata <- read.csv('./metadata.csv')
which(colnames(metadata)=="rel_g001")
which(colnames(metadata)=="rel_g249")
otu <- metadata[ , c(65:313)]
dis <- vegdist(otu, method = 'bray')

### Homogeneity of group dispersion test
## Between four early stages

metadata$Group <- factor(metadata$Group, levels = c('Control', 'RBD_FDR', 'RBD','Early_PD'))
bd_group <- betadisper(dis, metadata$Group)
bd_group
boxplot(bd_group)
disper_group <- data.frame(bd_group$distances, metadata$Group)

# Tests of normality and homogeneity of variance

leveneTest(disper_group$bd_group.distances, disper_group$metadata.Group)
res_aov <- aov(bd_group.distances ~ metadata.Group,
               data = disper_group)
hist(res_aov$residuals)
shapiro.test(res_aov$residuals)

# Kruskal-Wallis rank sum test with post-hoc analysis

disper_group %>% kruskal_test(bd_group.distances ~ metadata.Group) 
disper_group_pwc <- disper_group %>% 
  dunn_test(bd_group.distances ~ metadata.Group, p.adjust.method = "BH") # adjusted p values with Benjamini-Hochberg method
disper_group_pwc

## Between controls, RBD-FDR, RBD_lr, RBD_hr, and early PD
# RBD_lr and RBD_hr correspond to lower and higher risk of developing Parkinson's disease 

metadata$Group_riskPD <- factor(metadata$Group_riskPD, levels = c('Control','RBD_FDR','RBD_lr', 'RBD_hr','Early_PD')) 
bd_group_riskPD <- betadisper(dis, metadata$Group_riskPD)
bd_group_riskPD
boxplot(bd_group_riskPD)
disper_group_riskPD <- data.frame(bd_group_riskPD$distances, metadata$Group_riskPD)

leveneTest(disper_group_riskPD$bd_group_riskPD.distances, disper_group_riskPD$metadata.Group_riskPD)
res_aov <- aov(bd_group_riskPD.distances ~ metadata.Group_riskPD,
               data = disper_group_riskPD)
hist(res_aov$residuals)
shapiro.test(res_aov$residuals)

disper_group_riskPD %>% kruskal_test(bd_group_riskPD.distances ~ metadata.Group_riskPD)
disper_group_riskPD_pwc <- disper_group_riskPD %>% 
  dunn_test(bd_group_riskPD.distances ~ metadata.Group_riskPD, p.adjust.method = "BH") 
disper_group_riskPD_pwc

## Between controls, RBD-FDR, RBD_less_5yr, RBD_5yr, and early PD
# RBD_less_5yr and RBD_5yr correspond to disease duration of RBD <= 5 and > 5 years

metadata$Group_durationRBD <- factor(metadata$Group_durationRBD, levels = c('Control','RBD_FDR','RBD_less_5yr', 'RBD_5yr',"Early_PD"))
bd_group_durationRBD <- betadisper(dis, metadata$Group_durationRBD)
bd_group_durationRBD
boxplot(bd_group_durationRBD)
disper_group_durationRBD <- data.frame(bd_group_durationRBD$distances, metadata$Group_durationRBD)

leveneTest(disper_group_durationRBD$bd_group_durationRBD.distances, disper_group_durationRBD$metadata.Group_durationRBD)
res_aov <- aov(bd_group_durationRBD.distances ~ metadata.Group_durationRBD,
               data = disper_group_durationRBD)
hist(res_aov$residuals)
shapiro.test(res_aov$residuals)

disper_group_durationRBD %>% kruskal_test(bd_group_durationRBD.distances ~ metadata.Group_durationRBD)
disper_group_durationRBD_pwc <- disper_group_durationRBD %>% 
  dunn_test(bd_group_durationRBD.distances ~ metadata.Group_durationRBD, p.adjust.method = "BH") 
disper_group_durationRBD_pwc


## Overall permanova test

adonis_group <- adonis2(dis~Group+age+sex, metadata, permutations = 99999)
adonis_group

## Pairwise comparisons in permanova
# Control vs RBD-FDR

metadata_ctrl_fdr <- metadata %>% filter(Group == "Control" | Group == "RBD_FDR")
which(colnames(metadata_ctrl_fdr)=="rel_g001")
which(colnames(metadata_ctrl_fdr)=="rel_g249")
otu_ctrl_fdr <- metadata_ctrl_fdr[ , c(65:313)]
row.names(otu_ctrl_fdr) <- metadata_ctrl_fdr$Sample_ID
dis_ctrl_fdr <- vegdist(otu_ctrl_fdr, method = 'bray')
adonis_group_ctrl_fdr <- adonis2(dis_ctrl_fdr~Group+age+sex, metadata_ctrl_fdr, permutations = 99999)
adonis_group_ctrl_fdr
adonis_group_ctrl_fdr <- adonis_group_ctrl_fdr[1, ]
row.names(adonis_group_ctrl_fdr)[1] <- "Control_RBDFDR"

# pcoa
pcoa <- wcmdscale(dis_ctrl_fdr, eig = T, x.ret = T, add = "cailliez")
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_group <- data.frame({pcoa$point})[1:2]
names(sample_group)[1:2] <- c('PCoA1', 'PCoA2')
sample_group$Sample_ID <-row.names(sample_group)
metadata_ctrl_fdr <- merge(sample_group, metadata_ctrl_fdr, by = 'Sample_ID')

pcoa_ctrl_fdr <- ggplot(metadata_ctrl_fdr, aes(x = PCoA1, y = PCoA2, colour = Group)) +
  scale_shape_manual(values = c(16, 16)) +
  scale_color_manual(values = c('#ffc300', "#6fa8d6")) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_point(aes(color = Group, shape = Group), size = 4, alpha = 0.9) +
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  coord_fixed()+
  theme(legend.position="none")+
  stat_ellipse(level = 0.70, size = 1.5)
pcoa_ctrl_fdr

# Control vs RBD

metadata_ctrl_rbd <- metadata %>% filter(Group == "Control" | Group == "RBD")
which(colnames(metadata_ctrl_rbd)=="rel_g001")
which(colnames(metadata_ctrl_rbd)=="rel_g249")
otu_ctrl_rbd <- metadata_ctrl_rbd[ , c(28:276)]
row.names(otu_ctrl_rbd) <- metadata_ctrl_rbd$Sample_ID
dis_ctrl_rbd <- vegdist(otu_ctrl_rbd, method = 'bray')
adonis_group_ctrl_rbd <- adonis2(dis_ctrl_rbd~Group+age+sex, metadata_ctrl_rbd, permutations = 99999)
adonis_group_ctrl_rbd
adonis_group_ctrl_rbd <- adonis_group_ctrl_rbd[1, ]
row.names(adonis_group_ctrl_rbd)[1] <- "Control_RBD"

# pcoa
pcoa <- wcmdscale(dis_ctrl_rbd, eig = T, x.ret = T, add = "cailliez")
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_group <- data.frame({pcoa$point})[1:2]
names(sample_group)[1:2] <- c('PCoA1', 'PCoA2')
sample_group$Sample_ID <-row.names(sample_group)
metadata_ctrl_rbd <- merge(sample_group, metadata_ctrl_rbd, by = 'Sample_ID')

pcoa_ctrl_rbd <- ggplot(metadata_ctrl_rbd, aes(x = PCoA1, y = PCoA2, colour = Group)) +
  scale_shape_manual(values = c(16, 16)) +
  scale_color_manual(values = c('#ffc300', '#1f7d94')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_point(aes(color = Group, shape = Group), size = 4, alpha = 0.9) +
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  coord_fixed()+
  theme(legend.position="none")+
  stat_ellipse(level = 0.70, size = 1.5)+
  xlim(-0.54, 0.7)
pcoa_ctrl_rbd

# Control vs Early PD

metadata_ctrl_pd <- metadata %>% filter(Group == "Control" | Group == "Early_PD")
which(colnames(metadata_ctrl_pd)=="rel_g001")
which(colnames(metadata_ctrl_pd)=="rel_g249")
otu_ctrl_pd <- metadata_ctrl_pd[ , c(65:313)]
row.names(otu_ctrl_pd) <- metadata_ctrl_pd$Sample_ID
dis_ctrl_pd <- vegdist(otu_ctrl_pd, method = 'bray')
adonis_group_ctrl_pd <- adonis2(dis_ctrl_pd~Group+age+sex, metadata_ctrl_pd, permutations = 99999)
adonis_group_ctrl_pd
adonis_group_ctrl_pd <- adonis_group_ctrl_pd[1, ]
row.names(adonis_group_ctrl_pd)[1] <- "Control_PD"

# pcoa
pcoa <- wcmdscale(dis_ctrl_pd, eig = T, x.ret = T, add = "cailliez")
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_group <- data.frame({pcoa$point})[1:2]
names(sample_group)[1:2] <- c('PCoA1', 'PCoA2')
sample_group$Sample_ID <-row.names(sample_group)
metadata_ctrl_pd <- merge(sample_group, metadata_ctrl_pd, by = 'Sample_ID')

pcoa_ctrl_pd <- ggplot(metadata_ctrl_pd, aes(x = PCoA1, y = PCoA2, colour = Group)) +
  scale_shape_manual(values = c(16, 16)) +
  scale_color_manual(values = c('#ffc300','#013540')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_point(aes(color = Group, shape = Group), size = 4, alpha = 0.9) +
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  coord_fixed()+
  theme(legend.position="none")+
  stat_ellipse(level = 0.70, size = 1.5)+
  xlim(-0.47, 0.65)
pcoa_ctrl_pd


# RBD-FDR vs RBD

metadata_fdr_rbd <- metadata %>% filter(Group == "RBD_FDR" | Group == "RBD")
which(colnames(metadata_fdr_rbd)=="rel_g001")
which(colnames(metadata_fdr_rbd)=="rel_g249")
otu_fdr_rbd <- metadata_fdr_rbd[ , c(65:313)]
row.names(otu_fdr_rbd) <- metadata_fdr_rbd$Sample_ID
dis_fdr_rbd <- vegdist(otu_fdr_rbd, method = 'bray')
adonis_group_fdr_rbd <- adonis2(dis_fdr_rbd~Group+age+sex, metadata_fdr_rbd, permutations = 99999)
adonis_group_fdr_rbd
adonis_group_fdr_rbd <- adonis_group_fdr_rbd[1, ]
row.names(adonis_group_fdr_rbd)[1] <- "RBDFDR_RBD" 

# pcoa
pcoa <- wcmdscale(dis_fdr_rbd, eig = T, x.ret = T, add = "cailliez")
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_group <- data.frame({pcoa$point})[1:2]
names(sample_group)[1:2] <- c('PCoA1', 'PCoA2')
sample_group$Sample_ID <-row.names(sample_group)
metadata_fdr_rbd <- merge(sample_group, metadata_fdr_rbd, by = 'Sample_ID')

pcoa_fdr_rbd <- ggplot(metadata_fdr_rbd, aes(x = PCoA1, y = PCoA2, colour = Group)) +
  scale_shape_manual(values = c(16, 16)) +
  scale_color_manual(values = c("#8dadc7",'#1f7d94')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_point(aes(color = Group, shape = Group), size = 4, alpha = 0.9) +
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  coord_fixed()+
  theme(legend.position="none")+
  stat_ellipse(level = 0.70, size = 1.5)+
  xlim(-0.52, 0.69)
pcoa_fdr_rbd


# RBD-FDR vs Early PD

metadata_fdr_pd <- metadata %>% filter(Group == "RBD_FDR" | Group == "Early_PD")
which(colnames(metadata_fdr_pd)=="rel_g001")
which(colnames(metadata_fdr_pd)=="rel_g249")
otu_fdr_pd <- metadata_fdr_pd[ , c(65:313)]
row.names(otu_fdr_pd) <- metadata_fdr_pd$Sample_ID
dis_fdr_pd <- vegdist(otu_fdr_pd, method = 'bray')
adonis_group_fdr_pd <- adonis2(dis_fdr_pd~Group+age+sex, metadata_fdr_pd, permutations = 99999)
adonis_group_fdr_pd
adonis_group_fdr_pd <- adonis_group_fdr_pd[1, ]
row.names(adonis_group_fdr_pd)[1] <- "RBDFDR_PD"

# pcoa
pcoa <- wcmdscale(dis_fdr_pd, eig = T, x.ret = T, add = "cailliez")
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_group <- data.frame({pcoa$point})[1:2]
names(sample_group)[1:2] <- c('PCoA1', 'PCoA2')
sample_group$Sample_ID <-row.names(sample_group)
metadata_fdr_pd <- merge(sample_group, metadata_fdr_pd, by = 'Sample_ID')

pcoa_fdr_pd <- ggplot(metadata_fdr_pd, aes(x = PCoA1, y = PCoA2, colour = Group)) +
  scale_shape_manual(values = c(16, 16)) +
  scale_color_manual(values = c("#8dadc7",'#013540')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_point(aes(color = Group, shape = Group), size = 4, alpha = 0.9) +
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  coord_fixed()+
  theme(legend.position="none")+
  stat_ellipse(level = 0.70, size = 1.5)+
  xlim(-0.65, 0.44)
pcoa_fdr_pd

# RBD vs Early PD

metadata_rbd_pd <- metadata %>% filter(Group == "RBD" | Group == "Early_PD")
which(colnames(metadata_rbd_pd)=="rel_g001")
which(colnames(metadata_rbd_pd)=="rel_g249")
otu_rbd_pd <- metadata_rbd_pd[ , c(65:313)]
row.names(otu_rbd_pd) <- metadata_rbd_pd$Sample_ID
dis_rbd_pd <- vegdist(otu_rbd_pd, method = 'bray')
adonis_group_rbd_pd <- adonis2(dis_rbd_pd~Group+age+sex, metadata_rbd_pd, permutations = 99999)
adonis_group_rbd_pd
adonis_group_rbd_pd <- adonis_group_rbd_pd[1, ]
row.names(adonis_group_rbd_pd)[1] <- "RBD_PD"

# pcoa
pcoa <- wcmdscale(dis_rbd_pd, eig = T, x.ret = T, add = "cailliez")
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_group <- data.frame({pcoa$point})[1:2]
names(sample_group)[1:2] <- c('PCoA1', 'PCoA2')
sample_group$Sample_ID <-row.names(sample_group)
metadata_rbd_pd <- merge(sample_group, metadata_rbd_pd, by = 'Sample_ID')

pcoa_rbd_pd <- ggplot(metadata_rbd_pd, aes(x = PCoA1, y = PCoA2, colour = Group)) +
  scale_shape_manual(values = c(16, 16)) +
  scale_color_manual(values = c('#1f7d94','#013540')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_point(aes(color = Group, shape = Group), size = 4, alpha = 0.9) +
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  coord_fixed()+
  theme(legend.position="none")+
  stat_ellipse(level = 0.70, size = 1.5)
pcoa_rbd_pd

# Adjust p values for multiple comparisons in the permanova tests

p_group_pwc <- rbind(adonis_group_ctrl_fdr,adonis_group_ctrl_rbd,adonis_group_ctrl_pd,
                     adonis_group_fdr_rbd,adonis_group_fdr_pd,adonis_group_rbd_pd)
p.adjust(p_group_pwc$`Pr(>F)`, method = "BH", n = length(p_group_pwc$`Pr(>F)`))

### Principal coordinates analysis
metadata <- read.csv('./metadata.csv')

## Between four early stages
which(colnames(metadata)=="rel_g001")
which(colnames(metadata)=="rel_g249")
otu <- metadata[ , c(65:313)]
row.names(otu) <- metadata$Sample_ID
dis <- vegdist(otu, method = 'bray')
pcoa <- wcmdscale(dis, eig = T, x.ret = T, add = "cailliez")
summary(pcoa)
ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')

# Scree plot
var.per <- as.data.frame(round(pcoa$eig/sum(pcoa$eig)*100, 1))
colnames(var.per)<- "variance_explained"
row.names(var.per) <-paste0("", 1:441)
var.per$PC <- row.names(var.per)
var.per <- var.per[c(1:9),] # show first nine principal components

screeplot <- var.per %>%
  ggplot(aes(x=PC,y=variance_explained, group=1))+
  geom_point(size=4, shape = 1)+
  geom_line()+
  labs(title="Scree plot")+
  xlab("Principal component") + 
  ylab("% of variance explained")+
  theme(panel.grid = (element_line(colour = "grey")), panel.background = element_rect(color = 'black', fill = 'transparent'))
screeplot

# Add the first two principal components to the metadata
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_group <- data.frame({pcoa$point})[1:2]
names(sample_group)[1:2] <- c('PCoA1', 'PCoA2')
sample_group$Sample_ID <-row.names(sample_group)
metadata <- merge(sample_group, metadata, by = 'Sample_ID')

## PCOA plot - Four early stages (Figure 2)
# Calculate the centroid of each group

metadata$Group <- factor(metadata$Group, levels = c('Control', 'RBD_FDR', 'RBD','Early_PD'))
centroids <- aggregate(cbind(PCoA1, PCoA2)~Group, metadata, mean)
f <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
std_err <- aggregate(cbind(se.x=PCoA1,se.y=PCoA2)~Group, metadata, f)
centroids <- merge(centroids,std_err, by="Group")   

segs <- merge(metadata, setNames(centroids, c('Group','oPCoA1','oPCoA2')),
              by = 'Group', sort = FALSE)

# Main pcoa plot
pcoa_group <- ggplot(metadata, aes(x = PCoA1, y = PCoA2, colour = Group)) +
  geom_segment(data = segs,
               mapping = aes(xend = oPCoA1, yend = oPCoA2), size = 0.7) +
  scale_shape_manual(values = c(16, 16, 16, 16)) +
  scale_color_manual(values = c('#ffc300', "#8dadc7",'#1f7d94','#013540')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_point(aes(color = Group, shape = Group), size = 4, alpha = 0.9) +
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  coord_fixed()+
  geom_label(data = centroids, aes(x = PCoA1, y =PCoA2, label=Group), label.padding = unit(0.25, "lines"), color = "black",
             fill="white", size = 2)+
  theme(legend.position="none")+
  xlim(-0.78, 0.54)+
  ylim(-0.83, 0.501)
pcoa_group

# Add the boxplots of first two principal components
my_comparisons = list(c("Control","RBD_FDR"), c("Control", "RBD"), c("Control", "Early_PD"), c("RBD_FDR", "RBD"), c("RBD_FDR", "Early_PD"), c("RBD", "Early_PD"))
boxplot_pcoa2 <- ggplot(metadata, aes(x = Group, y = PCoA2, colour =Group))+
  geom_boxplot()+
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")+
  scale_color_manual(values = c('#ffc300', "#6fa8d6",'#1f7d94','#013540')) +
  xlab("") + ylab("")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")
boxplot_pcoa2

metadata$Group <- factor(metadata$Group, levels = c('Early_PD', 'RBD','RBD_FDR', 'Control'))
boxplot_pcoa1 <- ggplot(metadata, aes(x = Group, y = PCoA1, colour =Group))+
  geom_boxplot()+
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        legend.position = "none")+
  scale_color_manual(values = c('#013540','#1f7d94',"#6fa8d6",'#ffc300')) +
  xlab("") + ylab("")+
  coord_flip()+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")
boxplot_pcoa1

# Arrange pcoa and boxplots on the same page.
pcoa_boxplot_group <- ggarrange(pcoa_group, boxplot_pcoa2, boxplot_pcoa1, NULL, widths = c(5,1.2), heights = c(5,1.3), 
                vjust = 1.0)

pcoa_boxplot_group

# Group centroids and standard errors of the centroids
pcoa_group_centroid <- ggplot(data = metadata, aes(PCoA1, PCoA2, color = factor (Group))) +
  geom_point(data=centroids, size=10,aes(color = Group, shape = Group)) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_color_manual(values = c('#ffc300', "#6fa8d6",'#37abc8','#005f73')) +
  geom_errorbar(data=centroids, mapping = aes(ymin=PCoA2-se.y,ymax=PCoA2+se.y),width=0.008, size = 1)+
  geom_errorbarh(data=centroids,aes(xmin=PCoA1-se.x,xmax=PCoA1+se.x),height=0.008, size = 1)+
  labs(x = paste('PC1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PC2: ', round(100 * pcoa_eig[2], 2), '%')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.position = "none")

pcoa_group_centroid


## PCOA plot - five groups (Supplementary Figure 5A-B)
# Between controls, RBD-FDR, RBD_lr, RBD_hr, and early PD

centroids <- aggregate(cbind(PCoA1, PCoA2)~Group_riskPD, metadata, mean)
f <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
std_err <- aggregate(cbind(se.x=PCoA1,se.y=PCoA2)~Group_riskPD, metadata, f)
centroids <- merge(centroids,std_err, by="Group_riskPD")   
centroids$Group_riskPD <- factor(centroids$Group_riskPD, levels = c("Control","RBD_FDR",'RBD_lr', 'RBD_hr',"Early_PD"))

metadata$Group_riskPD <- factor(metadata$Group_riskPD, levels = c("Control","RBD_FDR",'RBD_lr', 'RBD_hr',"Early_PD"))
pcoa_group_riskPD <- ggplot(data = metadata, aes(PCoA1, PCoA2, color = factor (Group_riskPD))) +
  geom_point(data=centroids, size=8,aes(color = Group_riskPD, shape = Group_riskPD)) +
  scale_shape_manual(values = c(16, 17, 15, 15, 18)) +
  scale_color_manual(values = c('#ffc300', "#6aa0cc",'#858585','#005f73','#14213d')) +
  geom_errorbar(data=centroids, mapping = aes(ymin=PCoA2-se.y,ymax=PCoA2+se.y),width=0.01, size = 1)+
  geom_errorbarh(data=centroids,aes(xmin=PCoA1-se.x,xmax=PCoA1+se.x),height=0.01, size = 1)+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.position = "none")

pcoa_group_riskPD

# Between controls, RBD-FDR, RBD_less_5yr, RBD_5yr, and early PD

centroids <- aggregate(cbind(PCoA1, PCoA2)~Group_durationRBD, metadata, mean)
f <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
std_err <- aggregate(cbind(se.x=PCoA1,se.y=PCoA2)~Group_durationRBD, metadata, f)
centroids <- merge(centroids,std_err, by="Group_durationRBD")  
centroids$Group_durationRBD <- factor(centroids$Group_durationRBD, levels = c('Control','RBD_FDR','RBD_less_5yr', 'RBD_5yr',"Early_PD"))          

metadata$Group_durationRBD <- factor(metadata$Group_durationRBD, levels = c('Control','RBD_FDR','RBD_less_5yr', 'RBD_5yr',"Early_PD"))          
pcoa_group_durationRBD <- ggplot(data = metadata, aes(PCoA1, PCoA2, color = factor (Group_durationRBD))) +
  geom_point(data=centroids, size=8,aes(color = Group_durationRBD, shape = Group_durationRBD)) +
  scale_shape_manual(values = c(16, 17, 15, 15, 18)) +
  scale_color_manual(values = c('#ffc300', "#6aa0cc",'#858585','#005f73','#14213d')) +
  geom_errorbar(data=centroids, mapping = aes(ymin=PCoA2-se.y,ymax=PCoA2+se.y),width=0.01, size = 1)+
  geom_errorbarh(data=centroids,aes(xmin=PCoA1-se.x,xmax=PCoA1+se.x),height=0.01, size = 1)+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.position = "none")

pcoa_group_durationRBD         


# Correlation between genus faecalibacterium and disease progression (Figure 3D)

metadata <- read.csv('./metadata.csv')
          
which(colnames(metadata)=="rel_g001")
which(colnames(metadata)=="rel_g249")
which(colnames(metadata)=="Group_n") # control=1, RBD-FDR=2, RBD=3, Early PD=4
trait <- metadata[, 4]
otu <- metadata[ ,c(65:313)]
row.names(otu) <- metadata$Sample_ID
dis <- vegdist(otu, method = 'bray')
pcoa <- wcmdscale(dis, eig = T, x.ret = T, add = "cailliez")
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_group <- data.frame({pcoa$point})[1:2]
names(sample_group)[1:2] <- c('PCoA1', 'PCoA2')
sample_group$Sample_ID <-row.names(sample_group)
metadata <- merge(sample_group, metadata, by = 'Sample_ID')

trait <- scale(trait)
n <- nrow(trait)
points.stand <- scale(pcoa$points)
S <- cov(trait_df, points.stand)
colnames(S) <- colnames(pcoa$points)
pcoa$S <- S # Add values to pcoa
pcoa
pcoa$S[, 1:3]

# Pcoa plot with colors indicate the clr transformed abundance of faecalibacterium          
pcoa_faecalibacterium = ggplot(metadata, aes(PCoA1, PCoA2, color = factor (Faecalibacterium_10gp))) +
            geom_point(data=metadata, size=4, aes(color = Faecalibacterium_10gp, shape = Faecalibacterium_10gp)) +
            scale_shape_manual(values = c(16, 16, 16, 16, 16,16, 16, 16, 16, 16)) +
            scale_color_manual(values = c("#023459","#115f9a", "#1984c5", "#22a7f0", "#48b5c4", "#76c68f", "#a6d75b", "#c9e52f", "#d0ee11", "#d0f400")) +
            labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%')) +
            theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
            theme(legend.title = element_blank(), legend.background = element_rect(fill = "white"))
          
pcoa_faecalibacterium
          
data = as.data.frame(pcoa$S*3.5)
          names(data)[1:2] <- c('Axis1', 'Axis2')

pcoa_faecalibacterium + geom_segment(x = 0, y = 0, 
                         mapping = aes(xend = -0.7157209, yend =-0.7011793), # Add arrow head
                         colour = "#000000", arrow = arrow(length = unit(3, "mm")))
          


########################### Differential Abundance Analysis ###########################

## Microbiome multivariable associations with linear model (MaAsLin 2)
## The input data is centered log-ratio (clr) transformed abundance of differential 
## genera (n = 26) and family (n = 16) that identified by using Kruskal-Wallis test 
## between the four groups in SPSS (Benjamini-Hochberg method adjusted p-values < 0.05) 

library(Maaslin2)
input_metadata <- read.csv('./metadata.csv')

# Genus
which(colnames(input_metadata)=="Collinsella")
which(colnames(input_metadata)=="Akkermansia")
input_data <- input_metadata[ , c(314:339)]

fit_data_genus <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_group_genus_unadjusted', transform = "NONE",
  fixed_effects = c('Group'),
  random_effects = c('family_id'),
  reference = "Group,Control",
  normalization = 'NONE',
  correction = "BH",
  standardize = FALSE) 

fit_data_genus <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_group_genus_adjusted', transform = "NONE",
  fixed_effects = c('Group','sex','age'),
  random_effects = c('family_id'),
  reference = "Group,Control",
  normalization = 'NONE',
  correction = "BH",
  standardize = FALSE)

View(fit_data_genus$results)
res_genus <- as.data.frame(fit_data_genus$results)  

# Family
which(colnames(input_metadata)=="Coriobacteriaceae")
which(colnames(input_metadata)=="Sutterellaceae")
input_data <- input_metadata[, c(402:417)]

fit_data_family <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_group_family_unadjusted', transform = "NONE",
  fixed_effects = c('Group'),
  random_effects = c('family_id'),
  reference = "Group,Control",
  normalization = 'NONE',
  correction = "BH",
  standardize = FALSE)

fit_data_family <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_group_family_adjusted', transform = "NONE",
  fixed_effects = c('Group','sex','age'),
  random_effects = c('family_id'),
  reference = "Group,Control",
  normalization = 'NONE',
  correction = "BH",
  standardize = FALSE)

View(fit_data_family$results)
res_family <- as.data.frame(fit_data_family$results)  


########################### Analysis of Host-Microbiome Interactions ###########################

## The effect of host factors on the overall microbial compositions
# Inter-individual variation explained by each host factor

metadata <- read.csv('./metadata.csv')
which(colnames(metadata)=="rel_g001")
which(colnames(metadata)=="rel_g249")
otu <- metadata[ , c(65:313)]
dis <- vegdist(otu, method = 'bray')

adonis_age <- adonis2(dis~age, metadata, permutations = 99999)
adonis_age <-adonis_age[1,]
row.names(adonis_age)[1] <- "age"  
adonis_sex <- adonis2(dis~sex, metadata, permutations = 99999)
adonis_sex <-adonis_sex[1,]
row.names(adonis_sex)[1] <- "sex"  
adonis_bmf_score <- adonis2(dis~BMF_score, metadata, permutations = 99999)
adonis_bmf_score <-adonis_bmf_score[1,]
row.names(adonis_bmf_score)[1] <- "bmf_score" 
adonis_statin <- adonis2(dis~Statins_2gp, metadata, permutations = 99999)
adonis_statin <-adonis_statin[1,]
row.names(adonis_statin)[1] <- "statin" 
adonis_ppi <- adonis2(dis~Proton_pump_inhibitors_2gp, metadata, permutations = 99999)
adonis_ppi <-adonis_ppi[1,]
row.names(adonis_ppi)[1] <- "ppi" 
adonis_laxative <- adonis2(dis~Osmotic_laxative_2gp, metadata, permutations = 99999)
adonis_laxative <-adonis_laxative[1,]
row.names(adonis_laxative)[1] <- "laxative" 
adonis_antidepressant <- adonis2(dis~Antidepressants_2gp, metadata, permutations = 99999)
adonis_antidepressant <-adonis_antidepressant[1,]
row.names(adonis_antidepressant)[1] <- "antidepressant" 
adonis_benzodiazepine <- adonis2(dis~Benzodiazepines_2gp, metadata, permutations = 99999)
adonis_benzodiazepine <-adonis_benzodiazepine[1,]
row.names(adonis_benzodiazepine)[1] <- "benzodiazepines" 
adonis_carbidopa_levodopa <- adonis2(dis~Carbidopa_Levodopa_2gp, metadata, permutations = 99999)
adonis_carbidopa_levodopa <-adonis_carbidopa_levodopa[1,]
row.names(adonis_carbidopa_levodopa)[1] <- "carbidopa_levodopa" 
adonis_maob_inhibitor <- adonis2(dis~MAOB_inhibitor_2gp, metadata, permutations = 99999)
adonis_maob_inhibitor <-adonis_maob_inhibitor[1,]
row.names(adonis_maob_inhibitor)[1] <- "maob_inhibitor" 
adonis_dopamine_agonist <- adonis2(dis~Dopamine_agonist_2gp, metadata, permutations = 99999)
adonis_dopamine_agonist <-adonis_dopamine_agonist[1,]
row.names(adonis_dopamine_agonist)[1] <- "dopamine_agonist" 

p_sum <- rbind(adonis_age,adonis_sex,adonis_bmf_score,adonis_statin,adonis_ppi,adonis_laxative,
               adonis_antidepressant,adonis_benzodiazepine, adonis_carbidopa_levodopa,adonis_maob_inhibitor,
               adonis_dopamine_agonist)


## The effect of host factors in permanova pairwise comparisons 
# Control vs RBD-FDR
adonis_group_ctrl_fdr_cov <- adonis2(dis_ctrl_fdr~Group+age+sex+BMF_score+Statins_2gp+
                                   Proton_pump_inhibitors_2gp+Osmotic_laxative_2gp+
                                   Antidepressants_2gp+Benzodiazepines_2gp, metadata_ctrl_fdr, permutations = 99999)
adonis_group_ctrl_fdr_cov
adonis_group_ctrl_fdr_cov_padj <-p.adjust(adonis_group_ctrl_fdr_cov$`Pr(>F)`, method = "BH", n = length(adonis_group_ctrl_fdr_cov$`Pr(>F)`)-2)

# Control vs RBD
adonis_group_ctrl_rbd_cov <- adonis2(dis_ctrl_rbd~Group+age+sex+BMF_score+Statins_2gp+
                                       Proton_pump_inhibitors_2gp+Osmotic_laxative_2gp+
                                       Antidepressants_2gp+Benzodiazepines_2gp, metadata_ctrl_rbd, permutations = 99999)
adonis_group_ctrl_rbd_cov
adonis_group_ctrl_rbd_cov_padj <-p.adjust(adonis_group_ctrl_rbd_cov$`Pr(>F)`, method = "BH", n = length(adonis_group_ctrl_rbd_cov$`Pr(>F)`)-2)

# Control vs Early PD
adonis_group_ctrl_pd_cov <- adonis2(dis_ctrl_pd~Group+age+sex+BMF_score+Statins_2gp+
                                       Proton_pump_inhibitors_2gp+Osmotic_laxative_2gp+
                                       Antidepressants_2gp+Benzodiazepines_2gp+Carbidopa_Levodopa_2gp+
                                      MAOB_inhibitor_2gp+Dopamine_agonist_2gp, metadata_ctrl_pd, permutations = 99999)
adonis_group_ctrl_pd_cov
adonis_group_ctrl_pd_cov_padj <-p.adjust(adonis_group_ctrl_pd_cov$`Pr(>F)`, method = "BH", n = length(adonis_group_ctrl_pd_cov$`Pr(>F)`)-2)

# RBD-FDR vs RBD
adonis_group_fdr_rbd_cov <- adonis2(dis_fdr_rbd~Group+age+sex+BMF_score+Statins_2gp+
                                       Proton_pump_inhibitors_2gp+Osmotic_laxative_2gp+
                                       Antidepressants_2gp+Benzodiazepines_2gp, metadata_fdr_rbd, permutations = 99999)
adonis_group_fdr_rbd_cov
adonis_group_fdr_rbd_cov_padj <-p.adjust(adonis_group_fdr_rbd_cov$`Pr(>F)`, method = "BH", n = length(adonis_group_fdr_rbd_cov$`Pr(>F)`)-2)

# RBD-FDR vs Early PD
adonis_group_fdr_pd_cov <- adonis2(dis_fdr_pd~Group+age+sex+BMF_score+Statins_2gp+
                                      Proton_pump_inhibitors_2gp+Osmotic_laxative_2gp+
                                      Antidepressants_2gp+Benzodiazepines_2gp+Carbidopa_Levodopa_2gp+
                                      MAOB_inhibitor_2gp+Dopamine_agonist_2gp, metadata_fdr_pd, permutations = 99999)
adonis_group_fdr_pd_cov
adonis_group_fdr_pd_cov_padj <-p.adjust(adonis_group_fdr_pd_cov$`Pr(>F)`, method = "BH", n = length(adonis_group_fdr_pd_cov$`Pr(>F)`)-2)

# RBD vs Early PD
adonis_group_rbd_pd_cov <- adonis2(dis_rbd_pd~Group+age+sex+BMF_score+Statins_2gp+
                                     Proton_pump_inhibitors_2gp+Osmotic_laxative_2gp+
                                     Antidepressants_2gp+Benzodiazepines_2gp+Carbidopa_Levodopa_2gp+
                                     MAOB_inhibitor_2gp+Dopamine_agonist_2gp, metadata_rbd_pd, permutations = 99999)
adonis_group_rbd_pd_cov
adonis_group_rbd_pd_cov_padj <-p.adjust(adonis_group_rbd_pd_cov$`Pr(>F)`, method = "BH", n = length(adonis_group_rbd_pd_cov$`Pr(>F)`)-2)


## The effect of host factors on taxa abundance
# Genus abundance
input_metadata <- read.csv('./metadata.csv')
which(colnames(input_metadata)=="Collinsella")
which(colnames(input_metadata)=="Akkermansia")
input_data <- input_metadata[ , c(314:339)]

fit_data_genus_cov <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_group_genus_cov', transform = "NONE",
  fixed_effects = c('Group', "Statins_2gp",'sex','age','CAQ_constipation_frequency_ADJUSTED',"Proton_pump_inhibitors_2gp","Osmotic_laxative_2gp", 
                    "Antidepressants_2gp", "Benzodiazepines_2gp"), 
  random_effects = c('family_id'),
  reference = "Group,Control",
  normalization = 'NONE',
  correction = "BH",
  standardize = FALSE)

View(fit_data_genus_cov$results)
res_genus_cov <- as.data.frame(fit_data_genus_cov$results)  

# Family abundance
which(colnames(input_metadata)=="Coriobacteriaceae")
which(colnames(input_metadata)=="Sutterellaceae")
input_data <- input_metadata[, c(402:417)]

fit_data_family_cov <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_group_family_cov', transform = "NONE",
  fixed_effects = c('Group', "Statins_2gp",'sex','age','CAQ_constipation_frequency_ADJUSTED',"Proton_pump_inhibitors_2gp","Osmotic_laxative_2gp", 
                    "Antidepressants_2gp", "Benzodiazepines_2gp"), 
  random_effects = c('family_id'),
  reference = "Group, Control",
  normalization = 'NONE',
  correction = "BH",
  standardize = FALSE)

View(fit_data_family_cov$results)
res_family_cov <- as.data.frame(fit_data_family_cov$results)  

# Genus abundance (PD specific drugs)
input_metadata <- input_metadata %>% filter(Group == "Early_PD")
which(colnames(input_metadata)=="Collinsella")
which(colnames(input_metadata)=="Akkermansia")
input_data <- input_metadata[ , c(314:339)]

fit_data_genus_pd_drugs <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_genus_pd_drugs', transform = "NONE",
  fixed_effects = c("Carbidopa_Levodopa_2gp",'MAOB_inhibitor_2gp','Dopamine_agonist_2gp'),
  normalization = 'NONE',
  correction = "BH",
  standardize = FALSE)

View(fit_data_genus_pd_drugs$results)
res_genus_pd_drugs <- as.data.frame(fit_data_genus_pd_drugs$results)  

# Family abundance (PD specific drugs)
which(colnames(input_metadata)=="Coriobacteriaceae")
which(colnames(input_metadata)=="Sutterellaceae")
input_data <- input_metadata[ , c(402:417)]

fit_data_family_pd_drugs <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_family_pd_drugs', transform = "NONE",
  fixed_effects = c("Carbidopa_Levodopa_2gp",'MAOB_inhibitor_2gp','Dopamine_agonist_2gp'), 
  normalization = 'NONE',
  correction = "BH",
  standardize = FALSE)

View(fit_data_family_pd_drugs$results)
res_family_pd_drugs <- as.data.frame(fit_data_family_pd_drugs$results)  




########################### Random Forest classification ###########################

library(randomForest)
library(pROC)
library(caret)
library(vegan)
library(ggplot2)
library(plyr)

# The below code only demonstrate random forest model for classification of control and RBD 
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)

# Repeat 01 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
which(colnames(metadata)=="Group")
which(colnames(metadata)=="Collinsella")
which(colnames(metadata)=="Klebsiella")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_01 <- metadata[idx, ]
otu_test_0.8_01 <- metadata[-idx, ]
otu_train_m_0.8_01 <- otu_train_0.8_01[, c(3,314:401)]
otu_test_m_0.8_01 <- otu_test_0.8_01[, c(3,314:401)]
y <- otu_train_m_0.8_01$Group
x <- otu_train_m_0.8_01[, 2:89]
set.seed(123)
lmProfile_0.8_01<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_01
importance_otu_0.8_01 <- as.data.frame(lmProfile_0.8_01$fit$importance)

resample <- c("Resample 01")
otu_select <- rownames(importance_otu_0.8_01)[1:20]
otu_train_m_select_0.8_01 <- otu_train_m_0.8_01[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_01 <- otu_test_m_0.8_01[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_01$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_01$Group
predictions_final_0.8_01 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_01 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_01$fit, otu_test_m_select_0.8_01, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_01$Group
predictions_0.8_01 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_01 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 02 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_02 <- metadata[idx, ]
otu_test_0.8_02 <- metadata[-idx, ]
otu_train_m_0.8_02 <- otu_train_0.8_02[, c(3,314:401)]
otu_test_m_0.8_02 <- otu_test_0.8_02[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_02$Group
x <- otu_train_m_0.8_02[, 2:89]

lmProfile_0.8_02<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_02
importance_otu_0.8_02 <- as.data.frame(lmProfile_0.8_02$fit$importance)

resample <- c("Resample 02")
otu_select <- rownames(importance_otu_0.8_02)[1:19]
otu_train_m_select_0.8_02 <- otu_train_m_0.8_02[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_02 <- otu_test_m_0.8_02[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_02$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_02$Group
predictions_final_0.8_02 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_02 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_02$fit, otu_test_m_select_0.8_02, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_02$Group
predictions_0.8_02 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_02 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

                              
# Repeat 03 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_03 <- metadata[idx, ]
otu_test_0.8_03 <- metadata[-idx, ]
otu_train_m_0.8_03 <- otu_train_0.8_03[, c(3,314:401)]
otu_test_m_0.8_03 <- otu_test_0.8_03[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_03$Group
x <- otu_train_m_0.8_03[, 2:89]

lmProfile_0.8_03<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_03
importance_otu_0.8_03 <- as.data.frame(lmProfile_0.8_03$fit$importance)

resample <- c("Resample 03")
otu_select <- rownames(importance_otu_0.8_03)[1:8]
otu_train_m_select_0.8_03 <- otu_train_m_0.8_03[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_03 <- otu_test_m_0.8_03[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_03$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_03$Group
predictions_final_0.8_03 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_03 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_03$fit, otu_test_m_select_0.8_03, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_03$Group
predictions_0.8_03 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_03 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

                                     
# Repeat 04 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_04 <- metadata[idx, ]
otu_test_0.8_04 <- metadata[-idx, ]
otu_train_m_0.8_04 <- otu_train_0.8_04[, c(3,314:401)]
otu_test_m_0.8_04 <- otu_test_0.8_04[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_04$Group
x <- otu_train_m_0.8_04[, 2:89]

lmProfile_0.8_04<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_04
importance_otu_0.8_04 <- as.data.frame(lmProfile_0.8_04$fit$importance)

resample <- c("Resample 04")
otu_select <- rownames(importance_otu_0.8_04)[1:16] # based on the number of selected predictors
otu_train_m_select_0.8_04 <- otu_train_m_0.8_04[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_04 <- otu_test_m_0.8_04[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_04$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_04$Group
predictions_final_0.8_04 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_04 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_04$fit, otu_test_m_select_0.8_04, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_04$Group
predictions_0.8_04 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_04 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

                                     
# Repeat 05 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_05 <- metadata[idx, ]
otu_test_0.8_05 <- metadata[-idx, ]
otu_train_m_0.8_05 <- otu_train_0.8_05[, c(3,314:401)]
otu_test_m_0.8_05 <- otu_test_0.8_05[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_05$Group
x <- otu_train_m_0.8_05[, 2:89]

lmProfile_0.8_05<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_05
importance_otu_0.8_05 <- as.data.frame(lmProfile_0.8_05$fit$importance)
                                     
resample <- c("Resample 05")
otu_select <- rownames(importance_otu_0.8_05)[1:32]
otu_train_m_select_0.8_05 <- otu_train_m_0.8_05[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_05 <- otu_test_m_0.8_05[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_05$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_05$Group
predictions_final_0.8_05 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_05 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_05$fit, otu_test_m_select_0.8_05, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_05$Group
predictions_0.8_05 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_05 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

                                     
# Repeat 06 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_06 <- metadata[idx, ]
otu_test_0.8_06 <- metadata[-idx, ]
otu_train_m_0.8_06 <- otu_train_0.8_06[, c(3,314:401)]
otu_test_m_0.8_06 <- otu_test_0.8_06[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_06$Group
x <- otu_train_m_0.8_06[, 2:89]

lmProfile_0.8_06<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_06
importance_otu_0.8_06 <- as.data.frame(lmProfile_0.8_06$fit$importance)
                                     
resample <- c("Resample 06")
otu_select <- rownames(importance_otu_0.8_06)[1:14]
otu_train_m_select_0.8_06 <- otu_train_m_0.8_06[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_06 <- otu_test_m_0.8_06[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_06$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_06$Group
predictions_final_0.8_06 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_06 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_06$fit, otu_test_m_select_0.8_06, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_06$Group
predictions_0.8_06 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_06 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 07 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_07 <- metadata[idx, ]
otu_test_0.8_07 <- metadata[-idx, ]
otu_train_m_0.8_07 <- otu_train_0.8_07[, c(3,314:401)]
otu_test_m_0.8_07 <- otu_test_0.8_07[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_07$Group
x <- otu_train_m_0.8_07[, 2:89]

lmProfile_0.8_07<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_07
importance_otu_0.8_07 <- as.data.frame(lmProfile_0.8_07$fit$importance)
                                     
resample <- c("Resample 07")
otu_select <- rownames(importance_otu_0.8_07)[1:10]
otu_train_m_select_0.8_07 <- otu_train_m_0.8_07[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_07 <- otu_test_m_0.8_07[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_07$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_07$Group
predictions_final_0.8_07 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_07 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_07$fit, otu_test_m_select_0.8_07, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_07$Group
predictions_0.8_07 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_07 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 08 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_08 <- metadata[idx, ]
otu_test_0.8_08 <- metadata[-idx, ]
otu_train_m_0.8_08 <- otu_train_0.8_08[, c(3,314:401)]
otu_test_m_0.8_08 <- otu_test_0.8_08[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_08$Group
x <- otu_train_m_0.8_08[, 2:89]

lmProfile_0.8_08<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_08
importance_otu_0.8_08 <- as.data.frame(lmProfile_0.8_08$fit$importance)

resample <- c("Resample 08")
otu_select <- rownames(importance_otu_0.8_08)[1:13]
otu_train_m_select_0.8_08 <- otu_train_m_0.8_08[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_08 <- otu_test_m_0.8_08[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_08$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_08$Group
predictions_final_0.8_08 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_08 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_08$fit, otu_test_m_select_0.8_08, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_08$Group
predictions_0.8_08 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_08 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
                                     
                                     
# Repeat 09 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_09 <- metadata[idx, ]
otu_test_0.8_09 <- metadata[-idx, ]
otu_train_m_0.8_09 <- otu_train_0.8_09[, c(3,314:401)]
otu_test_m_0.8_09 <- otu_test_0.8_09[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_09$Group
x <- otu_train_m_0.8_09[, 2:89]

lmProfile_0.8_09<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_09
importance_otu_0.8_09 <- as.data.frame(lmProfile_0.8_09$fit$importance)
resample <- c("Resample 09")
otu_select <- rownames(importance_otu_0.8_09)[1:31]
otu_train_m_select_0.8_09 <- otu_train_m_0.8_09[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_09 <- otu_test_m_0.8_09[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_09$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_09$Group
predictions_final_0.8_09 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_09 <-cbind(sen,speci,resample) 
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_09$fit, otu_test_m_select_0.8_09, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_09$Group
predictions_0.8_09 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_09 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 10 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_10 <- metadata[idx, ]
otu_test_0.8_10 <- metadata[-idx, ]
otu_train_m_0.8_10 <- otu_train_0.8_10[, c(3,314:401)]
otu_test_m_0.8_10 <- otu_test_0.8_10[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_10$Group
x <- otu_train_m_0.8_10[, 2:89]

lmProfile_0.8_10<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_10
importance_otu_0.8_10 <- as.data.frame(lmProfile_0.8_10$fit$importance)
resample <- c("Resample 10")
otu_select <- rownames(importance_otu_0.8_10)[1:15]
otu_train_m_select_0.8_10 <- otu_train_m_0.8_10[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_10 <- otu_test_m_0.8_10[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_10$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_10$Group
predictions_final_0.8_10 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_10 <-cbind(sen,speci,resample)
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_10$fit, otu_test_m_select_0.8_10, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_10$Group
predictions_0.8_10 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_10 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 11 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_11 <- metadata[idx, ]
otu_test_0.8_11 <- metadata[-idx, ]
otu_train_m_0.8_11 <- otu_train_0.8_11[, c(3,277:364)]
otu_test_m_0.8_11 <- otu_test_0.8_11[, c(3,277:364)]
set.seed(123)
y <- otu_train_m_0.8_11$Group
x <- otu_train_m_0.8_11[, 2:89]

lmProfile_0.8_11<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_11
importance_otu_0.8_11 <- as.data.frame(lmProfile_0.8_11$fit$importance)
                                     
resample <- c("Resample 11")
otu_select <- rownames(importance_otu_0.8_11)[1:12]
otu_train_m_select_0.8_11 <- otu_train_m_0.8_11[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_11 <- otu_test_m_0.8_11[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_11$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_11$Group
predictions_final_0.8_11 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_11 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_11$fit, otu_test_m_select_0.8_11, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_11$Group
predictions_0.8_11 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_11 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

                                     
# Repeat 12 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_12 <- metadata[idx, ]
otu_test_0.8_12 <- metadata[-idx, ]
otu_train_m_0.8_12 <- otu_train_0.8_12[, c(3,314:401)]
otu_test_m_0.8_12 <- otu_test_0.8_12[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_12$Group
x <- otu_train_m_0.8_12[, 2:89]

lmProfile_0.8_12<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_12
importance_otu_0.8_12 <- as.data.frame(lmProfile_0.8_12$fit$importance)
                                     
resample <- c("Resample 12")
otu_select <- rownames(importance_otu_0.8_12)[1:14]
otu_train_m_select_0.8_12 <- otu_train_m_0.8_12[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_12 <- otu_test_m_0.8_12[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_12$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_12$Group
predictions_final_0.8_12 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_12 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_12$fit, otu_test_m_select_0.8_12, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_12$Group
predictions_0.8_12 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_12 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 13 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_13 <- metadata[idx, ]
otu_test_0.8_13 <- metadata[-idx, ]
otu_train_m_0.8_13 <- otu_train_0.8_13[, c(3,314:401)]
otu_test_m_0.8_13 <- otu_test_0.8_13[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_13$Group
x <- otu_train_m_0.8_13[, 2:89]

lmProfile_0.8_13<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_13
importance_otu_0.8_13 <- as.data.frame(lmProfile_0.8_13$fit$importance)
resample <- c("Resample 13")
otu_select <- rownames(importance_otu_0.8_13)[1:17]
otu_train_m_select_0.8_13 <- otu_train_m_0.8_13[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_13 <- otu_test_m_0.8_13[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_13$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_13$Group
predictions_final_0.8_13 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_13<-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_13$fit, otu_test_m_select_0.8_13, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_13$Group
predictions_0.8_13 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_13 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 14 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_14 <- metadata[idx, ]
otu_test_0.8_14 <- metadata[-idx, ]
otu_train_m_0.8_14 <- otu_train_0.8_14[, c(3,314:401)]
otu_test_m_0.8_14 <- otu_test_0.8_14[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_14$Group
x <- otu_train_m_0.8_14[, 2:89]

lmProfile_0.8_14<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_14
importance_otu_0.8_14 <- as.data.frame(lmProfile_0.8_14$fit$importance)
resample <- c("Resample 14")
otu_select <- rownames(importance_otu_0.8_14)[1:12]
otu_train_m_select_0.8_14 <- otu_train_m_0.8_14[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_14 <- otu_test_m_0.8_14[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_14$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_14$Group
predictions_final_0.8_14 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_14 <-cbind(sen,speci,resample)
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_14$fit, otu_test_m_select_0.8_14, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_14$Group
predictions_0.8_14 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_14 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 15 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_15 <- metadata[idx, ]
otu_test_0.8_15 <- metadata[-idx, ]
otu_train_m_0.8_15 <- otu_train_0.8_15[, c(3,314:401)]
otu_test_m_0.8_15 <- otu_test_0.8_15[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_15$Group
x <- otu_train_m_0.8_15[, 2:89]

lmProfile_0.8_15<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_15
importance_otu_0.8_15 <- as.data.frame(lmProfile_0.8_15$fit$importance)
resample <- c("Resample 15")
otu_select <- rownames(importance_otu_0.8_15)[1:17]
otu_train_m_select_0.8_15 <- otu_train_m_0.8_15[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_15 <- otu_test_m_0.8_15[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_15$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_15$Group
predictions_final_0.8_15 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_15 <-cbind(sen,speci,resample) 
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_15$fit, otu_test_m_select_0.8_15, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_15$Group
predictions_0.8_15 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_15 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 16 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_16 <- metadata[idx, ]
otu_test_0.8_16 <- metadata[-idx, ]
otu_train_m_0.8_16 <- otu_train_0.8_16[, c(3,314:401)]
otu_test_m_0.8_16 <- otu_test_0.8_16[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_16$Group
x <- otu_train_m_0.8_16[, 2:89]

lmProfile_0.8_16<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_16
importance_otu_0.8_16 <- as.data.frame(lmProfile_0.8_16$fit$importance)
resample <- c("Resample 16")
otu_select <- rownames(importance_otu_0.8_16)[1:13]
otu_train_m_select_0.8_16 <- otu_train_m_0.8_16[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_16 <- otu_test_m_0.8_16[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_16$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_16$Group
predictions_final_0.8_16 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_16 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_16$fit, otu_test_m_select_0.8_16, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_16$Group
predictions_0.8_16 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_16 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 17 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_17 <- metadata[idx, ]
otu_test_0.8_17 <- metadata[-idx, ]
otu_train_m_0.8_17 <- otu_train_0.8_17[, c(3,314:401))]
otu_test_m_0.8_17 <- otu_test_0.8_17[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_17$Group
x <- otu_train_m_0.8_17[, 2:89]

lmProfile_0.8_17<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_17
importance_otu_0.8_17 <- as.data.frame(lmProfile_0.8_17$fit$importance)
resample <- c("Resample 17")
otu_select <- rownames(importance_otu_0.8_17)[1:15]
otu_train_m_select_0.8_17 <- otu_train_m_0.8_17[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_17 <- otu_test_m_0.8_17[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_17$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_17$Group
predictions_final_0.8_17 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_17 <-cbind(sen,speci,resample) 
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_17$fit, otu_test_m_select_0.8_17, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_17$Group
predictions_0.8_17 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_17 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 18 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_18 <- metadata[idx, ]
otu_test_0.8_18 <- metadata[-idx, ]
otu_train_m_0.8_18 <- otu_train_0.8_18[, c(3,314:401)]
otu_test_m_0.8_18 <- otu_test_0.8_18[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_18$Group
x <- otu_train_m_0.8_18[, 2:89]

lmProfile_0.8_18<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_18
importance_otu_0.8_18 <- as.data.frame(lmProfile_0.8_18$fit$importance)


importance_otu_0.8_18 <- cbind(resample, importance_otu_0.8_18)
importance_otu_0.8_18$predictors <- row.names(importance_otu_0.8_18)
resample <- c("Resample 18")
otu_select <- rownames(importance_otu_0.8_18)[1:19]
otu_train_m_select_0.8_18 <- otu_train_m_0.8_18[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_18 <- otu_test_m_0.8_18[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_18$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_18$Group
predictions_final_0.8_18 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_18 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_18$fit, otu_test_m_select_0.8_18, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_18$Group
predictions_0.8_18 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_18 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 19 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_19 <- metadata[idx, ]
otu_test_0.8_19 <- metadata[-idx, ]
otu_train_m_0.8_19 <- otu_train_0.8_19[, c(3,314:401)]
otu_test_m_0.8_19 <- otu_test_0.8_19[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_19$Group
x <- otu_train_m_0.8_19[, 2:89]

lmProfile_0.8_19<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_19
importance_otu_0.8_19 <- as.data.frame(lmProfile_0.8_19$fit$importance)
resample <- c("Resample 19")
otu_select <- rownames(importance_otu_0.8_19)[1:28]
otu_train_m_select_0.8_19 <- otu_train_m_0.8_19[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_19 <- otu_test_m_0.8_19[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_19$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_19$Group
predictions_final_0.8_19 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_19 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_19$fit, otu_test_m_select_0.8_19, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_19$Group
predictions_0.8_19<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_19 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 20 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_20 <- metadata[idx, ]
otu_test_0.8_20 <- metadata[-idx, ]
otu_train_m_0.8_20 <- otu_train_0.8_20[, c(3,314:401)]
otu_test_m_0.8_20 <- otu_test_0.8_20[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_20$Group
x <- otu_train_m_0.8_20[, 2:89]

lmProfile_0.8_20<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_20
importance_otu_0.8_20 <- as.data.frame(lmProfile_0.8_20$fit$importance)
resample <- c("Resample 20")
otu_select <- rownames(importance_otu_0.8_20)[1:12]
otu_train_m_select_0.8_20 <- otu_train_m_0.8_20[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_20 <- otu_test_m_0.8_20[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_20$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_20$Group
predictions_final_0.8_20 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_20 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_20$fit, otu_test_m_select_0.8_20, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_20$Group
predictions_0.8_20<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_20 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 21 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_21 <- metadata[idx, ]
otu_test_0.8_21 <- metadata[-idx, ]
otu_train_m_0.8_21 <- otu_train_0.8_21[, c(3,314:401)]
otu_test_m_0.8_21 <- otu_test_0.8_21[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_21$Group
x <- otu_train_m_0.8_21[, 2:89]

lmProfile_0.8_21<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_21
importance_otu_0.8_21 <- as.data.frame(lmProfile_0.8_21$fit$importance)
resample <- c("Resample 21")
otu_select <- rownames(importance_otu_0.8_21)[1:23]
otu_train_m_select_0.8_21 <- otu_train_m_0.8_21[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_21 <- otu_test_m_0.8_21[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_21$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_21$Group
predictions_final_0.8_21 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_21 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_21$fit, otu_test_m_select_0.8_21, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_21$Group
predictions_0.8_21<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_21 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 22 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_22 <- metadata[idx, ]
otu_test_0.8_22 <- metadata[-idx, ]
otu_train_m_0.8_22 <- otu_train_0.8_22[, c(3,314:401)]
otu_test_m_0.8_22 <- otu_test_0.8_22[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_22$Group
x <- otu_train_m_0.8_22[, 2:89]

lmProfile_0.8_22<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_22
importance_otu_0.8_22 <- as.data.frame(lmProfile_0.8_22$fit$importance)
resample <- c("Resample 22")
otu_select <- rownames(importance_otu_0.8_22)[1:13]
otu_train_m_select_0.8_22 <- otu_train_m_0.8_22[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_22 <- otu_test_m_0.8_22[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_22$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_22$Group
predictions_final_0.8_22 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_22 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_22$fit, otu_test_m_select_0.8_22, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_22$Group
predictions_0.8_22<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_22 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 23 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_23 <- metadata[idx, ]
otu_test_0.8_23 <- metadata[-idx, ]
otu_train_m_0.8_23 <- otu_train_0.8_23[, c(3,314:401)]
otu_test_m_0.8_23 <- otu_test_0.8_23[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_23$Group
x <- otu_train_m_0.8_23[, 2:89]

lmProfile_0.8_23<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_23

importance_otu_0.8_23 <- as.data.frame(lmProfile_0.8_23$fit$importance)
resample <- c("Resample 23")
otu_select <- rownames(importance_otu_0.8_23)[1:27]
otu_train_m_select_0.8_23 <- otu_train_m_0.8_23[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_23 <- otu_test_m_0.8_23[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_23$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_23$Group
predictions_final_0.8_23 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_23 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_23$fit, otu_test_m_select_0.8_23, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_23$Group
predictions_0.8_23 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_23 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 24 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_24 <- metadata[idx, ]
otu_test_0.8_24 <- metadata[-idx, ]
otu_train_m_0.8_24 <- otu_train_0.8_24[, c(3,314:401)]
otu_test_m_0.8_24 <- otu_test_0.8_24[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_24$Group
x <- otu_train_m_0.8_24[, 2:89]

lmProfile_0.8_24<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_24
importance_otu_0.8_24 <- as.data.frame(lmProfile_0.8_24$fit$importance)
resample <- c("Resample 24")
otu_select <- rownames(importance_otu_0.8_24)[1:9]
otu_train_m_select_0.8_24 <- otu_train_m_0.8_24[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_24 <- otu_test_m_0.8_24[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_24$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_24$Group
predictions_final_0.8_24 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_24 <-cbind(sen,speci,resample)
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_24$fit, otu_test_m_select_0.8_24, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_24$Group
predictions_0.8_24<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_24 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))


# Repeat 25 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)

p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_25 <- metadata[idx, ]
otu_test_0.8_25 <- metadata[-idx, ]
otu_train_m_0.8_25 <- otu_train_0.8_25[, c(3,314:401)]
otu_test_m_0.8_25 <- otu_test_0.8_25[, c(3,314:401)]
set.seed(123)
y <- otu_train_m_0.8_25$Group
x <- otu_train_m_0.8_25[, 2:89]

lmProfile_0.8_25<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_25
importance_otu_0.8_25 <- as.data.frame(lmProfile_0.8_25$fit$importance)
resample <- c("Resample 25")
otu_select <- rownames(importance_otu_0.8_25)[1:44]
otu_train_m_select_0.8_25 <- otu_train_m_0.8_25[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_25 <- otu_test_m_0.8_25[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_25$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_25$Group
predictions_final_0.8_25<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_25 <-cbind(sen,speci,resample) 
                                     
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_25$fit, otu_test_m_select_0.8_25, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_25$Group
predictions_0.8_25<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_25 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))

## ROC curve and AUC
sen_speci_merged_final_control_RBD <- rbind(sen_speci_final_0.8_01, sen_speci_final_0.8_02, sen_speci_final_0.8_03, sen_speci_final_0.8_04
                                         , sen_speci_final_0.8_05, sen_speci_final_0.8_06, sen_speci_final_0.8_07, sen_speci_final_0.8_08
                                         , sen_speci_final_0.8_09, sen_speci_final_0.8_10, sen_speci_final_0.8_11, sen_speci_final_0.8_12
                                         , sen_speci_final_0.8_13, sen_speci_final_0.8_14, sen_speci_final_0.8_15, sen_speci_final_0.8_16
                                         , sen_speci_final_0.8_17, sen_speci_final_0.8_18, sen_speci_final_0.8_19, sen_speci_final_0.8_20
                                         , sen_speci_final_0.8_21, sen_speci_final_0.8_22, sen_speci_final_0.8_23, sen_speci_final_0.8_24, 
                                                sen_speci_final_0.8_25)

sen_speci_merged_control_RBD <- rbind(sen_speci_0.8_01, sen_speci_0.8_02, sen_speci_0.8_03, sen_speci_0.8_04
                                          , sen_speci_0.8_05, sen_speci_0.8_06, sen_speci_0.8_07, sen_speci_0.8_08
                                          , sen_speci_0.8_09, sen_speci_0.8_10, sen_speci_0.8_11, sen_speci_0.8_12
                                          , sen_speci_0.8_13, sen_speci_0.8_14, sen_speci_0.8_15, sen_speci_0.8_16
                                          , sen_speci_0.8_17, sen_speci_0.8_18, sen_speci_0.8_19, sen_speci_0.8_20
                                          , sen_speci_0.8_21, sen_speci_0.8_22, sen_speci_0.8_23, sen_speci_0.8_24, sen_speci_0.8_25)


roc_control_RBD_train_per_sample <- ggplot(sen_speci_merged_final_control_RBD, aes(x = 1-specificities, y = sensitivities)) + 
                                     geom_step() + geom_point() +
                                     theme(aspect.ratio = 1)
                                     
roc_control_RBD_test_per_sample <- ggplot(sen_speci_merged_test_control_RBD, aes(x = 1-specificities, y = sensitivities)) + 
                                     geom_step() + geom_point() +
                                     theme(aspect.ratio = 1)
library(pROC)
predictions_test_RBD_control_all <- rbind(predictions_0.8_01, predictions_0.8_02,predictions_0.8_03,
                             predictions_0.8_04,predictions_0.8_05,predictions_0.8_06,
                             predictions_0.8_07,predictions_0.8_08,predictions_0.8_09,
                             predictions_0.8_10,predictions_0.8_11,predictions_0.8_12,
                             predictions_0.8_13,predictions_0.8_14,predictions_0.8_15,
                             predictions_0.8_16,predictions_0.8_17,predictions_0.8_18,
                             predictions_0.8_19,predictions_0.8_20,predictions_0.8_21,
                             predictions_0.8_22,predictions_0.8_23,predictions_0.8_24,
                             predictions_0.8_25)
predictions_final_RBD_control_all <- rbind(predictions_final_0.8_01, predictions_final_0.8_02,predictions_final_0.8_03,
                                          predictions_final_0.8_04,predictions_final_0.8_05,predictions_final_0.8_06,
                                          predictions_final_0.8_07,predictions_final_0.8_08,predictions_final_0.8_09,
                                          predictions_final_0.8_10,predictions_final_0.8_11,predictions_final_0.8_12,
                                          predictions_final_0.8_13,predictions_final_0.8_14,predictions_final_0.8_15,
                                          predictions_final_0.8_16,predictions_final_0.8_17,predictions_final_0.8_18,
                                          predictions_final_0.8_19,predictions_final_0.8_20,predictions_final_0.8_21,
                                          predictions_final_0.8_22,predictions_final_0.8_23,predictions_final_0.8_24,
                                          predictions_final_0.8_25) 
                                     

roc_final_RBD_control_mean <-roc(ifelse(predictions_final_RBD_control_all$observed=="RBD", "RBD", "Control"), as.numeric(predictions_final_RBD_control_all$RBD), plot = TRUE,
                           direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
                           thresholds="best", print.thres="best")
auc(roc_final_RBD_control_mean)
ci.auc(roc_final_RBD_control_mean, conf.level=0.95, method=c("delong",
                                                      "bootstrap"), boot.n = 2000)
                                     
roc_test_RBD_control_mean <-roc(ifelse(predictions_test_RBD_control_all$observed=="RBD", "RBD", "Control"), as.numeric(predictions_test_RBD_control_all$RBD), plot = TRUE,
              direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
              thresholds="best", print.thres="best")
auc(roc_test_RBD_control_mean)
ci.auc(roc_test_RBD_control_mean, conf.level=0.95, method=c("delong",
                                         "bootstrap"), boot.n = 2000) 
                                     
########################### Mediation Analysis ###########################

library(mediation)
metadata <- read.csv('./metadata.csv')
which(colnames(metadata)=="rel_g001")
which(colnames(metadata)=="rel_g249")
otu <- metadata[ , c(28:276)]
row.names(otu) <- metadata$Sample_ID
dis <- vegdist(otu, method = 'bray')
pcoa <- wcmdscale(dis, eig = T, x.ret = T, add = "cailliez")
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_group <- data.frame({pcoa$point})[1:2]
names(sample_group)[1:2] <- c('PCoA1', 'PCoA2')
sample_group$Sample_ID <-row.names(sample_group)
metadata <- merge(sample_group, metadata, by = 'Sample_ID')


# In the mediation analysis, we used total likelihood ratio of prodromal PD (excluding items of RBD and constipation) 
# to represent the risk of developing PD, the value of first principal component (derived from the PCoA analysis) 
# and bowel movement frequency score to indicate the microbial composition and the severity of constipation, respectively.

# Here is the first mediation analysis

#                microbiota
#                (mediator)
#              /            \
#             /              \
# constipation ---------------> risk of PD
#  (exposure)                    (outcome)

fit.totaleffect=glm(LR_prodromal_pd_excluding_FC_RBD_LGtransition~BMF_score+age+gender, family =gaussian, metadata)
summary(fit.totaleffect) # total effect

fit.mediator=glm(PCoA1~BMF_score+age+gender, family=gaussian, metadata)
summary(fit.mediator) # the effect of exposure onto the mediator

fit.dv=glm(LR_prodromal_pd_excluding_FC_RBD_LGtransition~PCoA1+BMF_score+age+gender, family =gaussian, metadata)
summary(fit.dv) # the effect of the mediator on the outcome

results = mediate(fit.mediator, fit.dv, treat='BMF_score', mediator='PCoA1', sims = 10000, boot=T)
summary(results) # mediation analysis


# the second mediation analysis

#              constipation
#               (mediator)
#             /           \
#            /             \
# microbiota ----------------> risk of PD
# (exposure)                   (outcome)

fit.totaleffect=glm(LR_prodromal_pd_excluding_FC_RBD_LGtransition~PCoA1+age+gender, family =gaussian, metadata)
summary(fit.totaleffect)

fit.mediator=glm(BMF_score~PCoA1+age+gender,family=gaussian, metadata)
summary(fit.mediator)

fit.dv=glm(LR_prodromal_pd_excluding_FC_RBD_LGtransition~PCoA1+BMF_score+age+gender, family =gaussian, metadata)
summary(fit.dv)

results = mediate(fit.mediator, fit.dv, treat='PCoA1', mediator='BMF_score', sims = 10000, boot=T)
summary(results)


########################### pathway abundance analysis ###########################

## Microbiome multivariable associations with linear model (MaAsLin 2)
# The input data is centered log-ratio (clr) transformed abundance of differential pathways (n = 18) 
# that identified by using Kruskal-Wallis test in SPSS (adjusted p-values < 0.05) 

library(Maaslin2)

input_metadata <- read.csv('./metadata.csv')
which(colnames(input_metadata)=="COBALSYN.PWY")
which(colnames(input_metadata)=="TCA")
input_data <- input_metadata[ , c(381:398)]

# Correlations of pathway abundance with early stages of alpha-synucleinopathy
fit_data_pathway <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_group_pathway_unadjusted', transform = "NONE",
  fixed_effects = c('Group'),
  reference = "Group, Control",
  random_effects = c('family_id'),
  normalization = 'NONE',
  correction = "BH", # multiple comparison adjustment
  standardize = FALSE)
                                     
fit_data_pathway <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_group_pathway_adjusted', transform = "NONE",
  fixed_effects = c('Group','sex','age'),
  reference = "Group, Control",
  random_effects = c('family_id'),
  normalization = 'NONE',
  correction = "BH", # multiple comparison adjustment
  standardize = FALSE)

View(fit_data_pathway$results)
res_pathway <- as.data.frame(fit_data_pathway$results)  

# The effect of host factors on pathway abundance
fit_data_pathway_cov <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_group_pathway_cov', transform = "NONE",
  fixed_effects = c('Group', "Statins_2gp",'gender','age','BMF_score',"Proton_pump_inhibitors_2gp","Osmotic_laxative_2gp", 
                    "Antidepressants_2gp", "Benzodiazepines_2gp"), 
  random_effects = c('family_id'),
  reference = "Group, Control",
  normalization = 'NONE',
  correction = "BH", # multiple comparison adjustment
  standardize = FALSE)

View(fit_data_pathway_cov$results)
res_pathway_cov <- as.data.frame(fit_data_pathway_cov$results)

# Pathway abundance (PD specific drugs)
input_metadata <- read.csv('./metadata.csv')
input_metadata <- input_metadata %>% filter(Group == "Early_PD")
which(colnames(input_metadata)=="COBALSYN.PWY")
which(colnames(input_metadata)=="TCA")
input_data <- input_metadata[ , c(381:398)]


fit_data_pathway_pd_drugs <- Maaslin2(
  input_data, input_metadata, 'maaslin_output_group_pathway_pd_drugs', transform = "NONE",
  fixed_effects = c("Carbidopa_Levodopa_2gp",'MAOB_inhibitor_2gp','Dopamine_agonist_2gp'), 
  normalization = 'NONE',
  correction = "BH",
  standardize = FALSE)

View(fit_data_pathway_pd_drugs$results)
res_pathway_pd_drugs <- as.data.frame(fit_data_pathway_pd_drugs$results)  

