library(SCOPE)
major.marker.genes  <- Reduce(rbind,list(
  data.frame(Gene=c("CCL2", "COL14A1", "CXCL14", "CXCL1", "CXCL12", "CXCL2", "CXCL8", "CCL13",'IL6', 'IL11',
                    'CCL5','CCL22','CLEC3B','CCL17','LIF','HGF',
                    'APOD','IGF1','C3','C7','ITM2A','MGP','CCL11'),Cell='iCAF'),
  data.frame(Gene=c("AOC3", "COL12A1", "COL15A1", "COL1A2", "COL1A1", "COL8A1","FAP", "MYL9", "TAGLN","MYLK", "TPM1", "TPM2",
                    "POSTN","HOPX","PDGFA",'COL5A1','VIM','COL6A1','CTHRC1','THY1',
                    'CTGF','CTA2','FBLN1','SERPINF1','VCAN','COL11A1','COL10A1','ACTC2'),Cell='mCAF'),
  data.frame(Gene=c("CD74", "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "CD52","IGFBP3"),Cell='apCAF')))
load("CAF_markers.rda")
mCAF <- unique(c(CAF_markers$myCAF_Cell,CAF_markers$Myofibro,
                 CAF_markers$mCAF_MC,major.marker.genes[major.marker.genes$Cell=="mCAF",]$Gene))
iCAF <- unique(c(na.omit(CAF_markers$iCAF_Cell),na.omit(CAF_markers$iCAF_MC),
                 major.marker.genes[major.marker.genes$Cell=="iCAF",]$Gene))
apCAF <- unique(c(na.omit(CAF_markers$apCAF_Cell),
                  major.marker.genes[major.marker.genes$Cell=="apCAF",]$Gene))
CAFS <- list(mCAF,iCAF,apCAF)
scob <- AddModuleScore(sc, features = CAFS,name = c("mCAF", "iCAF", "apCAF"))
colnames(scob@meta.data)[37:39] <- c("mCAFs Signature","iCAFs Signature","apCAFs Signature")
table(scob$CAF.type)
mycolors=c("#1E90FF","#20B2AA","#FFA500","#40E0D0","#9370DB","#DC143C")
p1 <- VlnPlot(scob, features= "mCAFs Signature", cols = mycolors, pt.size = 0) +  theme_classic(base_line_size = 0.7) +
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = " mCAFs Signature", x="") + theme(legend.position="none") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")
p2 <- VlnPlot(scob, features= "iCAFs Signature", cols = mycolors, pt.size = 0) +  theme_classic(base_line_size = 0.7) +
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = " iCAFs Signature", x="") + theme(legend.position="none") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")
p3 <- VlnPlot(scob, features= "apCAFs Signature", cols = mycolors, pt.size = 0) +  theme_classic(base_line_size = 0.7) +
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = " apCAFs Signature", x="") + theme(legend.position="bottom") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")
library(patchwork)
p1 + p2 + p3 + plot_layout(nrow = 3, byrow = FALSE)