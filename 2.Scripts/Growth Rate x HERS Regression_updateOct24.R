# Regression between HERS and growth rates
rm(list = ls())
graphics.off()

setwd("~/Documents/GitHub/SIBER")
#growth<-read.csv("Clam growth rates and HERS.csv", header=TRUE)

# new dataset with mean growth rates extracted from Tan et al 2022
#growth<-read.csv("Clam growth rates and HERS Tan et al 2022 data.csv", header=TRUE)

# new dataset 18-Oct-2024 with new HERS scores with SD
setwd("~/Documents/GitHub/SIBER")
growth<-read.csv("Clam growth rates Tan et al 2022 and new HERS scores CI_updatedOct24.csv", header=TRUE)

head(growth)


growth.hers.reg <- lm(Growth.rate ~ HERS, data=growth)
summary(growth.hers.reg)

model.growth.HERS = lm(growth$Growth.rate ~ growth$HERS)
summary(model.growth.HERS)


# HERS x Growth Rate Regression Figure

library(ggplot2)
#library(grid)
#library(scales)
#library(gridExtra)
#library(dplyr)

axis.ticks.length=unit(-0.5, units = "cm")
panel.border = element_rect(colour = "black", fill=NA, size=0.2)

# T. maximum = undetermined
HERS.reg.plot<-ggplot(growth, aes(x=HERS, y=Growth.rate))+
  geom_segment(aes(x=0.15, y=1.685, xend=0.5713, yend=10), color="grey60", size = 1)+
  geom_point(aes(color=Trophic.1, shape=Trophic.1), size=5)+
  scale_color_manual(values=c("#39BEB1", "#7DB0DD", "grey50"))+
  theme_classic()+
  theme(text=element_text(size=14), 
        plot.margin = margin(0.2, 0.1, 0.1, 0.2, "cm"), 
        legend.title=element_blank(), 
        legend.position = c(0.12,0.93), legend.text=element_text(size=11), 
        legend.background=element_rect(color=NA, fill=NA),
        legend.margin = margin(0,2,1,0.5), 
        panel.border = element_rect(fill = NA, size = 1))+
  scale_y_continuous(limits=c(0,10))+
#  scale_x_continuous(limits=c(0,0.6))+
  ylab(expression("Growth rate (mm "~ month^{-1}~")"))+
  xlab("HERS Score")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  
  annotate("text", x=0.61, y=1.3, label="y = 19.735x-1.2748",size=4,hjust=1) +
  annotate("text", x=0.61, y=0.8, parse=TRUE, label=as.character(expression(paste(italic(R^2),"= 0.97"))), size=4,hjust=1) +
  annotate("text", x=0.61, y=0.2, parse=TRUE, label=as.character(expression(paste(italic(p), " = 0.0018"))),size=4,hjust=1)

HERS.reg.plot


# Colored according to conservation status
# Updated to reflect Tan et al 2022 growth rates
# Changed from SD to 0.89 CI 
HERS.reg.plot<-ggplot(growth, aes(x=HERS, y=Growth.rate, label=Species))+
  geom_segment(aes(x=0.05, y=1.3389, xend=0.57, yend=7.34802), color="grey70", size = 1)+ #for the graph to have 0 sligthly on the right of x axis
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=0.2, alpha=0.6)+
  geom_point(aes(color=Conservation.status), shape=19, size=3)+
  geom_text(nudge_x = -0.03, nudge_y = 0.3, check_overlap = T, size=12/.pt, fontface = 'italic') +
  labs(shape='Trophic strategy') +
  scale_color_manual("Conservation status", values=c("#56B4E9", "#D55E00"))+
  theme_classic()+
  theme(text=element_text(size=14), legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        plot.margin = margin(0.2, 0.1, 0.1, 0.2, "cm"), 
#        legend.title=element_blank(), 
        legend.position = c(0.18,0.9), legend.text=element_text(size=12), 
        legend.background=element_rect(color=NA, fill=NA),
        legend.margin = margin(0,2,1,0.5), 
        panel.border = element_rect(fill = NA, size = 0.5))+
  scale_y_continuous(limits=c(1,8))+
  #  scale_x_continuous(limits=c(0,0.6))+
  ylab(expression("Growth rate (mm "~ month^{-1}~")"))+
  xlab("Trophic niche score")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  
  annotate("text", x=0.61, y=1.7, label="y = -0.2775+12.985x",size=10/.pt,hjust=1) +
  annotate("text", x=0.61, y=1.4, parse=TRUE, label=as.character(expression(paste(italic(R^2),"= 0.58"))), size=10/.pt,hjust=1) +
  annotate("text", x=0.61, y=1.0, parse=TRUE, label=as.character(expression(paste(italic(p), " = 0.04719"))),size=10/.pt,hjust=1) 
#  guides(shape = guide_legend(order = 1, override.aes = list(size=2)),
#         col = guide_legend(order = 2, override.aes = list(size=2))) 
#  guides(shape = guide_legend(override.aes = list(size=1))) +
#  guides(colour = guide_legend(override.aes = list(size=1))) 

HERS.reg.plot

#ggsave(path="Plots",filename = "GrowthRate_TNS.pdf",device = "pdf")

#####test changing geam_segment
growth.hers.reg <- lm(Growth.rate ~ HERS, data=growth)
summary(growth.hers.reg)

# Intecept
intercept <- growth.hers.reg$coefficients[1]
slope <- growth.hers.reg$coefficients[2]

# Get the range of HERS
x_start <- min(growth$HERS)
x_end <- max(growth$HERS)

# Calculate the corresponding y values
y_start <- intercept + slope * x_start
y_end <- intercept + slope * x_end

#Another way to get the values
summary(growth.hers.reg)$r.squared
summary(growth.hers.reg)$coefficients[2,4]

#plot again
HERS.reg.plot<-ggplot(growth, aes(x=HERS, y=Growth.rate, label=Species))+
  geom_segment(aes(x=x_start, y=y_start, xend=x_end, yend=y_end), color="grey70", size = 1)+ #for the graph to have 0 sligthly on the right of x axis
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=0.2, alpha=0.6)+
  geom_point(aes(color=Conservation.status), shape=19, size=3)+
  geom_text(nudge_x = -0.05, nudge_y = 0.3, check_overlap = T, size=10/.pt, fontface = 'italic') +
  labs(shape='Trophic strategy') +
  scale_color_manual("Conservation status", values=c("#56B4E9", "#D55E00"))+
  theme_classic()+
  theme(text=element_text(size=14), legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        plot.margin = margin(0.2, 0.1, 0.1, 0.2, "cm"), 
        #        legend.title=element_blank(), 
        legend.position = c(0.18,0.9), legend.text=element_text(size=12), 
        legend.background=element_rect(color=NA, fill=NA),
        legend.margin = margin(0,2,1,0.5), 
        panel.border = element_rect(fill = NA, size = 0.5))+
  scale_y_continuous(limits=c(1,8))+
  #  scale_x_continuous(limits=c(0,0.6))+
  ylab(expression("Growth rate (mm "~ month^{-1}~")"))+
  xlab("Trophic niche score")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  
  annotate("text", x=0.61, y=1.7, label="y = -0.2775+12.985x",size=10/.pt,hjust=1) +
  annotate("text", x=0.61, y=1.4, parse=TRUE, label=as.character(expression(paste(italic(R^2),"= 0.68"))), size=10/.pt,hjust=1) +
  annotate("text", x=0.61, y=1.0, parse=TRUE, label=as.character(expression(paste(italic(p), " = 0.04719"))),size=10/.pt,hjust=1) 
#  guides(shape = guide_legend(order = 1, override.aes = list(size=2)),
#         col = guide_legend(order = 2, override.aes = list(size=2))) 
#  guides(shape = guide_legend(override.aes = list(size=1))) +
#  guides(colour = guide_legend(override.aes = list(size=1))) 



####
####
#FInal used with geom smooth even if no differences
HERS.reg.plot <- ggplot(growth, aes(x=HERS, y=Growth.rate)) +
  geom_smooth(method = "lm", se = FALSE, color="grey70", size = 1) + 
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=0.2, alpha=0.6) +
  geom_point(aes(color=Conservation.status), shape=19, size=3) +
  geom_text(aes(label=Species), nudge_x = 0.01, nudge_y = 0.2, check_overlap = T, size=10/.pt, fontface = 'italic') +
  expand_limits(y = max(growth$Growth.rate) + 1) + # increase y-axis limit
  labs(shape='Trophic strategy') +
  scale_color_manual("Conservation status", values=c("#56B4E9", "#D55E00"))+
  theme_classic()+
  theme(text=element_text(size=14), legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        plot.margin = margin(0.3, 0.2, 0.2, 0.3, "cm"), 
        #        legend.title=element_blank(), 
        legend.position = c(0.18,0.9), legend.text=element_text(size=12), 
        legend.background=element_rect(color=NA, fill=NA),
        legend.margin = margin(0,2,1,0.5), 
        panel.border = element_rect(fill = NA, size = 0.5))+
  scale_y_continuous(limits=c(1,8))+
   scale_x_continuous(limits=c(0.2,0.6))+
  ylab(expression("Growth rate (mm "~ month^{-1}~")"))+
  xlab("Trophic niche score")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  
  annotate("text", x=0.61, y=1.7, label="y = -0.2775+12.985x",size=10/.pt,hjust=1) +
  annotate("text", x=0.61, y=1.4, parse=TRUE, label=as.character(expression(paste(italic(R^2),"= 0.68"))), size=10/.pt,hjust=1) +
  annotate("text", x=0.61, y=1.0, parse=TRUE, label=as.character(expression(paste(italic(p), " = 0.04719"))),size=10/.pt,hjust=1) 
#  guides(shape = guide_legend(order = 1, override.aes = list(size=2)),
#         col = guide_legend(order = 2, override.aes = list(size=2))) 
#  guides(shape = guide_legend(override.aes = list(size=1))) +
#  guides(colour = guide_legend(override.aes = list(size=1))) 
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) + # adjust plot margins
  print(HERS.reg.plot)
  ggsave(path="Plots",filename = "GrowthRate_TNS_updatedgeomsmooth.pdf",device = "pdf")


#without t.maxima if considered as outlier
#remove T.max
growth2 <-growth[-5, ]

#lm
growth.hers.reg <- lm(Growth.rate ~ HERS, data=growth2)
summary(growth.hers.reg)

model.growth.HERS = lm(growth2$Growth.rate ~ growth2$HERS)
summary(model.growth.HERS)

# no T.max
# Colored according to conservation status
# Updated to reflect Tan et al 2022 growth rates
# Changed from SD to 0.89 CI 
#geom_segment not changed here, would need to that if figure needed
HERS.reg.plot<-ggplot(growth2, aes(x=HERS, y=Growth.rate, label=Species))+
  geom_segment(aes(x=0.05, y=1.3389, xend=0.57, yend=7.34802), color="grey70", size = 1)+ #for the graph to have 0 sligthly on the right of x axis
  geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=0.2, alpha=0.6)+
  geom_point(aes(color=Conservation.status), shape=19, size=3)+
  geom_text(nudge_x = -0.03, nudge_y = 0.3, check_overlap = T, size=12/.pt, fontface = 'italic') +
  labs(shape='Trophic strategy') +
  scale_color_manual("Conservation status", values=c("#56B4E9", "#D55E00"))+
  theme_classic()+
  theme(text=element_text(size=14), legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        plot.margin = margin(0.2, 0.1, 0.1, 0.2, "cm"), 
        #        legend.title=element_blank(), 
        legend.position = c(0.18,0.9), legend.text=element_text(size=12), 
        legend.background=element_rect(color=NA, fill=NA),
        legend.margin = margin(0,2,1,0.5), 
        panel.border = element_rect(fill = NA, size = 0.5))+
  scale_y_continuous(limits=c(1,8))+
  #  scale_x_continuous(limits=c(0,0.6))+
  ylab(expression("Growth rate (mm "~ month^{-1}~")"))+
  xlab("Trophic niche score")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  
  annotate("text", x=0.61, y=1.7, label="y = 0.6450+11.2537x",size=10/.pt,hjust=1) +
  annotate("text", x=0.61, y=1.4, parse=TRUE, label=as.character(expression(paste(italic(R^2),"= 0.85"))), size=10/.pt,hjust=1) +
  annotate("text", x=0.61, y=1.0, parse=TRUE, label=as.character(expression(paste(italic(p), " = 0.0262"))),size=10/.pt,hjust=1) 
#  guides(shape = guide_legend(order = 1, override.aes = list(size=2)),
#         col = guide_legend(order = 2, override.aes = list(size=2))) 
#  guides(shape = guide_legend(override.aes = list(size=1))) +
#  guides(colour = guide_legend(override.aes = list(size=1))) 

HERS.reg.plot

#ggsave(path="Plots",filename = "GrowthRate_TNS_noTmax.pdf",device = "pdf")
