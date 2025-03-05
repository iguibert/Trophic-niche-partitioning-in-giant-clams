####################################################################################
## Phyllogeny : phyllogenetic tree & comparison with giant clams trait            ####
## Data: DNA data from Tan.al, traits from the same publication                   ####
## Data: HERS results                                                             ####
## Leonard Pons                                                                   ####
##################################################################################

# Packages #####

#A list of all packages used for this script
#Clean our global environment
rm(list=ls())
graphics.off()
Sys.setenv(LANG = "en")

library(data.table)
library(tidyverse)
library(dplyr)
library(readxl)
library(cluster)
library(ggdendro)
library(dendextend)
library(phylotools)
library(phylobase)
library(scales)
library(rphylopic)
library(picante)
library(ggtree)


# Data ####
# set work directory
setwd("~/Travail/Hong Kong/Clams/SIA data/Script/Phylosignal analyses")

# giant clam isotopic data with metadata


LeoNeo_traits <- read_excel("LeoNeo_traits_new_2024.xlsx",
                            col_types = c("text", rep("numeric",8)))


traits <- list()
it <- 1
for (i in 2:9) {
  
  
  
  traits[[it]] <- aggregate(as.formula(paste(names(LeoNeo_traits[,i]), "~ species")),
                            LeoNeo_traits, mean)
  names(traits)[[it]] <- names(LeoNeo_traits[i])
  it <- it + 1
  
}
traits_df <- traits[[1]]
for (i in 2:8) {
  
  traits_df <- inner_join(traits_df,traits[[i]])
  
}


row.names(traits_df) <- traits_df$species
traits_df <- traits_df[,-1]

traits_df_HERS <- data.frame("HERS" = traits_df[,"HERS_leo_average"],
                             'Species' = row.names(traits_df))

relatives_traits <- list()
it <- 1
for (i in 1:ncol(traits_df)) {
  
  relatives_traits[[it]] <- scales::rescale(c(0,traits_df[,i]), to = c(0,1))[-1]
  names(relatives_traits)[[it]] <- names(traits_df[i])
  
  it <- it + 1
}
relatives_traits <- as.data.frame(relatives_traits)
row.names(relatives_traits) <- row.names(traits_df)


phylofinal <- read.tree("Sixgiantclamsandoutgroup")

#cool isn't it ?
plot(phylofinal)
#remame tip.label
#to be sure of the order, display the tips label name, it will show you the right order
phylofinal$tip.label
#create a list of new name you want to change. BE CAREFULL THE ORDER IS IMPORTANT !!
phylofinal$tip.label <- c("Fulvia mutica","Hippopus hippopus","Hippopus porcellanus","Tridacna derasa",'Tridacna squamosa','Tridacna maxima',"Tridacna gigas")
plot(phylofinal)
phylofinal <- dendextend::rotate(phylofinal,10)
plot(phylofinal)

tree <- drop.tip(phylofinal,1)
plot(tree)

# PhyloLoop ####


rbind(traits_df_HERS, data.frame( "HERS" = sample(seq(from = 0, to = 1, by = 0.001),6),
                                  'Species' = traits_df_HERS$Species))
#  significant or close to significant patern partern (or pvalue within 0.1)

multiphylo_list <- list()
it <- 1
for (i in 1:1000) {
  
  multiphylo_list[[it]] <-multiPhylosignal(x = traits_df, phy = tree, reps = 720)
  message((it/1000)*100,"%")
  it <- it + 1
  
}

multiphylo_df <- multiphylo_list[[1]]
multiphylo_df$traitsmultiphylo_df <- row.names(multiphylo_list[[1]])

for (i in 2:length(multiphylo_list)) {
  
  tmp_df <- multiphylo_list[[i]]
  tmp_df$traitsmultiphylo_df <- row.names(multiphylo_list[[i]])
  
  multiphylo_df <- rbind(multiphylo_df,tmp_df)
  
}

row.names(multiphylo_df) <-  NULL
couleursrdm <- c()
for (i in seq(from=1, to = 24, by = 6)) {
  
  R <- sample((0:255),1)
  G <- sample((0:255),1)
  B <- sample((0:255),1)
  
  revR <- abs(R - 255)
  revG <- abs(G - 255)
  revB <- abs(B - 255)
  
  couleursrdm[i] <- rgb(R,G,B,maxColorValue = 255)
  couleursrdm[i+1] <- rgb(revR,revG,revB,maxColorValue = 255)
  couleursrdm[i+2] <- rgb(B,R,G,maxColorValue = 255)
  couleursrdm[i+3] <- rgb(revB,revR,revG,maxColorValue = 255)
  couleursrdm[i+4] <- rgb(G,B,R,maxColorValue = 255)
  couleursrdm[i+5] <- rgb(revG,revB,revR,maxColorValue = 255)
}


multiphylo_plot <- multiphylo_df
labels <- unique(multiphylo_plot$traitsmultiphylo_df)

multiphylo_plot[multiphylo_plot$traitsmultiphylo_df == labels[1], "traitsmultiphylo_df"] <- "DEC"
multiphylo_plot[multiphylo_plot$traitsmultiphylo_df == labels[2], "traitsmultiphylo_df"] <- "HERS score average"
multiphylo_plot[multiphylo_plot$traitsmultiphylo_df == labels[3], "traitsmultiphylo_df"] <- "Average D15N"
multiphylo_plot[multiphylo_plot$traitsmultiphylo_df == labels[4], "traitsmultiphylo_df"] <- "Average D13C"
multiphylo_plot[multiphylo_plot$traitsmultiphylo_df == labels[5], "traitsmultiphylo_df"] <- "Provinces"
multiphylo_plot[multiphylo_plot$traitsmultiphylo_df == labels[6], "traitsmultiphylo_df"] <- "Growth rate average"
multiphylo_plot[multiphylo_plot$traitsmultiphylo_df == labels[7], "traitsmultiphylo_df"] <- "Max depth"
multiphylo_plot[multiphylo_plot$traitsmultiphylo_df == labels[8], "traitsmultiphylo_df"] <- "Shell length max."

ggplot(multiphylo_plot,aes(group = traitsmultiphylo_df, x = multiphylo_plot$traitsmultiphylo_df,
                           y = PIC.variance.P, colour = traitsmultiphylo_df))+
  geom_jitter(width=0.25, alpha = 0.35, show.legend = FALSE)+
  geom_boxplot(alpha=0.25, outlier.alpha=0, colour = "black",show.legend = FALSE) +
  scale_colour_manual(values = couleursrdm)+
  geom_hline(yintercept = 0.05, size = 1, col = "red", alpha = 0.5)+
  geom_hline(yintercept = 0.1, size = 1, linetype = 2, alpha = 0.5)+
  theme_bw(base_size = 22, base_family = "Helvetica", base_line_size = 22/22, base_rect_size = 22/22) %+replace% 
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", size = 0.4),
        panel.grid.minor.x = element_line(colour = "grey80", size = 0.4),
        panel.grid.major.y = element_line(colour = "grey60", size = 0.5),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.line = element_line(colour = "black",size = rel(1)),
        legend.key = element_blank(), complete = TRUE,
        strip.background = element_rect(fill = "white",colour = "black", size = rel(2)))+
  ylim(0,1)+
  xlab("")+
  ylab("P.value")+
  labs(colour = " ")





multiphylo_df_agregated <- aggregate(multiphylo_df[-6] ,by = list(multiphylo_df$traitsmultiphylo_df), FUN = "mean")

signif_multiphylo_df_agregated <- filter(multiphylo_df_agregated, multiphylo_df_agregated$PIC.variance.P<=0.1)

signif_multiphylo_df_agregated$Group.1[1] <- "Growth rate average"
signif_multiphylo_df_agregated$Group.1[2] <- "Max depth"
signif_multiphylo_df_agregated$Group.1[3] <- "HERS score average"
signif_multiphylo_df_agregated$Group.1[4] <- "Shell length max."
signif_multiphylo_df_agregated

colnames(relatives_traits)[6] <- "Growth rate average"
colnames(relatives_traits)[7] <- "Max depth"
colnames(relatives_traits)[2] <- "HERS score average"
colnames(relatives_traits)[8] <- "Shell length max."
colnames(relatives_traits)

relatives_traits_signif <- select( relatives_traits , signif_multiphylo_df_agregated$Group.1)
relatives_traits_signif
relatives_traits_signif[7,] <- rep(NA,ncol(relatives_traits_signif))
rownames(relatives_traits_signif)[7] <- "Fulvia mutica"

tree_fig <- phylo4d(phylofinal, relatives_traits_signif, missing.data="warn")
tree_fig@data
tree_fig_metadata <- data.frame(x = seq(0.75, 2.05, length.out = length(relatives_traits_signif)),
                                lab = colnames(relatives_traits_signif),
                                p_value= format(round(signif_multiphylo_df_agregated$PIC.variance.P,3),nsmall=3),
                                K = format(round(signif_multiphylo_df_agregated$K,3),nsamll=3),
                                mark = format(c("*","*","*","*")))


ggtree(tree_fig,size = 1.3) + #loag the tree into ggtree (package which works like ggplot2)
  geom_tiplab(align=TRUE, linesize=0, font = 3 ,offset=0, linetype = "blank", size = 5, alpha = 1) + # parameters for little ligne
  geom_tiplab(align=TRUE, linesize=1, font = 3 ,offset=0, linetype = "dotted", size = 5, alpha = 0.4) + # parameters for little ligne
  geom_tippoint(aes(size = seq(from = 0, to = 0.8, length.out= 12)) , shape  = 21, colour = "black", x=tree_fig_metadata$x[4], show.legend =  TRUE, stroke = 1.15)+
  geom_tippoint(aes(size = Growth.rate.average) , shape  = 21, colour = "black",fill = couleursrdm[1], x=tree_fig_metadata$x[1], show.legend =  FALSE, stroke = 1.15)+
  geom_tippoint(aes(size = Max.depth) , shape  = 21, colour = "black", fill = couleursrdm[2], x=tree_fig_metadata$x[2], show.legend =  FALSE, stroke = 1.15)+
  geom_tippoint(aes(size = HERS.score.average) , shape  = 21, colour = "black", fill = couleursrdm[3], x=tree_fig_metadata$x[3], show.legend =  FALSE, stroke = 1.15)+
  geom_tippoint(aes(size = Shell.length.max.) , shape  = 21, colour = "black", fill = couleursrdm[4], x=tree_fig_metadata$x[4], show.legend =  FALSE, stroke = 1.15)+
  geom_text(aes(x = x, y = 7.5, label = lab),size = 5.5, data = tree_fig_metadata, angle = -45, hjust = 1, family = "Helvetica", font = 2) +
  geom_text(aes(x = x, y =0.5 , label = mark, family = "Helvetica", font = 2), size = 9, data = tree_fig_metadata,  angle = 0, hjust = 0.5) +
  geom_text(aes(x = x, y = 0.3 , label = p_value, family = "Helvetica", font = 2), size = 5.5, data = tree_fig_metadata, angle = 0, hjust = 0.5) +
  geom_text(aes(x = 0.35, y =0.3 , label = "p-value"), size = 5, angle = 0, hjust = 0.5) +
  geom_text(aes(x = x, y =0 , label = K, family = "Helvetica", font = 2), size = 5.5, data = tree_fig_metadata, angle = 0, hjust = 0.5) +
  geom_text(aes(x = 0.35, y =0 , label = "K"), size = 6, angle = 0, hjust = 0.5) +
  scale_size_continuous(range = c(0,15), name="") + 
  xlim(0,2.7)+ # let some space on the right to see labels text
  ylim(0,9.5)+
  geom_treescale(x = 0 , y= 0.5, offset = 0.05, linesize = 0.7, fontsize = 5, family = "Helvetica") # add a little scale at the bottom of the tree

