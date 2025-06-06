#######################################################################################
## HERS score script                                                                   ####
## Data: giant clams isotope data - Bakers Lab - HongKongUni - SWIRE                   ####
## Script Author:Leonard Pons (Not optimized)                                          ####
## Script DOI:
## Manuscript: Guibert et al. Niche partitioning in giant clams
##
## If you use this script, please cite the manuscript and script DOI. Thank you! ####
#################################################################################################
# Packages #####

#A list of all packages used for this script

#Clean our global environment
rm(list=ls())
graphics.off()
Sys.setenv(LANG = "en")

library(data.table)
library(dplyr)
library(readr)
library(SIBER)
library(tidyverse)
library(rjags)
library(xlsx)

# Data ####

# your dataset should contain
# d13C
# d15N
# TYPE (Host/symb)
# Species names or site (species or site)
Yourdataset <- read.csv("~/Travail/Korea/Travail/Publication/HERS - Methodology publication - Inga/Script/GiantClams - Isotopes - HKU - Isis Guibert.csv", sep=";")

Yourdataset <- Yourdataset[,c(1,2,6,7)]
Yourdataset$Tissues <- gsub("[^A-Z]","",Yourdataset$Tissues)

# Run HERS 20 times and plot ellipses ######


# R doesn't provide a statistical mode function
# Here is a ok way to compute it
Mode <- function(x){
  
  uni.vect <- unique(x)
  
  digit.max.min <- c(8,3)
  while ((length(uni.vect)) > 0.1*length(x) & digit.max.min[1] > digit.max.min[2]) {
    
    x <- round(x, digits = digit.max.min[1])
    uni.vect <- unique(x)
    digit.max.min[1] <- digit.max.min[1] - 1
    
  }
  
  count.values <- tabulate(match(x, uni.vect))
  most.occur.ind <- sort(count.values, decreasing = TRUE)[1:5]
  most.occur.val <- c(uni.vect[count.values == most.occur.ind[1]],
                      uni.vect[count.values == most.occur.ind[2]],
                      uni.vect[count.values == most.occur.ind[3]],
                      uni.vect[count.values == most.occur.ind[4]],
                      uni.vect[count.values == most.occur.ind[5]])
  mean(most.occur.val, na.rm = TRUE)
  
}


# This loop is a bit long but will do everything: (8 sec * iteration * number of synthetic datasets)
# from creating a SIBER object
# to plot biplot with ellipses including relevant information

# How many times SIBER should run over the same dataset 
# The goal is to take into account the inherent variability of Bayesian statistics
# According to stability test 20 repetition is enough
iteration <- 20

# Split your data frame with all the species into a list of one data frame per species 
Yourdataset_splited <- split(Yourdataset, Yourdataset$Species)

# Set up list to retreive raw results
Master_list <- list()

# Get all SIBER metrics needed to compute HERS score
for (i in 1:length(Yourdataset_splited) #For all data frame created (1 per specie)
) {
  
  # You might want to change this line !!! You just have to select:
  # d13C
  # d15N
  # TYPE (Host/symb)
  # Species names (species)
  # in this specific order
  df_siber <- data.frame(Yourdataset_splited[[i]])[,c(1,2,3,4)]
  
  
  # Create objects wich will be used in the script
  Master_list[[i]] <- list(
    Species = c(), 
    # Standard metrics (40% computed with Bayesian model)
    SEAH = c(), # Standard Ellipse Area - Host
    SEAS = c(), # Standard Ellipse Area - Symbiont
    SEAO = c(), # Standard Ellipse Area - Overlap
    SEABH = c(), # SEAO / SEAH (mode based on 50 posterior distribution)
    SEABS = c(), # SEAO / SEAS (mode based on 50 posterior distribution)
    # Majors metrics (95% computed with Bayesian model)
    MEAH = c(), # Major Ellipse Area - Host
    MEAS = c(), # Major Ellipse Area - Symbiont
    MEAO = c(), # Major Ellipse Area - Overlap
    MEABH = c(), # MEAO / MEAH (mode based on 50 posterior distribution)
    MEABS = c() # MEAO / MEAS (mode based on 50 posterior distribution)
  )
  
  # Get the species name
  names(Master_list)[i] <- paste0(df_siber$Species[1])
  Master_list[[i]][["Species"]][1] <- df_siber$Species[1]
  
  # Artificialy create a community metric for SIBER
  df_siber$Species <- 1
  # be sure your (Host/Symb) colum is named => TYPE
  df_siber$Tissues[df_siber$Tissues == "H"] <- 1
  df_siber$Tissues[df_siber$Tissues == "S"] <- 2
  # Rename all colums to be SIBER friendly
  df_siber <- rename(df_siber, 
                     community = Species,
                     iso1 = d13C,
                     iso2 = d15N,
                     group = Tissues)
  
  #Create SIBER object
  df_siberized <- createSiberObject(df_siber)
  
  for (it in 1:iteration) {
    
    
    Master_list[[i]][["Species"]][it] <- Master_list[[i]][["Species"]][1]
    
    # Message
    message("Iteration: ", it,"/",iteration,
            " Total: ", round(((((i-1)*iteration) + it) / (length(Yourdataset_splited)*iteration)) * 100, digit = 3),
            "%"," (Specie ", Master_list[[i]][["Species"]][it], ")")
    
    
    
    # FROM HERE IT JUST SIBER ANALYSIS
    
    # Parameter: change them if you want
    parms <- list()
    parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
    parms$n.burnin <- 1 * 10^3 # discard the first set of values
    parms$n.thin <- 10     # thin the posterior by this many
    parms$n.chains <- 2        # run this many chains
    
    # define the priors: don't change it
    priors <- list()
    priors$R <- 1 * diag(2)
    priors$k <- 2
    priors$tau.mu <- 1.0E-3
    
    
    # fit the ellipses which uses an Inverse Wishart prior
    # on the covariance matrix Sigma, and a vague normal prior on the 
    # means. Fitting is via the JAGS method.
    ellipses.posterior_df <- siberMVN(df_siberized, parms, priors)
    
    #"1.2" mean community 1 group 2 
    ellipse1 <- "1.1" # Remember group 1 -> Host tissue in our dataframe 
    ellipse2 <- "1.2" # Remember group 2 -> Symbionts tissue in our dataframe
    
    
    # and the corresponding Bayesian estimates for the overlap between the 
    # 40% ellipses sizes is given by:
    bayes40.overlap_df <- bayesianOverlap(ellipse1, ellipse2, ellipses.posterior_df,
                                          draws = 50, p.interval = 0.40, n = 100)
    # with 
    # draws = number of ellipses taken in account (ellipses from posterior fit)
    # p.interval = percentage of points used by SIBER used to build the SEAc
    # n = number of point to create a smooth ellipse (mainly graphical)
    
    
    # get the basic mode of symb area / host area / overlap area
    # mode host
    Master_list[[i]][["SEAH"]][it] <- Mode(bayes40.overlap_df[,1])
    
    # mode symb
    Master_list[[i]][["SEAS"]][it] <- Mode(bayes40.overlap_df[,2])
    
    # mode overlap
    Master_list[[i]][["SEAO"]][it] <- Mode(bayes40.overlap_df[,3])
    
    
    # Proportion Overlap : compare the area of overlap to the area of Symb
    # Get the mode of our percentage overlap over the total area of ellipses
    Master_list[[i]][["SEABS"]][it]<- Mode(bayes40.overlap_df[,3] / bayes40.overlap_df[,2])
    
    # Proportion Overlap : compare the area of overlap to the area of Host
    # Get the mode of our percentage overlap over the total area of ellipses
    Master_list[[i]][["SEABH"]][it] <- Mode(bayes40.overlap_df[,3] / bayes40.overlap_df[,1])
    
    
    
    
    # Same but for 95% ellipses size
    bayes95.overlap_df <- bayesianOverlap(ellipse1, ellipse2, ellipses.posterior_df,
                                          draws = 50, p.interval = 0.95, n = 100)
    
    # get the basic mode of symb area / host area / overlap area
    # mode host
    Master_list[[i]][["MEAH"]][it] <- Mode(bayes95.overlap_df[,1])
    
    # mode symb
    Master_list[[i]][["MEAS"]][it] <- Mode(bayes95.overlap_df[,2])
    
    # mode overlap
    Master_list[[i]][["MEAO"]][it] <- Mode(bayes95.overlap_df[,3])
    
    # Proportion Overlap : compare the area of overlap to the area of Symb
    # Get the mode of our percentage overlap over the total area of ellipses
    Master_list[[i]][["MEABS"]][it] <- Mode(bayes95.overlap_df[,3] / bayes95.overlap_df[,2]) 
    
    #Proportion Overlap : compare the area of overlap to the area of HOST
    # Get the mode of our percentage overlap over the total area of ellipses
    Master_list[[i]][["MEABH"]][it] <- Mode(bayes95.overlap_df[,3] / bayes95.overlap_df[,1])
    
    
    
    # Compute HERS score
    
    # MEAB => Score Major Ellipse Area (95% ellipses size)
    # Computed from:
    # MEABH => Major Ellipse Area overlap / Hosts Major Ellipse Area 
    # MEABS => Major Ellipse Area overlap / Symbiont Major Ellipse Area 
    Master_list[[i]][["MEAB"]][it] <- Master_list[[i]][["MEABH"]][it]^exp(-Master_list[[i]][["MEABS"]][it])
    
    # SEAB => Score Standard Ellipse Area (40% ellipses size)
    # Computed from:
    # SEABH => Standard Ellipse Area overlap / Hosts Standard Ellipse Area
    # SEABS => Standard Ellipse Area overlap / Symbiont Standard Ellipse Area
    Master_list[[i]][["SEAB"]][it] <- Master_list[[i]][["SEABH"]][it]^exp(-Master_list[[i]][["SEABS"]][it])
    
    
    # Zscore => HERS score
    # Computed fom the two previous parts
    # SEAB => Score Standard Ellipse Area (40% ellipses size)
    # MEAB => Score Major Ellipse Area (95% ellipses size)
    Master_list[[i]][["HERS_score"]][it] <- (Master_list[[i]][["MEAB"]][it] + Master_list[[i]][["SEAB"]][it]) / 2
    
  }
  
  # Extract HERS score - Mode
  HERSScoreMode <- round(mean(Master_list[[i]][["HERS_score"]]), digit = 2)
  
  
  # Plot
  community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
  group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
  group.hull.args      <- list(lty = 2, col = "grey20")
  palette(c("blue4", "chartreuse4"))
  
  # Parameters
  par(mar=c(t = 3,r = 4, b = 1, l =1),
      oma=c(t = 1,r = 1, b = 1, l = 1),
      mgp=c(2,0.7,0))
  
  # Plot ellipses
  plotSiberObject(df_siberized,
                  ax.pad = 0.5, 
                  hulls = F, community.hulls.args, 
                  ellipses = T, group.ellipses.args,
                  group.hulls = T, group.hull.args,
                  bty = "L",
                  iso.order = c(1,2),
                  ylab = expression({delta}^15*N~'\u2030'), 
                  xlab = expression({delta}^13*C~'(\u2030)'),
                  x.limits = c(min(df_siber$iso1)-4,max(df_siber$iso1)+4), 
                  y.limits = c(min(df_siber$iso2)-2,max(df_siber$iso2)+2), cex=2)
  
  # replot some ellipses to get some shade
  car::ellipse(center = c(df_siberized[["ML.mu"]][["1"]][1,1,1],df_siberized[["ML.mu"]][["1"]][1,2,1]), # Centroid
               shape = df_siberized[["ML.cov"]][[1]][,,1], # Covariance matrix
               radius = sqrt(qchisq(0.95, df=2)), # Ellipse size
               draw=T, col = "blue4", center.pch = NULL, lwd= 0.8, segment = 100,fill=TRUE,fill.alpha=0.1, lty = 0) # Aesthetic parameters
  
  car::ellipse(center = c(df_siberized[["ML.mu"]][["1"]][1,1,2],df_siberized[["ML.mu"]][["1"]][1,2,2]), # Centroid
               shape = df_siberized[["ML.cov"]][[1]][,,2], # Covariance matrix
               radius = sqrt(qchisq(.95, df=2)), # Ellipse size
               draw=T, col = "chartreuse4", center.pch = NULL, lwd= 0.8, segment = 100,fill=TRUE, fill.alpha=0.1, lty = 0) # Aesthetic parameters
  
  
  plotGroupEllipses(df_siberized, n = 100, p.interval = 0.4,
                    lty = 2, lwd = 2)
  
  # replot some ellipses to get even more shade
  car::ellipse(center = c(df_siberized[["ML.mu"]][["1"]][1,1,1],df_siberized[["ML.mu"]][["1"]][1,2,1]), # Centroid
               shape = df_siberized[["ML.cov"]][[1]][,,1], # Covariance matrix
               radius = sqrt(qchisq(0.40, df=2)), # Ellipse size
               draw=T, col = "blue4", center.pch = NULL, lwd= 0.8, segment = 100,fill=TRUE,fill.alpha=0.2, lty = 0) # Aesthetic parameters
  
  car::ellipse(center = c(df_siberized[["ML.mu"]][["1"]][1,1,2],df_siberized[["ML.mu"]][["1"]][1,2,2]), # Centroid
               shape = df_siberized[["ML.cov"]][[1]][,,2], # Covariance matrix
               radius = sqrt(qchisq(.40, df=2)), # Ellipse size
               draw=T, col = "chartreuse4", center.pch = NULL, lwd= 0.8, segment = 100,fill=TRUE, fill.alpha=0.2, lty = 0) # Aesthetic parameters
  
  # add the legend
  legend(min(df_siber$iso1)-3.5,
         max(df_siber$iso2)+1.5,
         horiz = F ,legend=c(expression(paste(italic("Host"))),
                             expression(paste(italic("Symbionts")))),
         bty='n', pch=19,cex = 1.5,
         col= c("blue4", "chartreuse4"))
  
  title(paste0(Master_list[[i]][["Species"]][1], " - HERS score:", HERSScoreMode))
  
}

# Have a look to our raw results
Master_list

NewFile <- createWorkbook()

for (i in 1:length(Master_list)) {
  
  OneSpecies <- as.data.frame(Master_list[[i]])
  
  sheet <- createSheet(NewFile,paste(OneSpecies$Species[1]))
  
  addDataFrame(OneSpecies,
               sheet = sheet,
               startRow = 1,
               row.names = FALSE)
  
}

saveWorkbook(NewFile, file = "GiantClamsHERSResults.xlsx")
