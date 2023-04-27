library(readxl)
library(writexl)
library(dplyr) # Df tools 
library(ggplot2)
library(patchwork) # For adding plots together
library(viridis) # Colormaps
library(colorspace) # Color manipulation (e.g., desaturation colors)
library(comprehenr) # For python comprehensions

df <- read.csv("data/SEER RPD 17 Nov 2021_firsts.csv")
age_bins <- sort(unique(df$`Age.recode.with..1.year.olds`)) # Get ordered labels
df$`Age.recode.with..1.year.olds` <- factor(df$`Age.recode.with..1.year.olds`, order=TRUE, levels=age_bins) # Assign labels to ordering to make ordinal (as.factor() doesn't have this option, is just a wrapper around factor() hence use the latter instead)
df$`Age` <- as.integer(df$`Age.recode.with..1.year.olds`) # Convert to integer since R seems to use dummy variable coding even for ordinal 



df$`GBM` <- as.logical(as.factor(df$`Non.first.GBM`))
df$`Sex` <- as.factor(df$`Sex`)
df$`Site` <- as.factor(df$`Site.recode.ICD.O.3.WHO.2008`)
df$`Type` <- as.factor(df$`ICD.O.3.Hist.behav`)
df$`Radiation` <- as.factor(df$`Radiation.recode`)
df$`Chemotherapy` <- as.factor(df$`Chemotherapy.recode..yes..no.unk.`)




lm_base <- glm(`GBM` ~
                    `Age` +
                    `Sex` +
                    `Survival.months` +
                    `Radiation` +
                    `Chemotherapy`,
              family="binomial",
              data=df
                    )
summary(lm_base)

lm_site <- glm(`GBM` ~
                    `Age` +
                    `Sex` +
                    `Survival.months` +
                    `Site` +
                    `Radiation` +
                    `Chemotherapy`,
              family="binomial",
              data=df
                    )
summary(lm_site)

type_dist <- tabulate(df$Type)
most_common100 <- tail(sort(type_dist), 100)[1] # Get 100th most common value
thresh <- most_common100/sum(type_dist)

condenseFactor <- function(vector, threshold = 0.02, newName = "Other") {
  toCondense <- names(which(prop.table(table(vector)) <= threshold))
  vector[vector %in% toCondense] <- newName
  vector
}

condenseFactor2 <- function(vector, thresh_count, newName = "Other") { # Use count instead of proportion
  toCondense <- names(which(table(vector) <= thresh_count))
  vector[vector %in% toCondense] <- newName
  vector
}

df$Type_comp <- as.factor(condenseFactor(df$`ICD.O.3.Hist.behav`, thresh))


lm_type <- glm(`GBM` ~
                    `Age` +
                    `Sex` +
                    `Survival.months` +
                    `Type_comp` +
                    `Radiation` +
                    `Chemotherapy`,
              family="binomial",
              data=df
                    )

summary(lm_type)

glm(`GBM` ~ 1, family="binomial", data=df)
forward_type <- step(glm(`GBM` ~ 1, family="binomial", data=df), direction="forward", scope=formula(lm_type))
summary(forward_type)
