
#### Meta-analysis 2022 ####


#Set WD
  
setwd("~/PhD work/meta-analysis/Code")

dev.off()

### Load packages


library(metafor)
library(tidyverse)
library(orchaRd)
library(janitor)
library(rotl)
library(ape)
library(patchwork)
library(ggtext)
library(knitr)
library(kableExtra)

## Load data

ef <- read.csv("effect_sizes.csv") 

## calculate absolute ef

yi_ab <- abs(ef$yi)
ef_ab <- ef
ef_ab$yi <- yi_ab  #converts effect sizes to absolute values

#################################################

# construct phylogentic tree matrix for use as random factor

species <- unique(ef$species) # list of unique species in meta-analysis
species <- as.character(species) # change to character object
taxa <- tnrs_match_names(species)
tree <- rotl::tol_induced_subtree(taxa$ott_id)
tree$tip.label <-
  strip_ott_ids(tree$tip.label, remove_underscores = TRUE) # change ids to the names from dataset

# calculate correlations between all species = cor

tree2 <- compute.brlen(tree)
cor <- vcv(tree2, cor = T)


##############################################################################
    
#                     M E T A - A N A L Y S E S                              #

##############################################################################


# Meta-analyses for 8 leaf traits
# random factors = ACC (study), experiment (diversity study site), species  R = list(species = cor) (phylogeny)
# and individual effect size (ef)
# both standard and absolute meta-analyses conducted 
# added 'control = list(optimizer = "optim")' argument to all analysis to prevent 
# "Optimizer (nlminb) did not achieve convergence (convergence = 1)" error


#--------------------------------#
       ### thickness ###
#--------------------------------#


# thickness

thick<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                trait_category=="thickness"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), control = list(optimizer = "optim"))

summary(thick)


# thickness absolute

thick_ab<- rma.mv(yi,vi,data=ef_ab, subset=(structure == "leaf" & diversity_type == "species" &
                         trait_category=="thickness"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"))

summary(thick_ab)


#--------------------------------#
       ### toughness ###
#--------------------------------#

# toughness

tough<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="toughness"),method= "REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"))

summary(tough)



# toughness absolute

tough_ab<- rma.mv(yi,vi,data=ef_ab, subset=(structure == "leaf" & diversity_type == "species" &
                                              trait_category=="toughness"),method="REML",
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                  R = list(species = cor), control = list(optimizer = "optim"))

summary(tough_ab)


#--------------------------------#
        ### terpenes ###
#--------------------------------#

# terpenes 

terp<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                     trait_category=="terpenes"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"))

summary(terp)

# terpenes absolute

terp_ab<- rma.mv(yi,vi,data=ef_ab, subset=(structure == "leaf" & diversity_type == "species" &
                                              trait_category=="terpenes"),method="REML",
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                  R = list(species = cor), control = list(optimizer = "optim"))

summary(terp_ab)



terp_o
terp_ab_o



#--------------------------------#
           ### LDMC ###
#--------------------------------#

# LDMC

LDMC<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                   trait_category=="LDMC"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"))

summary(LDMC)


# LDMC absolute

LDMC_ab<- rma.mv(yi,vi,data=ef_ab, subset=(structure == "leaf" & diversity_type == "species" &
                                             trait_category=="LDMC"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), control = list(optimizer = "optim"))

summary(LDMC_ab)


#--------------------------------#
           ### SLA ###
#--------------------------------#


# SLA

SLA<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                 trait_category=="SLA"),method="REML",
             random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"))

summary(SLA)


# SLA absolute

SLA_ab<- rma.mv(yi,vi,data=ef_ab, subset=(structure == "leaf" & diversity_type == "species" &
                                             trait_category=="SLA"),method="REML",
                random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), control = list(optimizer = "optim"))

summary(SLA_ab)


#--------------------------------#
           ### C ###
#--------------------------------#

# C

C<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="C"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"))

summary(C)


# C absolute

C_ab<- rma.mv(yi,vi,data=ef_ab, subset=(structure == "leaf" & diversity_type == "species" &
                                             trait_category=="C"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), control = list(optimizer = "optim"))

summary(C_ab)



#--------------------------------#
            ### N ###
#--------------------------------#

# N

N<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="N"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"))

summary(N)



# N absolute

N_ab<- rma.mv(yi,vi,data=ef_ab, subset=(structure == "leaf" & diversity_type == "species" &
                                             trait_category=="N"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), control = list(optimizer = "optim"))

summary(N_ab)


#--------------------------------#
        ### phenolics ###
#--------------------------------#

# phenolics

phenolics<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                               trait_category=="phenolics"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"))


summary(phenolics)


phenolics_o


phenolics_o

# phenolics absolute

phenolics_ab<- rma.mv(yi,vi,data=ef_ab, subset=(structure == "leaf" & diversity_type == "species" &
                                             trait_category=="phenolics"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), control = list(optimizer = "optim"))

summary(phenolics_ab)


# phenolics by type



phenolics_type <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                       trait_category=="phenolics" ), method="REML",
                         random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                         R = list(species = cor), control = list(optimizer = "optim"),
                         mods = ~ defence_type_specific-1)

orchard_plot(phenolics_type, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "defence_type_specific")


### prediction intervals ##


predict(tough)
predict(thick)
predict(LDMC)
predict(SLA)
predict(terp)
predict(phenolics)
predict(N)
predict(C)

predict(tough_ab)
predict(thick_ab)
predict(LDMC_ab)
predict(SLA_ab)
predict(terp_ab)
predict(phenolics_ab)
predict(N_ab)
predict(C_ab)

### orchard plots ###

# thickness

thick_o <- orchard_plot(thick, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#E69F00") + scale_colour_manual(values = "#E69F00") +
  geom_richtext(x = 3,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 20(3)", size = 6) +
  geom_richtext(x = -4, y = 1.3,label.color = NA, label = "Thickness", size = 6) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# toughness

tough_o <- orchard_plot(tough, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#56B4E9") + scale_colour_manual(values = "#56B4E9") +
  geom_richtext(x = 3,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 20(3)", size = 6) +
  geom_richtext(x = -4, y = 1.3,label.color = NA, label = "Toughness", size = 6) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# thickness

terp_o <- orchard_plot(terp, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#CC79A7") + scale_colour_manual(values = "#CC79A7") +
  geom_richtext(x = 3,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 24(6)", size = 6) +
  geom_richtext(x = -4, y = 1.3,label.color = NA, label = "Terpenoids", size = 6) +
    xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# LDMC

LDMC_o <- orchard_plot(LDMC, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#F0E442") + scale_colour_manual(values = "#F0E442") +
  geom_richtext(x = 3,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 119(9)", size = 6) +
  geom_richtext(x = -4, y = 1.3,label.color = NA, label = "LDMC", size = 6) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# SLA

SLA_o <- orchard_plot(SLA, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#0072B2") + scale_colour_manual(values = "#0072B2") +
  geom_richtext(x = 3,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 251(17)", size = 6) +
  geom_richtext(x = -4, y = 1.3,label.color = NA, label = "SLA", size = 6) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# C

C_o <- orchard_plot(C, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#999999" ) + scale_colour_manual(values = "#999999" ) +
  geom_richtext(x = 3,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 138(11)", size = 6) +
  geom_richtext(x = -4, y = 1.3,label.color = NA, label = "Carbon", size = 6) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
# N

N_o <- orchard_plot(N, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#009E73" ) + scale_colour_manual(values = "#009E73" ) +
  geom_richtext(x = 3,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 206(27)", size = 6) +
  geom_richtext(x = -4, y = 1.3,label.color = NA, label = "Nitrogen", size = 6) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# phenolics

phenolics_o <- orchard_plot(phenolics, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#D55E00" ) + scale_colour_manual(values = "#D55E00") +
  geom_richtext(x = 3,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 228(15)", size = 6) +
  geom_richtext(x = -4, y = 1.3,label.color = NA, label = "Phenolics", size = 6) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


###  Absolute orchard plots ###

# thickness absolute

thick_ab_o <- orchard_plot(thick_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#E69F00") + scale_colour_manual(values = "#E69F00") +
  geom_richtext(x = 4,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 20(3)", size = 6) +
  geom_richtext(x = -1, y = 1.3,label.color = NA, label = "Thickness", size = 6) +
  xlim(-2, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# toughness absolute

tough_ab_o <- orchard_plot(tough_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#56B4E9") + scale_colour_manual(values = "#56B4E9") +
  geom_richtext(x = 4,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 20(3)", size = 6) +
  geom_richtext(x = -1, y = 1.3,label.color = NA, label = "Toughness", size = 6) +
  xlim(-2, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# terpenes absolute

terp_ab_o <- orchard_plot(terp_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45) +
  scale_fill_manual(values = "#CC79A7") + scale_colour_manual(values = "#CC79A7") +
  geom_richtext(x = 4,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 24(6)", size = 6) +
  geom_richtext(x = -1, y = 1.3,label.color = NA, label = "Terpenoids", size = 6) +
  xlim(-2, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# LDMC absolute

LDMC_ab_o <- orchard_plot(LDMC_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45) +
  scale_fill_manual(values = "#F0E442") + scale_colour_manual(values = "#F0E442") +
  geom_richtext(x = 4,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 119(9)", size = 6) +
  geom_richtext(x = -1, y = 1.3,label.color = NA, label = "LDMC", size = 6) +
  xlim(-2, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


# SLA absolute

SLA_ab_o <- orchard_plot(SLA_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45) +
  scale_fill_manual(values = "#0072B2") + scale_colour_manual(values = "#0072B2") +
  geom_richtext(x = 4,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 251(17)", size = 6) +
  geom_richtext(x = -1, y = 1.3,label.color = NA, label = "SLA", size = 6) +
  xlim(-2, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# C absolute

C_ab_o <- orchard_plot(C_ab, xlab = "Absolute standardised mean difference", transfm = "none", angle = 45) +
  scale_fill_manual(values = "#999999" ) + scale_colour_manual(values = "#999999" ) +
  geom_richtext(x = 4,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 138(11)", size = 6) +
  geom_richtext(x = -1, y = 1.3,label.color = NA, label = "Carbon", size = 6) +
  xlim(-2, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# N absolute

N_ab_o <- orchard_plot(N_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45) +
  scale_fill_manual(values = "#009E73" ) + scale_colour_manual(values = "#009E73" ) +
  geom_richtext(x = 4,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 206(27)", size = 6) +
  geom_richtext(x = -1, y = 1.3,label.color = NA, label = "Nitrogen", size = 6) +
  xlim(-2, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# phenolics absolute 

phenolics_ab_o <- orchard_plot(phenolics_ab, xlab = "Standardised mean difference", transfm = "none", angle = 45, k = FALSE) +
  scale_fill_manual(values = "#D55E00" ) + scale_colour_manual(values = "#D55E00") +
  geom_richtext(x = 4,  y = 1.3, label.color = NA, label = "<i>k</i>(N) = 228(15)", size = 6) +
  geom_richtext(x = -1, y = 1.3,label.color = NA, label = "Phenolics", size = 6) +
  xlim(-2, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


##############################################################################

#                      M E T A - R E G R E S S I O N S                      #

##############################################################################


# meta regressions only performed for SLA, LDMC, C, N and phenolics
# terpenes, thickness and toughness had too few effect sizes to be considered
# SR, PD are continous variables
# experiment type, plant age, N-fixing neighbours = categorical mods



### C A T E G O R I C A L  M O D E R A T O R S ###


#--------------------------------#
###    N-fixing neighbours     ###
#--------------------------------#

options(na.action = "na.omit")

# SLA

SLA_N <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                             trait_category=="SLA"), method="REML",
                 mods = ~ Nfixing-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), control = list(optimizer = "optim"))

summary(SLA_N)


SLA_N_o <- orchard_plot(SLA_N, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "Nfixing") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("SLA"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
SLA_N_o

SLA_N_R2 <- r2_ml(SLA_N)
100*SLA_N_R2[1]


# Phenolics

phenolics_N <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="phenolics"),method="REML",
                 mods = ~ Nfixing-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), control = list(optimizer = "optim"))

summary(phenolics_N)

phenolics_N_o <- orchard_plot(phenolics_N, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "Nfixing") +
  annotate(geom = "text", x = -4, y = 2.3, label = paste0("Phenolics"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

phenolics_N_o 

phenolics_N_R2 <- R2 <- r2_ml(phenolics_N)
100*phenolics_N_R2[1]


# LDMC

LDMC_N <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                                trait_category=="LDMC"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                       R = list(species = cor), control = list(optimizer = "optim"))

summary(LDMC_N)

LDMC_N_o <- orchard_plot(LDMC_N, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "Nfixing") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("LDMC"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

LDMC_N_o 

LDMC_N_R2 <- R2 <- r2_ml(LDMC_N)
100*LDMC_N_R2[1]


#Nitrogen (with phylogentic matrix argument)

N_N <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                          trait_category=="N"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                       R = list(species = cor), control = list(optimizer = "optim"))

summary(N_N)


N_N_o <- orchard_plot(N_N, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "Nfixing") +
  annotate(geom = "text", x = -4.5, y = 2.4, label = paste0("A"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

N_N_o
  
# Nitrogen (without phylogenetic argument)  
  
N_N2 <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                          trait_category=="N"),method="REML",
                       mods = ~ Nfixing-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                        control = list(optimizer = "optim"))

summary(N_N2)

N_N2_o <- orchard_plot(N_N2, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "Nfixing") +
    annotate(geom = "text", x = -4.5, y = 2.4, label = paste0("B"), color = "black", parse = TRUE, size = 8,) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

N_N2_o  

N_N_R2 <-r2_ml(N_N)
100*N_N_R2[1]

# carbon

C_N <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                                trait_category=="C"),method="REML",
                       mods = ~ broadleaf_conifer_mix-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                       R = list(species = cor), control = list(optimizer = "optim"))

summary(C_N)

C_N_o <- orchard_plot(C_N, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "Nfixing") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("C"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

C_N_o

C_N_R2 <-r2_ml(C_N)
100*C_N_R2[1]


#--------------------------------#
###         Study type         ###
#--------------------------------#


#SLA

SLA_ST<- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="SLA"),method="REML",
                 mods = ~ study_type-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), control = list(optimizer = "optim"))

summary(SLA_ST)

SLA_ST_o <- orchard_plot(SLA_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "study_type") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("SLA"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5)  +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

SLA_ST_o

SLA_ST_R2 <- r2_ml(SLA_ST)
100*SLA_ST_R2[1]

#phenolics

phenolics_ST <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                                trait_category=="phenolics"),method="REML",
                       mods = ~ study_type-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                       R = list(species = cor), control = list(optimizer = "optim"))

summary(phenolics_ST)

phenolics_ST_R2 <- r2_ml(phenolics_ST)
100*phenolics_ST_R2[1]

phenolics_ST_o <- orchard_plot(phenolics_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "study_type") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("phenolics"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

phenolics_ST_o

# LDMC

LDMC_ST <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="LDMC"),method="REML",
                  mods = ~ study_type-1,
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                  R = list(species = cor), control = list(optimizer = "optim"))

summary(LDMC_ST)

LDMC_ST_o <- orchard_plot(LDMC_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "study_type") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("LDMC"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

LDMC_ST_o 

LDMC_ST_R2 <- r2_ml(LDMC_ST)
100*LDMC_ST_R2[1]


#Nitrogen

N_ST <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="N"),method="REML",
               mods = ~ study_type-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"))

summary(N_ST)

N_ST_o <- orchard_plot(N_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "study_type") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("N"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

N_ST_o 

N_ST_R2 <- r2_ml(N_ST)
100*N_ST_R2[1]


#carbon

C_ST <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="C"),method="REML",
               mods = ~ study_type-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"))

summary(C_ST)

C_ST_o <- orchard_plot(C_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "study_type") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("C"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

C_ST_o 

C_ST_R2 <- r2_ml(C_ST)
100*C_ST_R2[1]


#--------------------------------#
###          plant age         ###
#--------------------------------#

#SLA

SLA_LS <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                 trait_category=="SLA" & study_type == "experimental"),method="REML",
                 mods = ~ plant_lifestage_2-1,
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), control = list(optimizer = "optim"))

summary(SLA_LS)

SLA_LS_o <- orchard_plot(SLA_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "plant_lifestage_2") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("SLA"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

SLA_LS_o

SLA_LS_R2 <- r2_ml(SLA_LS)
100*SLA_LS_R2[1]

#phenolics

phenolics_LS <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                   trait_category=="phenolics"),method="REML",
                       mods = ~ plant_lifestage_2-1,
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                       R = list(species = cor), control = list(optimizer = "optim"))

summary(phenolics_LS)

phenolics_LS_o <- orchard_plot(phenolics_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "plant_lifestage_2") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("phenolics"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

phenolics_LS_o

phenolics_LS_R2 <- r2_ml(phenolics_LS)
100*phenolics_LS_R2[1]


phenolics_LS_o 

# LDMC

LDMC_LS <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="LDMC"),method="REML",
                  mods = ~ plant_lifestage_2-1,
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                  R = list(species = cor), control = list(optimizer = "optim"))

summary(LDMC_LS)

LDMC_LS_o <- orchard_plot(LDMC_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "plant_lifestage_2") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("LDMC"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

LDMC_LS_o  

LDMC_LS_R2 <- r2_ml(LDMC_LS)
100*LDMC_LS_R2[1]

# Nitrogen

N_LS <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="N"),method="REML",
               mods = ~ plant_lifestage_2-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"))

summary(N_LS)

N_LS_o <- orchard_plot(N_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "plant_lifestage_2") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("N"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

N_LS_o

N_LS_R2 <- r2_ml(N_LS)
100*N_LS_R2[1]

#carbon

C_LS <- rma.mv (yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="C"),method="REML",
               mods = ~ plant_lifestage_2-1,
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"))

summary(C_LS)

C_LS_o <- orchard_plot(C_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0, mod = "plant_lifestage_2") +
  annotate(geom = "text", x = -4, y = 1.5, label = paste0("C"), color = "black", parse = TRUE, size = 8) +
  xlim(-5, 5) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

C_LS_o 

C_LS_R2 <- r2_ml(C_LS)
100*C_LS_R2[1]


### C O N T I N I O U S   M O D E R A T O R S ###

#--------------------------------#
###          density           ###
#--------------------------------#

#SLA


SLA_density <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="SLA" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), mods = ~ density, control = list(optimizer = "optim"))

summary(SLA_density)

SLA_density_r <- regplot(SLA_density, mod="density", pi=TRUE, refline=0, 
                    xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "SLA", ylim =c(-5,5))

SLA_density_R2 <- r2_ml(SLA_density)
100*SLA_density_R2[1]


# phenolics

phenolics_density <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                               trait_category=="phenolics"),method="REML",
                      random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                      R = list(species = cor), mods = ~ density,  control = list(optimizer = "optim"))

summary(phenolics_density)

phenolics_density_r <- regplot(phenolics_density, mod="density", pi=TRUE, refline=0, 
                         xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "phenolics", ylim =c(-5,5))

phenolics_density_R2 <- r2_ml(phenolics_density)
100*phenolics_density_R2[1]


#LDMC

LDMC_density <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                                     trait_category=="LDMC"),method="REML",
                            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                            R = list(species = cor), mods = ~ density, control = list(optimizer = "optim"))
summary(LDMC_density)

LDMC_density_r <- regplot(LDMC_density, mod="density", pi=TRUE, refline=0, 
                               xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "LDMC", ylim =c(-5,5))


LDMC_density_R2 <- r2_ml(LDMC_density)
100*LDMC_density_R2[1]

# N

### excluded experiments with much higher density -- result not changed

N_density <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                                     trait_category=="N" & experiment != "IDENT Cloquet"),method="REML",
                            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                            R = list(species = cor), mods = ~ density, control = list(optimizer = "optim"))

summary(N_density)

N_density_r <- regplot(N_density, mod="density", pi=TRUE, refline=0, 
                               xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "N", ylim =c(-5,5))

N_density_R2 <- r2_ml(N_density)
100*N_density_R2[1]



# C

C_density <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                                     trait_category=="C" & experiment != "IDENT Cloquet"),method="REML",
                            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                            R = list(species = cor), mods = ~ density, control = list(optimizer = "optim"))

C_density_r <- regplot(C_density, mod="density", pi=TRUE, refline=0, 
                               xlab="trees per m2", ylab="SMD", bg= "#009e73", main = "C", ylim =c(-5,5))

summary(C_density)

C_density_R2 <- r2_ml(C_density)
100*C_density_R2[1]

#--------------------------------#
###      species richness      ###
#--------------------------------#

# SLA

SLA_SR <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                      trait_category=="SLA"),method="REML",
                   random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                   R = list(species = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(SLA_SR)

SLA_SR_r <- regplot(SLA_SR, mod="SR", pi=TRUE, refline=0, 
                    xlab="Species richness", ylab="Hedges' d", bg= "#009e73", main = "SLA", ylim =c(-5,5))

SLA_SR_R2 <- r2_ml(SLA_SR)
100*SLA_SR_R2[1]


# phenolics

phenolics_SR <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="phenolics"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(phenolics_SR)

phenolics_SR_r <- regplot(phenolics_SR, mod="SR", pi=TRUE, refline=0, 
                    xlab="Species richness", ylab="SMD", bg= "#009e73", main = "phenolics", ylim =c(-5,5))

phenolics_SR_R2 <- r2_ml(phenolics_SR)
100*phenolics_SR_R2[1]


#LDMC

LDMC_SR <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                                trait_category=="LDMC"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                       R = list(species = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(LDMC_SR)

LDMC_SR_r <- regplot(LDMC_SR, mod="SR", pi=TRUE, refline=0, 
                          xlab="Species richness", ylab="SMD", bg= "#009e73", main = "LDMC", ylim =c(-5,5))

LDMC_SR_R2 <- r2_ml(LDMC_SR)
100*LDMC_SR_R2[1]

#N

N_SR <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                                trait_category=="N"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                       R = list(species = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(N_SR)

N_SR_r <- regplot(N_SR, mod="SR", pi=TRUE, refline=0, 
                          xlab="Species richness", ylab="SMD", bg= "#009e73", main = "N", ylim =c(-5,5))

N_SR_R2 <- r2_ml(N_SR)
100*N_SR_R2[1]

#C

C_SR <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                                trait_category=="C"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                       R = list(species = cor), mods = ~ SR, control = list(optimizer = "optim"))

summary(C_SR)

C_SR_r <- regplot(C_SR, mod="SR", pi=TRUE, refline=0, 
                          xlab="Species richness", ylab="SMD", bg= "#009e73", main = "C", ylim =c(-5,5))

C_SR_R2 <- r2_ml(C_SR)
100*C_SR_R2[1]

#--------------------------------#
###   phylogentic diversity    ###
#--------------------------------#

#SLA

SLA_PD2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="SLA"),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                 R = list(species = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(SLA_PD2)

SLA_PD2_r <- regplot(SLA_PD2, mod="PD2", pi=TRUE, refline=0, 
                    xlab="Phylognetic diversity", ylab="SMD", bg= "#009e73", main = "SLA", ylim =c(-5,5))

SLA_PD2_R2 <- r2_ml(SLA_PD2)
100*SLA_PD2_R2[1]

#phenolics

phenolics_PD2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                                trait_category=="phenolics"),method="REML",
                       random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                       R = list(species = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(SLA_PD2)

phenolics_PD2_r <- regplot(phenolics_PD2, mod="PD2", pi=TRUE, refline=0, 
                          xlab="Phylognetic diversity", ylab="SMD", bg= "#009e73", main = "phenolics", ylim =c(-5,5))


phenolics_PD2_R2 <- r2_ml(phenolics_PD2)
100*phenolics_PD2_R2[1]

#LDMC

LDMC_PD2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="LDMC"),method="REML",
                  random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
                  R = list(species = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(LDMC_PD2)

LDMC_PD2_r <- regplot(LDMC_PD2, mod="PD2", pi=TRUE, refline=0, 
                     xlab="Phylognetic diversity", ylab="SMD", bg= "#009e73", main = "LDMC", ylim =c(-5,5))

LDMC_PD2_R2 <- r2_ml(LDMC_PD2)
100*LDMC_PD2_R2[1]


# N

N_PD2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
             trait_category=="N" ) ,method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(N_PD2)

N_PD2_r <- regplot(N_PD2, mod="PD2", pi=TRUE, refline=0, 
                  xlab="Phylognetic diversity", ylab="SMD", bg= "#009e73", main = "N", ylim =c(-5,5))

N_PD2_R2 <- r2_ml(N_PD2)
100*N_PD2_R2[1]  

# C

C_PD2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="C"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), mods = ~ PD2, control = list(optimizer = "optim"))

summary(C_PD2)

C_PD2_r <- regplot(C_PD2, mod="PD2", pi=TRUE, refline=0, 
                  xlab="Phylognetic diversity", ylab="SMD", bg= "#009e73", main = "C", ylim =c(-5,5))


C_PD2_R2 <- r2_ml(C_PD2)
100*C_PD2_R2[1]




##############################################################################

#                             SENSITIVITY ANALYSIS                           #
 
##############################################################################


# What species are over-represented ?
  
# C

 # 15.2 % Betula pendula
 # 9.4 % Fagus sylvetcia
 # 5.1 % Quercus robur
 # All other species < 5 %

# LDMC

 # 10.1 % Betula pendula
 # 8.4 % Fagus sylvetica
 # 5 % Quercus robur
 # All other species < 5 % 

# N

 # 10.1 % Betula pendula
 # 6.2 % Fagus sylvetica
 # All other species < 5 % 


# phenolics

 # 67.1 % Betula pendula
 # 7.1 % Fagus slyvetica

# SLA

 # 8.7 % Fagus sylvetica
 # 7.1 % Betula pendula


# Terpenes 

 # 24 % Plantago lancelota
 # 24 % Pinus halpensis
 # 24 % Rosmarinus officentalis
 # 12 % Cistus albidus


# Thickness

 # 80 % Betula Pendula
 # 15 % Quercus robur
 # 5 % Pinus pinaster 

# Toughness

 # 80 % Betula pendula
 # 15 % Quercus robur
 # 5 % Shorea leprosula



##### Re-analysis with/without Betula pendula #####


# Bpy = Betula pendula yes
# Bpn = Betula pendula no

# need to remove phylogeny random factor for bp only to do this (would equal 0 anyway)


### C +/- Betula Pendula ###



C_Bpy<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
            trait_category=="C" & species == "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"))
summary(C_Bpy)

orchard_plot(C_Bpy, xlab = "Standardised mean difference", transfm = "none", angle = 0)


C_Bpn<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="C" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"))


summary(C_Bpn)

orchard_plot(C_Bpn, xlab = "Standardised mean difference", transfm = "none", angle = 0)


## meta-regressions ##

## Phylogenetic diversity ##


C_Bpy_PD2<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="C" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ PD2)


summary(C_Bpy_PD2)

C_Bpy_PD2_r <- regplot(C_Bpy_PD2, mod="PD2", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


C_Bpn_PD2<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(C_Bpn_PD2)

C_Bpn_PD2_r <- regplot(C_Bpn_PD2, mod="PD2", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


# species richness


C_Bpy_SR<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="C" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ SR)


summary(C_Bpy_SR)

C_Bpy_SR_r <- regplot(C_Bpy_SR, mod="SR", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


C_Bpn_SR<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(C_Bpn_SR)

C_Bpn_SR_r <- regplot(C_Bpn_SR, mod="SR", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


# density


C_Bpy_density<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="C" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ density)


summary(C_Bpy_density)

C_Bpy_density_r <- regplot(C_Bpy_density, mod="density", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


C_Bpn_density<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(C_Bpn_density)

C_Bpn_density_r <- regplot(C_Bpn_density, mod="density", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


## study type

C_Bpy_ST<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="C" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ study_type-1)


summary(C_Bpy_ST)

orchard_plot(C_Bpy_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0)


# no observational birch studies


C_Bpn_ST<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(C_Bpn_ST)

orchard_plot(C_Bpn_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0)



### Nfixing

C_Bpy_N<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="C" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ Nfixing-1)


summary(C_Bpy_N)

orchard_plot(C_Bpy_N, xlab = "Standardised mean difference", transfm = "none", angle = 0)




C_Bpn_N<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(C_Bpn_N)

orchard_plot(C_Bpn_N, xlab = "Standardised mean difference", transfm = "none", angle = 0)


# all Nfixing = birch 


## age 


C_Bpy_LS<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="C" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)


summary(C_Bpy_LS)

orchard_plot(C_Bpy_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0)


C_Bpn_LS<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="C" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(C_Bpn_LS)

orchard_plot(C_Bpn_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0)




### LDMC + Betula pendula ###


LDMC_Bpy<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
            trait_category=="LDMC" & species == "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"))

summary(LDMC_Bpy)

orchard_plot(LDMC_Bpy, xlab = "Standardised mean difference", transfm = "none", angle = 0)


LDMC_Bpn<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="LDMC" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  R = list(species = cor))

summary(LDMC_Bpn)

orchard_plot(LDMC_Bpn, xlab = "Standardised mean difference", transfm = "none", angle = 0)

## meta-regressions ##

## Phylogenetic diversity ##


LDMC_Bpy_PD2<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="LDMC" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ PD2)


summary(LDMC_Bpy_PD2)

LDMC_Bpy_PD2_r <- regplot(C_Bpy_PD2, mod="PD2", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "LDMC", ylim =c(-5,5))


LDMC_Bpn_PD2<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(LDMC_Bpn_PD2)

LDMC_Bpn_PD2_r <- regplot(C_Bpn_PD2, mod="PD2", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "LDMC", ylim =c(-5,5))


# species richness


LDMC_Bpy_SR<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="LDMC" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ SR)


summary(LDMC_Bpy_SR)

LDMC_Bpy_SR_r <- regplot(C_Bpy_SR, mod="SR", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "LDMC", ylim =c(-5,5))


LDMC_Bpn_SR<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(LDMC_Bpn_SR)

LDMC_Bpn_SR_r <- regplot(C_Bpn_SR, mod="SR", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


# density


LDMC_Bpy_density<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                   trait_category=="LDMC" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ density)


summary(LDMC_Bpy_density)


LDMC_Bpy_density_r <- regplot(LDMC_Bpy_density, mod="density", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


LDMC_Bpn_density<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(C_Bpn_density)

LDMC_Bpn_density_r <- regplot(C_Bpn_density, mod="density", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))



## study type


LDMC_Bpy_ST<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="LDMC" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ study_type-1)


summary(LDMC_Bpy_ST)

orchard_plot(LDMC_Bpy_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0)


# only 2 observational birch studies


LDMC_Bpn_ST<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(LDMC_Bpn_ST)

orchard_plot(LDMC_Bpn_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0)



### Nfixing


LDMC_Bpy_N<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="LDMC" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ Nfixing-1)


summary(LDMC_Bpy_N)

orchard_plot(LDMC_Bpy_N, xlab = "Standardised mean difference", transfm = "none", angle = 0)


LDMC_Bpn_N<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor),  control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(LDMC_Bpn_N)

orchard_plot(LDMC_Bpn_N, xlab = "Standardised mean difference", transfm = "none", angle = 0)


# all Nfixing = birch 



## age 

LDMC_Bpy_LS<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="LDMC" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)


summary(LDMC_Bpy_LS)

orchard_plot(LDMC_Bpy_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0)


LDMC_Bpn_LS<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="LDMC" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(LDMC_Bpn_LS)

orchard_plot(LDMC_Bpn_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0)



### N + Betula Pendula ###


N_Bpy<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
            trait_category=="N" & species == "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"))

summary(N_Bpy)

orchard_plot(N_Bpy, xlab = "Standardised mean difference", transfm = "none", angle = 0)



N_Bpn<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="N" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  R = list(species = cor))

summary(N_Bpn)

orchard_plot(N_Bpn, xlab = "Standardised mean difference", transfm = "none", angle = 0)


## meta-regressions ##

## Phylogenetic diversity ##


N_Bpy_PD2<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="N" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ PD2)


summary(N_Bpy_PD2)

N_Bpy_PD2_r <- regplot(N_Bpy_PD2, mod="PD2", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "N", ylim =c(-5,5))


N_Bpn_PD2<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor),  control = list(optimizer = "optim"),  mods = ~ PD2)

summary(N_Bpn_PD2)

N_Bpn_PD2_r <- regplot(N_Bpn_PD2, mod="PD2", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "N", ylim =c(-5,5))



# species richness

N_Bpy_SR<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="N" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ SR)


summary(N_Bpy_SR)

N_Bpy_SR_r <- regplot(N_Bpy_SR, mod="SR", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "N", ylim =c(-5,5))


N_Bpn_SR<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ SR)

summary(N_Bpn_SR)

N_Bpn_SR_r <- regplot(N_Bpn_SR, mod="SR", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "N", ylim =c(-5,5))



# density

N_Bpy_density<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                   trait_category=="N" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ density)


summary(N_Bpy_density)


N_Bpy_density_r <- regplot(N_Bpy_density, mod="density", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "N", ylim =c(-5,5))


N_Bpn_density<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(N_Bpn_density)

N_Bpn_density_r <- regplot(N_Bpn_density, mod="density", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


## study type

N_Bpy_ST<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="N" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ study_type-1)


summary(N_Bpy_ST)

orchard_plot(N_Bpy_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0)


# only observational birch studies


N_Bpn_ST<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor),  control = list(optimizer = "optim"),  mods = ~ study_type-1, )

summary(N_Bpn_ST)

orchard_plot(N_Bpn_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0)


### Nfixing


N_Bpy_N<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="N" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ Nfixing-1)


summary(N_Bpy_N)

orchard_plot(N_Bpy_N, xlab = "Standardised mean difference", transfm = "none", angle = 0)


N_Bpn_N<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1, R = list(species = cor))

summary(N_Bpn_N)

orchard_plot(N_Bpn_N, xlab = "Standardised mean difference", transfm = "none", angle = 0)


## age 

N_Bpy_LS<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="N" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)


summary(N_Bpy_LS)

orchard_plot(N_Bpy_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0)


# no differences


N_Bpn_LS<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="N" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1, R = list(species = cor))

summary(N_Bpn_LS)

orchard_plot(N_Bpn_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0)



### phenolics + Betula pendula ###

phenolics_Bpy<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
            trait_category=="phenolics" & species == "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"))

summary(phenolics_Bpy)

orchard_plot(phenolics_Bpy, xlab = "Standardised mean difference", transfm = "none", angle = 0)


phenolics_Bpn<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="phenolics" & species != "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"))

summary(phenolics_Bpn)

orchard_plot(phenolics_Bpn, xlab = "Standardised mean difference", transfm = "none", angle = 0)



## meta-regressions ##

## Phylogenetic diversity ##


phenolics_Bpy_PD2<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="phenolics" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ PD2)


summary(phenolics_Bpy_PD2)

phenolics_Bpy_PD2_r <- regplot(C_Bpy_PD2, mod="PD2", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "phenolics", ylim =c(-5,5))


phenolics_Bpn_PD2<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ PD2)

summary(phenolics_Bpn_PD2)

phenolics_Bpn_PD2_r <- regplot(phenolics_Bpn_PD2, mod="PD2", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "phenolics", ylim =c(-5,5))


###

# species richness


phenolics_Bpy_SR<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="phenolics" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ SR)


summary(phenolics_Bpy_SR)

phenolics_Bpy_SR_r <- regplot(phenolics_Bpy_SR, mod="SR", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "phenolics", ylim =c(-5,5))


phenolics_Bpn_SR<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
               R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ SR, R = list(species = cor))

summary(phenolics_Bpn_SR)

phenolics_Bpn_SR_r <- regplot(phenolics_Bpn_SR, mod="SR", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


# density


phenolics_Bpy_density<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                   trait_category=="phenolics" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ density)


summary(phenolics_Bpy_density)


phenolics_Bpy_density_r <- regplot(phenolics_Bpy_density, mod="density", pi=TRUE, refline=0, 
                   xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


phenolics_Bpn_density<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ density)

summary(phenolics_Bpn_density)

phenolics_Bpn_density_r <- regplot(phenolics_Bpn_density, mod="density", pi=TRUE, refline=0, 
                       xlab="Phylognetic diversity", ylab="SMD", bg= "#56b4e9", main = "C", ylim =c(-5,5))


# no differences

## study type


phenolics_Bpy_ST<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="phenolics" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ study_type-1)


summary(phenolics_Bpy_ST)

orchard_plot(phenolics_Bpy_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0)


# no observational birch studies


phenolics_Bpn_ST<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ study_type-1)

summary(phenolics_Bpn_ST)

orchard_plot(phenolics_Bpn_ST, xlab = "Standardised mean difference", transfm = "none", angle = 0)


### Nfixing


phenolics_Bpy_N<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="phenolics" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ Nfixing-1)


summary(phenolics_Bpy_N)

orchard_plot(phenolics_Bpy_N, xlab = "Standardised mean difference", transfm = "none", angle = 0)


# phenolics not reduced when only BP


phenolics_Bpn_N<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),  mods = ~ Nfixing-1)

summary(phenolics_Bpn_N)

orchard_plot(phenolics_Bpn_N, xlab = "Standardised mean difference", transfm = "none", angle = 0)


# phenolics reduced with N when BP not present



## age 


phenolics_Bpy_LS<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="phenolics" & species == "Betula pendula"),method="REML",
               random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)


summary(phenolics_Bpy_LS)

orchard_plot(phenolics_Bpy_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0)


# no differences


phenolics_Bpn_LS<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                           trait_category=="phenolics" & species != "Betula pendula" ),method="REML",
                 random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
             R = list(species = cor),  control = list(optimizer = "optim"),  mods = ~ plant_lifestage_2-1)

summary(phenolics_Bpn_LS)

orchard_plot(phenolics_Bpn_LS, xlab = "Standardised mean difference", transfm = "none", angle = 0)

## phenolics reduced in mature trees when no BP


### Thickness ##


thick_Bpy<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
            trait_category=="thickness" & species == "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"))

summary(thick_Bpy)

orchard_plot(thick_Bpy, xlab = "Standardised mean difference", transfm = "none", angle = 0)



thick_Bpn<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
            trait_category=="thickness" & species != "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
             R = list(species = cor),  control = list(optimizer = "optim"))

summary(thick_Bpn)

orchard_plot(thick_Bpn, xlab = "Standardised mean difference", transfm = "none", angle = 0)


### Toughness ###



tough_Bpy<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
            trait_category=="toughness" & species == "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              control = list(optimizer = "optim"))

summary(tough_Bpy)

orchard_plot(tough_Bpy, xlab = "Standardised mean difference", transfm = "none", angle = 0)


tough_Bpn<- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
            trait_category=="toughness" & species != "Betula pendula"),method="REML",
            random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"))

summary(tough_Bpn)

orchard_plot(tough_Bpn, xlab = "Standardised mean difference", transfm = "none", angle = 0)

############################

#-------------------#
#   funnel plots    #
#-------------------#


# xlim = c(-5, 5), ylim= c(2, 0))

f_thick <- funnel(thick, xlim = c(-5, 5), ylim= c(2, 0), main = "Thickness")

f_tough <- funnel(tough, xlim = c(-5, 5), ylim= c(2, 0), main = "Toughness")

f_LDMC <- funnel(LDMC, xlim = c(-5, 5), ylim= c(2, 0), main = "LDMC")

f_SLA <- funnel(SLA, xlim = c(-5, 5), ylim= c(2, 0), main = "SLA")

f_terp <- funnel(terp, xlim = c(-5, 5), ylim= c(2, 0), main = "Terpenoids")

f_phenolics <- funnel(phenolics, xlim = c(-5, 5), ylim= c(2, 0), main = "Phenolics")

f_N <- funnel(N, xlim = c(-5, 5), ylim= c(2, 0), main = "Nitrogen")

f_C <- funnel(C, xlim = c(-5, 5), ylim= c(2, 0), main = "Carbon")


#outliers removed from new csv file

ef3 <- read.csv("ef_removed_outliers.csv")

SLA_ro<- rma.mv(yi,vi,data=ef3, subset=(structure == "leaf" & diversity_type == "species" &
                                      trait_category=="SLA"),method="REML",
             random = list( ~ 1 | experiment, ~ 1 | ACC))
funnel(SLA_ro)

N_ro<- rma.mv(yi,vi,data=ef3, subset=(structure == "leaf" & diversity_type == "species" &
                                        trait_category=="N"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC))
funnel(N_ro)


(f_tough + f_thick )/
  (f_LDMC + f_SLA)

(f_thick | f_tough) /
  (f_LDMC | f_SLA) /
  (f_terp | f_phenolics) /
  (f_N| f_C)

f_LDMC + f_SLA



### trim and fill ###

SLA_taf <- trimfill(SLA)
funnel(SLA, legend=TRUE)


###### Publication bias tests #####


# 1) sample size bias 
#run meta-regressions with standard error as moderator 

ef$sei <- sqrt(ef$vi) ### calculate sampling error

## SLA

### removed regression data as this had huge varience

SLA_pb1 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
              trait_category=="SLA" & ACC != "18" & ACC != "57" & EF != "ef877"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods= sei)

summary(SLA_pb1)

SLA_pb1_r <- regplot(SLA_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "SLA")

# phenolics

phenolics_pb1 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="phenolics"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(phenolics_pb1)


phenolics_pb1_r <- regplot(phenolics_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "Phenolics")


# LDMC

LDMC_pb1 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="LDMC"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(LDMC_pb1)


LDMC_pb1_r <- regplot(LDMC_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "LDMC")


# N

N_pb1 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                  trait_category=="N" & ACC != c("18", "57") & EF != c("ef1215", "ef901",	"ef1235")), method="REML", 
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ sei)

#here i have remvoved studies 18 and 57 which had massive variences 

N_pb1 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                  trait_category=="N" & sei < 2), method="REML", 
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ sei)

# removed outlier here with argument  # sei < 2 #

summary(N_pb1)

N_pb1_r <- regplot(N_pb1, mod ="sei", pi=TRUE, refline=0, 
    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "Nitrogen")
  
  
# C
  

C_pb1 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="C" & sei < 2),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(C_pb1)


C_pb1_r <- regplot(C_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "Carbon")






SLA_pb1_r + phenolics_pb1_r





# terpenes



terpenes_pb1 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="terpenes"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(terpenes_pb1)


terpenes_pb1_r <- regplot(terpenes_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "terpenes")


## tough



toughness_pb1 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="toughness"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(toughness_pb1)


toughness_pb1_r <- regplot(toughness_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "toughness")


### thick



thickness_pb1 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="thickness"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ sei)


summary(thickness_pb1)


thickness_pb1_r <- regplot(thickness_pb1, pi=TRUE, refline=0, 
                    xlab="sampling error", ylab="SMD", bg= "#009e73", main = "thickness")

## no convergance 

#######################################################

# 2) publication bias

# run meta-regressions with year as moderator 


## SLA


SLA_pb2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
              trait_category=="SLA"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods= ~ year)

summary(SLA_pb2)

SLA_pb2_r <- regplot(SLA_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "SLA")


# phenolics

phenolics_pb2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="phenolics"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(phenolics_pb2)


phenolics_pb2_r <- regplot(phenolics_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "Phenolics")


# LDMC

LDMC_pb2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="LDMC"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(LDMC_pb2)


LDMC_pb2_r <- regplot(LDMC_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "LDMC")


# N


N_pb2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="N"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(N_pb2)

N_pb2_r <- regplot(N_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "Nitrogen")
  
  
# C
  

C_pb2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="C"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(C_pb2)


C_pb2_r <- regplot(C_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "Carbon")


# terpenes

terpenes_pb2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="terpenes"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(terpenes_pb2)


terpenes_pb2_r <- regplot(terpenes_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "terpenes")


## tough


oughness_pb2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="toughness"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(toughness_pb2)


toughness_pb2_r <- regplot(toughness_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "toughness")


# thickness

thickness_pb2 <- rma.mv(yi,vi,data=ef, subset=(structure == "leaf" & diversity_type == "species" &
                                          trait_category=="thickness"),method="REML",
              random = list( ~ 1 | experiment, ~ 1 | ACC, ~1 | species, ~1 | EF),
              R = list(species = cor), control = list(optimizer = "optim"),
             mods = ~ year)


summary(thickness_pb2)


thickness_pb2_r <- regplot(thickness_pb2, pi=TRUE, refline=0, 
                    xlab="year", ylab="SMD", bg= "#E69F00", main = "thickness")



### not significant effect apart from terpnes 
### too few effect sizes to be accurate

