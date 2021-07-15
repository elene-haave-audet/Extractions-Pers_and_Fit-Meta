# Data analysis for:
# Differences in resource acquisition, not allocation, mediate the relationship between behaviour and fitness: A systematic review and meta-analysis 
# Haave-Audet, E., Besson, A.A., Nakagawa, S., & Mathot, K.J.

# Download script and data files onto personal computer, keeping folder structure the same to run code

# Libraries----

# install the following packages to get orchaRd and metaAidR:
#install.packages("devtools")
#install.packages("tidyverse")
#install.packages("metafor")
#install.packages("patchwork")
#install.packages("R.rsp")

#devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, build_vignettes = TRUE)

library(orchaRd)

#install_github("daniel1noble/metaAidR")
# run line 20 is metaAidR is not installed on computer
library(metaAidR)

library(devtools)
library(tidyverse)
library(metafor) #meta-regression
library(rotl) #for phylogeny
library(ape) # for phylogeny
library(here) #to set working directory
library(matrixcalc) #for variance-covariance matrix
library(viridis) #ggplot color palatte
library(ggpubr)
library(patchwork)
library(lmodel2)

# Working Directory----
dir<-setwd(here::here())

# Functions----
  
  #* r to Zr#### 

# By S. Nakagawa for Mathot et al 2018
r.to.Zr<-function(r){
  Zr<-0.5*(log(1+r)-log(1-r))}

  #* Zr to r####

# By S. Nakagawa for Mathot et al 2018
Zr.to.r<-function(Zr){
  r<-(exp(2*Zr)-1)/(exp(2*Zr)+1)}

  #* Estimates from rma objects#### 

# From Hayward, Poulin & Nakagawa 2019
get_est <- function(model, mod = " ") {
  
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  estimate <- as.numeric(model$beta)
  lowerCL <- model$ci.lb
  upperCL <- model$ci.ub
  
  table <- tibble(name = name, estimate = estimate, lowerCL = lowerCL, upperCL = upperCL)
}

  #* Prediction intervals for rma objects#### 

# From Hayward, Poulin & Nakagawa 2019
get_pred <- function(model, mod = " ") {
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  len <- length(name)
  
  if (len != 1) {
    newdata <- matrix(NA, ncol = len, nrow = len)
    for (i in 1:len) {
      # getting the position of unique case from X (design matrix)
      pos <- which(model$X[, i] == 1)[[1]]
      newdata[, i] <- model$X[pos, ]
    }
    pred <- predict.rma(model, newmods = newdata)
  } else {
    pred <- predict.rma(model)
  }
  lowerPR <- pred$cr.lb
  upperPR <- pred$cr.ub
  
  table <- tibble(name = name, lowerPR = lowerPR, upperPR = upperPR)
}

  #* Function for VCV####

#' @title Covariance and correlation matrix function basing on shared level ID
#' @description Function for generating simple covariance and correlation matrices 
#' @param data Dataframe object containing effect sizes, their variance, unique IDs and clustering variable
#' @param V Name of the variable (as a string – e.g, "V1") containing effect size variances variances
#' @param cluster Name of the variable (as a string – e.g, "V1") indicating which effects belong to the same cluster. Same value of 'cluster' are assumed to be nonindependent (correlated).
#' @param obs Name of the variable (as a string – e.g, "V1") containing individual IDs for each value in the V (Vector of variances). If this parameter is missing, label will be labelled with consecutive integers starting from 1.
#' @param rho Known or assumed correlation value among effect sizes sharing same 'cluster' value. Default value is 0.5.
#' @param type Optional logical parameter indicating whether a full variance-covariance matrix (default or "vcv") is needed or a correlation matrix ("cor") for the non-independent blocks of variance values.
#' @export

make_VCV_matrix <- function(data, V, cluster, obs, type=c("vcv", "cor"), rho=0.5){
  type <- match.arg(type)
  if (missing(data)) {
    stop("Must specify dataframe via 'data' argument.")
  }
  if (missing(V)) {
    stop("Must specify name of the variance variable via 'V' argument.")
  }
  if (missing(cluster)) {
    stop("Must specify name of the clustering variable via 'cluster' argument.")
  }
  if (missing(obs)) {
    obs <- 1:length(V)   
  }
  if (missing(type)) {
    type <- "vcv" 
  }
  
  new_matrix <- matrix(0,nrow = dim(data)[1],ncol = dim(data)[1]) #make empty matrix of the same size as data length
  rownames(new_matrix) <- data[ ,obs]
  colnames(new_matrix) <- data[ ,obs]
  # find start and end coordinates for the subsets
  shared_coord <- which(data[ ,cluster] %in% data[duplicated(data[ ,cluster]), cluster]==TRUE)
  # matrix of combinations of coordinates for each experiment with shared control
  combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,cluster], function(x) t(utils::combn(x,2))))
  
  if(type == "vcv"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho * sqrt(data[p1,V]) * sqrt(data[p2,V])
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- data[ ,V]   #add the diagonal
  }
  
  if(type == "cor"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- 1   #add the diagonal of 1
  }
  
  return(new_matrix)
}

# Data----

  #* Load====
data<-read.csv("Analyses/HaaveAudet_meta_extractions_20210617.csv")

  #* Zr to r====

# Convert Zr estimates from Moiron et al 2019 to r
data$r<-ifelse(is.na(data$r),Zr.to.r(data$Zr),(data$r))

  #* Get Zr====

# Zr
data$Zr<-r.to.Zr(data$r)

# Variance (VZr)
VZr<-1/(data$NInd-3)
data$VZr<- VZr

  #* Observation level random effects====
data$obs<- as.factor(1:dim(data)[1])

  #* Names for phylogeny====

# Convert taxon names to lower case
data$LatinName<-tolower(data$LatinName)

# replace borken names in dataset with names that will work & change common name to match
data<-mutate_if(data, is.character, 
                str_replace_all, patter="aquarius remigis", 
                replacement = "gerris remigis") #different water strider
data<-mutate_if(data, is.character, 
                str_replace_all, patter="stream water strider", 
                replacement = "water strider")
data<-mutate_if(data, is.character, 
                str_replace_all, patter="zootoca vivipara", 
                replacement = "lacerta vivipara") #have same common name, no change needed
data<-mutate_if(data, is.character, 
                str_replace_all, patter="pomacentrus wardi", 
                replacement = "pomacentrus moluccensis")
data<-mutate_if(data, is.character, 
                str_replace_all, patter="ward's damsel", 
                replacement = "lemon damselfish")
data<-mutate_if(data, is.character, 
                str_replace_all, patter="poecilia reticulata", 
                replacement = "gambusia geiseri")
data<-mutate_if(data, is.character, 
                str_replace_all, patter="guppy", 
                replacement = "largespring mosquitofish")
data<-mutate_if(data, is.character, 
                str_replace_all, patter="pan troglodytes", 
                replacement = "pan troglodytes troglodytes")

# 1) Phylogenetic Tree (Figure 3)----
phylo.2<-tnrs_match_names(unique(data$LatinName), context="Animals")
head(phylo.2) # Check that number of species matches number of distinct LatinName in data

# Create a named vector that maps the names I have for each species
  ## to the names Open Tree uses:
phylo.map.2<-structure(phylo.2$search_string, names=phylo.2$unique_name)

# Get a tree 
tree.2<-tol_induced_subtree(ott_id(phylo.2)[is_in_tree(ott_id(phylo.2))])
  # Removed [is_in_tree(ott_id(phylo))] so that all taxa would be included

plot(tree.2, show.tip.label = FALSE)

otl_tips.2 <- strip_ott_ids(tree.2$tip.label, remove_underscores = TRUE)
tree.2$tip.label <- phylo.map.2[ otl_tips.2 ]

# from RMarkdown for Mathot et al 2019
#tree.2<-drop.tip(tree.2,tree.2$tip.label[-match(data.names$LatinName, tree.2$tip.label)])

tree.2<-compute.brlen(tree.2)
is.ultrametric(tree.2)
cor.2<-vcv(tree.2,cor=T)

  #* Tree by Class====

# Keep one representative species from each class to plot reduced tree
to.keep<-c("larinioides sclopetarius", "cherax destructor", "ischnura genei", "euprymna tasmanica", "chlorostoma funebralis", "gadus morhua", "lithobates sylvaticus", "lacerta vivipara", "columba livia", "felis catus")

plot(keep.tip(tree.2, to.keep ),show.tip.label=FALSE)

  #* Effect of Class====

# Extract mean effect sizes by class to plot with phylogenetic tree
# VCV matrix
matrix.phylo <- make_VCV_matrix(data, V = "VZr", cluster = "GroupID_global", obs = "obs")

# meta-regression with class as a moderator
m.class<-rma.mv(yi = Zr, V = matrix.phylo, mod= ~Class -1, random = list(~1|LatinName,~1|RecordNo,~1|obs),R = list(LatinName = cor.2),data = data)

  #* Figure 3====

# Get estimates
res_class <- get_est(m.class, mod = "Class")

# Add sample size (k) for each category
k_class<-data %>% group_by(Class) %>% count()
n_class<-data %>% group_by(Class) %>% summarise(N=n_distinct(RecordNo))

k_class<-left_join(k_class,n_class, by="Class")

# Get estimates and predictions
pred_class<-get_pred(m.class, mod = "Class")
res_class<-left_join(res_class, k_class, by=c("name"="Class")) %>% left_join(pred_class)

# Order Latin Names to match phylogenetic tree branches
res_class$name<-factor(res_class$name, levels=c("Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Gastropoda", "Cephalopoda", "Insecta", "Malacostraca","Arachnida"))
res_class<-as.data.frame(arrange(res_class, name))

# Plot
plot.class<-ggplot(data=res_class, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1, 1), breaks = seq(-1, 1, by = 0.5) )+
  #geom_point(data = pheno.s, 
  #aes(x= Zr, y = LatinName, size = ((1/VZr) + 3)), alpha=0.3)+
  #scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  #geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, colour = "gray", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(fill = "black", size = 2, shape = 21)+
  #scale_color_viridis(discrete = TRUE)+
  annotate('text', x = 0.75, y= 1:10, label= res_class$n, parse = TRUE, hjust = "center", size = 5)+
  annotate('text', x = 1, y= 1:10, label= res_class$N, parse = TRUE, size = 5)+
  labs(x = "Zr", y = "") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5),
        axis.text.x = element_text(size = 11))
plot.class

# 2) Allocation vs Acquisition----
  #* Match paired estimates====

# Subset data by fitness proxy and merge by study and individual ID
surv<- subset(data, FitnessProxy=="survival")

repro<- subset(data, FitnessProxy=="repro")

fitness_wide=inner_join(surv, repro, suffix=c("_surv", "_repro"), 
                        by=c("RecordNo", "Author", "Year", "LatinName", "CommonName", "Class", "StudyPop", "Age", "Sex", "Origin", "GroupID_global", "Behaviour", "LevelBehav", "Conditions"))

# Give each row a unique pair ID
fitness_wide$Fit_GroupID<-as.factor(1:dim(fitness_wide)[1]) 

# Give each estimate in a pair its own row (i.e. re-split surv & repro)
surv_paired<-fitness_wide %>% 
  select(RecordNo, Author, Year, Class, CommonName, LatinName, GroupID_global, Fit_GroupID, FitnessProxy_surv, Behaviour, NInd_surv, r_surv, Zr_surv, VZr_surv, sign_surv, Conditions) %>% 
  rename(Zr=Zr_surv, VZr=VZr_surv, NInd=NInd_surv, FitnessProxy=FitnessProxy_surv, r=r_surv, sign=sign_surv)

repro_paired<-fitness_wide %>% 
  select(RecordNo, Author, Year, Class, CommonName, LatinName, GroupID_global, Fit_GroupID, FitnessProxy_repro, Behaviour, NInd_repro, r_repro, Zr_repro, VZr_repro, sign_repro, Conditions) %>% 
  rename(Zr=Zr_repro, VZr=VZr_repro, NInd=NInd_repro, FitnessProxy=FitnessProxy_repro, r=r_repro, sign=sign_repro)

# Merge single observations of paired estimates into same df
fitness_long<-bind_rows(surv_paired, repro_paired)

  #* Among-individuals====

# Subset level=among
amg<-subset(fitness_wide, LevelBehav=="among")

# Correlation b/w estimates using surv and repro
cor.test(amg$Zr_surv, amg$VZr_repro) #r=0.206, -0.288-0.614

# Extract regression line for plotting
mod.amg<-lmodel2(Zr_repro~Zr_surv, data=amg, "interval", "interval", 99)
reg.amg<-mod.amg$regression.results
# rma= row 4, slope= column 3, intercept= column 2
reg.amg<-subset(reg.amg, Method=="RMA")

# Plot correlation
amg.cor<-ggplot(data = amg, aes(Zr_surv, Zr_repro, size=(1/VZr_surv)+3))+
  geom_abline(data=reg.amg, aes(intercept=Intercept, slope= Slope), linetype="dashed")+
  geom_point(shape= 21, fill="grey90")+
  xlab(label="Zr survival")+
  ylab(label="Zr reproduction")+
  xlim(-0.5,1)+
  ylim(-0.5,1)+
  labs(size= expression(paste(italic(n),"(# of individuals)")))+
  guides(size=guide_legend(override.aes = list(linetype=c(0,0,0))))+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size=14))+
  theme(legend.position=c(1,0), legend.justification = c(1,0))+
  theme(legend.direction = "horizontal")
amg.cor

  #* Within-individuals (phenotypic)====

# Subset level=phenotypic
pheno<-subset(fitness_wide, LevelBehav=="phenotypic")

# Get correlation b/w estimates of surv and repro
cor.test(pheno$Zr_surv, pheno$Zr_repro) #r=0.387, 0.154-0.579

# Get regression line for plot
mod.pheno<-lmodel2(Zr_repro~ Zr_surv, data=pheno, "interval", "interval", 99)
reg.pheno<-mod.pheno$regression.results
# rma= row 4, slope= column 3, intercept= column 2
reg.pheno<-subset(reg.pheno, Method=="RMA")

# Plot correlation
pheno.cor<-ggplot(data = pheno, aes(Zr_surv, Zr_repro, size=(1/VZr_surv)+3))+
  geom_abline(intercept = reg.pheno$Intercept, slope = reg.pheno$Slope, linetype="dashed")+
  geom_point(shape= 21, fill="grey90")+
  xlab(label="Zr survival")+
  ylab(label="Zr reproduction")+
  xlim(-0.5, 1)+
  ylim(-0.5, 1)+
  labs(size= expression(paste(italic(n),"(# of individuals)")))+
  guides(size=guide_legend(override.aes = list(linetype=c(0,0,0))))+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size=14))+
  theme(legend.position=c(1,0), legend.justification = c(1,0))+
  theme(legend.direction = "horizontal")
pheno.cor

  #* Figure 4====
fit.plots<-ggarrange(amg.cor, pheno.cor, 
                     labels="auto",
                     font.label = list(size = 12, face = "bold"))
fit.plots

# 3) Effect of behaviour and condition----

  #* a) Phenotypic-Survival====

# Subset data: level=phenotypic, fitness=survival
pheno.s<- subset(data, LevelBehav=="phenotypic" & FitnessProxy=="survival")

# Var-covar matrix
matrix.ps <- make_VCV_matrix(pheno.s, V = "VZr", cluster = "GroupID_global", obs = "obs")

    #** Phylogeny####
phylo.ps<-tnrs_match_names(unique(pheno.s$LatinName), context="Animals")
head(phylo.ps) 

# Create a named vector that maps the names I have for each species
  ## to the names Open Tree uses:
phylo.map.ps<-structure(phylo.ps$search_string, names=phylo.ps$unique_name)

# Get a tree 
tree.ps<-tol_induced_subtree(ott_id(phylo.ps)[is_in_tree(ott_id(phylo.ps))])
### removed [is_in_tree(ott_id(phylo))] so that all taxa would be included

plot(tree.ps, show.tip.label = FALSE)

otl_tips.ps <- strip_ott_ids(tree.ps$tip.label, remove_underscores = TRUE)
tree.ps$tip.label <- phylo.map.ps[ otl_tips.ps ]

# from RMarkdown for Mathot et al 2019
tree.ps<-drop.tip(tree.ps,tree.ps$tip.label[-match(pheno.s$LatinName, tree.ps$tip.label)])

tree.ps<-compute.brlen(tree.ps)
is.ultrametric(tree.ps)
phylo.cor.ps<-vcv(tree.ps,cor=T)

    #** Random effects####

# Full model (all random effects) 
#ps.full<-rma.mv(yi = Zr, V = matrix.ps, random = list( ~1|LatinName, ~1|CommonName,~1|RecordNo,~1|GroupID_global, ~1|obs),R = list(LatinName = phylo.cor.ps), data = pheno.s)
#i2_ml(ps.full)
  # Group ID has very small effect, so will be removed

# Full model w/out Group ID (TO USE)
ps.full1<-rma.mv(yi = Zr, V = matrix.ps, random = list(~1|LatinName, ~1|CommonName, ~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.ps), data = pheno.s)
summary(ps.full1)
i2_ml(ps.full1)
Zr_to_r(ps.full1$b) #0.04779

    #** Behaviour####
ps.behav<-rma.mv(yi = Zr, V = matrix.ps, mod= ~Behaviour -1, random = list(~1|LatinName, ~1|CommonName,~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.ps),data = pheno.s)
summary(ps.behav)
r2_ml(ps.behav) 
i2_ml(ps.behav)

    #** Fitness Condition####
ps.fit<-rma.mv(yi = Zr, V = matrix.ps, mod= ~FitnessMeasured -1, random = list(~1|LatinName, ~1|CommonName,~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.ps),data = pheno.s)
summary(ps.fit)
r2_ml(ps.fit) #r2m=0.0118
i2_ml(ps.fit)

  #* b) Phenotypic-Reproduction====

# Subset level=phenotypic, fitness=repro
pheno.r<- subset(data, LevelBehav=="phenotypic" & FitnessProxy=="repro")

# Var-Covar matrix

# Create square matrix matching the number of Zr, filled with zeros
matrix.pr <- matrix(0,nrow = dim(pheno.r)[1],ncol = dim(pheno.r)[1])

rownames(matrix.pr) <- pheno.r$obs
colnames(matrix.pr) <- pheno.r$obs 

# Find start and end coordinates for the subsets
shared_coord_pr <- which(pheno.r$GroupID_global %in% pheno.r$GroupID_global[duplicated(pheno.r$GroupID_global)]==TRUE)
# Matrix of combinations of coordinates for each experiment with shared control
combinations_pr <- do.call("rbind", tapply(shared_coord_pr, data[shared_coord_pr,"GroupID_global"], function(x) t(combn(x,2))))

# Calculate covariance values between  values at the positions in shared_list and place them on the matrix

for (i in 1:dim(combinations_pr)[1]){
  p1.pr <- combinations_pr[i,1]
  p2.pr <- combinations_pr[i,2]
  p1_p2_cov_pr <- 0.5*sqrt(pheno.r[p1.pr,"VZr"]) * sqrt(pheno.r[p2.pr,"VZr"])
  matrix.pr[p1.pr,p2.pr] <- p1_p2_cov_pr
  matrix.pr[p2.pr,p1.pr] <- p1_p2_cov_pr
}

# Add the diagonal - use "varlnAA"
# Create variance-covariance matrix
diag(matrix.pr) <- pheno.r$VZr
is.positive.definite(matrix.pr)

    #** Phylogeny####
phylo.pr<-tnrs_match_names(unique(pheno.r$LatinName), context="Animals")
head(phylo.pr)

# Create a named vector that maps the names I have for each species
  ## to the names Open Tree uses:
phylo.map.pr<-structure(phylo.pr$search_string, names=phylo.pr$unique_name)

# Get a tree 
tree.pr<-tol_induced_subtree(ott_id(phylo.pr)[is_in_tree(ott_id(phylo.pr))])
  ### removed [is_in_tree(ott_id(phylo))] so that all taxa would be included

plot(tree.pr, show.tip.label = FALSE)

otl_tips.pr <- strip_ott_ids(tree.pr$tip.label, remove_underscores = TRUE)
tree.pr$tip.label <- phylo.map.pr[ otl_tips.pr ]

# From RMarkdown for Mathot et al 2019
tree.pr<-drop.tip(tree.pr,tree.pr$tip.label[-match(pheno.r$LatinName, tree.pr$tip.label)])

tree.pr<-compute.brlen(tree.pr)
is.ultrametric(tree.pr)
phylo.cor.pr<-vcv(tree.pr,cor=T)

    #** Random effects####

# Full model with all random effect
#pr.full<-rma.mv(yi = Zr, V = matrix.pr, random = list(~1|CommonName,~1|LatinName, ~1|GroupID_global, ~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.pr), data = pheno.r)
#summary(pr.full)
# check contribution to heterogeneity to select random effects
#i2_ml(pr.full)

# Model w/out Group ID (TO USE)
pr.full1<-rma.mv(yi = Zr, V = matrix.pr, random = list(~1|CommonName,~1|LatinName, ~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.pr), data = pheno.r)
i2_ml(pr.full1)
summary(pr.full1)
Zr_to_r(pr.full1$b) #0.057479

    #** Behaviour####
pr.behav<-rma.mv(yi = Zr, V = matrix.pr, mod= ~Behaviour -1, random = list(~1|LatinName,~1|CommonName,~1|RecordNo, ~1|obs), R=list(LatinName=phylo.cor.pr), data = pheno.r)
summary(pr.behav)
r2_ml(pr.behav) #R2m=0.0907
i2_ml(pr.behav)

    #** Fitness Condition####
pr.fit<-rma.mv(yi = Zr, V = matrix.pr, mod= ~FitnessMeasured -1, random = list(~1|LatinName,~1|CommonName,~1|RecordNo,~1|obs), R=list(LatinName=phylo.cor.pr), data = pheno.r)
summary(pr.fit)
r2_ml(pr.fit) #r2m=0.00041
i2_ml(pr.fit)

  #* c) Among-Survival====

# Subset level=among, fitness=survival
amg.s<- subset(data, LevelBehav=="among" & FitnessProxy=="survival")

# Var-covar matrix
# Create square matrix matching the number of Zr, filled with zeros
matrix.as <- matrix(0,nrow = dim(amg.s)[1],ncol = dim(amg.s)[1])

rownames(matrix.as) <- amg.s$obs
colnames(matrix.as) <- amg.s$obs 

# Find start and end coordinates for the subsets
shared_coord_as <- which(amg.s$GroupID_global %in% amg.s$GroupID_global[duplicated(amg.s$GroupID_global)]==TRUE)
# Matrix of combinations of coordinates for each experiment with shared control
combinations_as <- do.call("rbind", tapply(shared_coord_as, data[shared_coord_as,"GroupID_global"], function(x) t(combn(x,2))))

# Calculate covariance values between  values at the positions in shared_list and place them on the matrix

for (i in 1:dim(combinations_as)[1]){
  p1.as <- combinations_as[i,1]
  p2.as <- combinations_as[i,2]
  p1_p2_cov_as <- 0.5*sqrt(amg.s[p1.as,"VZr"]) * sqrt(amg.s[p2.as,"VZr"])
  matrix.as[p1.as,p2.as] <- p1_p2_cov_as
  matrix.as[p2.as,p1.as] <- p1_p2_cov_as
}

# Add the diagonal - use "varlnAA"
# Create variance-covariance matrix
diag(matrix.as) <- amg.s$VZr
is.positive.definite(matrix.as)

    #** Phylogeny####
phylo.as<-tnrs_match_names(unique(amg.s$LatinName), context="Animals")
head(phylo.as) #31 unique spp names

# Create a named vector that maps the names I have for each species
  ## to the names Open Tree uses:
phylo.map.as<-structure(phylo.as$search_string, names=phylo.as$unique_name)

# Get a tree 
tree.as<-tol_induced_subtree(ott_id(phylo.as)[is_in_tree(ott_id(phylo.as))])
### removed [is_in_tree(ott_id(phylo))] so that all taxa would be included

plot(tree.as, show.tip.label = FALSE)

otl_tips.as <- strip_ott_ids(tree.as$tip.label, remove_underscores = TRUE)
tree.as$tip.label <- phylo.map.as[ otl_tips.as ]

# From RMarkdown for Mathot et al 2019
tree.as<-drop.tip(tree.as,tree.as$tip.label[-match(amg.s$LatinName, tree.as$tip.label)])

tree.as<-compute.brlen(tree.as)
is.ultrametric(tree.as)
phylo.cor.as<-vcv(tree.as,cor=T)

    #** Random effects####

# Full model with all random effects
#as.full<-rma.mv(yi = Zr, V = matrix.as, random = list(~1|CommonName, ~1|GroupID_global, ~1|LatinName,~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.as), data = amg.s)
#summary(as.full)

# check contribution to heterogeneity to select random effects
#i2_ml(as.full)
#total:954764351
#CommonName: 0.083497677 
#GroupID: 0.004236229 
#ObservationID:0.638749884 
#RecordID: 0.135683334  
#Phylo: 0.092597226

# Model w/out Group ID (TO USE)
as.reduced<-rma.mv(yi = Zr, V = matrix.as, random = list(~1|CommonName,  ~1|LatinName, ~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.as), data = amg.s)
i2_ml(as.reduced) #used as null model
Zr_to_r(as.reduced$b) #-0.0213
summary(as.reduced)

    #** Behaviour####
as.behav<-rma.mv(yi = Zr, V = matrix.as, mod= ~Behaviour -1, random = list(~1|LatinName, ~1|CommonName,~1|RecordNo, ~1|obs),R = list(LatinName = phylo.cor.as), data = amg.s)
summary(as.behav)
r2_ml(as.behav) #R2m=0.062
i2_ml(as.behav)

    #** Fitness Condition####
as.fit<-rma.mv(yi = Zr, V = matrix.as, mod= ~FitnessMeasured -1, random = list( ~1|LatinName,~1|CommonName,~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.as), data = amg.s)
summary(as.fit)
r2_ml(as.fit) #R2m=0.081
i2_ml(as.fit)

  #* d) Among-Reproduction====

# Subset level=among, fitness=repro
amg.r<- subset(data, LevelBehav=="among" & FitnessProxy=="repro")

# Var-Covar matrix
matrix.ar <- matrix(0,nrow = dim(amg.r)[1],ncol = dim(amg.r)[1])

rownames(matrix.ar) <- amg.r$obs
colnames(matrix.ar) <- amg.r$obs 

# Find start and end coordinates for the subsets
shared_coord_ar <- which(amg.r$GroupID_global %in% amg.r$GroupID_global[duplicated(amg.r$GroupID_global)]==TRUE)
# Matrix of combinations of coordinates for each experiment with shared control
combinations_ar <- do.call("rbind", tapply(shared_coord_ar, data[shared_coord_ar,"GroupID_global"], function(x) t(combn(x,2))))

# Calculate covariance values between  values at the positions in shared_list and place them on the matrix

for (i in 1:dim(combinations_ar)[1]){
  p1.ar <- combinations_ar[i,1]
  p2.ar <- combinations_ar[i,2]
  p1_p2_cov_ar <- 0.5*sqrt(amg.r[p1.ar,"VZr"]) * sqrt(amg.r[p2.ar,"VZr"])
  matrix.ar[p1.ar,p2.ar] <- p1_p2_cov_ar
  matrix.ar[p2.ar,p1.ar] <- p1_p2_cov_ar
}

# Add the diagonal - use "varlnAA"
# Create variance-covariance matrix
diag(matrix.ar) <- amg.r$VZr
is.positive.definite(matrix.ar)

    #** Phylogeny####
phylo.ar<-tnrs_match_names(unique(amg.r$LatinName), context="Animals")
head(phylo.ar) 

# Create a named vector that maps the names I have for each species
  ## to the names Open Tree uses:
phylo.map.ar<-structure(phylo.ar$search_string, names=phylo.ar$unique_name)

# Get a tree 
tree.ar<-tol_induced_subtree(ott_id(phylo.ar)[is_in_tree(ott_id(phylo.ar))])
### removed [is_in_tree(ott_id(phylo))] so that all taxa would be included

plot(tree.ar, show.tip.label = FALSE)

otl_tips.ar <- strip_ott_ids(tree.ar$tip.label, remove_underscores = TRUE)
tree.ar$tip.label <- phylo.map.ar[ otl_tips.ar ]

# From RMarkdown for Mathot et al 2019
tree.ar<-drop.tip(tree.ar,tree.ar$tip.label[-match(amg.r$LatinName, tree.ar$tip.label)])

tree.ar<-compute.brlen(tree.ar)
is.ultrametric(tree.ar)
phylo.cor.ar<-vcv(tree.ar,cor=T)

    #** Random effects####

# Full model with all random effects
#ar.full<-rma.mv(yi = Zr, V = matrix.ar, random = list(~1|CommonName, ~1|GroupID_global, ~1|LatinName, ~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.ar), data = amg.r)
#summary(ar.full)

# check contribution to heterogeneity to select random effects
#i2_ml(ar.full)
#total:0.9158521
#CommonName: 0.1822306 
#GroupID: 8.241572e-11 
#ObservationID:0.7336214 
#RecordID: 3.144623e-11   
#Phylo: 6.588322e-11

# Model w/out Group ID (TO USE)
ar.null<-rma.mv(yi = Zr, V = matrix.ar, random = list(~1|CommonName, ~1|LatinName, ~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.ar), data = amg.r)
i2_ml(ar.null)
summary(ar.null)
Zr_to_r(ar.null$b) #-0.0136

    #** Behaviour####
ar.behav<-rma.mv(yi = Zr, V = matrix.ar, mod= ~Behaviour -1, random = list(~1|CommonName, ~1|LatinName, ~1|RecordNo, ~1|obs), R=list(LatinName=phylo.cor.ar), data = amg.r)
summary(ar.behav)
r2_ml(ar.behav) #R2m=0.1168
i2_ml(ar.behav)
    
    #** Fitness Condition####
ar.fit<-rma.mv(yi = Zr, V = matrix.ar, mod= ~FitnessMeasured -1, random = list( ~1|CommonName, ~1|LatinName, ~1|RecordNo, ~1|obs), R=list(LatinName=phylo.cor.ar), data = amg.r)
summary(ar.fit)
r2_ml(ar.fit) #R2m=0.00013
i2_ml(ar.fit)

  #* e) Within-Reproduction====

# Subset level=within, fitness=repro
wthn.r<- subset(data, LevelBehav=="within" & FitnessProxy=="repro")

# Var-covar matrix

# Create square matrix matching the number of Zr, filled with zeros
matrix.wr <- matrix(0,nrow = dim(wthn.r)[1],ncol = dim(wthn.r)[1])

rownames(matrix.wr) <- wthn.r$obs
colnames(matrix.wr) <- wthn.r$obs 

# Find start and end coordinates for the subsets
shared_coord_wr <- which(wthn.r$GroupID_global %in% wthn.r$GroupID_global[duplicated(wthn.r$GroupID_global)]==TRUE)
# Matrix of combinations of coordinates for each experiment with shared control
combinations_wr <- do.call("rbind", tapply(shared_coord_wr, data[shared_coord_wr,"GroupID_global"], function(x) t(combn(x,2))))

# Calculate covariance values between  values at the positions in shared_list and place them on the matrix

for (i in 1:dim(combinations_wr)[1]){
  p1.wr <- combinations_wr[i,1]
  p2.wr <- combinations_wr[i,2]
  p1_p2_cov_wr <- 0.5*sqrt(wthn.r[p1.wr,"VZr"]) * sqrt(wthn.r[p2.wr,"VZr"])
  matrix.wr[p1.wr,p2.wr] <- p1_p2_cov_wr
  matrix.wr[p2.wr,p1.wr] <- p1_p2_cov_wr
}

# Add the diagonal - use "varlnAA"
# Create variance-covariance matrix
diag(matrix.wr) <- wthn.r$VZr
is.positive.definite(matrix.wr)

    #** Phylogeny####
phylo.wr<-tnrs_match_names(unique(wthn.r$LatinName), context="Animals")
head(phylo.wr) #7 unique spp names

# Create a named vector that maps the names I have for each species
  ## to the names Open Tree uses:
phylo.map.wr<-structure(phylo.wr$search_string, names=phylo.wr$unique_name)

# Get a tree 
tree.wr<-tol_induced_subtree(ott_id(phylo.wr)[is_in_tree(ott_id(phylo.wr))])
### removed [is_in_tree(ott_id(phylo))] so that all taxa would be included

plot(tree.wr, show.tip.label = FALSE)

otl_tips.wr <- strip_ott_ids(tree.wr$tip.label, remove_underscores = TRUE)
tree.wr$tip.label <- phylo.map.wr[ otl_tips.wr ]

# From RMarkdown for Mathot et al 2019
tree.wr<-drop.tip(tree.wr,tree.wr$tip.label[-match(wthn.r$LatinName, tree.wr$tip.label)])

tree.wr<-compute.brlen(tree.wr)
is.ultrametric(tree.wr)
phylo.cor.wr<-vcv(tree.wr,cor=T)

    #** Random effects####
wr.null<-rma.mv(yi = Zr, V = matrix.wr, random = list(~1|CommonName, ~1|LatinName, ~1|RecordNo,~1|obs),R = list(LatinName = phylo.cor.wr), data = wthn.r)
i2_ml(wr.null)
summary(wr.null)
Zr_to_r(wr.null$b) #-0.00068

    #** Behaviour####
wr.behav<-rma.mv(yi = Zr, V = matrix.wr, mod= ~Behaviour -1, random = list(~1|CommonName, ~1|LatinName, ~1|RecordNo, ~1|obs), R=list(LatinName=phylo.cor.wr), data = wthn.r)
summary(wr.behav)
r2_ml(wr.behav) #R2m=0.5732
i2_ml(wr.behav)

    #** Fitness Condition####
wr.fit<-rma.mv(yi = Zr, V = matrix.wr, mod= ~FitnessMeasured -1, random = list( ~1|CommonName, ~1|LatinName, ~1|RecordNo, ~1|obs), R=list(LatinName=phylo.cor.wr), data = wthn.r)
summary(wr.fit)
r2_ml(wr.fit) #R2m=0.135
i2_ml(ar.fit)

  #* Figure 5 & 6====

# Phenotypic: survival
# Behaviour

## get estimates
res_ps.b <- get_est(ps.behav, mod = "Behaviour")

## prepare for plotting
### add sample size (k) for each category
k_ps.b<-pheno.s %>% group_by(Behaviour) %>% count()

### get estimates and predictions
pred_ps.b<-get_pred(ps.behav, mod = "Behaviour")
res_ps.b<-left_join(res_ps.b, k_ps.b, by=c("name"="Behaviour")) %>% left_join(pred_ps.b)

res_ps.b<-add_row(.data=res_ps.b,.after=3, name="courtship",n=0)

#create colour scale for plotting
my_colours<-c("#440154" ,"#472D7B","#3B528B","#2C728E","#21908C","#27AD81","#5DC863","#AADC32","#FDE725")
names(my_colours)<-levels(factor(c(levels(res_ps.b$name))))
my_scale<-scale_fill_manual(name="name", values=my_colours)

### plot
pheno.surv.behav<-ggplot(data=res_ps.b, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  geom_point(data = pheno.s, 
             aes(x= Zr, y = Behaviour, size = ((1/VZr) + 3), colour = Behaviour), alpha=0.3)+
  scale_fill_manual(values=my_colours)+
  #scale_fill_manual(values = c("#440154" ,"#472D7B","#3B528B","#2C728E","#21908C","#27AD81","#5DC863","#AADC32","#FDE725"))+
  #scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(aes(fill = name), size = 3, shape = 21)+
  #scale_colour_manual(values=c("#440154" ,"#472D7B","#3B528B","#2C728E","#21908C","#27AD81","#5DC863","#AADC32","#FDE725"))+
  scale_color_viridis(discrete = TRUE)+
  annotate('text', x = 3, y= 1:9, label= paste("italic(k)==", res_ps.b$n), parse = TRUE, hjust = "right", size = 3.5)+
  labs(x = "Zr", y = "") +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5))
pheno.surv.behav

# Overall
## get estimates
res_ps.full <- get_est(ps.full1)

## prepare for plotting

### get estimates and predictions
pred_ps.full<-get_pred(ps.full1)
res_ps.full<-full_join(res_ps.full, pred_ps.full)

### plot
pheno.surv.full<-ggplot(data=res_ps.full, aes(x=estimate, y="Phenotypic-Survival"))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  scale_y_discrete()+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(size = 3, shape = 21, fill="black")+
  annotate("text", x=3, y=1, label= "k= 218", size=3.5, hjust="right")+
  labs(x = "", y = "") +
  ggtitle("a")+
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=12, colour = "black"),
        plot.title = element_text(face="bold", size = 10, hjust = -0.1, vjust = -6))
pheno.surv.full

# Fitness

## get estimates
res_ps.f <- get_est(ps.fit, mod = "FitnessMeasured")

## prepare for plotting
### add sample size (k) for each category
k_ps.f<-pheno.s %>% group_by(FitnessMeasured) %>% count()

### get estimates and predictions
pred_ps.f<-get_pred(ps.fit, mod = "FitnessMeasured")
res_ps.f<-left_join(res_ps.f, k_ps.f, by=c("name"="FitnessMeasured")) %>% left_join(pred_ps.f)

### plot
pheno.surv.fit<-ggplot(data=res_ps.f, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  scale_y_discrete(labels=c("lab", "wild"))+
  geom_point(data = pheno.s, 
             aes(x= Zr, y = FitnessMeasured, size = ((1/VZr) + 3), colour = FitnessMeasured), alpha=0.2)+
  scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(aes(fill = name), size = 3, shape = 21)+
  scale_color_viridis(discrete = TRUE)+
  annotate('text', x = 3, y= 1:2, label= paste("italic(k)==", res_ps.f$n), parse = TRUE, hjust = "right", size = 3.5)+
  labs(x = "Zr", y = "") +
  labs(title = "Phenotypic-Survival", subtitle = "a")+
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5))+
  theme(plot.title = element_text(size=12, hjust = 0), plot.subtitle = element_text(size=10, hjust = -0.05, vjust = -6, face = "bold"))
pheno.surv.fit

# Phenotypic: reproduction
# Behaviour

## get estimates
res_pr.b <- get_est(pr.behav, mod = "Behaviour")

## prepare for plotting
### add sample size (k) for each category
k_pr.b<-pheno.r %>% group_by(Behaviour) %>% count()

### get estimates and predictions
pred_pr.b<-get_pred(pr.behav, mod = "Behaviour")
res_pr.b<-left_join(res_pr.b, k_pr.b, by=c("name"="Behaviour")) %>% left_join(pred_pr.b)

### plot
pheno.repro.behav<-ggplot(data=res_pr.b, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  geom_point(data = pheno.r, 
             aes(x= Zr, y = Behaviour, size = ((1/VZr) + 3), colour = Behaviour), alpha=0.3)+
  scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(aes(fill = name), size = 3, shape = 21)+
  scale_color_viridis(discrete = TRUE)+
  annotate('text', x = 3, y= 1:9, label= paste("italic(k)==", res_pr.b$n), parse = TRUE, hjust = "right", size = 3.5)+
  labs(x = "Zr", y = "") +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5))
pheno.repro.behav

# Overall
## get estimates
res_pr.full <- get_est(pr.full1)

## prepare for plotting

### get estimates and predictions
pred_pr.full<-get_pred(pr.full1)
res_pr.full<-full_join(res_pr.full, pred_pr.full)

### plot
pheno.repro.full<-ggplot(data=res_pr.full, aes(x=estimate, y="Phenotypic-Reproduction"))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  scale_y_discrete()+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(size = 3, shape = 21, fill="black")+
  annotate("text", x=3, y=1, label="k= 239", size=3.5, hjust="right")+
  labs(x = "", y = "") +
  ggtitle("b")+
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=12, colour = "black"),
        plot.title = element_text(face="bold", size = 10, hjust = -0.1, vjust = -6))
pheno.repro.full

# Fitness

## get estimates
res_pr.f <- get_est(pr.fit, mod = "FitnessMeasured")

## prepare for plotting
### add sample size (k) for each category
k_pr.f<-pheno.r %>% group_by(FitnessMeasured) %>% count()

### get estimates and predictions
pred_pr.f<-get_pred(pr.fit, mod = "FitnessMeasured")
res_pr.f<-left_join(res_pr.f, k_pr.f, by=c("name"="FitnessMeasured")) %>% left_join(pred_pr.f)

### plot
pheno.repro.fit<-ggplot(data=res_pr.f, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  scale_y_discrete(labels=c("lab", "wild"))+
  geom_point(data = pheno.r, 
             aes(x= Zr, y = FitnessMeasured, size = ((1/VZr) + 3), colour = FitnessMeasured), alpha=0.2)+
  scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(aes(fill = name), size = 3, shape = 21)+
  scale_color_viridis(discrete = TRUE)+
  annotate('text', x = 3, y= 1:2, label= paste("italic(k)==", res_pr.f$n), parse = TRUE, hjust = "right", size = 3.5)+
  labs(x = "Zr", y = "") +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  labs(title="Phenotypic-Reproduction", subtitle = "b")+
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5))+
  theme(plot.title = element_text(size =12, hjust = 0), plot.subtitle = element_text(size=10, hjust=-0.05, vjust=-6, face="bold"))
pheno.repro.fit

# Among: survival
# Behaviour

## get estimates
res_as.b <- get_est(as.behav, mod = "Behaviour")

## prepare for plotting
### add sample size (k) for each category
k_as.b<-amg.s %>% group_by(Behaviour) %>% count()

### get estimates and predictions
pred_as.b<-get_pred(as.behav, mod = "Behaviour")
res_as.b<-left_join(res_as.b, k_as.b, by=c("name"="Behaviour")) %>% left_join(pred_as.b)

### plot
amg.surv.behav<-ggplot(data=res_as.b, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  geom_point(data = amg.s, 
             aes(x= Zr, y = Behaviour, size = ((1/VZr) + 3), colour = Behaviour), alpha=0.3)+
  scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(aes(fill = name), size = 2.7, shape = 21)+
  scale_color_viridis(discrete = TRUE)+
  annotate('text', x = 3, y= 1:9, label= paste("italic(k)==", res_as.b$n), parse = TRUE, hjust = "right", size = 3.5)+
  labs(x = "Zr", y = "") +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5))
amg.surv.behav

# Overall
## get estimates
res_as.full <- get_est(as.reduced)

## prepare for plotting

### get estimates and predictions
pred_as.full<-get_pred(as.reduced)
res_as.full<-full_join(res_as.full, pred_as.full)

### plot
amg.surv.full<-ggplot(data=res_as.full, aes(x=estimate, y="Among-Survival"))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  scale_y_discrete()+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(size = 3, shape = 21, fill="black")+
  annotate("text", x=3, y=1, label= "k= 142", size=3.5, hjust="right")+
  labs(x = "", y = "") +
  ggtitle("c")+
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=12, colour = "black"),
        plot.title = element_text(face="bold", size = 10, hjust = -0.1, vjust = -6))
amg.surv.full

# Fitness

## get estimates
res_as.f <- get_est(as.fit, mod = "FitnessMeasured")

## prepare for plotting
### add sample size (k) for each category
k_as.f<-amg.s %>% group_by(FitnessMeasured) %>% count()

### get estimates and predictions
pred_as.f<-get_pred(as.fit, mod = "FitnessMeasured")
res_as.f<-left_join(res_as.f, k_as.f, by=c("name"="FitnessMeasured")) %>% left_join(pred_as.f)

### plot
amg.surv.fit<-ggplot(data=res_as.f, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  scale_y_discrete(labels=c("lab", "wild"))+
  geom_point(data = amg.s, 
             aes(x= Zr, y = FitnessMeasured, size = ((1/VZr) + 3), colour = FitnessMeasured), alpha=0.2)+
  scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(aes(fill = name), size = 3, shape = 21)+
  scale_color_viridis(discrete = TRUE)+
  annotate('text', x = 3, y= 1:2, label= paste("italic(k)==", res_as.f$n), parse = TRUE, hjust = "right", size = 3.5)+
  labs(x = "Zr", y = "") +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  labs(title="Among-Survival", subtitle="c")+
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5))+
  theme(plot.title = element_text(size = 12, hjust = 0), plot.subtitle = element_text(size=10, hjust=-0.05, vjust=-6, face="bold"))
amg.surv.fit

# Among: repro
# Behaviour

## get estimates
res_ar.b <- get_est(ar.behav, mod = "Behaviour")

## prepare for plotting
### add sample size (k) for each category
k_ar.b<-amg.r %>% group_by(Behaviour) %>% count()

### get estimates and predictions
pred_ar.b<-get_pred(ar.behav, mod = "Behaviour")
res_ar.b<-left_join(res_ar.b, k_ar.b, by=c("name"="Behaviour")) %>% left_join(pred_ar.b)

### plot
amg.repro.behav<-ggplot(data=res_ar.b, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  geom_point(data = amg.r, 
             aes(x= Zr, y = Behaviour, size = ((1/VZr) + 3), colour = Behaviour), alpha=0.3)+
  scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(aes(fill = name), size = 3, shape = 21)+
  scale_color_viridis(discrete = TRUE)+
  annotate('text', x =3, y= 1:9, label= paste("italic(k)==", res_ar.b$n), parse = TRUE, hjust = "right", size = 3.5)+
  labs(x = "Zr", y = "") +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5))
amg.repro.behav

# Overall
## get estimates
res_ar.full <- get_est(ar.null)

## prepare for plotting

### get estimates and predictions
pred_ar.full<-get_pred(ar.null)
res_ar.full<-full_join(res_ar.full, pred_ar.full)

### plot
amg.repro.full<-ggplot(data=res_ar.full, aes(x=estimate, y="Among-Reproduction"))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  scale_y_discrete()+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(size = 3, shape = 21, fill="black")+
  annotate("text", x=3, y=1, label= "k= 144", size=3.5, hjust="right")+
  labs(x = "", y = "") +
  ggtitle("d")+
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(face="bold", size = 10, hjust = -0.1, vjust = -6))
amg.repro.full

# Fitness

## get estimates
res_ar.f <- get_est(ar.fit, mod = "FitnessMeasured")

## prepare for plotting
### add sample size (k) for each category
k_ar.f<-amg.r %>% group_by(FitnessMeasured) %>% count()

### get estimates and predictions
pred_ar.f<-get_pred(ar.fit, mod = "FitnessMeasured")
res_ar.f<-left_join(res_ar.f, k_ar.f, by=c("name"="FitnessMeasured")) %>% left_join(pred_ar.f)

### plot
amg.repro.fit<-ggplot(data=res_ar.f, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  scale_y_discrete(labels= c("lab", "wild"))+
  geom_point(data = amg.r, 
             aes(x= Zr, y = FitnessMeasured, size = ((1/VZr) + 3), colour = FitnessMeasured), alpha=0.2)+
  scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(aes(fill = name), size = 3, shape = 21)+
  scale_color_viridis(discrete = TRUE)+
  annotate('text', x = 3, y= 1:2, label= paste("italic(k)==", res_ar.f$n), parse = TRUE, hjust = "right", size = 3.5)+
  labs(x = "Zr", y = "") +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  labs(title="Among-Reproduction", subtitle="d")+
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5))+
  theme(plot.title = element_text(size=12, hjust = 0), plot.subtitle = element_text(size=10, hjust = -0.05, vjust=-6, face="bold"))
amg.repro.fit

# Within: repro
# Behaviour

## get estimates
res_wr.b <- get_est(wr.behav, mod = "Behaviour")

## prepare for plotting
### add sample size (k) for each category
k_wr.b<-wthn.r %>% group_by(Behaviour) %>% count()

### get estimates and predictions
pred_wr.b<-get_pred(wr.behav, mod = "Behaviour")
res_wr.b<-left_join(res_wr.b, k_wr.b, by=c("name"="Behaviour")) %>% left_join(pred_wr.b)

res_wr.b<-add_row(.data=res_wr.b,.before=1, name="activity",n=0)
res_wr.b<-add_row(.data=res_wr.b,.after=1, name="aggression",n=0)
res_wr.b<-add_row(.data=res_wr.b,.after=5, name="foraging",n=0)
res_wr.b<-add_row(.data=res_wr.b,.after=7, name="social behaviour",n=0)

### plot
wthn.repro.behav<-ggplot(data=res_wr.b, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  geom_point(data = wthn.r, 
             aes(x= Zr, y = Behaviour, size = ((1/VZr) + 3), colour = Behaviour), alpha=0.3)+
  scale_fill_manual(values=my_colours)+
  #scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(aes(fill = name), size = 3, shape = 21)+
  scale_color_manual(values = my_colours)+
  annotate('text', x =3, y= 1:9, label= paste("italic(k)==", res_wr.b$n), parse = TRUE, hjust = "right", size = 3.5)+
  labs(x = "Zr", y = "") +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5))
wthn.repro.behav

# Fitness

## get estimates
res_wr.f <- get_est(wr.fit, mod = "FitnessMeasured")

## prepare for plotting
### add sample size (k) for each category
k_wr.f<-wthn.r %>% group_by(FitnessMeasured) %>% count()

### get estimates and predictions
pred_wr.f<-get_pred(wr.fit, mod = "FitnessMeasured")
res_wr.f<-left_join(res_wr.f, k_wr.f, by=c("name"="FitnessMeasured")) %>% left_join(pred_wr.f)

### plot
wthn.repro.fit<-ggplot(data=res_wr.f, aes(x=estimate, y=name))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  scale_y_discrete(labels= c("lab", "wild"))+
  geom_point(data = wthn.r, 
             aes(x= Zr, y = FitnessMeasured, size = ((1/VZr) + 3), colour = FitnessMeasured), alpha=0.2)+
  scale_fill_viridis(discrete = TRUE)+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(aes(fill = name), size = 3, shape = 21)+
  scale_color_viridis(discrete = TRUE)+
  annotate('text', x = 3, y= 1:2, label= paste("italic(k)==", res_wr.f$n), parse = TRUE, hjust = "right", size = 3.5)+
  labs(x = "Zr", y = "") +
  guides(fill = "none", colour = "none") +
  theme_classic() +
  labs(title="Within-Reproduction", subtitle="e")+
  theme(legend.position="none") +
  theme(axis.text.y = element_text(size = 12, colour ="black",hjust = 0.5))+
  theme(plot.title = element_text(size=12, hjust = 0), plot.subtitle = element_text(size=10, hjust=-0.05, vjust = -6, face="bold"))
wthn.repro.fit

# Overall
## get estimates
res_wr.full <- get_est(wr.null)

## prepare for plotting

### get estimates and predictions
pred_wr.full<-get_pred(wr.null)
res_wr.full<-full_join(res_wr.full, pred_wr.full)

### plot
wthn.repro.full<-ggplot(data=res_wr.full, aes(x=estimate, y="Within-Reproduction"))+
  scale_x_continuous(limits=c(-1.5, 3), breaks = seq(-1.5, 3, by = 0.5) )+
  scale_y_discrete()+
  # 95 %predition interval (PI)
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = F, size = 0.5, alpha = 0.6) +
  # 95 %CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = F, size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.3)+
  # creating dots and different size (bee-swarm and bubbles)
  geom_point(size = 3, shape = 21, fill="black")+
  annotate("text", x=3, y=1, label= "k= 17", size=3.5, hjust="right")+
  labs(x = "", y = "") +
  ggtitle("e")+
  guides(fill = "none", colour = "none") +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(face="bold", size = 10, hjust = -0.1, vjust = -6))
wthn.repro.full

      #*** Figure 5####
fig.behav<-(pheno.surv.full/pheno.surv.behav/amg.surv.full/amg.surv.behav/plot_spacer()+
              plot_layout(heights = c(1.6, 3.8, 1.6, 3.8, 5.4)))|
  (pheno.repro.full/ pheno.repro.behav/amg.repro.full/amg.repro.behav/wthn.repro.full/wthn.repro.behav + 
     plot_layout(heights=c(1.6, 3.8, 1.6, 3.8,1.6,3.8)))
fig.behav

      #*** Figure 6####
fit.plots<-(pheno.surv.fit/amg.surv.fit/plot_spacer()+
              plot_layout(heights = c(5.4,5.4, 5.4)))|
  (pheno.repro.fit/amg.repro.fit/wthn.repro.fit + 
     plot_layout(heights=c(5.4,5.4,5.4)))
fit.plots

# 4) Publication Bias----
  #* Funnel Plot====
res_funnel_plot <- rma.mv(yi = Zr, V = VZr, mods = ~Behaviour + 
                            BehavCondition + FitnessProxy+ LevelBehav, random = list(~1|CommonName,~1|RecordNo, ~1|obs) , data = data)

funnel<-funnel(res_funnel_plot, yaxis = "seinv", level = c(90, 95, 99), 
               shade = c("white","gray55", "gray75"), refline = 0, legend = FALSE)

  #* Egger Regression====
egger.test <-rma.mv(yi = Zr, V = VZr, mods = ~ sqrt(VZr) + Behaviour + FitnessProxy + LevelBehav, random = list(~1|CommonName, ~1|RecordNo,~1|obs), data = data)
summary(egger.test)

pred_egger_regression_uni <- predict.rma(egger.test)
pred_egger_regression_uni

fit_egger_regression_uni <- data %>% mutate(ymin = pred_egger_regression_uni$ci.lb, 
                                            ymax = pred_egger_regression_uni$ci.ub, ymin2 = pred_egger_regression_uni$cr.lb, 
                                            ymax2 = pred_egger_regression_uni$cr.ub, pred = pred_egger_regression_uni$pred) %>% 
  ggplot(aes(x = sqrt(VZr), y = Zr, size = (1/VZr) + 3)) + 
  geom_point(shape = 21,fill = "grey90") +
  geom_smooth(aes(y = ymin2), method = "loess", se = FALSE, lty = "dotted", lwd = 0.25, 
              colour = "#0072B2") + 
  geom_smooth(aes(y = ymax2), method = "loess", se = FALSE, 
              lty = "dotted", lwd = 0.25, colour = "#0072B2") + 
  geom_smooth(aes(y = ymin), method = "loess", se = FALSE, lty = "dotted", lwd = 0.25, colour = "#D55E00") + 
  geom_smooth(aes(y = ymax), method = "loess", se = FALSE, lty = "dotted", 
              lwd = 0.25, colour = "#D55E00") + 
  geom_smooth(aes(y = pred), method = "loess", se = FALSE, lty = "dashed", lwd = 0.5, colour = "black") + 
  labs(x = "sqrt(sampling variance)", y = expression(paste(italic(Zr), " (effect size)")), 
       size = expression(paste(italic(n), " (# of individuals)"))) + 
  guides(fill = "none", colour = "none")+
  theme_classic() + 
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) + 
  theme(legend.direction = "horizontal", legend.background = element_blank())

fit_egger_regression_uni

  #* Time Lag Effect====
time_lag_effect_uni <- rma.mv(yi = Zr, V = VZr, mods = ~Year + Behaviour + FitnessProxy + LevelBehav, random = list(~1|CommonName,~1 | 
                                                                                                                      RecordNo, ~1|obs), data = data)
summary(time_lag_effect_uni)

pred_time_lag_effect_uni <- predict.rma(time_lag_effect_uni)

fit_time_lag_effect <- data %>% mutate(ymin = pred_time_lag_effect_uni$ci.lb, 
                                       ymax = pred_time_lag_effect_uni$ci.ub, ymin2 = pred_time_lag_effect_uni$cr.lb, 
                                       ymax2 = pred_time_lag_effect_uni$cr.ub, pred = pred_time_lag_effect_uni$pred) %>% 
  ggplot(aes(x = Year, y = Zr, size = (1/VZr) + 3)) + geom_point(shape = 21, 
                                                                 fill = "grey90")+
  geom_smooth(aes(y = ymin2), method = "loess", se = FALSE, lty = "dotted", lwd = 0.25, 
              colour = "#0072B2") + geom_smooth(aes(y = ymax2), method = "loess", se = FALSE, 
                                                lty = "dotted", lwd = 0.25, colour = "#0072B2") + geom_smooth(aes(y = ymin), 
                                                                                                              method = "loess", se = FALSE, lty = "dotted", lwd = 0.25, colour = "#D55E00") + 
  geom_smooth(aes(y = ymax), method = "loess", se = FALSE, lty = "dotted", 
              lwd = 0.25, colour = "#D55E00") + geom_smooth(aes(y = pred), method = "loess", 
                                                            se = FALSE, lty = "dashed", lwd = 0.5, colour = "black") + ylim(-1, 2) + 
  xlim(1984, 2019) + scale_x_continuous(breaks = c(1985,1990,1995, 2000, 2005, 2010, 
                                                   2015, 2020)) +
  labs(x = "Year", y = expression(paste(italic(Zr), " (effect size)")), size = expression(paste(italic(n), 
                                                                                                " (# of indiviudals)"))) + guides(fill = "none", colour = "none") + # themses
  theme_classic() + theme(legend.position = c(1, 1), legend.justification = c(1, 1)) + 
  theme(legend.direction = "horizontal") +
  theme(legend.background = element_blank()) + theme(axis.text.y = element_text(size = 10, 
                                                                                colour = "black", hjust = 0.5, angle = 90))

fit_time_lag_effect

  #* Figure 7====
pub.bias<-ggarrange(fit_egger_regression_uni,fit_time_lag_effect)
pub.bias

# 5) Partitioning (supplementary material)----

## "Level_GroupID" indicates whether estimate has a matching 
## partitioned/unpartitioned estimate
### remove estimates where Level_GroupID=0 (no matching estimate)

data_level<-subset(data, Level_GroupID!=0)

## for estimates to match, the following variables must be the same:
### Record ID
### GroupID_global
### FitnessProxy
### Behaviour
### Level_GroupID

## "LevelBehav" gives level of partitioning
### To make comparisons b/w phenotypic and among-individual estimates
### use LevelBehav ="phenotypic" and "among"
### within-individual estimates coded as LevelBehav=="within"

phenotypic<- data_level %>% 
  subset(LevelBehav=="phenotypic") %>% 
  select(RecordNo, Year, Author, GroupID_global, Class, LatinName,CommonName, Age, Sex, FitnessProxy, Behaviour, BehavDetails,Level_GroupID, Conditions,sign, SignReversed, Zr, VZr, NInd)

among<-data_level %>% 
  subset(LevelBehav=="among") %>% 
  select(RecordNo, Year, Author, GroupID_global, Class, LatinName,CommonName, FitnessProxy, Behaviour,Level_GroupID, Conditions,sign, SignReversed, Zr, VZr, NInd)

wthn<-data_level %>% 
  subset(LevelBehav=="within") %>% 
  select(RecordNo, Year, Author, GroupID_global, Class, LatinName,CommonName, FitnessProxy, Behaviour,Level_GroupID, Conditions,sign,SignReversed, Zr, VZr, NInd)

# Merge levels into a wide format
pa_wide<-inner_join(phenotypic, among, by=c("RecordNo", "Year", "Author", "GroupID_global", "Class", "CommonName","LatinName", "FitnessProxy", "Behaviour","Level_GroupID", "Conditions"), suffix=c(".pheno", ".amg"))
pw_wide<- inner_join(phenotypic, wthn, by=c("RecordNo", "Year", "Author", "GroupID_global", "Class", "CommonName","LatinName", "FitnessProxy", "Behaviour","Level_GroupID", "Conditions"), suffix=c(".pheno", ".wth"))

# Phenotypic-Among correlation
cor.test(pa_wide$Zr.pheno, pa_wide$Zr.amg) #r=0.228, 0.042-0.399

pa.cor<-ggplot(data = pa_wide, aes(Zr.pheno, Zr.amg, size= (1/VZr.pheno)+3))+
  geom_smooth(method = "lm", se=FALSE, lty="dashed", colour="black")+
  geom_point(shape= 21, fill ="grey90")+
  xlab(label="Zr phenotypic")+
  ylab(label="Zr among-individual")+
  xlim(-0.5,0.75)+
  ylim(-1.25,2)+
  labs(size= expression(paste(italic(n),"(# of individuals)")))+
  guides(size=guide_legend(override.aes = list(linetype=c(0,0,0))))+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size=14))+
  theme(legend.position=c(1,0), legend.justification = c(1,0))+
  theme(legend.direction = "horizontal")
pa.cor

# Phenotypic-within correlation
cor.test(pw_wide$Zr.pheno, pw_wide$Zr.wth) #r=0.84, 0.538-0.951

pw.cor<-ggplot(data = pw_wide, aes(Zr.pheno, Zr.wth, size= (1/VZr.pheno)+3))+
  geom_smooth(method = "lm", se=FALSE, color="black", lty="dashed")+
  geom_point(shape= 21, fill="grey90")+
  xlab(label="Zr phenotypic")+
  ylab(label="Zr within-individual")+
  xlim(-0.5,0.75)+
  ylim(-1.25,2)+
  labs(size= expression(paste(italic(n),"(# of individuals)")))+
  guides(size=guide_legend(override.aes = list(linetype=c(0,0,0))))+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size=14))+
  theme(legend.position=c(1,0), legend.justification = c(1,0))+
  theme(legend.direction = "horizontal")
pw.cor
  #* Figure S1====
part.plots<-ggarrange(pa.cor, pw.cor,
                      labels="auto", label.x = 0.04, label.y = 1,
                      font.label = list(size = 12, face = "bold"))
part.plots
