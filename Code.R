#Load packages for analysis
library(glmmTMB) #for glmm models
library(MuMIn) #for some R2 computations and AIC  
library(bbmle) #for AICs 
library(performance) #for assumptions and R2 values 
library(gllvm) #for gllvm analysis
library(FD) #for FRic, FEve and FDiv

#Load packages for plots
library(ggplot2) 
library(patchwork) 
library(tidyr)
library(dplyr)
library(gridExtra)
library(lattice)


## Generalized Linear Mixed Models (GLMMs) ##

# Load dataset (change the path as needed)
Sunbridge <- read.csv("C:/Path/To/Your/Folder/Sunbridge.csv")

#Transformation of categorical data into factors
Sunbridge$fPLOT <- factor(Sunbridge$PLOT) 
Sunbridge$fTRAP <- factor(Sunbridge$TRAP) 
Sunbridge$fC <- factor(Sunbridge$C) 
Sunbridge$fR <- factor(Sunbridge$R) 
Sunbridge$fBLOCK <- factor(Sunbridge$BLOCK) 
Sunbridge$fSHR <- factor(Sunbridge$SHR) 
Sunbridge$fGWF <- factor(Sunbridge$GWF) 
Sunbridge$fUWF <- factor(Sunbridge$UWF) 
Sunbridge$fTRE <- factor(Sunbridge$TRE) 
Sunbridge$fGRA <- factor(Sunbridge$GRA) 
Sunbridge$fMUL <- factor(Sunbridge$MUL) 
Sunbridge$fSEASON <- factor(Sunbridge$SEASON) 

#Check quality dataset 
str(Sunbridge) 

#Scale variables in different units to enable comparison in the best-fitting model
Sunbridge$sBGT <- scale(Sunbridge$BGT) 
Sunbridge$sVWC <- scale(Sunbridge$VWC) 
Sunbridge$sTCTN <- scale(Sunbridge$TCTN) 
Sunbridge$spH <- scale(Sunbridge$pH) 
Sunbridge$sBD <- scale(Sunbridge$BD) 
Sunbridge$sPOR <- scale(Sunbridge$POR) 
Sunbridge$sPRCP <- scale(Sunbridge$PRCP) 
Sunbridge$sTMIN <- scale(Sunbridge$TMIN) 
Sunbridge$sTMAX <- scale(Sunbridge$TMAX) 
Sunbridge$sTAVG <- scale(Sunbridge$TAVG) 
Sunbridge$sOM <- scale(Sunbridge$OM) 

#SPECIES RICHNESS 

#(1|fSEASON) for random effect for season 
#(1|fBLOCK/fPLOT/fTRAP) nested random effect 

null<- glmmTMB(SR ~ 1, data = Sunbridge, family=poisson) # Null 

OM_veg <- glmmTMB(SR ~ OM + fR + fSHR + fGWF + fUWF + fTRE + fGRA +(1|fSEASON) + (1|fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Vegetation only

OM_GC1 <- glmmTMB(SR ~ OM + fR  + fGWF + BGT + VWC + I(VWC^2) + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Ground cover 1: GWF + BGT + VWC 

OM_GC2 <- glmmTMB(SR ~ OM + fR  +  fMUL + BGT + VWC + I(VWC^2) + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Ground cover 2: MUL + BGT + VWC 

OM_GC3 <- glmmTMB(SR ~ OM + fR  + fGWF + fMUL + BGT + VWC + I(VWC^2) +  (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Ground cover 3: GWF + MUL + BGT + VWC 

OM_phys <- glmmTMB(SR ~ OM + fR + BD  + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Physical soil properties: BD 

OM_chem <- glmmTMB(SR ~ OM + fR + TCTN + pH + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Chemical soil properties: TC-TN + pH 

OM_GC_soil <- glmmTMB(SR ~ OM + fR  + fGWF + fMUL + BGT + VWC + I(VWC^2) + BD + TCTN + pH + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Ground cover, physical and chemical soil properties 

OM_weather <- glmmTMB(SR ~ OM + fR  + PRCP + TMIN + TMAX + (1 | fSEASON) + (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Weather 

OM_full <- glmmTMB(SR ~ OM + fR  + fSHR + fGWF + fUWF + fTRE + fGRA + fMUL + BGT + VWC + I(VWC^2) + PRCP + TAVG + + BD + TCTN + pH + (1 | fSEASON) + (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Full model 

BD_veg <- glmmTMB(SR ~ BD + fR  + fSHR + fGWF + fUWF + fTRE + fGRA +  (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Vegetation only 

BD_GC1 <- glmmTMB(SR ~ sBD + fR  + fGWF + sBGT + sVWC + I(sVWC^2) + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Ground cover 1: GWF + BGT + VWC 

BD_GC2 <- glmmTMB(SR ~ BD + fR  +  fMUL + BGT + VWC + I(VWC^2) + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Ground cover 2: MUL + BGT + VWC 

BD_GC3 <- glmmTMB(SR ~ BD + fR  + fGWF + fMUL + BGT + VWC + I(VWC^2) +  (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Ground cover 3: GWF + MUL + BGT + VWC 

BD_chem <- glmmTMB(SR ~ BD + fR + TCTN + pH + OM + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Chemical soil properties: TC-TN + pH + OM 

BD_GC_soil <- glmmTMB(SR ~ BD + fR  + fGWF + fMUL + BGT + VWC + I(VWC^2)  + TCTN + pH + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Ground cover, physical and chemical soil properties 

BD_weather <- glmmTMB(SR ~ BD + fR  + PRCP + TMIN + TMAX + (1 | fSEASON) + (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Weather 

BD_full <- glmmTMB(SR ~ BD + fR   + fSHR + fGWF + fUWF + fTRE + fGRA + fMUL + BGT + VWC + I(VWC^2) + PRCP + TAVG + TCTN + pH + (1 | fSEASON) + (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge, family=poisson) # Full model 

#AIC model selection 
AICctab(null, OM_veg, OM_GC1, OM_GC2, OM_GC3, OM_phys, OM_chem, OM_GC_soil, OM_weather, OM_full, BD_veg, BD_GC1, BD_GC2, BD_GC3, BD_chem, BD_GC_soil, BD_weather, BD_full, weights=T, base = T) 

summary() 
check_model() 
r.squaredGLMM() 

#BEETLE ABUNDANCE 

#The same set of hypotheses used to construct the GLMMs for SR were applied to the models for abundance

#(1|fSEASON) for random effect for season 
#(1|fBLOCK/fPLOT/fTRAP) nested random effect 

Sunbridge$logCOUNT<- log(COUNT+1) 

null<- glmmTMB(logCOUNT ~ 1, data = Sunbridge) 

OM_veg <- glmmTMB(logCOUNT ~ OM + fR  + fSHR + fGWF + fUWF + fTRE + fGRA +  (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

OM_GC1 <- glmmTMB(logCOUNT ~ OM + fR  + fGWF + BGT + VWC + I(VWC^2) + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

OM_GC2 <- glmmTMB(logCOUNT ~ OM + fR  +  fMUL + BGT + VWC + I(VWC^2) + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

OM_GC3 <- glmmTMB(logCOUNT ~ OM + fR  + fGWF + fMUL + BGT + VWC + I(VWC^2) +  (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

OM_phys <- glmmTMB(logCOUNT ~ OM + fR + BD  + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

OM_chem <- glmmTMB(logCOUNT ~ OM + fR + TCTN + pH + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge) 

OM_GC_soil <- glmmTMB(logCOUNT ~ OM + fR  + fGWF + fMUL + BGT + VWC + I(VWC^2) + BD + TCTN + pH + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge) 

OM_weather <- glmmTMB(logCOUNT ~ OM + fR  + PRCP + TMIN + TMAX + (1 | fSEASON) + (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

OM_full <- glmmTMB(logCOUNT ~ OM + fR   + fSHR + fGWF + fUWF + fTRE + fGRA + fMUL + BGT + VWC + I(VWC^2) + PRCP + TAVG + BD + TCTN + pH + (1 | fSEASON) + (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

BD_veg <- glmmTMB(logCOUNT ~ BD + fR  + fSHR + fGWF + fUWF + fTRE + fGRA +  (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

BD_GC1 <- glmmTMB(logCOUNT ~ BD + fR  + fGWF + BGT + VWC + I(VWC^2) + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge) 

BD_GC2 <- glmmTMB(logCOUNT ~ BD + fR  +  fMUL + BGT + VWC + I(VWC^2) + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge) 

BD_GC3 <- glmmTMB(logCOUNT ~ BD + fR  + fGWF + fMUL + BGT + VWC + I(VWC^2) +  (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

BD_chem <- glmmTMB(logCOUNT ~ BD + fR + TCTN + pH + OM + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

BD_GC_soil <- glmmTMB(logCOUNT ~ BD + fR  + fGWF + fMUL + BGT + VWC + I(VWC^2)  + TCTN + pH + (1 | fSEASON)+ (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

BD_weather <- glmmTMB(logCOUNT ~ BD + fR  + PRCP + TMIN + TMAX + (1 | fSEASON) + (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge)  

BD_full <- glmmTMB(logCOUNT ~ BD + fR   + fSHR + fGWF + fUWF + fTRE + fGRA + fMUL + BGT + VWC + I(VWC^2) + PRCP + TAVG + TCTN + pH + (1 | fSEASON) + (1 |fBLOCK/fPLOT/fTRAP), data=Sunbridge) 

#AIC model selection 
AICctab(null, OM_veg, OM_GC1, OM_GC2, OM_GC3, OM_phys, OM_chem, OM_GC_soil, OM_weather, OM_full, BD_veg, BD_GC1, BD_GC2, BD_GC3, BD_chem, BD_GC_soil, BD_weather, BD_full, weights=T, base = T) 

summary() 
check_model() 
r.squaredGLMM() 

#Plot data

#Set up the 2x3 plot grid
par(mfrow = c(2, 3))  

#Plot for BGT and species richness
par(mar = c(5, 5, 4, 2))  
plot(Sunbridge$BGT, Sunbridge$SR, 
     main = "", 
     xlab = "BGT", 
     ylab = "Species Richness",
     pch = 19, col = "slategray",
     cex.axis = 1.2,  
     cex.lab = 1.4)   
abline(lm(SR ~ BGT, data = Sunbridge), col = "darkred", lwd = 2)

#Plot for VWC (with both linear and quadratic regression) and species richness
plot(Sunbridge$VWC, Sunbridge$SR, 
     main = "", 
     xlab = "VWC", 
     ylab = "Species Richness",
     pch = 19, col = "slategray",
     cex.axis = 1.2,  
     cex.lab = 1.4)  
#Quadratic regression line
lines(sort(Sunbridge$VWC), predict(lm(SR ~ poly(VWC, 2), data = Sunbridge), 
                                   newdata = data.frame(VWC = sort(Sunbridge$VWC))), 
      col = "darkred", lwd = 2)

#Boxplot for GWF and species richness
boxplot(SR ~ GWF, data = Sunbridge, 
        main = "", 
        xlab = "GWF", 
        ylab = "Species Richness",
        cex.axis = 1.2,  
        cex.lab = 1.4,   
        col = "white")

#Plot for BGT and abundance
plot(Sunbridge$BGT, Sunbridge$COUNT, 
     main = "", 
     xlab = "BGT", 
     ylab = "Abundance",
     pch = 19, col = "slategray",
     cex.axis = 1.2,  
     cex.lab = 1.4)   
abline(lm(COUNT ~ BGT, data = Sunbridge), col = "darkred", lwd = 2)

#Plot for VWC (with both linear and quadratic regression) and abundance
plot(Sunbridge$VWC, Sunbridge$COUNT, 
     main = "", 
     xlab = "VWC", 
     ylab = "Abundance",
     pch = 19, col = "slategray",
     cex.axis = 1.2, 
     cex.lab = 1.4)  
#Quadratic regression line
lines(sort(Sunbridge$VWC), predict(lm(COUNT ~ poly(VWC, 2), data = Sunbridge), 
                                   newdata = data.frame(VWC = sort(Sunbridge$VWC))), 
      col = "darkred", lwd = 2)

#Boxplot for GWF and abundance
boxplot(COUNT ~ MUL, data = Sunbridge, 
        main = "", 
        xlab = "MUL", 
        ylab = "Abundance",
        cex.axis = 1.2,  
        cex.lab = 1.4,    
        col = "white")


## Generalized linear latent variable model (GLLVM) ##

#Following steps from Jenni Niku 2024 gllvm 1.4.7 

#Load dataset (change the path as needed)
sp <- read.csv("C:/Path/To/Your/Folder/sp.csv")
env <- read.csv("C:/Path/To/Your/Folder/env.csv")
traits <- read.csv("C:/Path/To/Your/Folder/traits.csv")

library(gllvm)
library(lattice)

#rownames(traits) <- traits[, 1] 
#traits <- traits[, -1] # Remove the first column after setting row names 

#Test which family to use 

y <- as.matrix(sp)  
X <- as.matrix(env)  
TR <- traits  

#Plot residuals for the poisson model 
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1)) 
plot(fitp, var.colors = 1) 

#Plot residuals for the negative binomial model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(fit_ord, var.colors = 1)

#Construct an ordination as a scatter plot of the predicted latent variables via the ordiplot function
par(mfrow = c(1, 2))

ordiplot(fit_ord, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3),
         main = "Biplot")

ordiplot(fit_ord, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3),
         main = "Ordination plot", predict.region = TRUE)

#Model selection for num.lv
fitx1 <- gllvm(y = sp, X = env, family = "negative.binomial", num.lv = 1)
fitx2 <- gllvm(y = sp, X = env, family = "negative.binomial", num.lv = 2)
fitx3 <- gllvm(y = sp, X = env, family = "negative.binomial", num.lv = 3)

AIC(fitx1)
AIC(fitx2)
AIC(fitx3) 

#Residual Analysis for poisson and negative binomial
par(mfrow = c(1,2))
plot(fitx1, which = 1:2)

fitp <- gllvm(y, family = poisson())  
fitp    

fit_ord <- gllvm(y, family = "negative.binomial")   
fit_ord   

#Plot residuals for the negative.binomial 
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))  
plot(fit_ord, var.colors = 1) 

#Following residual analysis, a negative binomial distribution with one latent variable was chosen for all GLLVMs 

#GLLVM model selection

#Ground cover classes were converted into numerical values. Other categorical variables were converted into factors
substitutions <- c("1" = 0, "2" = 3, "3" = 15.5, "4" = 38, "5" = 63, "6" = 85.5, "7" = 97.5) 
env$GWF <- substitutions[env$GWF] 
env$MUL <- substitutions[env$MUL] 
env$SEASON <- as.factor(env$SEASON) 
substitutions <- c("1" = "SP", "2" = "SU", "3" = "FA") 
env$SEASON <- substitutions[env$SEASON] 
traits$Diet <- as.factor(traits$Diet ) 
traits$Wings <- as.factor(traits$Wings ) 
traits$Larval_Substrate <- as.factor(traits$Larval_Substrate ) 
traits$Body_Size <- as.numeric(traits$Body_Size)  

#Check quality of dataset 
str(traits)  
str(env)  

y <- as.matrix(sp) 
X <- as.matrix(env) 
TR <- traits 

#Alternative GLLVMs were informed by GLMM model results and varied in their inclusion of environmental predictors and functional traits 

#Separate environmental and trait models were built to test whether including the interaction between environmental predictors and traits improves model performance
gllvm1a <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ BGT, seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm1b <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,
                 formula = y ~ VWC, seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm1c <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ GWF, seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm1d <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ MUL, seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm2a <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ (BGT) : (Diet + Wings + Body_Size + Larval_Substrate), seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm2b <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,
                 formula = y ~ (VWC) : (Diet + Wings + Body_Size + Larval_Substrate), seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm2c <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ (GWF) : (Diet + Wings + Body_Size + Larval_Substrate), seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm2d <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ (MUL) : (Diet + Wings + Body_Size + Larval_Substrate), seed = 123,
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

#Likelihood ratio test: test whether including traits improves the model fit
anova(gllvm1a, gllvm2a) 
anova(gllvm1b, gllvm2b) 
anova(gllvm1c, gllvm2c) 
anova(gllvm1d, gllvm2d) 

#Since the traits improved the model, I run model with interaction env/traits, and run AIC model selection 
gllvm01 <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ (GWF + MUL) : (Diet + Wings + Body_Size + Larval_Substrate), seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm02 <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ (BGT + VWC) : (Diet + Wings + Body_Size + Larval_Substrate), seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm03 <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ (GWF + BGT) : (Diet + Wings + Body_Size + Larval_Substrate), seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm04 <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ (GWF + VWC) : (Diet + Wings + Body_Size + Larval_Substrate), seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON)

gllvm05 <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,
                 formula = y ~ (MUL + BGT) : (Diet + Wings + Body_Size + Larval_Substrate), seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm06 <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,
                 formula = y ~ (MUL + VWC) : (Diet + Wings + Body_Size + Larval_Substrate), seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm07 <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ (GWF + BGT + VWC) : (Diet + Wings + Body_Size + Larval_Substrate), 
                 seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm08 <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ (MUL + BGT + VWC) : (Diet + Wings + Body_Size + Larval_Substrate), 
                 seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

gllvm09 <- gllvm(y, env, TR, family = "negative.binomial", num.lv = 1,  
                 formula = y ~ (GWF + MUL + BGT + VWC) : (Diet + Wings + Body_Size + Larval_Substrate), 
                 seed = 123, 
                 row.eff = "random",  
                 control.start = list(n.init = 3, jitter.var = 0.01), 
                 randomX = ~ SEASON) 

#AIC model selection 
AIC_values <- AIC(gllvm01, gllvm02, gllvm03, gllvm04, gllvm05, gllvm06, gllvm07, gllvm08, gllvm09) 

AIC_values[order(AIC_values$AIC), ]  # Sort by AIC in ascending order 

summary () 


#Plot for full model 
coefplot.gllvm(gllvm06, y.label = TRUE, mar = c(4,15, 1, 1), cex.ylab = 0.8) 

fourth <- gllvm09$fourth.corner
a <- 4
b <- 4
colort <- colorRampPalette(c("steelblue", "white", "salmon")) #try steelblue and salmon
plot.5th <- levelplot((as.matrix(fourth)), xlab = "Environmental Variables", 
                      ylab = "Species traits", col.regions = colort(100), cex.lab = 1.3, 
                      at = seq(-a, b, length = 100), scales = list(x = list(rot = 45)))
plot.5th    

#Plot data for functional traits abundances

#Note: add site column (TRAP1-TRAP309) in sp.csv, and 'Species' as column name in traits.csv

#Add Site column or use rownames
species <- sp %>% 
  mutate(TRAP = rownames(.))

#Keep only numeric species columns
abund_long <- species %>%
  select(where(is.numeric), TRAP) %>%
  pivot_longer(cols = -TRAP, names_to = "Species", values_to = "Abundance") %>%
  filter(Abundance > 0)

combined <- abund_long %>%
  left_join(traits, by = "Species")

p1 <-ggplot(combined, aes(x = Body_Size, weight = Abundance)) +
  geom_density(fill = "slategray", alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "",
       x = "Body Size",
       y = "Density")

p2<- ggplot(combined, aes(x = Larval_Substrate, y = Abundance)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, color = "slategray") +
  theme_minimal() +
  labs(title = "",
       x = "Larval Substrate",
       y = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3<-ggplot(combined, aes(x = Diet, y = Abundance)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, color = "slategray") +
  theme_minimal() +
  labs(title = "",
       x = "Diet",
       y = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p4<-ggplot(combined, aes(x = Wings, y = Abundance)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, color = "slategray") +
  theme_minimal() +
  labs(title = "",
       x = "Wings",
       y = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(p1, p4, p2, p3, nrow = 2, ncol = 2)


## Functional Diversity (FRic, FEve, FDiv) ##

#rownames(traits) <- traits[, 1] 
#traits <- traits[, -1]  #Remove the first column after setting row names 

#Obtain values FRic, FEve, FDiv 
library(FD) 

ex1 <- dbFD(traits, sp) 
ex1 

#ex1 can be downloaded as csv file 

#Combine parts into a single data frame 
df_export <- data.frame(FRic = ex1$FRic, FEve = ex1$FEve, FDiv = ex1$FDiv) 
write.csv(df_export, "FD.csv", row.names = FALSE) 

#GLMM selection on FRic, FEve and FDiv 

#Transformation of categorical variables into factors
Sunbridge$fPLOT <- factor(Sunbridge$PLOT) 
Sunbridge$fTRAP <- factor(Sunbridge$TRAP) 
Sunbridge$fC <- factor(Sunbridge$C) 
Sunbridge$fR <- factor(Sunbridge$R) 
Sunbridge$fBLOCK <- factor(Sunbridge$BLOCK) 
Sunbridge$fSHR <- factor(Sunbridge$SHR) 
Sunbridge$fGWF <- factor(Sunbridge$GWF) 
Sunbridge$fUWF <- factor(Sunbridge$UWF) 
Sunbridge$fTRE <- factor(Sunbridge$TRE) 
Sunbridge$fGRA <- factor(Sunbridge$GRA) 
Sunbridge$fMUL <- factor(Sunbridge$MUL) 
Sunbridge$fSEASON <- factor(Sunbridge$SEASON) 

#The same set of hypotheses used to construct the GLMMs for species richness and beetle abundance were applied here

#Substitute response variable with FRic, FEve and FDiv

null <- glmmTMB(FRic ~ 1, data = Sunbridge, family = beta_family) 

OM_veg <- glmmTMB(FRic ~ OM + fR + fSHR + fGWF + fUWF + fTRE + fGRA + (1 | fSEASON), data = Sunbridge, family = beta_family) 

OM_GC1 <- glmmTMB(FRic ~ OM + fR + fGWF + BGT + VWC + I(VWC^2) + (1 | fSEASON), data = Sunbridge, family = beta_family) 

OM_GC2 <- glmmTMB(FRic ~ OM + fR + fMUL + BGT + VWC + I(VWC^2) + (1 | fSEASON), data = Sunbridge, family = beta_family) 

OM_GC3 <- glmmTMB(FRic ~ OM + fR + fGWF + fMUL + BGT + VWC + I(VWC^2) + (1 | fSEASON), data = Sunbridge, family = beta_family)

OM_phys <- glmmTMB(FRic ~ OM + fR + BD + (1 | fSEASON), data = Sunbridge, family = beta_family) 

OM_chem <- glmmTMB(FRic ~ OM + fR + TCTN + pH + (1 | fSEASON), data = Sunbridge, family = beta_family) 

OM_GC_soil <- glmmTMB(FRic ~ OM + fR + fGWF + fMUL + BGT + VWC + I(VWC^2) + BD + TCTN + pH + (1 | fSEASON), data = Sunbridge, family = beta_family) 

OM_weather <- glmmTMB(FRic ~ OM + fR + PRCP + TMIN + TMAX + (1 | fSEASON), data = Sunbridge, family = beta_family)

OM_full <- glmmTMB(FRic ~ OM + fR + fSHR + fGWF + fUWF + fTRE + fGRA + fMUL + BGT + VWC + I(VWC^2) + PRCP + TAVG + BD + TCTN + pH + (1 | fSEASON), data = Sunbridge, family = beta_family) 

BD_veg <- glmmTMB(FRic ~ BD + fR + fSHR + fGWF + fUWF + fTRE + fGRA + (1 | fSEASON), data = Sunbridge, family = beta_family) 

BD_GC1 <- glmmTMB(FRic ~ BD + fR + fGWF + BGT + VWC + I(VWC^2) + (1 | fSEASON), data = Sunbridge, family = beta_family) 

BD_GC2 <- glmmTMB(FRic ~ BD + fR + fMUL + BGT + VWC + I(VWC^2) + (1 | fSEASON), data = Sunbridge, family = beta_family) 

BD_GC3 <- glmmTMB(FRic ~ BD + fR + fGWF + fMUL + BGT + VWC + I(VWC^2) + (1 | fSEASON), data = Sunbridge, family = beta_family) 

BD_chem <- glmmTMB(FRic ~ BD + fR + TCTN + pH + OM + (1 | fSEASON), data = Sunbridge, family = beta_family) 

BD_GC_soil <- glmmTMB(FRic ~ BD + fR + fGWF + fMUL + BGT + VWC + I(VWC^2) + TCTN + pH + (1 | fSEASON), data = Sunbridge, family = beta_family) 

BD_weather <- glmmTMB(FRic ~ BD + fR + PRCP + TMIN + TMAX + (1 | fSEASON), data = Sunbridge, family = beta_family) 

BD_full <- glmmTMB(FRic ~ BD + fR + fSHR + fGWF + fUWF + fTRE + fGRA + fMUL + BGT + VWC + I(VWC^2) + PRCP + TAVG + TCTN + pH + (1 | fSEASON), data = Sunbridge, family = beta_family) 

#AIC model selection
AICctab(null, OM_veg, OM_GC1, OM_GC2, OM_GC3, OM_phys, OM_chem, OM_GC_soil, OM_weather, OM_full, BD_veg, BD_GC1, BD_GC2, BD_GC3, BD_chem, BD_GC_soil, BD_weather, BD_full, 
        weights = TRUE, base = TRUE)  

summary() 
check_model() 
r.squaredGLMM() 

#Plot Species Richness (SR) and Functional diversity metrics

#Create the data frames
df1 <- data.frame(SR = SR, FRic = FRic)
df2 <- data.frame(SR = SR, FEve = FEve)
df3 <- data.frame(SR = SR, FDiv = FDiv)

#Create individual plots
p1 <- ggplot(df1, aes(x = SR, y = FRic)) +
  geom_point(size = 2, color = "slategray") +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed", color = "darkred") +
  labs(title = "",
       x = "Species Richness",
       y = "Functional Richness (FRic)") +
  theme_minimal()

p2 <- ggplot(df2, aes(x = SR, y = FEve)) +
  geom_point(size = 2, color = "slategray") +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed", color = "darkred") +
  labs(title = "",
       x = "Species Richness",
       y = "Functional Evenness (FEve)") +
  theme_minimal()

p3 <- ggplot(df3, aes(x = SR, y = FDiv)) +
  geom_point(size = 2, color = "slategray") +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed", color = "darkred") +
  labs(title = "",
       x = "Species Richness",
       y = "Functional Divergence (FDiv)") +
  theme_minimal()

p1 + p2 + p3