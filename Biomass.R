#### Load data ####
biomass <- read.csv("~/Desktop/Projects/QUB/Data/Decomposition/Biomass.csv")

require(tidyverse)
biomass <- biomass %>% mutate(BB = Blotted/Buoyant) # compute blotted to buoyant ratio
biomass <- biomass %>% mutate(LB = Lyophilised/Blotted) # compute lyophilised to blotted ratio

ggplot(biomass, aes(Buoyant, Blotted)) + 
  geom_point() + 
  facet_grid(~Species) # looks linear

ggplot(biomass, aes(Buoyant, BB)) + 
  geom_point() + 
  facet_grid(~Species) 
# remove outlier for S. latissima where blotted-buoyant ratio > 1 (logically impossible)

biomass[biomass$BB > 1,] # row 165
biomass <- biomass[-c(165),] # remove row 165
rownames(biomass) <-  NULL # reset rownames

ggplot(biomass, aes(Blotted, Lyophilised)) + 
  geom_point() + 
  facet_grid(~Species) # looks non-linear

ggplot(biomass, aes(Blotted, LB)) + 
  geom_point() + 
  facet_grid(~Species) # looks non-linear

# no other clear outliers

#### Buoyant vs. Blotted ####
m1 <- lm(Blotted ~ 0 + Buoyant * Species,
         data = biomass) # intercept set to zero because that is the only logical intercept
boxplot(resid(m1) ~ biomass[!is.na(biomass$Buoyant),]$Species) # homogenous
plot(resid(m1) ~ biomass[!is.na(biomass$Buoyant),]$Buoyant)
abline(0,0) # nonlinearity
hist(resid(m1)) # somewhat left-skewed

biomass$Species <- factor(biomass$Species)
m2 <- nls(Blotted ~ k[Species]*Buoyant^n[Species], 
          start = list(k = c(0.8, 0.8, 0.8, 0.8),
                       n = c(1, 1, 1, 1)),
          data = biomass)
boxplot(resid(m2) ~ biomass[!is.na(biomass$Buoyant),]$Species) # more heterogenous
plot(resid(m2) ~ biomass[!is.na(biomass$Buoyant),]$Buoyant) 
abline(0,0) # nonlinearity solved
hist(resid(m2)) # more left-skewed

# try solving heterogeneity
summary(m2) # get new starting values
require(nlme)
m3 <- gnls(Blotted ~ k*Buoyant^n,
           start = list(k = c(0.798306, 0.732038, 0.847876, 0.816139),
                        n = c(1.044268, 1.067880, 1.035529, 1.060429)),
           params = list(k ~ Species, n ~ Species),
           weights = varConstPower(),
           data = biomass[!is.na(biomass$Buoyant),1:4])

boxplot(resid(m3,  type = "normalized") ~ biomass[!is.na(biomass$Buoyant),]$Species) # homogeneity improved
plot(resid(m3,  type = "normalized") ~ biomass[!is.na(biomass$Buoyant),]$Buoyant) 
abline(0,0) # homogeneity improved
hist(resid(m3,  type = "normalized")) # almost perfect normality
# m3 is chosen as the optimal model

summary(m3)
#                                      Value  Std.Error  t-value p-value
# k.(Intercept)                    0.7490144 0.01134767 66.00601  0.0000
# k.SpeciesLaminaria hyperborea   -0.0416514 0.01705768 -2.44180  0.0153
# k.SpeciesSaccharina latissima    0.0835219 0.01629920  5.12429  0.0000
# k.SpeciesSaccorhiza polyschides  0.1010537 0.01655213  6.10518  0.0000
# n.(Intercept)                    1.0919480 0.01143070 95.52765  0.0000
# n.SpeciesLaminaria hyperborea   -0.0029331 0.01615849 -0.18152  0.8561
# n.SpeciesSaccharina latissima   -0.0414320 0.01570883 -2.63750  0.0088
# n.SpeciesSaccorhiza polyschides -0.0582685 0.01572485 -3.70550  0.0003

coef(m3)
# y = 0.749014421*x^1.091948012

biomass$Species <- factor(biomass$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                                      "Laminaria digitata", "Saccharina latissima"))
m3 <- gnls(Blotted ~ k*Buoyant^n,
           start = list(k = c(0.732038, 0.816139, 0.798306, 0.847959),
                        n = c(1.067880, 1.060429, 1.044268, 1.035479)),
           params = list(k ~ Species, n ~ Species),
           weights = varConstPower(),
           data = biomass[!is.na(biomass$Buoyant),1:4])

summary(m3)
#                                      Value  Std.Error  t-value p-value
# k.(Intercept)                    0.7073603 0.01273550 55.54240  0.0000
# k.SpeciesSaccorhiza polyschides  0.1427087 0.01753263  8.13961  0.0000
# k.SpeciesLaminaria digitata      0.0416510 0.01705755  2.44179  0.0153
# k.SpeciesSaccharina latissima    0.1251753 0.01729406  7.23805  0.0000
# n.(Intercept)                    1.0890180 0.01142096 95.35254  0.0000
# n.SpeciesSaccorhiza polyschides -0.0553394 0.01571777 -3.52082  0.0005
# n.SpeciesLaminaria digitata      0.0029352 0.01615863  0.18165  0.8560
# n.SpeciesSaccharina latissima   -0.0385007 0.01570178 -2.45199  0.0148

coef(m3)
# y = 0.707360349*x^1.089017974

biomass$Species <- factor(biomass$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                                      "Saccharina latissima", "Laminaria hyperborea"))
m3 <- gnls(Blotted ~ k*Buoyant^n,
           start = list(k = c(0.816139, 0.798306, 0.847959, 0.732038),
                        n = c(1.060429, 1.044268, 1.035479, 1.067880)),
           params = list(k ~ Species, n ~ Species),
           weights = varConstPower(),
           data = biomass[!is.na(biomass$Buoyant),1:4])

summary(m3)
#                                    Value  Std.Error  t-value p-value
# k.(Intercept)                  0.8500695 0.01204984 70.54611  0.0000
# k.SpeciesLaminaria digitata   -0.1010596 0.01655187 -6.10563  0.0000
# k.SpeciesSaccharina latissima -0.0175341 0.01679550 -1.04398  0.2974
# k.SpeciesLaminaria hyperborea -0.1427104 0.01753257 -8.13973  0.0000
# n.(Intercept)                  1.0336782 0.01079860 95.72333  0.0000
# n.SpeciesLaminaria digitata    0.0582774 0.01572493  3.70605  0.0003
# n.SpeciesSaccharina latissima  0.0168398 0.01525508  1.10388  0.2706
# n.SpeciesLaminaria hyperborea  0.0553412 0.01571781  3.52092  0.0005

coef(m3)
# y = 0.85006949*x^1.03367817

biomass$Species <- factor(biomass$Species, levels = c("Saccharina latissima", "Laminaria digitata", 
                                                      "Laminaria hyperborea", "Saccorhiza polyschides"))
m3 <- gnls(Blotted ~ k*Buoyant^n,
           start = list(k = c(0.847959, 0.798306, 0.732038, 0.816139),
                        n = c(1.035479, 1.044268, 1.067880, 1.060429)),
           params = list(k ~ Species, n ~ Species),
           weights = varConstPower(),
           data = biomass[!is.na(biomass$Buoyant),1:4])

summary(m3)
#                                      Value  Std.Error  t-value p-value
# k.(Intercept)                    0.8325358 0.01170008 71.15642  0.0000
# k.SpeciesLaminaria digitata     -0.0835239 0.01629905 -5.12447  0.0000
# k.SpeciesLaminaria hyperborea   -0.1251750 0.01729408 -7.23803  0.0000
# k.SpeciesSaccorhiza polyschides  0.0175332 0.01679561  1.04391  0.2975
# n.(Intercept)                    1.0505171 0.01077531 97.49302  0.0000
# n.SpeciesLaminaria digitata      0.0414352 0.01570890  2.63769  0.0088
# n.SpeciesLaminaria hyperborea    0.0385004 0.01570176  2.45198  0.0148
# n.SpeciesSaccorhiza polyschides -0.0168384 0.01525507 -1.10379  0.2707

coef(m3)
# y = 0.83253578*x^1.05051710

biomass$Species <- factor(biomass$Species, levels = c("Laminaria digitata", "Laminaria hyperborea", 
                                                      "Saccharina latissima", "Saccorhiza polyschides"))
m3 <- gnls(Blotted ~ k*Buoyant^n,
           start = list(k = c(0.798306, 0.732038, 0.847959, 0.816139),
                        n = c(1.044268, 1.067880, 1.035479, 1.060429)),
           params = list(k ~ Species, n ~ Species),
           weights = varConstPower(),
           data = biomass[!is.na(biomass$Buoyant),1:4])

# the maximal buoyant mass measured in the decomposition experiment was 47.744 g
# although m3 fits best in the specified range due to nonlinearity at low mass,
# the data only cover 0.051 to 16.661 g buoyant mass
# hence, the robustness of m3's predictions need to be tested outside this range 
new <- data.frame(Buoyant = rep(seq(0, 48, 0.1), 4), # create prediction dataframe
                  Species = c(rep("Laminaria digitata", 481),
                              rep("Laminaria hyperborea", 481),
                              rep("Saccharina latissima", 481),
                              rep("Saccorhiza polyschides", 481))) 

new$fit <- predict(m3, newdata = new)

ggplot() +
  geom_abline(intercept = 0, slope = 1, colour = "#c9d2d7") +
  geom_line(data = new, aes(Buoyant, fit, colour = Species)) +
  ylab("Blotted") +
  theme_minimal()
# clearly m3 is not suitable for prediction beyond the data range since 
# L. digitata and S. latissima buoyant mass is predicted to be greater than
# blotted mass (a logical impossibility) above ~20 g and ~30 g buoyant mass

# hence a linear model is required
# perhaps m1 can be optimised
m1 <- gls(Blotted ~ 0 + Buoyant * Species,
          weights = varExp(form = ~Buoyant|Species),
          data = biomass[!is.na(biomass$Buoyant),1:4]) # intercept set to zero because that is the only logical intercept
boxplot(resid(m1, type = "normalized") ~ biomass[!is.na(biomass$Buoyant),1:4]$Species) # heterogeneity introduced
plot(resid(m1, type = "normalized") ~ biomass[!is.na(biomass$Buoyant),1:4]$Buoyant)
abline(0,0) # homogeneity improved at lower end but fit is less good for large values, which is contrary to the goal
hist(resid(m1)) # normality improved

# the simple linear model is better for the purpose of prediction beyond the data range
m1 <- lm(Blotted ~ 0 + Buoyant * Species,
         data = biomass)

coef(m1) # L. digitata: y = 0.89487955*x
biomass$Species <- fct_shift(biomass$Species)
m1 <- lm(Blotted ~ 0 + Buoyant * Species,
         data = biomass)
coef(m1) # L. hyperborea: y = 0.87910760*x
biomass$Species <- fct_shift(biomass$Species)
m1 <- lm(Blotted ~ 0 + Buoyant * Species,
         data = biomass)
coef(m1) # S. latissima: y = 0.923208693*x
biomass$Species <- fct_shift(biomass$Species)
m1 <- lm(Blotted ~ 0 + Buoyant * Species,
         data = biomass)
coef(m1) # S. polyschides: y = 0.932246089*x
biomass$Species <- fct_shift(biomass$Species)
m1 <- lm(Blotted ~ 0 + Buoyant * Species,
         data = biomass)

#### Blotted vs. Lyophilised ####
m4 <- lm(Lyophilised ~ 0 + Blotted * Species,
         data = biomass)
boxplot(resid(m4) ~ biomass[!is.na(biomass$Lyophilised),]$Species) # fairly homogenous but unbalanced
plot(resid(m4) ~ biomass[!is.na(biomass$Lyophilised),]$Blotted) 
abline(0,0) # clear nonlinearity and heterogenous
hist(resid(m4)) # somewhat right-skewed

m5 <- nls(Lyophilised ~ k[Species]*Blotted^n[Species], 
          start = list(k = c(0.2, 0.2, 0.2, 0.2),
                       n = c(1, 1, 1, 1)),
          data = biomass)
boxplot(resid(m5) ~ biomass[!is.na(biomass$Lyophilised),]$Species) # more balanced
plot(resid(m5) ~ biomass[!is.na(biomass$Lyophilised),]$Blotted) 
abline(0,0) # nonlinearity solved but still heterogenous
hist(resid(m5)) # more normal

# try solving heterogeneity
summary(m5) # get new starting values
m6 <- gnls(Lyophilised ~ k*Blotted^n,
           start = list(k = c(0.176742, 0.224231, 0.039560, 0.067314),
                        n = c(0.958811, 1.047887, 1.510932, 1.136141)),
           params = list(k ~ Species, n ~ Species),
           weights = varConstPower(),
           data = biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)])

boxplot(resid(m6, type = "normalized") ~ biomass[!is.na(biomass$Lyophilised),]$Species) # less homogenous
plot(resid(m6, type = "normalized") ~ biomass[!is.na(biomass$Lyophilised),]$Blotted) 
abline(0,0) # less balanced
hist(resid(m6, type = "normalized")) # more right-skew

# try exponential fit
m7 <- nls(Lyophilised ~ I[Species]*exp(k[Species]*Blotted),
          start = list(I = c(0.15, 0.15, 0.15, 0.15),
                       k = c(0.1, 0.1, 0.1, 0.1)),
          data = biomass)
boxplot(resid(m7) ~ biomass[!is.na(biomass$Lyophilised),]$Species) # more homogenous
plot(resid(m7) ~ biomass[!is.na(biomass$Lyophilised),]$Blotted) 
abline(0,0) # nonlinearity introduced
hist(resid(m7)) # less normal

summary(m7) # get new starting values
m8 <- gnls(Lyophilised ~ I*exp(k*Blotted),
           start = list(I = c(0.554338, 0.825246, 0.461785, 0.347630),
                        k = c(0.066531, 0.074699, 0.089432, 0.063225)),
           params = list(I ~ Species, k ~ Species),
           weights = varConstPower(),
           data = biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)])

boxplot(resid(m8, type = "normalized") ~ biomass[!is.na(biomass$Lyophilised),]$Species) # more homogenous
plot(resid(m8, type = "normalized") ~ biomass[!is.na(biomass$Lyophilised),]$Blotted) 
abline(0,0) # nonlinearity worse
hist(resid(m8, type = "normalized")) # less normal
# the exponential equation does not fit the data

# m5 is chosen as the optimal model but needs to be rerun as gnls object to compute 95% CI
m5 <- gnls(Lyophilised ~ k*Blotted^n,
           start = list(k = c(0.176742, 0.224231, 0.039560, 0.067314),
                        n = c(0.958811, 1.047887, 1.510932, 1.136141)),
           params = list(k ~ Species, n ~ Species),
           data = biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)])

boxplot(resid(m5, type = "normalized") ~ biomass[!is.na(biomass$Lyophilised),]$Species)
plot(resid(m5, type = "normalized") ~ biomass[!is.na(biomass$Lyophilised),]$Blotted) 
abline(0,0)
hist(resid(m5, type = "normalized"))

summary(m5)
#                                      Value  Std.Error   t-value p-value
# k.(Intercept)                    0.1767421 0.02384433  7.412333  0.0000
# k.SpeciesLaminaria hyperborea    0.0474887 0.03147943  1.508562  0.1326
# k.SpeciesSaccharina latissima   -0.1371816 0.02515582 -5.453276  0.0000
# k.SpeciesSaccorhiza polyschides -0.1094283 0.02795846 -3.913961  0.0001
# n.(Intercept)                    0.9588112 0.04198502 22.836984  0.0000
# n.SpeciesLaminaria hyperborea    0.0890753 0.05052455  1.763009  0.0790
# n.SpeciesSaccharina latissima    0.5521396 0.07583770  7.280542  0.0000
# n.SpeciesSaccorhiza polyschides  0.1773302 0.07527813  2.355667  0.0192

coef(m5)
# y = 0.17674214*x^0.95881125

biomass$Species <- factor(biomass$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                                      "Laminaria digitata", "Saccharina latissima"))
m5 <- gnls(Lyophilised ~ k*Blotted^n,
           start = list(k = c(0.224231, 0.067314, 0.176742, 0.039560),
                        n = c(1.047887, 1.136141, 0.958811, 1.510932)),
           params = list(k ~ Species, n ~ Species),
           data = biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)])

summary(m5)
#                                      Value  Std.Error  t-value p-value
# k.(Intercept)                    0.2242308 0.02055242 10.91019  0.0000
# k.SpeciesSaccorhiza polyschides -0.1569170 0.02520963 -6.22449  0.0000
# k.SpeciesLaminaria digitata     -0.0474887 0.03147943 -1.50856  0.1326
# k.SpeciesSaccharina latissima   -0.1846703 0.02206049 -8.37109  0.0000
# n.(Intercept)                    1.0478865 0.02810674 37.28240  0.0000
# n.SpeciesSaccorhiza polyschides  0.0882549 0.06851309  1.28815  0.1988
# n.SpeciesLaminaria digitata     -0.0890753 0.05052455 -1.76301  0.0790
# n.SpeciesSaccharina latissima    0.4630641 0.06912743  6.69870  0.0000

coef(m5)
# y = 0.22423082*x^1.04788650

biomass$Species <- factor(biomass$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                                      "Saccharina latissima", "Laminaria hyperborea"))
m5 <- gnls(Lyophilised ~ k*Blotted^n,
           start = list(k = c(0.067314, 0.176742, 0.039560, 0.224231),
                        n = c(1.136141, 0.958811, 1.510932, 1.047887)),
           params = list(k ~ Species, n ~ Species),
           data = biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)])

summary(m5)
#                                    Value  Std.Error   t-value p-value
# k.(Intercept)                  0.0673138 0.01459874  4.610932  0.0000
# k.SpeciesLaminaria digitata    0.1094283 0.02795846  3.913961  0.0001
# k.SpeciesSaccharina latissima -0.0277533 0.01665492 -1.666369  0.0968
# k.SpeciesLaminaria hyperborea  0.1569170 0.02520963  6.224486  0.0000
# n.(Intercept)                  1.1361414 0.06248244 18.183371  0.0000
# n.SpeciesLaminaria digitata   -0.1773302 0.07527813 -2.355667  0.0192
# n.SpeciesSaccharina latissima  0.3748089 0.08884069  4.218888  0.0000
# n.SpeciesLaminaria hyperborea -0.0882549 0.06851310 -1.288147  0.1988

coef(m5)
# y = 0.06731382*x^1.13614145

biomass$Species <- factor(biomass$Species, levels = c("Saccharina latissima", "Laminaria digitata", 
                                                      "Laminaria hyperborea", "Saccorhiza polyschides"))
m5 <- gnls(Lyophilised ~ k*Blotted^n,
           start = list(k = c(0.039560, 0.176742, 0.224231, 0.067314),
                        n = c(1.510932, 0.958811, 1.047887, 1.136141)),
           params = list(k ~ Species, n ~ Species),
           data = biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)])

summary(m5)
#                                      Value  Std.Error   t-value p-value
# k.(Intercept)                    0.0395612 0.00801653  4.934955  0.0000
# k.SpeciesLaminaria digitata      0.1371803 0.02515581  5.453225  0.0000
# k.SpeciesLaminaria hyperborea    0.1846692 0.02206051  8.371033  0.0000
# k.SpeciesSaccorhiza polyschides  0.0277528 0.01665499  1.666338  0.0968
# n.(Intercept)                    1.5109455 0.06315525 23.924305  0.0000
# n.SpeciesLaminaria digitata     -0.5521331 0.07583755 -7.280472  0.0000
# n.SpeciesLaminaria hyperborea   -0.4630585 0.06912724 -6.698639  0.0000
# n.SpeciesSaccorhiza polyschides -0.3748049 0.08884049 -4.218853  0.0000

coef(m5)
# y = 0.03956119*x^1.51094550

biomass$Species <- factor(biomass$Species, levels = c("Laminaria digitata", "Laminaria hyperborea", 
                                                      "Saccharina latissima", "Saccorhiza polyschides"))
m5 <- gnls(Lyophilised ~ k*Blotted^n,
           start = list(k = c(0.176742, 0.224231, 0.039560, 0.067314),
                        n = c(0.958811, 1.047887, 1.510932, 1.136141)),
           params = list(k ~ Species, n ~ Species),
           data = biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)])


# the maximal buoyant mass measured in the decomposition experiment was 47.744 g
# although m5 fits best in the specified range, the data only cover 0.026 to 44.27 g blotted mass
# hence, the robustness of m5's predictions need to be tested outside this range
# blotted mass is on average 90.73605% of buoyant mass, so the maximum mass to be tested is
47.744*0.9073605 # 43.32102 g

new <- data.frame(Blotted = rep(seq(0, 45, 0.1), 4), # create prediction dataframe
                  Species = c(rep("Laminaria digitata", 451),
                              rep("Laminaria hyperborea", 451),
                              rep("Saccharina latissima", 451),
                              rep("Saccorhiza polyschides", 451))) 

new$fit <- predict(m5, newdata = new)

ggplot() +
  geom_abline(intercept = 0, slope = 1, colour = "#c9d2d7") +
  geom_point(data = biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)], 
             aes(Blotted, Lyophilised, colour = Species)) +
  geom_line(data = new, aes(Blotted, fit, colour = Species)) +
  ylab("Lyophilised") +
  theme_minimal()

# m5 fits well in the predicted range, plus species-specific data tend to be
# in the species-specific prediction ranges
# this model can safely be used to convert blotted to lyophilised biomass

#### Data visualisation ####
#### Buoyant vs. Blotted ####
aggregate(Buoyant ~ Species, min, data = biomass)
aggregate(Buoyant ~ Species, max, data = biomass)

BBnew <- data.frame(Buoyant = c(seq(0.109, 13.150, 0.01), # generate new data
                                seq(0.051, 16.661, 0.01),
                                seq(0.101, 11.032, 0.01),
                                seq(0.066, 11.748, 0.01)),
                    Species = c(rep("Laminaria digitata", 1305),
                                rep("Laminaria hyperborea", 1662),
                                rep("Saccharina latissima", 1094),
                                rep("Saccorhiza polyschides", 1169)))

BBnew$fit <- predict(m3, newdata = BBnew)

# bootstrap confidence interval
bootfun <- function(newdata) {
  start <- coef(m3)
  boot <- biomass[!is.na(biomass$Buoyant),1:4][sample(nrow(biomass[!is.na(biomass$Buoyant),1:4]), 
                  size = nrow(biomass[!is.na(biomass$Buoyant),1:4]), replace = TRUE),]
  bootfit <- try(update(m3,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(BBnew))
BBnew$lwr <- apply(bmat, 1, quantile, 0.005, na.rm = TRUE)
BBnew$upr <- apply(bmat, 1, quantile, 0.995, na.rm = TRUE)

#### Blotted vs. Lyophilised ####
aggregate(Blotted ~ Species, min, data = biomass)
aggregate(Blotted ~ Species, max, data = biomass)

LBnew <- data.frame(Blotted = c(seq(0.026, 37.562, 0.01), # generate new data
                                seq(0.036, 34.008, 0.01),
                                seq(0.036, 34.784, 0.01),
                                seq(0.034, 44.270, 0.01)),
                    Species = c(rep("Laminaria digitata", 3754),
                                rep("Laminaria hyperborea", 3398),
                                rep("Saccharina latissima", 3475),
                                rep("Saccorhiza polyschides", 4424)))

LBnew$fit <- predict(m5, newdata = LBnew)

# bootstrap confidence interval
bootfun <- function(newdata) {
  start <- coef(m5)
  boot <- biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)][sample(nrow(biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)]), 
                  size = nrow(biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)]), replace = TRUE),]
  bootfit <- try(update(m5,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(LBnew))
LBnew$lwr <- apply(bmat, 1, quantile, 0.005, na.rm = TRUE)
LBnew$upr <- apply(bmat, 1, quantile, 0.995, na.rm = TRUE)



BBann <- aggregate(Buoyant ~ Species, length, data = biomass)
colnames(BBann)[2] <- "n"
BBann$n <- c(expression(italic("n ")*"= 70"),
             expression(italic("n ")*"= 70"),
             expression(italic("n ")*"= 70"),
             expression(italic("n ")*"= 70"))
BBann$n <- as.character(BBann$n)
BBann$equation <- c(expression(italic(y)*" = 0.75"*italic(x)^1.09),
                    expression(italic(y)*" = 0.71"*italic(x)^1.09),
                    expression(italic(y)*" = 0.83"*italic(x)^1.05),
                    expression(italic(y)*" = 0.85"*italic(x)^1.03))
BBann$equation <- as.character(BBann$equation)


LBann <- aggregate(Lyophilised ~ Species, length, data = biomass)
colnames(LBann)[2] <- "n"
LBann$n <- c(expression(italic("n ")*"= 70"),
             expression(italic("n ")*"= 70"),
             expression(italic("n ")*"= 70"),
             expression(italic("n ")*"= 70"))
LBann$n <- as.character(LBann$n)
LBann$equation <- c(expression(italic(y)*" = 0.18"*italic(x)^0.96),
                    expression(italic(y)*" = 0.22"*italic(x)^1.05),
                    expression(italic(y)*" = 0.04"*italic(x)^1.51),
                    expression(italic(y)*" = 0.07"*italic(x)^1.14))
LBann$equation <- as.character(LBann$equation)


mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .3, .2, .2),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15, hjust = 0),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.key.size = unit(.3, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.spacing.x = unit(.1, "cm"),
                 legend.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title.align = 0,
                 legend.title = element_text(size = 12, face = "bold"),
                 text = element_text(family = "Helvetica Neue"))

BBp <- ggplot(data = BBnew) +
          geom_abline(intercept = 0, slope = 1, colour = "#c9d2d7") +
          geom_line(aes(Buoyant, fit, colour = Species)) +
          geom_ribbon(aes(Buoyant, ymin = lwr, ymax = upr, fill = Species), alpha = 0.5) +
          geom_point(data = biomass[!is.na(biomass$Buoyant),1:4], aes(Buoyant, Blotted, colour = Species),
                     size = 2, shape = 16, alpha = 0.3) +
          scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                              guide = "none") +
          scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                            guide = "none") +
          facet_grid(~Species) +
          geom_text(data = BBann, aes(0.58, 15.5, label = n),
                    family = "Helvetica Neue", parse = T, size = 4.2, hjust = 0) +
          geom_text(data = BBann, aes(0.58, 14, label = equation),
                    family = "Helvetica Neue", parse = T, size = 4.2, hjust = 0) +
          labs(y = "Blotted mass (g)",
               x = "Buoyant mass (g)") +
          scale_x_continuous(breaks = seq(0, 18, by = 6)) +
          coord_cartesian(xlim = c(0, 18), ylim = c(0, 16), expand = FALSE) +
          mytheme +
          theme(strip.background = element_blank(),
                strip.text = element_text(size = 12, hjust = 0,
                                          face = "italic"),
                panel.spacing = unit(.5, "cm"))

BBp # dimensions: 3 x 9 in

aggregate(LB ~ Species, mean, data = biomass)

LBp <- ggplot(data = LBnew) +
          geom_abline(intercept = 0, slope = 1, colour = "#c9d2d7") +
          geom_line(aes(Blotted, fit, colour = Species)) +
          geom_ribbon(aes(Blotted, ymin = lwr, ymax = upr, fill = Species), alpha = 0.5) +
          geom_point(data = biomass[!is.na(biomass$Lyophilised),c(1:2, 4:5)], aes(Blotted, Lyophilised, colour = Species),
                     size = 2, shape = 16, alpha = 0.3) +
          scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                              guide = "none") +
          scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                            guide = "none") +
          facet_grid(~Species) +
          geom_text(data = LBann, aes(23, 15.5, label = n),
                    family = "Helvetica Neue", parse = T, size = 4.2, hjust = 0) +
          geom_text(data = LBann, aes(23, 14, label = equation),
                    family = "Helvetica Neue", parse = T, size = 4.2, hjust = 0) +
          labs(y = "Lyophilised mass (g)",
               x = "Blotted mass (g)") +
          scale_x_continuous(breaks = seq(0, 45, by = 15)) +
          coord_cartesian(xlim = c(0, 45), ylim = c(0, 16), expand = FALSE) +
          mytheme +
          theme(strip.background = element_blank(),
                strip.text = element_text(size = 12, hjust = 0,
                                          face = "italic", colour = NA),
                panel.spacing = unit(.5, "cm"))

LBp # dimensions: 3 x 9 in

require(cowplot)
BBLB <- plot_grid(BBp, LBp, nrow = 2, labels = "auto", 
                  label_size = 15, label_fontfamily = "Helvetica Neue")
BBLB # dimensions: 6 x 9 in

#### Load data ####
deco <- read.csv("~/Desktop/Projects/QUB/Data/Pigments/Decomposition.csv")

#### Build function ####
butoly <- function(buoyant, species){
  if(species == "Laminaria digitata"){
    0.17674214*(0.89487955*buoyant)^0.95881125
  } else if(species == "Laminaria hyperborea"){
    0.22423082*(0.87910760*buoyant)^1.04788650
  } else if(species == "Saccharina latissima"){
    0.03956119*(0.923208693*buoyant)^1.51094550
  } else if(species == "Saccorhiza polyschides"){
    0.06731382*(0.932246089*buoyant)^1.13614145
  } else{
    'Not a valid species'
  }
}
butoly <- Vectorize(butoly)

deco$Decomposition <- with(deco, (butoly(Initial, Species) - butoly(Final, Species))/Age)
s <- deco %>% group_by(Species) %>%
  summarise(mean = mean(Decomposition*1000),
            sd = sd(Decomposition*1000),
            n = length(Decomposition),
            se = sd/sqrt(n))
s

# unlike the decomposition experiment, the grazing experiment used tissue samples 
# of up to only 8.666 g buoyant mass, so a conversion function incorporating the 
# better fitting buoyant-blotted power equation is more appropriate
butolyg <- function(buoyant, species){
  if(species == "Laminaria digitata"){
    0.17674214*(0.749014421*buoyant^1.091948012)^0.95881125
  } else if(species == "Laminaria hyperborea"){
    0.22423082*(0.707360349*buoyant^1.089017974)^1.04788650
  } else if(species == "Saccharina latissima"){
    0.03956119*(0.83253578*buoyant^1.05051710)^1.51094550
  } else if(species == "Saccorhiza polyschides"){
    0.06731382*(0.85006949*buoyant^1.03367817)^1.13614145
  } else{
    'Not a valid species'
  }
}
butolyg <- Vectorize(butolyg)



