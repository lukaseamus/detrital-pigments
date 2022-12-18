###################################################################
##### Project: Genus-specific response of kelp photosynthetic #####
#####          pigments to decomposition                      #####
##### Script purpose: Data analysis and visualisation         #####
##### Author: Luka Seamus Wright                              #####
###################################################################

#### 1. Averaged data preparation ####
#### 1.1 Load averaged data ####
pig <- read.csv("~/Desktop/Projects/QUB/Data/Pigments/Pigments_averaged.csv")

#### 1.2 Rename variables ####
sp <- pig$Species
age <- pig$Age # Detrital age (d)
mesh <- factor(pig$Diameter) # Mesh diameter (mm)
chla <- pig$Chlorophyll.a # Chlorophyll a (µg g-1)
chlc <- pig$Chlorophyll.c # Chlorophyll c (µg g-1)
fuco <- pig$Fucoxanthin # Fucoxanthin (µg g-1)
ant <- pig$Antenna # Antenna size (antenna pigment to chlorophyll a ratio)
tot <- pig$Total # Total pigment concentration including minor pigments (µg g-1)

#### 1.3 Visual data exploration ####
require(ggplot2)
ggplot(pig, aes(age, chla, colour = sp)) +
  geom_point() +
  geom_smooth(span = 10, se = F) +
  theme_minimal()

ggplot(pig, aes(age, chlc, colour = sp)) +
  geom_point() +
  geom_smooth(span = 10, se = F) +
  theme_minimal()

ggplot(pig, aes(age, fuco, colour = sp)) +
  geom_point() +
  geom_smooth(span = 10, se = F) +
  theme_minimal()

ggplot(pig, aes(age, ant, colour = sp)) +
  geom_point() +
  geom_smooth(span = 10, se = F) +
  theme_minimal()

ggplot(pig, aes(age, tot, colour = sp)) +
  geom_point() +
  geom_smooth(span = 10, se = F) +
  theme_minimal()

#### 2. Averaged data analysis ####
#### 2.1 Mesh size effect ####
m1 <- lm(chla[-c(1:100)] ~ mesh[-c(1:100)])
boxplot(resid(m1) ~ mesh[-c(1:100)])
hist(resid(m1)) # right-skewed

m2 <- glm(chla[-c(1:100)] ~ mesh[-c(1:100)],
          family = Gamma(link = "log"))
boxplot(resid(m2) ~ mesh[-c(1:100)])
hist(resid(m2)) # normal

require(car)
Anova(m2)
# Response: chla[-c(1:100)]
#                 LR Chisq Df Pr(>Chisq)
# mesh[-c(1:100)]  0.13526  1      0.713

m3 <- lm(chlc[-c(1:100)] ~ mesh[-c(1:100)])
boxplot(resid(m3) ~ mesh[-c(1:100)])
hist(resid(m3)) # right-skewed

m4 <- glm(chlc[-c(1:100)] ~ mesh[-c(1:100)],
          family = Gamma(link = "log"))
boxplot(resid(m4) ~ mesh[-c(1:100)])
hist(resid(m4)) # normal

Anova(m4)
# Response: chlc[-c(1:100)]
#                 LR Chisq Df Pr(>Chisq)
# mesh[-c(1:100)] 0.001131  1     0.9732

m5 <- lm(fuco[-c(1:100)] ~ mesh[-c(1:100)])
boxplot(resid(m5) ~ mesh[-c(1:100)])
hist(resid(m5)) # right-skewed

m6 <- glm(fuco[-c(1:100)] ~ mesh[-c(1:100)],
          family = Gamma(link = "log"))
boxplot(resid(m6) ~ mesh[-c(1:100)])
hist(resid(m6)) # normal

Anova(m6)
# Response: fuco[-c(1:100)]
#                 LR Chisq Df Pr(>Chisq)
# mesh[-c(1:100)]  0.61443  1     0.4331

m7 <- lm(ant[-c(1:100)] ~ mesh[-c(1:100)])
boxplot(resid(m7) ~ mesh[-c(1:100)])
hist(resid(m7)) # normal

Anova(m7)
# Response: ant[-c(1:100)]
#                  Sum Sq Df F value Pr(>F)
# mesh[-c(1:100)] 0.01789  1  2.2451 0.1381
# Residuals       0.62139 78   

m8 <- lm(tot[-c(1:100)] ~ mesh[-c(1:100)])
boxplot(resid(m8) ~ mesh[-c(1:100)])
hist(resid(m8)) # right-skewed

m9 <- glm(tot[-c(1:100)] ~ mesh[-c(1:100)],
          family = Gamma(link = "log"))
boxplot(resid(m9) ~ mesh[-c(1:100)])
hist(resid(m9)) # normal

Anova(m9)
# Response: tot[-c(1:100)]
#                 LR Chisq Df Pr(>Chisq)
# mesh[-c(1:100)]  0.34509  1     0.5569 


#### 2.2 Chlorophyll a ####
pig2 <- pig[pig$Species != "Saccharina latissima",] # generate dataframe excluding S. latissima

m10 <- lm(Chlorophyll.a ~ Age * Species, 
          data = pig2)
boxplot(resid(m10) ~ sp[sp != "Saccharina latissima"]) 
plot(resid(m10) ~ age[sp != "Saccharina latissima"]) # heterogenous
hist(resid(m10)) # normal

require(nlme)
m11 <- gls(Chlorophyll.a ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)
boxplot(resid(m11, type = "normalized") ~ sp[sp != "Saccharina latissima"]) 
plot(resid(m11, type = "normalized") ~ age[sp != "Saccharina latissima"]) # homogenous
hist(resid(m11, type = "normalized")) # normal

Anova(m11, type = 3)
# Response: Chlorophyll.a
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 122.3609  < 2.2e-16 ***
# Age          1   3.0844   0.079044 .  
# Species      2  64.9833  7.746e-15 ***
# Age:Species  2  12.3013   0.002132 ** 

coef(m11)/1000 # convert from µg g-1 to mg g-1
# y = -0.004247106x + 1.012052579

pig2$Species <- factor(pig2$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                                "Laminaria digitata"))
m11 <- gls(Chlorophyll.a ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

Anova(m11, type = 3)
# Response: Chlorophyll.a
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 58.5780  1.954e-14 ***
# Age          1  0.3031   0.581959    
# Species      2 64.9833  7.746e-15 ***
# Age:Species  2 12.3013   0.002132 ** 

coef(m11)/1000
# y = -0.001032310x + 0.582271796

pig2$Species <- factor(pig2$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                                "Laminaria hyperborea"))
m11 <- gls(Chlorophyll.a ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

Anova(m11, type = 3)
# Response: Chlorophyll.a
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 242.493  < 2.2e-16 ***
# Age          1  16.316  5.360e-05 ***
# Species      2  64.983  7.747e-15 ***
# Age:Species  2  12.301   0.002132 ** 

coef(m11)/1000
# y = -0.01778542x + 1.61679195

pig2$Species <- factor(pig2$Species, levels = c("Laminaria digitata", "Laminaria hyperborea",
                                                "Saccorhiza polyschides"))
m11 <- gls(Chlorophyll.a ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

pig3 <- pig[pig$Species == "Saccharina latissima",]

m12 <- lm(Chlorophyll.a ~ poly(Age, 2, raw = TRUE), 
          data = pig3)
plot(resid(m12) ~ age[sp == "Saccharina latissima"]) # homogenous
hist(resid(m12)) # normal

Anova(m12, type = 2)
# Response: Chlorophyll.a
#                           Sum Sq Df F value    Pr(>F)    
# poly(Age, 2, raw = TRUE) 6234018  2  49.552 1.958e-11 ***
# Residuals                2453248 39  

coef(m12)/1000
# y = 0.108502135x - 0.002585103x^2 + 0.655552309


#### 2.3 Chlorophyll c ####
m13 <- lm(Chlorophyll.c ~ Age * Species, 
          data = pig2)
boxplot(resid(m13) ~ sp[sp != "Saccharina latissima"]) 
plot(resid(m13) ~ age[sp != "Saccharina latissima"]) # heterogenous
hist(resid(m13)) # normal

m14 <- gls(Chlorophyll.c ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)
boxplot(resid(m14, type = "normalized") ~ sp[sp != "Saccharina latissima"]) 
plot(resid(m14, type = "normalized") ~ age[sp != "Saccharina latissima"]) # homogenous
hist(resid(m14, type = "normalized")) # normal

Anova(m14, type = 3)
# Response: Chlorophyll.c
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 22.6003  1.995e-06 ***
# Age          1  0.0106    0.91781    
# Species      2 60.0893  8.949e-14 ***
# Age:Species  2  6.9548    0.03089 *
  
coef(m14)/1000
# y = -1.741735e-05x + 3.241781e-02

pig2$Species <- factor(pig2$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                                "Laminaria digitata"))
m14 <- gls(Chlorophyll.c ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

Anova(m14, type = 3)
# Response: Chlorophyll.c
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 12.9170  0.0003256 ***
# Age          1  1.5607  0.2115633    
# Species      2 60.0893  8.949e-14 ***
# Age:Species  2  6.9548  0.0308879 * 

coef(m14)/1000
# y = 0.0001882696x + 0.0222851260

pig2$Species <- factor(pig2$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                                "Laminaria hyperborea"))
m14 <- gls(Chlorophyll.c ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

Anova(m14, type = 3)
# Response: Chlorophyll.c
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 129.9728  < 2.2e-16 ***
# Age          1   5.5024    0.01899 *  
# Species      2  60.0891   8.95e-14 ***
# Age:Species  2   6.9548    0.03089 *  

coef(m14)/1000
# y = -0.001212425x + 0.107082467

pig2$Species <- factor(pig2$Species, levels = c("Laminaria digitata", "Laminaria hyperborea",
                                                "Saccorhiza polyschides"))
m14 <- gls(Chlorophyll.c ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

m15 <- lm(Chlorophyll.c ~ poly(Age, 2, raw = TRUE), 
          data = pig3)
plot(resid(m15) ~ age[sp == "Saccharina latissima"]) # homogenous
hist(resid(m15)) # normal

Anova(m15, type = 2)
# Response: Chlorophyll.c
#                          Sum Sq Df F value  Pr(>F)   
# poly(Age, 2, raw = TRUE)  18594  2  8.1155 0.00113 **
# Residuals                 44679 39

coef(m15)/1000
# y = 0.0059025067x - 0.0001403812x^2 + 0.0543497540


#### 2.4 Fucoxanthin ####
m16 <- lm(Fucoxanthin ~ Age * Species, 
          data = pig2)
boxplot(resid(m16) ~ sp[sp != "Saccharina latissima"]) 
plot(resid(m16) ~ age[sp != "Saccharina latissima"]) # heterogenous
hist(resid(m16)) # normal

m17 <- gls(Fucoxanthin ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)
boxplot(resid(m17, type = "normalized") ~ sp[sp != "Saccharina latissima"]) 
plot(resid(m17, type = "normalized") ~ age[sp != "Saccharina latissima"]) # homogenous
hist(resid(m17, type = "normalized")) # normal

Anova(m17, type = 3)
# Response: Fucoxanthin
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 121.4658  < 2.2e-16 ***
# Age          1   0.0153    0.90171    
# Species      2  45.6963  1.194e-10 ***
# Age:Species  2   7.1336    0.02825 * 

coef(m17)/1000
# y = -0.0002450320x + 0.8184203415

pig2$Species <- factor(pig2$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                                "Laminaria digitata"))
m17 <- gls(Fucoxanthin ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

Anova(m17, type = 3)
# Response: Fucoxanthin
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 49.6111  1.875e-12 ***
# Age          1  0.2126    0.64475    
# Species      2 45.6963  1.194e-10 ***
# Age:Species  2  7.1336    0.02825 * 

coef(m17)/1000
# y = -0.0008367267x + 0.4974846836

pig2$Species <- factor(pig2$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                                "Laminaria hyperborea"))
m17 <- gls(Fucoxanthin ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

Anova(m17, type = 3)
# Response: Fucoxanthin
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 219.6666  < 2.2e-16 ***
# Age          1   8.7940   0.003022 ** 
# Species      2  45.6962  1.195e-10 ***
# Age:Species  2   7.1335   0.028247 * 

coef(m17)/1000
# y = -0.01141812x + 1.23614320

pig2$Species <- factor(pig2$Species, levels = c("Laminaria digitata", "Laminaria hyperborea",
                                                "Saccorhiza polyschides"))
m17 <- gls(Fucoxanthin ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

m18 <- lm(Fucoxanthin ~ poly(Age, 2, raw = TRUE), 
          data = pig3)
plot(resid(m18) ~ age[sp == "Saccharina latissima"]) # homogenous
hist(resid(m18)) # normal

Anova(m18, type = 2)
# Response: Fucoxanthin
#                           Sum Sq Df F value    Pr(>F)    
# poly(Age, 2, raw = TRUE) 4779016  2  61.393 8.941e-13 ***
# Residuals                1517931 39  

coef(m18)/1000

#### 2.5 Antenna size ####
m19 <- lm(Antenna ~ Age * Species, 
          data = pig2)
boxplot(resid(m19) ~ sp[sp != "Saccharina latissima"]) # homogenous
plot(resid(m19) ~ age[sp != "Saccharina latissima"]) # homogenous
hist(resid(m19)) # normal

Anova(m19, type = 3)
# Response: Antenna
#              Sum Sq  Df   F value    Pr(>F)    
# (Intercept) 17.6111   1 2865.3538 < 2.2e-16 ***
# Age          0.3721   1   60.5411 1.839e-12 ***
# Species      0.0559   2    4.5488   0.01229 *  
# Age:Species  0.1653   2   13.4473 4.838e-06 ***
# Residuals    0.8113 132 

coef(m19)
# y = 0.004655092x + 0.819808655

pig2$Species <- factor(pig2$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                                "Laminaria digitata"))
m19 <- lm(Antenna ~ Age * Species, 
          data = pig2)

Anova(m19, type = 3)
# Response: Antenna
#              Sum Sq  Df   F value    Pr(>F)    
# (Intercept) 20.3696   1 3314.1642 < 2.2e-16 ***
# Age          0.0015   1    0.2455   0.62109    
# Species      0.0559   2    4.5488   0.01229 *  
# Age:Species  0.1653   2   13.4473 4.838e-06 ***
# Residuals    0.8113 132

coef(m19)
# y = 0.0002943062x + 0.8819743964

pig2$Species <- factor(pig2$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                                "Laminaria hyperborea"))
m19 <- lm(Antenna ~ Age * Species, 
          data = pig2)

Anova(m19, type = 3)
# Response: Antenna
#              Sum Sq  Df   F value    Pr(>F)    
# (Intercept) 17.8000   1 2896.0844 < 2.2e-16 ***
# Age          0.0421   1    6.8576  0.009859 ** 
# Species      0.0559   2    4.5488  0.012291 *  
# Age:Species  0.1653   2   13.4473 4.838e-06 ***
# Residuals    0.8113 132

coef(m19)
# y = 0.002118332x + 0.833333020

pig2$Species <- factor(pig2$Species, levels = c("Laminaria digitata", "Laminaria hyperborea",
                                                "Saccorhiza polyschides"))
m19 <- lm(Antenna ~ Age * Species, 
          data = pig2)


m20 <- lm(Antenna ~ poly(Age, 2, raw = TRUE), 
          data = pig3)
plot(resid(m20) ~ age[sp == "Saccharina latissima"]) # homogeneity cannot be improved
hist(resid(m20)) # normal

Anova(m20, type = 2)
# Response: Antenna
#                           Sum Sq Df F value Pr(>F)
# poly(Age, 2, raw = TRUE) 0.02243  2  0.8186 0.4485
# Residuals                0.53435 39 

coef(m20)
# y = -0.0060350406x + 0.0001871924x^2 + 0.9685480480

#### 2.6 Total pigment ####
m21 <- lm(Total ~ Age * Species, 
          data = pig2)
boxplot(resid(m21) ~ sp[sp != "Saccharina latissima"]) 
plot(resid(m21) ~ age[sp != "Saccharina latissima"]) # heterogenous
hist(resid(m21)) # normal

m22 <- gls(Total ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)
boxplot(resid(m22, type = "normalized") ~ sp[sp != "Saccharina latissima"]) 
plot(resid(m22, type = "normalized") ~ age[sp != "Saccharina latissima"]) # homogenous
hist(resid(m22, type = "normalized")) # normal

Anova(m22, type = 3)
# Response: Total
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 117.2621  < 2.2e-16 ***
# Age          1   0.9178    0.33804    
# Species      2  62.6861  2.443e-14 ***
# Age:Species  2   8.8743    0.01183 *  

coef(m22)/1000
# y = -0.004357507x + 1.873352248

pig2$Species <- factor(pig2$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                                "Laminaria digitata"))
m22 <- gls(Total ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

Anova(m22, type = 3)
# Response: Total
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 51.0558  8.978e-13 ***
# Age          1  0.0763    0.78237    
# Species      2 62.6861  2.443e-14 ***
# Age:Species  2  8.8743    0.01183 *  

coef(m22)/1000
# y = -0.001069048x + 1.106535990

pig2$Species <- factor(pig2$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                                "Laminaria hyperborea"))
m22 <- gls(Total ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

Anova(m22, type = 3)
# Response: Total
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 243.5867  < 2.2e-16 ***
# Age          1  11.1099  0.0008587 ***
# Species      2  62.6867  2.442e-14 ***
# Age:Species  2   8.8743  0.0118298 * 

coef(m22)/1000
# y = -0.03108069x + 3.10096599

pig2$Species <- factor(pig2$Species, levels = c("Laminaria digitata", "Laminaria hyperborea",
                                                "Saccorhiza polyschides"))
m22 <- gls(Total ~ Age * Species, 
           weights = varExp(form = ~Age|Species),
           data = pig2)

m23 <- lm(Total ~ poly(Age, 2, raw = TRUE), 
          data = pig3)
plot(resid(m23) ~ age[sp == "Saccharina latissima"]) # homogenous
hist(resid(m23)) # normal

Anova(m23, type = 2)
# Response: Total
#                            Sum Sq Df F value    Pr(>F)    
# poly(Age, 2, raw = TRUE) 21483587  2  40.717 2.826e-10 ***
# Residuals                10288864 39 

coef(m23)/1000
# y = 0.197591285x - 0.004667492x^2 + 1.421183826


new <- data.frame(Age = c(rep(seq(0, 46, by = 0.1), 2),
                          rep(seq(0, 39, by = 0.1), 2)),
                  Species = c(rep("Laminaria digitata", 461),
                              rep("Laminaria hyperborea", 461),
                              rep("Saccharina latissima", 391),
                              rep("Saccorhiza polyschides", 391)))

fit1 <- predict(m11, newdata = new[new$Species != "Saccharina latissima",])
fit2 <- predict(m12, newdata = new[new$Species == "Saccharina latissima",])

new$Chla.fit <- c(fit1[1:922], fit2, fit1[923:1313])

modmat <-  model.matrix(formula(m11)[-2], new[new$Species != "Saccharina latissima",])
int <- diag(modmat %*% vcov(m11) %*% t(modmat))
lwr1 <- fit1 - qnorm(0.975)*sqrt(int)
upr1 <- fit1 + qnorm(0.975)*sqrt(int)

lwr2 <- predict(m12, interval = "confidence", 
                newdata = new[new$Species == "Saccharina latissima",])[,2]
upr2 <- predict(m12, interval = "confidence", 
                newdata = new[new$Species == "Saccharina latissima",])[,3]

new$Chla.lwr <- c(lwr1[1:922], lwr2, lwr1[923:1313])
new$Chla.upr <- c(upr1[1:922], upr2, upr1[923:1313])


fit1 <- predict(m14, newdata = new[new$Species != "Saccharina latissima",])
fit2 <- predict(m15, newdata = new[new$Species == "Saccharina latissima",])

new$Chlc.fit <- c(fit1[1:922], fit2, fit1[923:1313])

modmat <-  model.matrix(formula(m14)[-2], new[new$Species != "Saccharina latissima",])
int <- diag(modmat %*% vcov(m14) %*% t(modmat))
lwr1 <- fit1 - qnorm(0.975)*sqrt(int)
upr1 <- fit1 + qnorm(0.975)*sqrt(int)

lwr2 <- predict(m15, interval = "confidence", 
                newdata = new[new$Species == "Saccharina latissima",])[,2]
upr2 <- predict(m15, interval = "confidence", 
                newdata = new[new$Species == "Saccharina latissima",])[,3]

new$Chlc.lwr <- c(lwr1[1:922], lwr2, lwr1[923:1313])
new$Chlc.upr <- c(upr1[1:922], upr2, upr1[923:1313])


fit1 <- predict(m17, newdata = new[new$Species != "Saccharina latissima",])
fit2 <- predict(m18, newdata = new[new$Species == "Saccharina latissima",])

new$Fuco.fit <- c(fit1[1:922], fit2, fit1[923:1313])

modmat <-  model.matrix(formula(m17)[-2], new[new$Species != "Saccharina latissima",])
int <- diag(modmat %*% vcov(m17) %*% t(modmat))
lwr1 <- fit1 - qnorm(0.975)*sqrt(int)
upr1 <- fit1 + qnorm(0.975)*sqrt(int)

lwr2 <- predict(m18, interval = "confidence", 
                newdata = new[new$Species == "Saccharina latissima",])[,2]
upr2 <- predict(m18, interval = "confidence", 
                newdata = new[new$Species == "Saccharina latissima",])[,3]

new$Fuco.lwr <- c(lwr1[1:922], lwr2, lwr1[923:1313])
new$Fuco.upr <- c(upr1[1:922], upr2, upr1[923:1313])


fit1 <- predict(m19, newdata = new[new$Species != "Saccharina latissima",])
fit2 <- predict(m20, newdata = new[new$Species == "Saccharina latissima",])

new$Ant.fit <- c(fit1[1:922], fit2, fit1[923:1313])

lwr1 <- predict(m19, interval = "confidence", 
                newdata = new[new$Species != "Saccharina latissima",])[,2]
upr1 <- predict(m19, interval = "confidence", 
                newdata = new[new$Species != "Saccharina latissima",])[,3]


lwr2 <- predict(m20, interval = "confidence", 
                newdata = new[new$Species == "Saccharina latissima",])[,2]
upr2 <- predict(m20, interval = "confidence", 
                newdata = new[new$Species == "Saccharina latissima",])[,3]

new$Ant.lwr <- c(lwr1[1:922], lwr2, lwr1[923:1313])
new$Ant.upr <- c(upr1[1:922], upr2, upr1[923:1313])


fit1 <- predict(m22, newdata = new[new$Species != "Saccharina latissima",])
fit2 <- predict(m23, newdata = new[new$Species == "Saccharina latissima",])

new$Tot.fit <- c(fit1[1:922], fit2, fit1[923:1313])

modmat <-  model.matrix(formula(m22)[-2], new[new$Species != "Saccharina latissima",])
int <- diag(modmat %*% vcov(m22) %*% t(modmat))
lwr1 <- fit1 - qnorm(0.975)*sqrt(int)
upr1 <- fit1 + qnorm(0.975)*sqrt(int)

lwr2 <- predict(m23, interval = "confidence", 
                newdata = new[new$Species == "Saccharina latissima",])[,2]
upr2 <- predict(m23, interval = "confidence", 
                newdata = new[new$Species == "Saccharina latissima",])[,3]

new$Tot.lwr <- c(lwr1[1:922], lwr2, lwr1[923:1313])
new$Tot.upr <- c(upr1[1:922], upr2, upr1[923:1313])

require(dplyr)
new <- new %>% 
  mutate(Group = ifelse(Species == "Laminaria digitata" | Species == "Laminaria hyperborea", 
                        "Laminaria spp.", 
                        "Other genera"))
new$Group <- factor(new$Group)

pig <- pig %>% 
  mutate(Group = ifelse(Species == "Laminaria digitata" | Species == "Laminaria hyperborea", 
                        "Laminaria spp.", 
                        "Other genera"))
pig$Group <- factor(pig$Group)

require(psych)
Chla.stat <- describeBy(chla, list(age, sp), mat = T)
Chla.stat$group1 <- as.integer(Chla.stat$group1)
Chla.stat <- Chla.stat %>% 
  mutate(Group = ifelse(group2 == "Laminaria digitata" | group2 == "Laminaria hyperborea", 
                        "Laminaria spp.", 
                        "Other genera"))
Chla.stat$Group <- factor(Chla.stat$Group)

Chlc.stat <- describeBy(chlc, list(age, sp), mat = T)
Chlc.stat$group1 <- as.integer(Chlc.stat$group1)
Chlc.stat <- Chlc.stat %>% 
  mutate(Group = ifelse(group2 == "Laminaria digitata" | group2 == "Laminaria hyperborea", 
                        "Laminaria spp.", 
                        "Other genera"))
Chlc.stat$Group <- factor(Chlc.stat$Group)

Fuco.stat <- describeBy(fuco, list(age, sp), mat = T)
Fuco.stat$group1 <- as.integer(Fuco.stat$group1)
Fuco.stat <- Fuco.stat %>% 
  mutate(Group = ifelse(group2 == "Laminaria digitata" | group2 == "Laminaria hyperborea", 
                        "Laminaria spp.", 
                        "Other genera"))
Fuco.stat$Group <- factor(Fuco.stat$Group)

Ant.stat <- describeBy(ant, list(age, sp), mat = T)
Ant.stat$group1 <- as.integer(Ant.stat$group1)
Ant.stat <- Ant.stat %>% 
  mutate(Group = ifelse(group2 == "Laminaria digitata" | group2 == "Laminaria hyperborea", 
                        "Laminaria spp.", 
                        "Other genera"))
Ant.stat$Group <- factor(Ant.stat$Group)

Tot.stat <- describeBy(tot, list(age, sp), mat = T)
Tot.stat$group1 <- as.integer(Tot.stat$group1)
Tot.stat <- Tot.stat %>% 
  mutate(Group = ifelse(group2 == "Laminaria digitata" | group2 == "Laminaria hyperborea", 
                        "Laminaria spp.", 
                        "Other genera"))
Tot.stat$Group <- factor(Tot.stat$Group)

require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .3, .2, .5),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15, hjust = 0),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.width = unit(.45, "cm"),
                 legend.spacing.x = unit(.1, "cm"),
                 legend.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title.align = 0,
                 legend.title = element_blank(),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 15, hjust = 0),
                 panel.spacing = unit(.5, "cm"),
                 text = element_text(family = "Helvetica Neue"))

chlap <- ggplot() +
  geom_point(data = pig, mapping = aes(Age, Chlorophyll.a/1000, colour = Species),
             position = position_dodge(width = 3.5), size = 3, shape = 16, alpha = 0.2) +
  geom_line(data = new, mapping = aes(Age, Chla.fit/1000, colour = Species,
                                      lty = Species)) +
  geom_ribbon(data = new, mapping = aes(Age, ymin = Chla.lwr/1000, ymax = Chla.upr/1000,
                                        fill = Species), alpha = 0.4) +
  geom_pointrange(data = Chla.stat,
                  mapping = aes(group1, mean/1000, ymin = mean/1000 - (qnorm(0.975)*se)/1000,
                  ymax = mean/1000 + (qnorm(0.975)*se)/1000, colour = group2),
                  position = position_dodge(width = 3.5)) +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                    guide = guide_legend(ncol = 2)) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                      guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c(5, 5, 1, 1),
                        guide = "none") +
  facet_grid(~Group) +
  coord_cartesian(xlim = c(-2.5, 50), ylim = c(0, 4), expand = F) +
  labs(x = "Detrital age (d)",
       y = expression("Chlorophyll "*italic(a)*" (mg g"^-1*")")) +
  mytheme +
  theme(legend.position = c(0.41, 0.925))

chlap

chlcp <- ggplot() +
  geom_point(data = pig, mapping = aes(Age, Chlorophyll.c/1000, colour = Species),
             position = position_dodge(width = 3.5), size = 3, shape = 16, alpha = 0.2) +
  geom_line(data = new, mapping = aes(Age, Chlc.fit/1000, colour = Species,
                                      lty = Species)) +
  geom_ribbon(data = new, mapping = aes(Age, ymin = Chlc.lwr/1000, ymax = Chlc.upr/1000,
                                        fill = Species), alpha = 0.4) +
  geom_pointrange(data = Chlc.stat,
                  mapping = aes(group1, mean/1000, ymin = mean/1000 - (qnorm(0.975)*se)/1000,
                  ymax = mean/1000 + (qnorm(0.975)*se)/1000, colour = group2),
                  position = position_dodge(width = 3.5)) +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                    guide = guide_legend(ncol = 2)) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                      guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c(5, 5, 1, 1),
                        guide = "none") +
  facet_grid(~Group) +
  coord_cartesian(xlim = c(-2.5, 50), ylim = c(0, 0.3), expand = F) +
  scale_y_continuous(breaks = seq(0, 0.3, by = 0.1)) +
  labs(x = "Detrital age (d)",
       y = expression("Chlorophyll "*italic(c)*" (mg g"^-1*")")) +
  mytheme +
  theme(legend.position = c(0.41, 0.925))

chlcp

fucop <- ggplot() +
  geom_point(data = pig, mapping = aes(Age, Fucoxanthin/1000, colour = Species),
             position = position_dodge(width = 3.5), size = 3, shape = 16, alpha = 0.2) +
  geom_line(data = new, mapping = aes(Age, Fuco.fit/1000, colour = Species,
                                      lty = Species)) +
  geom_ribbon(data = new, mapping = aes(Age, ymin = Fuco.lwr/1000, ymax = Fuco.upr/1000,
                                        fill = Species), alpha = 0.4) +
  geom_pointrange(data = Fuco.stat,
                  mapping = aes(group1, mean/1000, ymin = mean/1000 - (qnorm(0.975)*se)/1000,
                  ymax = mean/1000 + (qnorm(0.975)*se)/1000, colour = group2),
                  position = position_dodge(width = 3.5)) +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                    guide = guide_legend(ncol = 2)) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                      guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c(5, 5, 1, 1),
                        guide = "none") +
  facet_grid(~Group) +
  coord_cartesian(xlim = c(-2.5, 50), ylim = c(0, 3), expand = F) +
  labs(x = "Detrital age (d)",
       y = expression("Fucoxanthin (mg g"^-1*")")) +
  mytheme +
  theme(legend.position = c(0.41, 0.925))

fucop

antp <- ggplot() +
  geom_hline(intercept = 1) +
  geom_point(data = pig, mapping = aes(Age, Antenna, colour = Species),
             position = position_dodge(width = 3.5), size = 3, shape = 16, alpha = 0.2) +
  geom_line(data = new, mapping = aes(Age, Ant.fit, colour = Species,
                                      lty = Species)) +
  geom_ribbon(data = new, mapping = aes(Age, ymin = Ant.lwr, ymax = Ant.upr,
                                        fill = Species), alpha = 0.4) +
  geom_pointrange(data = Ant.stat,
                  mapping = aes(group1, mean, ymin = mean - (qnorm(0.975)*se),
                  ymax = mean + (qnorm(0.975)*se), colour = group2),
                  position = position_dodge(width = 3.5)) +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                    guide = guide_legend(ncol = 2)) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                      guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c(1, 5, 5, 1),
                        guide = "none") +
  facet_grid(~Group) +
  coord_cartesian(xlim = c(-2.5, 50), ylim = c(0.5, 1.25), expand = F) +
  scale_y_continuous(breaks = seq(0.5, 1.25, by = 0.25)) +
  labs(x = "Detrital age (d)",
       y = "Relative antenna size") +
  mytheme +
  theme(legend.position = c(0.41, 0.925))

antp

totp <- ggplot() +
  geom_point(data = pig, mapping = aes(Age, Total/1000, colour = Species),
             position = position_dodge(width = 3.5), size = 3, shape = 16, alpha = 0.2) +
  geom_line(data = new, mapping = aes(Age, Tot.fit/1000, colour = Species,
                                      lty = Species)) +
  geom_ribbon(data = new, mapping = aes(Age, ymin = Tot.lwr/1000, ymax = Tot.upr/1000,
                                        fill = Species), alpha = 0.4) +
  geom_pointrange(data = Tot.stat,
                  mapping = aes(group1, mean/1000, ymin = mean/1000 - (qnorm(0.975)*se)/1000,
                  ymax = mean/1000 + (qnorm(0.975)*se)/1000, colour = group2),
                  position = position_dodge(width = 3.5)) +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                    guide = guide_legend(ncol = 2)) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                      guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c(5, 5, 1, 1),
                        guide = "none") +
  facet_grid(~Group) +
  coord_cartesian(xlim = c(-2.5, 50), ylim = c(0, 6), expand = F) +
  scale_y_continuous(breaks = seq(0, 6, by = 2)) +
  labs(x = "Detrital age (d)",
       y = expression("Total pigment (mg g"^-1*")")) +
  mytheme +
  theme(legend.position = c(0.41, 0.925))

totp

require(cowplot)
chlap <- chlap + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank())
chlcp <- chlcp + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       strip.text = element_blank(),
                       legend.position = "none")
fucop <- fucop + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       strip.text = element_blank(),
                       legend.position = "none")
antp <- antp + theme(strip.text = element_blank(),
                     legend.position = "none")

comb <- plot_grid(chlap, chlcp, fucop, antp, labels = "auto", label_size = 15, label_fontfamily = "Helvetica Neue",
                  ncol = 1, rel_heights = c(0.98, 0.865, 0.865, 1), align = "v", hjust = 0)
comb
# dimensions: 10 x 6 in



#### 3. Data preparation ####
p <- read.csv("~/Desktop/Projects/QUB/Data/Pigments/Pigments.csv")
p2 <- p[p$Species != "Saccharina latissima",] # generate dataframe excluding S. latissima
p3 <- p[p$Species == "Saccharina latissima",] # generate dataframe excluding other species

#### 4. Data analysis ####
#### 4.1 Mesh size effect ####
m24 <- lme(Chlorophyll.a ~ Diameter,
           random = ~1|ID,
           data = p[-c(1:300),])
boxplot(resid(m24) ~ Diameter, data = p[-c(1:300),])
hist(resid(m24)) 

Anova(m24)
# Response: Chlorophyll.a
#           Chisq Df Pr(>Chisq)
# Diameter 0.121  1      0.728

m25 <- lme(Chlorophyll.c ~ Diameter,
           random = ~1|ID,
           data = p[-c(1:300),])
boxplot(resid(m25) ~ Diameter, data = p[-c(1:300),])
hist(resid(m25)) 

Anova(m25)
# Response: Chlorophyll.c
#           Chisq Df Pr(>Chisq)
# Diameter 0.0043  1     0.9474

m26 <- lme(Fucoxanthin ~ Diameter,
           random = ~1|ID,
           data = p[-c(1:300),])
boxplot(resid(m26) ~ Diameter, data = p[-c(1:300),])
hist(resid(m26)) 

Anova(m26)
# Response: Fucoxanthin
#           Chisq Df Pr(>Chisq)
# Diameter 0.5883  1     0.4431

m27 <- lme(Antenna ~ Diameter,
           random = ~1|ID,
           data = p[-c(1:300),])
boxplot(resid(m27) ~ Diameter, data = p[-c(1:300),])
hist(resid(m27)) 

Anova(m27)
# Response: Antenna
#           Chisq Df Pr(>Chisq)
# Diameter 2.4555  1     0.1171

m28 <- lme(Total ~ Diameter,
           random = ~1|ID,
           data = p[-c(1:300),])
boxplot(resid(m28) ~ Diameter, data = p[-c(1:300),])
hist(resid(m28)) 

Anova(m28)
# Response: Total
#           Chisq Df Pr(>Chisq)
# Diameter 0.3282  1     0.5667



#### 4.2 Chlorophyll a ####
m29 <- lme(Chlorophyll.a ~ Age * Species,
           random = ~1|ID,
           data = p2)
boxplot(resid(m29) ~ Species, data = p2) 
plot(resid(m29) ~ Age, data = p2) # heterogenous
hist(resid(m29)) # normal

m30 <- lme(Chlorophyll.a ~ Age * Species,
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)
boxplot(resid(m30, type = "normalized") ~ Species, data = p2) 
plot(resid(m30, type = "normalized") ~ Age, data = p2) # homogenous
hist(resid(m11, type = "normalized")) # normal

Anova(m30, type = 3)
# Response: Chlorophyll.a
#                Chisq Df Pr(>Chisq)    
# (Intercept) 155.8211  1  < 2.2e-16 ***
# Age           2.6078  1    0.10634    
# Species      67.5810  2  2.113e-15 ***
# Age:Species   8.7402  2    0.01265 * 

summary(m30) # convert from µg g-1 to mg g-1
# y = -0.0052434x + 1.037932

p2$Species <- factor(p2$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                            "Laminaria digitata"))
m30 <- lme(Chlorophyll.a ~ Age * Species,
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)

Anova(m30, type = 3)
# Response: Chlorophyll.a
#               Chisq Df Pr(>Chisq)    
# (Intercept) 64.1235  1  1.169e-15 ***
# Age          1.1846  1    0.27643    
# Species     67.5810  2  2.113e-15 ***
# Age:Species  8.7402  2    0.01265 * 

summary(m30)
# y = -0.0035077x + 0.6660352

p2$Species <- factor(p2$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                            "Laminaria hyperborea"))
m30 <- lme(Chlorophyll.a ~ Age * Species,
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)

Anova(m30, type = 3)
# Response: Chlorophyll.a
#                Chisq Df Pr(>Chisq)    
# (Intercept) 376.1624  1  < 2.2e-16 ***
# Age          18.5008  1  1.698e-05 ***
# Species      67.5810  2  2.113e-15 ***
# Age:Species   8.7402  2    0.01265 *

summary(m30)
# y = -0.0188828x + 1.630569

p2$Species <- factor(p2$Species, levels = c("Laminaria digitata", "Laminaria hyperborea",
                                            "Saccorhiza polyschides"))
m30 <- lme(Chlorophyll.a ~ Age * Species,
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)



m31 <- lme(Chlorophyll.a ~ poly(Age, 2, raw = TRUE),
           random = ~1|ID,
           data = p3)
plot(resid(m31) ~ Age, data = p3) # homogenous
hist(resid(m31)) # normal

Anova(m31)
# Response: Chlorophyll.a
#                           Chisq Df Pr(>Chisq)    
# poly(Age, 2, raw = TRUE) 98.358  2  < 2.2e-16 *** 

summary(m31)
# y = 0.1082898x - 0.0025802x^2 + 0.6556015



#### 4.3 Chlorophyll c ####
m32 <- lme(Chlorophyll.c ~ Age * Species, 
           random = ~1|ID,
           data = p2)
boxplot(resid(m32) ~ Species, data = p2) 
plot(resid(m32) ~ Age, data = p2) # heterogenous
hist(resid(m32)) # normal

m33 <- lme(Chlorophyll.c ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)
boxplot(resid(m33, type = "normalized") ~ Species, data = p2) 
plot(resid(m33, type = "normalized") ~ Age, data = p2) # homogenous
hist(resid(m33, type = "normalized")) # normal

Anova(m33, type = 3)
# Response: Chlorophyll.c
#               Chisq Df Pr(>Chisq)    
# (Intercept) 20.2559  1  6.774e-06 ***
# Age          0.0771  1    0.78120    
# Species     68.0483  2  1.673e-15 ***
# Age:Species  6.8782  2    0.03209 * 

summary(m33)
# y = -8.311e-05x + 0.03456425

p2$Species <- factor(p2$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                            "Laminaria digitata"))
m33 <- lme(Chlorophyll.c ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)

Anova(m33, type = 3)
# Response: Chlorophyll.c
#               Chisq Df Pr(>Chisq)    
# (Intercept)  9.6274  1   0.001917 ** 
# Age          0.2231  1   0.636657    
# Species     68.0483  2  1.673e-15 ***
# Age:Species  6.8782  2   0.032093 * 

summary(m33)
# y = 0.00014028x + 0.02383584

p2$Species <- factor(p2$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                            "Laminaria hyperborea"))
m33 <- lme(Chlorophyll.c ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)

Anova(m33, type = 3)
# Response: Chlorophyll.c
#                Chisq Df Pr(>Chisq)    
# (Intercept) 189.1370  1  < 2.2e-16 ***
# Age           7.9500  1   0.004809 ** 
# Species      68.0483  2  1.673e-15 ***
# Age:Species   6.8782  2   0.032093 *  

summary(m33)
# y = -0.00114899x + 0.1068248

p2$Species <- factor(p2$Species, levels = c("Laminaria digitata", "Laminaria hyperborea",
                                            "Saccorhiza polyschides"))
m33 <- lme(Chlorophyll.c ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)



m34 <- lme(Chlorophyll.c ~ poly(Age, 2, raw = TRUE), 
           random = ~1|ID,
           data = p3)
plot(resid(m34) ~ Age, data = p3) # homogenous
hist(resid(m34)) # normal

Anova(m34)
# Response: Chlorophyll.c
#                           Chisq Df Pr(>Chisq)    
# poly(Age, 2, raw = TRUE) 16.122  2  0.0003156 ***

summary(m34)
# y = 0.00589108x - 0.00014012x^2 + 0.05435242


#### 2.4 Fucoxanthin ####
m35 <- lme(Fucoxanthin ~ Age * Species,
           random = ~1|ID,
           data = p2)
boxplot(resid(m35) ~ Species, data = p2) 
plot(resid(m35) ~ Age, data = p2) # heterogenous
hist(resid(m35)) # normal

m36 <- lme(Fucoxanthin ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)
boxplot(resid(m36, type = "normalized") ~ Species, data = p2) 
plot(resid(m36, type = "normalized") ~ Age, data = p2) # homogenous
hist(resid(m36, type = "normalized")) # normal

Anova(m36, type = 3)
# Response: Fucoxanthin
#                Chisq Df Pr(>Chisq)    
# (Intercept) 154.0034  1  < 2.2e-16 ***
# Age           0.1411  1     0.7072    
# Species      49.4911  2  1.791e-11 ***
# Age:Species   6.2112  2     0.0448 *  

summary(m36)
# y = -0.0009898x + 0.8373773

p2$Species <- factor(p2$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                            "Laminaria digitata"))
m36 <- lme(Fucoxanthin ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)

Anova(m36, type = 3)
# Response: Fucoxanthin
#               Chisq Df Pr(>Chisq)    
# (Intercept) 71.4575  1  < 2.2e-16 ***
# Age          1.4817  1     0.2235    
# Species     49.4912  2  1.791e-11 ***
# Age:Species  6.2112  2     0.0448 * 

summary(m36)
# y = -0.003184x + 0.5705735

p2$Species <- factor(p2$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                            "Laminaria hyperborea"))
m36 <- lme(Fucoxanthin ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)

Anova(m36, type = 3)
# Response: Fucoxanthin
#                Chisq Df Pr(>Chisq)    
# (Intercept) 331.1070  1  < 2.2e-16 ***
# Age          10.9886  1  0.0009167 ***
# Species      49.4911  2  1.791e-11 ***
# Age:Species   6.2112  2  0.0447976 * 

summary(m36)
# y = -0.0118439x + 1.241534

p2$Species <- factor(p2$Species, levels = c("Laminaria digitata", "Laminaria hyperborea",
                                            "Saccorhiza polyschides"))
m36 <- lme(Fucoxanthin ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)


m37 <- lme(Fucoxanthin ~ poly(Age, 2, raw = TRUE), 
           random = ~1|ID,
           data = p3)
plot(resid(m37) ~ Age, data = p3) # heterogenous
hist(resid(m37)) # normal

m38 <- lme(Fucoxanthin ~ poly(Age, 2, raw = TRUE), 
           random = ~1|ID,
           weights = varExp(form = ~Age),
           data = p3)
plot(resid(m38, type = "normalized") ~ Age, data = p3) # homogenous
hist(resid(m38)) # normal

Anova(m38)
# Response: Fucoxanthin
#                           Chisq Df Pr(>Chisq)    
# poly(Age, 2, raw = TRUE) 122.39  2  < 2.2e-16 ***  

summary(m38)
# y = 0.088107x - 0.0020334x^2 + 0.5716063



#### 4.5 Antenna size ####
m39 <- lme(Antenna ~ Age * Species, 
           random = ~1|ID,
           data = p2)
boxplot(resid(m39) ~ Species, data = p2)
plot(resid(m39) ~ Age, data = p2) # heterogenous
hist(resid(m39)) # normal

m40 <- lme(Antenna ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)
boxplot(resid(m40, type = "normalized") ~ Species, data = p2)
plot(resid(m40, type = "normalized") ~ Age, data = p2) # homogenous
hist(resid(m40, type = "normalized")) # normal

Anova(m40, type = 3)
# Response: Antenna
#                 Chisq Df Pr(>Chisq)    
# (Intercept) 2969.0135  1  < 2.2e-16 ***
# Age           61.7860  1  3.829e-15 ***
# Species        9.5451  2   0.008459 ** 
# Age:Species   28.0439  2  8.135e-07 ***

summary(m40)
# y = 0.0045716x + 0.8223343

p2$Species <- factor(p2$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                            "Laminaria digitata"))
m40 <- lme(Antenna ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)

Anova(m40, type = 3)
# Response: Antenna
#                 Chisq Df Pr(>Chisq)    
# (Intercept) 3430.8634  1  < 2.2e-16 ***
# Age            0.1642  1   0.685289    
# Species        9.5451  2   0.008459 ** 
# Age:Species   28.0439  2  8.135e-07 ***

summary(m40)
# y = 0.0002341x + 0.8843242

p2$Species <- factor(p2$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                            "Laminaria hyperborea"))
m40 <- lme(Antenna ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)

Anova(m40, type = 3)
# Response: Antenna
#                 Chisq Df Pr(>Chisq)    
# (Intercept) 2977.0099  1  < 2.2e-16 ***
# Age            7.6226  1   0.005764 ** 
# Species        9.5451  2   0.008459 ** 
# Age:Species   28.0439  2  8.135e-07 ***

summary(m40)
# y = 0.0022109x + 0.8336589

p2$Species <- factor(p2$Species, levels = c("Laminaria digitata", "Laminaria hyperborea",
                                            "Saccorhiza polyschides"))
m40 <- lme(Antenna ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)



m41 <- lme(Antenna ~ poly(Age, 2, raw = TRUE),
           random = ~1|ID,
           data = p3)
plot(resid(m41) ~ Age, data = p3) # homogeneity cannot be improved
hist(resid(m41)) # normal


Anova(m41)
# Response: Antenna
#                           Chisq Df Pr(>Chisq)
# poly(Age, 2, raw = TRUE) 1.6523  2     0.4377

summary(m41)
# y = -0.0060991x + 0.0001887x^2 + 0.9685628




#### 4.6 Total pigment ####
m42 <- lme(Total ~ Age * Species,
           random = ~1|ID,
           data = p2)
boxplot(resid(m42) ~ Species, data = p2) 
plot(resid(m42) ~ Age, data = p2) # heterogenous
hist(resid(m42)) # normal

m43 <- lme(Total ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)
boxplot(resid(m43, type = "normalized") ~ Species, data = p2) 
plot(resid(m43, type = "normalized") ~ Age, data = p2) # homogenous
hist(resid(m43, type = "normalized")) # normal

Anova(m43, type = 3)
# Response: Total
#                Chisq Df Pr(>Chisq)    
# (Intercept) 142.2283  1  < 2.2e-16 ***
# Age           0.9295  1    0.33500    
# Species      66.5432  2  3.551e-15 ***
# Age:Species   7.3529  2    0.02531 * 

summary(m43)
# y = -0.0060562x + 1.917934

p2$Species <- factor(p2$Species, levels = c("Laminaria hyperborea", "Saccorhiza polyschides",
                                            "Laminaria digitata"))
m43 <- lme(Total ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)

Anova(m43, type = 3)
# Response: Total
#               Chisq Df Pr(>Chisq)    
# (Intercept) 62.6549  1  2.463e-15 ***
# Age          0.9685  1    0.32505    
# Species     66.5432  2  3.551e-15 ***
# Age:Species  7.3529  2    0.02531 * 

summary(m43)
# y = -0.0061361x + 1.273365

p2$Species <- factor(p2$Species, levels = c("Saccorhiza polyschides", "Laminaria digitata", 
                                            "Laminaria hyperborea"))
m43 <- lme(Total ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)

Anova(m43, type = 3)
# Response: Total
#                Chisq Df Pr(>Chisq)    
# (Intercept) 366.6740  1  < 2.2e-16 ***
# Age          14.2556  1  0.0001596 ***
# Species      66.5432  2  3.551e-15 ***
# Age:Species   7.3529  2  0.0253127 * 

summary(m43)
# y = -0.0320614x + 3.113639

p2$Species <- factor(p2$Species, levels = c("Laminaria digitata", "Laminaria hyperborea",
                                            "Saccorhiza polyschides"))
m43 <- lme(Total ~ Age * Species, 
           random = ~1|ID,
           weights = varExp(form = ~Age|Species),
           data = p2)



m44 <- lme(Total ~ poly(Age, 2, raw = TRUE),
           random = ~1|ID,
           data = p3)
plot(resid(m44) ~ Age, data = p3) # homogenous
hist(resid(m44)) # normal

Anova(m44)
# Response: Total
#                           Chisq Df Pr(>Chisq)    
# poly(Age, 2, raw = TRUE) 80.797  2  < 2.2e-16 ***

summary(m44)
# y = 0.1971408x - 0.0046571x^2 + 1.421288


#### 5. Data visualisation ####
#### 5.1 Chlorophyll a estimates and uncertainty ####
new <- data.frame(Age = c(rep(seq(0, 46, by = 0.1), 2),
                          rep(seq(0, 39, by = 0.1), 2)),
                  Species = c(rep("Laminaria digitata", 461),
                              rep("Laminaria hyperborea", 461),
                              rep("Saccharina latissima", 391),
                              rep("Saccorhiza polyschides", 391)))
new2 <- new[new$Species != "Saccharina latissima",]
new3 <- new[new$Species == "Saccharina latissima",]

fit1 <- predict(m30, level = 0, newdata = new2)
fit2 <- predict(m31, level = 0, newdata = new3)

new$Chla.fit <- c(fit1[1:922], fit2, fit1[923:1313]) # estimate

modmat <- model.matrix(formula(m30)[-2], new2)
fvar <- diag(modmat %*% vcov(m30) %*% t(modmat)) # fixed variance
tvar <- fvar + as.numeric(VarCorr(m30)[1]) # fixed + random variance

lwr1 <- fit1 - qnorm(0.995)*sqrt(fvar) 
upr1 <- fit1 + qnorm(0.995)*sqrt(fvar)

modmat <- model.matrix(formula(m31)[-2], new3)
fvar <- diag(modmat %*% vcov(m31) %*% t(modmat)) # fixed variance
tvar <- fvar + as.numeric(VarCorr(m31)[1]) # fixed + random variance

lwr2 <- fit2 - qnorm(0.995)*sqrt(fvar)
upr2 <- fit2 + qnorm(0.995)*sqrt(fvar)

new$Chla.lwr <- c(lwr1[1:922], lwr2, lwr1[923:1313]) # 99% confidence intervals
new$Chla.upr <- c(upr1[1:922], upr2, upr1[923:1313])

#### 5.2 Chlorophyll c estimates and uncertainty ####
fit1 <- predict(m33, level = 0, newdata = new2)
fit2 <- predict(m34, level = 0, newdata = new3)

new$Chlc.fit <- c(fit1[1:922], fit2, fit1[923:1313]) # estimate

modmat <- model.matrix(formula(m33)[-2], new2)
fvar <- diag(modmat %*% vcov(m33) %*% t(modmat)) # fixed variance
tvar <- fvar + as.numeric(VarCorr(m33)[1]) # fixed + random variance

lwr1 <- fit1 - qnorm(0.995)*sqrt(fvar) 
upr1 <- fit1 + qnorm(0.995)*sqrt(fvar)

modmat <- model.matrix(formula(m34)[-2], new3)
fvar <- diag(modmat %*% vcov(m34) %*% t(modmat)) # fixed variance
tvar <- fvar + as.numeric(VarCorr(m34)[1]) # fixed + random variance

lwr2 <- fit2 - qnorm(0.995)*sqrt(fvar)
upr2 <- fit2 + qnorm(0.995)*sqrt(fvar)

new$Chlc.lwr <- c(lwr1[1:922], lwr2, lwr1[923:1313]) # 99% confidence intervals
new$Chlc.upr <- c(upr1[1:922], upr2, upr1[923:1313])

#### 5.3 Fucoxanthin estimates and uncertainty ####
fit1 <- predict(m36, level = 0, newdata = new2)
fit2 <- predict(m38, level = 0, newdata = new3)

new$Fuco.fit <- c(fit1[1:922], fit2, fit1[923:1313]) # estimate

modmat <- model.matrix(formula(m36)[-2], new2)
fvar <- diag(modmat %*% vcov(m36) %*% t(modmat)) # fixed variance
tvar <- fvar + as.numeric(VarCorr(m36)[1]) # fixed + random variance

lwr1 <- fit1 - qnorm(0.995)*sqrt(fvar) 
upr1 <- fit1 + qnorm(0.995)*sqrt(fvar)

modmat <- model.matrix(formula(m38)[-2], new3)
fvar <- diag(modmat %*% vcov(m38) %*% t(modmat)) # fixed variance
tvar <- fvar + as.numeric(VarCorr(m38)[1]) # fixed + random variance

lwr2 <- fit2 - qnorm(0.995)*sqrt(fvar)
upr2 <- fit2 + qnorm(0.995)*sqrt(fvar)

new$Fuco.lwr <- c(lwr1[1:922], lwr2, lwr1[923:1313]) # 99% confidence intervals
new$Fuco.upr <- c(upr1[1:922], upr2, upr1[923:1313])

#### 5.4 Antenna estimates and uncertainty ####
fit1 <- predict(m40, level = 0, newdata = new2)
fit2 <- predict(m41, level = 0, newdata = new3)

new$Ant.fit <- c(fit1[1:922], fit2, fit1[923:1313]) # estimate

modmat <- model.matrix(formula(m40)[-2], new2)
fvar <- diag(modmat %*% vcov(m40) %*% t(modmat)) # fixed variance
tvar <- fvar + as.numeric(VarCorr(m40)[1]) # fixed + random variance

lwr1 <- fit1 - qnorm(0.995)*sqrt(fvar) 
upr1 <- fit1 + qnorm(0.995)*sqrt(fvar)

modmat <- model.matrix(formula(m41)[-2], new3)
fvar <- diag(modmat %*% vcov(m41) %*% t(modmat)) # fixed variance
tvar <- fvar + as.numeric(VarCorr(m41)[1]) # fixed + random variance

lwr2 <- fit2 - qnorm(0.995)*sqrt(fvar)
upr2 <- fit2 + qnorm(0.995)*sqrt(fvar)

new$Ant.lwr <- c(lwr1[1:922], lwr2, lwr1[923:1313]) # 99% confidence intervals
new$Ant.upr <- c(upr1[1:922], upr2, upr1[923:1313])

#### 5.5 Total pigment estimates and uncertainty ####
fit1 <- predict(m43, level = 0, newdata = new2)
fit2 <- predict(m44, level = 0, newdata = new3)

new$Tot.fit <- c(fit1[1:922], fit2, fit1[923:1313]) # estimate

modmat <- model.matrix(formula(m43)[-2], new2)
fvar <- diag(modmat %*% vcov(m43) %*% t(modmat)) # fixed variance
tvar <- fvar + as.numeric(VarCorr(m43)[1]) # fixed + random variance

lwr1 <- fit1 - qnorm(0.995)*sqrt(fvar) 
upr1 <- fit1 + qnorm(0.995)*sqrt(fvar)

modmat <- model.matrix(formula(m44)[-2], new3)
fvar <- diag(modmat %*% vcov(m44) %*% t(modmat)) # fixed variance
tvar <- fvar + as.numeric(VarCorr(m44)[1]) # fixed + random variance

lwr2 <- fit2 - qnorm(0.995)*sqrt(fvar)
upr2 <- fit2 + qnorm(0.995)*sqrt(fvar)

new$Tot.lwr <- c(lwr1[1:922], lwr2, lwr1[923:1313]) # 99% confidence intervals
new$Tot.upr <- c(upr1[1:922], upr2, upr1[923:1313])

#### 5.6 Add group to separate plot panels by ####
new <- new %>% 
  mutate(Group = ifelse(Species == "Laminaria digitata" | Species == "Laminaria hyperborea", 
                        "Laminaria spp.", 
                        "Other genera"))
new$Group <- factor(new$Group)

# reduce the data.frame length for graphical reasons (dashed geom_line cannot handle highly resolved data)
new.filt <- new %>%
  filter(as.numeric(row.names(new))%%2 == 1) # delete every other row 

p <- p %>% 
  mutate(Group = ifelse(Species == "Laminaria digitata" | Species == "Laminaria hyperborea", 
                        "Laminaria spp.", 
                        "Other genera"))
p$Group <- factor(p$Group)

#### 5.7 Descriptive statistics ####
pavg <- p %>% 
  group_by(ID) %>% 
  summarise(Species = sample(Species, size = 1),
            Group = sample(Group, size = 1),
            Age = mean(Age),
            Chla.se = sd(Chlorophyll.a)/sqrt(length(Chlorophyll.a)),
            Chlorophyll.a = mean(Chlorophyll.a),
            Chlc.se = sd(Chlorophyll.c)/sqrt(length(Chlorophyll.c)),
            Chlorophyll.c = mean(Chlorophyll.c),
            Fuco.se = sd(Fucoxanthin)/sqrt(length(Fucoxanthin)),
            Fucoxanthin = mean(Fucoxanthin),
            Ant.se = sd(Antenna)/sqrt(length(Antenna)),
            Antenna = mean(Antenna),
            Tot.se = sd(Total)/sqrt(length(Total)),
            Total = mean(Total),
            n = length(ID))

tse <- pavg %>% 
  summarise(Chla.se = mean(Chla.se),
            Chlc.se = mean(Chlc.se),
            Fuco.se = mean(Fuco.se),
            Ant.se = mean(Ant.se),
            Tot.se = mean(Tot.se))

# # alternatively use
# with(pavg, lapply(list(Chla.se, Chlc.se, Fuco.se, Ant.se, Tot.se), mean))

Chla.stat <- pavg %>%
  group_by(Age, Species) %>%
  summarise(Species = sample(Species, size = 1),
            Group = sample(Group, size = 1),
            Age = mean(Age),
            mean = mean(Chlorophyll.a),
            sd = sd(Chlorophyll.a),
            se = sd(Chlorophyll.a)/sqrt(length(Chlorophyll.a)),
            n = length(Chlorophyll.a))

Chlc.stat <- pavg %>%
  group_by(Age, Species) %>%
  summarise(Species = sample(Species, size = 1),
            Group = sample(Group, size = 1),
            Age = mean(Age),
            mean = mean(Chlorophyll.c),
            se = sd(Chlorophyll.c)/sqrt(length(Chlorophyll.c)),
            n = length(Chlorophyll.c))

Fuco.stat <- pavg %>%
  group_by(Age, Species) %>%
  summarise(Species = sample(Species, size = 1),
            Group = sample(Group, size = 1),
            Age = mean(Age),
            mean = mean(Fucoxanthin),
            se = sd(Fucoxanthin)/sqrt(length(Fucoxanthin)),
            n = length(Fucoxanthin))


Ant.stat <- pavg %>%
  group_by(Age, Species) %>%
  summarise(Species = sample(Species, size = 1),
            Group = sample(Group, size = 1),
            Age = mean(Age),
            mean = mean(Antenna),
            se = sd(Antenna)/sqrt(length(Antenna)),
            n = length(Antenna))


Tot.stat <- pavg %>%
  group_by(Age, Species) %>%
  summarise(Species = sample(Species, size = 1),
            Group = sample(Group, size = 1),
            Age = mean(Age),
            mean = mean(Total),
            se = sd(Total)/sqrt(length(Total)),
            n = length(Total))



#### 5.8 Plots ####
chlap <- ggplot() +
  geom_point(data = p, mapping = aes(Age, Chlorophyll.a/1000, colour = Species),
             position = position_jitterdodge(dodge.width = 4, jitter.width = 1.5), 
             size = 3, shape = 16, alpha = 0.1) +
  geom_line(data = new.filt, mapping = aes(Age, Chla.fit/1000, colour = Species,
                                      lty = Species)) +
  geom_ribbon(data = new.filt, mapping = aes(Age, ymin = Chla.lwr/1000, ymax = Chla.upr/1000,
                                        fill = Species), alpha = 0.4) +
  geom_pointrange(data = Chla.stat,
                  mapping = aes(Age, mean/1000, ymin = mean/1000 - (qnorm(0.995)*se)/1000,
                                ymax = mean/1000 + (qnorm(0.995)*se)/1000, colour = Species),
                  position = position_dodge(width = 4)) +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                    guide = guide_legend(ncol = 2)) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      labels = c(expression(italic("Laminaria digitata")),
                                 expression(italic("Laminaria hyperborea")),
                                 expression(italic("Saccharina latissima")),
                                 expression(italic("Saccorhiza polyschides"))),
                      guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c(5, 5, 1, 1),
                        guide = "none") +
  facet_grid(~Group) +
  coord_cartesian(xlim = c(-3, 50), ylim = c(0, 4), expand = F) +
  labs(x = "Detrital age (d)",
       y = expression("Chlorophyll "*italic(a)*" (mg g"^-1*")")) +
  mytheme +
  theme(legend.position = c(0.41, 0.925))

chlap

chlcp <- ggplot() +
  geom_point(data = p, mapping = aes(Age, Chlorophyll.c/1000, colour = Species),
             position = position_jitterdodge(dodge.width = 4, jitter.width = 1.5), 
             size = 3, shape = 16, alpha = 0.1) +
  geom_line(data = new.filt, mapping = aes(Age, Chlc.fit/1000, colour = Species,
                                      lty = Species)) +
  geom_ribbon(data = new.filt, mapping = aes(Age, ymin = Chlc.lwr/1000, ymax = Chlc.upr/1000,
                                        fill = Species), alpha = 0.4) +
  geom_pointrange(data = Chlc.stat,
                  mapping = aes(Age, mean/1000, ymin = mean/1000 - (qnorm(0.995)*se)/1000,
                                ymax = mean/1000 + (qnorm(0.995)*se)/1000, colour = Species),
                  position = position_dodge(width = 4)) +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                    guide = guide_legend(ncol = 2)) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      labels = c(expression(italic("Laminaria digitata")),
                                 expression(italic("Laminaria hyperborea")),
                                 expression(italic("Saccharina latissima")),
                                 expression(italic("Saccorhiza polyschides"))),
                      guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c(5, 5, 1, 1),
                        guide = "none") +
  facet_grid(~Group) +
  coord_cartesian(xlim = c(-3, 50), ylim = c(0, 0.3), expand = F) +
  scale_y_continuous(breaks = seq(0, 0.3, by = 0.1)) +
  labs(x = "Detrital age (d)",
       y = expression("Chlorophyll "*italic(c)*" (mg g"^-1*")")) +
  mytheme +
  theme(legend.position = c(0.41, 0.925))

chlcp

fucop <- ggplot() +
  geom_point(data = p, mapping = aes(Age, Fucoxanthin/1000, colour = Species),
             position = position_jitterdodge(dodge.width = 4, jitter.width = 1.5), 
             size = 3, shape = 16, alpha = 0.1) +
  geom_line(data = new.filt, mapping = aes(Age, Fuco.fit/1000, colour = Species,
                                      lty = Species)) +
  geom_ribbon(data = new.filt, mapping = aes(Age, ymin = Fuco.lwr/1000, ymax = Fuco.upr/1000,
                                        fill = Species), alpha = 0.4) +
  geom_pointrange(data = Fuco.stat,
                  mapping = aes(Age, mean/1000, ymin = mean/1000 - (qnorm(0.995)*se)/1000,
                                ymax = mean/1000 + (qnorm(0.995)*se)/1000, colour = Species),
                  position = position_dodge(width = 4)) +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                    guide = guide_legend(ncol = 2)) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      labels = c(expression(italic("Laminaria digitata")),
                                 expression(italic("Laminaria hyperborea")),
                                 expression(italic("Saccharina latissima")),
                                 expression(italic("Saccorhiza polyschides"))),
                      guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c(5, 5, 1, 1),
                        guide = "none") +
  facet_grid(~Group) +
  coord_cartesian(xlim = c(-3, 50), ylim = c(0, 3), expand = F) +
  labs(x = "Detrital age (d)",
       y = expression("Fucoxanthin (mg g"^-1*")")) +
  mytheme +
  theme(legend.position = c(0.41, 0.925))

fucop

antp <- ggplot() +
  geom_hline(yintercept = 1, colour = "#c9d2d7") +
  geom_point(data = p, mapping = aes(Age, Antenna, colour = Species),
             position = position_jitterdodge(dodge.width = 4, jitter.width = 1.5), 
             size = 3, shape = 16, alpha = 0.1) +
  geom_line(data = new.filt, mapping = aes(Age, Ant.fit, colour = Species,
                                      lty = Species)) +
  geom_ribbon(data = new.filt, mapping = aes(Age, ymin = Ant.lwr, ymax = Ant.upr,
                                        fill = Species), alpha = 0.4) +
  geom_pointrange(data = Ant.stat,
                  mapping = aes(Age, mean, ymin = mean - (qnorm(0.995)*se),
                                ymax = mean + (qnorm(0.995)*se), colour = Species),
                  position = position_dodge(width = 4)) +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                    guide = guide_legend(ncol = 2)) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      labels = c(expression(italic("Laminaria digitata")),
                                 expression(italic("Laminaria hyperborea")),
                                 expression(italic("Saccharina latissima")),
                                 expression(italic("Saccorhiza polyschides"))),
                      guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c(1, 5, 5, 1),
                        guide = "none") +
  facet_grid(~Group) +
  coord_cartesian(xlim = c(-3, 50), ylim = c(0.5, 1.25), expand = F) +
  scale_y_continuous(breaks = seq(0.5, 1.25, by = 0.25)) +
  labs(x = "Detrital age (d)",
       y = "Relative antenna size") +
  mytheme +
  theme(legend.position = c(0.41, 0.925))

antp

totp <- ggplot() +
  geom_point(data = p, mapping = aes(Age, Total/1000, colour = Species),
             position = position_jitterdodge(dodge.width = 4, jitter.width = 1.5), 
             size = 3, shape = 16, alpha = 0.1) +
  geom_line(data = new.filt, mapping = aes(Age, Tot.fit/1000, colour = Species,
                                      lty = Species)) +
  geom_ribbon(data = new.filt, mapping = aes(Age, ymin = Tot.lwr/1000, ymax = Tot.upr/1000,
                                        fill = Species), alpha = 0.4) +
  geom_pointrange(data = Tot.stat,
                  mapping = aes(Age, mean/1000, ymin = mean/1000 - (qnorm(0.995)*se)/1000,
                                ymax = mean/1000 + (qnorm(0.995)*se)/1000, colour = Species),
                  position = position_dodge(width = 4)) +
  scale_fill_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                    labels = c(expression(italic("Laminaria digitata")),
                               expression(italic("Laminaria hyperborea")),
                               expression(italic("Saccharina latissima")),
                               expression(italic("Saccorhiza polyschides"))),
                    guide = guide_legend(ncol = 2)) +
  scale_colour_manual(values = c("#333b08", "#627d0e", "#b6b12d", "#5e5003"),
                      labels = c(expression(italic("Laminaria digitata")),
                                 expression(italic("Laminaria hyperborea")),
                                 expression(italic("Saccharina latissima")),
                                 expression(italic("Saccorhiza polyschides"))),
                      guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c(5, 5, 1, 1),
                        guide = "none") +
  facet_grid(~Group) +
  coord_cartesian(xlim = c(-3, 50), ylim = c(0, 6), expand = F) +
  scale_y_continuous(breaks = seq(0, 6, by = 2)) +
  labs(x = "Detrital age (d)",
       y = expression("Total pigment (mg g"^-1*")")) +
  mytheme +
  theme(legend.position = c(0.41, 0.925))

totp

chlap <- chlap + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank())
chlcp <- chlcp + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       strip.text = element_blank(),
                       legend.position = "none")
fucop <- fucop + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       strip.text = element_blank(),
                       legend.position = "none")
antp <- antp + theme(strip.text = element_blank(),
                     legend.position = "none")

comb <- plot_grid(chlap, chlcp, fucop, antp, labels = "auto", label_size = 15, label_fontfamily = "Helvetica Neue",
                  ncol = 1, rel_heights = c(0.98, 0.865, 0.865, 1), align = "v", hjust = 0)
comb
# dimensions: 10 x 6 in


