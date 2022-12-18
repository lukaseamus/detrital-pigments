# load data
path <- read.csv("~/Desktop/Projects/QUB/Data/Pigments/Pigment deconvolution/pathlength.csv")
eva <- read.csv("~/Desktop/Projects/QUB/Data/Pigments/Pigment deconvolution/evaporation.csv")
evam <- read.csv("~/Desktop/Projects/QUB/Data/Pigments/Pigment deconvolution/evaporation.manual.csv")

# calculate pathlength (cm) with water-peak pathlength equation
cm <- with(path, (A977-A900)/0.18)

# visualise
par(mfrow = c(1,2))
boxplot(cm)
plot(cm)
par(mfrow = c(1,1))

# remove 5 outliers
cm <- cm[-c(20, 21, 40, 50, 262)]

# visualise
par(mfrow = c(1,2))
boxplot(cm)
plot(cm)
par(mfrow = c(1,1))

require(psych)
as.data.frame(describe(cm))

# calculate acetone evaporation relationship
require(lme4)
m1 <- lm(data = eva, mm ~ min)
m2 <- lmer(data = eva, mm ~ min + (1|ID))
m3 <- lmer(data = eva, mm ~ min + (min|ID))
anova(m3, m2, m1) # m2 fits best

plot(m2)
plot(resid(m2) ~ eva$min) # homogenous

par(mfrow = c(1,3))
boxplot(resid(m2), horizontal = T)
hist(resid(m2))
qqnorm(resid(m2))
qqline(resid(m2))
par(mfrow = c(1,1)) # normal

summary(m2)
# y = -0.146786x + 9.823214
# 0.15 mm per min evaporation

eva$fit <- predict(m2, re.form = NA)

modmat <-  model.matrix(terms(m2), data = eva)
int <- diag(modmat %*% tcrossprod(vcov(m2), modmat))
eva$lwr <- with(eva, fit - qnorm(0.975)*sqrt(int))
eva$upr <- with(eva, fit + qnorm(0.975)*sqrt(int))

require(ggplot2)
ggplot(data = eva) + 
  geom_point(aes(min, mm)) + 
  geom_line(aes(min, fit)) +
  geom_ribbon(aes(min, ymin = lwr, ymax = upr), alpha = 0.5) +
  theme_minimal()

# cross-validate slope with image data
mean(evam$mm.min) # 0.13 mm per min evaporation
# very similar so the fluorescence measurement worked

# calculate median time after which half the wells in a plate are filled
t <- c(rep(1.468394, 3), rep(1.824369, 3), 1.846617, rep(1.868866, 2), 2.113598,
       2.269337, rep(2.959037, 4), 3.804476, 5.228374) + 0.6666667 # 0.67 min is the plate loading time

par(mfrow = c(1,2))
boxplot(t)
plot(t)
par(mfrow = c(1,1))

median(t)/2*0.146786 # 0.1860904 = the length (mm) by which the pathlength of the middle well
# is reduced after the median pipetting time (min)

mean(t)/2*0.146786 # 0.2247031 = the length (mm) by which the pathlength of the middle well
# is reduced after the mean pipetting time (min)

0.6144484-(0.2247031/10) # 0.5919781 cm = pathlength adjusted for acetone evaporation 
