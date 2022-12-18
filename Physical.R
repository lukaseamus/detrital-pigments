phys <- read.csv("~/Desktop/Projects/QUB/Data/Pigments/Physical.csv")

d <- factor(phys$D) # Categorical time (d)
day <- phys$Day # Continuous time (d)
PAR <- phys$PAR # Photon flux density (µmol photons m-2 s-1)
temp <- phys$Temperature # Temperature (°C)

require(psych)
# Calculate 24-h descriptive statistics
PAR24h.stats <- describeBy(PAR, d, mat = TRUE)
PAR24h.stats$group1 <- as.integer(PAR24h.stats$group1)

T.stats <- describeBy(temp, d, mat = TRUE)
T.stats$group1 <- as.integer(T.stats$group1)

# Calculate daytime descriptive statistics
PARD.stats <- describeBy(PAR[PAR != 0], d[PAR != 0], mat = TRUE)
PARD.stats$group1 <- as.integer(PARD.stats$group1)

# Combine
PAR.stats <- data.frame(rbind(PAR24h.stats, PARD.stats), 
                        time = c(rep("Day and night", 46),
                                 rep("Day", 46)))

require(ggplot2)
require(ggalt)
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
                 legend.key.size = unit(.3, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.spacing.x = unit(.1, "cm"),
                 legend.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title.align = 0,
                 legend.title = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

light <- ggplot() +
  # geom_point(phys, mapping = aes(day, PAR), shape = 16, alpha = 0.005) +
  geom_xspline(PAR.stats, mapping = aes(group1, mean, colour = time), spline_shape = -0.5,
               size = 0.5)  +
  geom_pointrange(PAR.stats, mapping = aes(group1, mean, ymin = mean - qnorm(0.995)*se, 
                                       ymax = mean + qnorm(0.995)*se, colour = time), 
                  fatten = 0.3, size = 0.5) +
  scale_color_manual(values = c("#f7b060", "#2e4a5b")) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 1500), expand = F, clip = "off") +
  labs(y = expression("PAR ("*mu*"mol photons m"^-2*" s"^-1*")"), 
       x = "Time (d)") +
  mytheme +
  theme(legend.position = c(0.85, 0.95))

light

Tplot <- ggplot() +
  # geom_point(phys, mapping = aes(day, temp), shape = 16, alpha = 0.005) +
  geom_xspline(T.stats, mapping = aes(group1, mean), spline_shape = 0.4,
               size = 0.5)  +
  geom_pointrange(T.stats, mapping = aes(group1, mean, ymin = mean - qnorm(0.995)*se, 
                                         ymax = mean + qnorm(0.995)*se), 
                  fatten = 0.3, size = 0.5) +
  coord_cartesian(xlim = c(0, 50), ylim = c(11, 13), expand = F, clip = "off") +
  labs(y = "Temperature (°C)", 
       x = "Time (d)") +
  mytheme

Tplot

require(cowplot)
comb <- plot_grid(light, Tplot, labels = "auto", label_size = 15, 
                  label_fontfamily = "Helvetica Neue", nrow = 1, align = "v", hjust = 0)
comb
# dimensions: 3.2 x 9 in