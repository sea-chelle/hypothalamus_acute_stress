# load required packages
library(here)
library(ggplot2)
library(car)
library(emmeans)
library(grid)

# set working directory
setwd(here("plots"))

# extract grey60 module eigengene expression for each sample
traits$MEgrey60 <- bwnet$MEs[, "MEgrey60"]

# convert experimental variables to factors
traits$sex <- factor(traits$sex)
traits$condition <- factor(traits$condition,
                           levels = c("basal", "10", "60"))

# set sum-to-zero contrasts for Type III ANOVA
contrasts(traits$sex) <- contr.sum(2)
contrasts(traits$condition) <- contr.sum(3)

# test effects of sex, ARS duration and their interaction
grey60_lm <- lm(MEgrey60 ~ sex * condition, data = traits)
Anova(grey60_lm, type = 3)

# compare ARS conditions averaged across sexes
emmeans(grey60_lm, pairwise ~ condition, adjust = "tukey")

# compare ARS conditions separately within each sex
emmeans(grey60_lm, pairwise ~ condition | sex, adjust = "tukey")
traits$MEgrey60 <- bwnet$MEs[, "MEgrey60"]

# visualise sex-specific changes in module eigengene expression across ARS
ggplot(traits, aes(x = factor(condition, levels = c("basal", "10", "60")), y = MEgrey60, colour = sex, group = sex)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12, linewidth = 0.5) +
  scale_colour_manual(values = c("m" = "#78A8D1", "f" = "#D78FA9"), labels = c("m" = "Male", "f" = "Female")) +
  labs(x = "Acute restraint stress (min)", y = "Grey60 module\neigengene expression", colour = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"),
        panel.spacing = unit(1, "lines")) +
  guides(colour = guide_legend(override.aes = list(linewidth = .5, size = 2)))

ggsave("fig4.10d_gre90_sex.tiff", width = 3, height = 2, dpi = 600, compression = "lzw")
