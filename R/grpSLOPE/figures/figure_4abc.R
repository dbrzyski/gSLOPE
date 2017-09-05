#---------------------------------------------------------------------------
# This code reproduces figures analogous to the ones from
# simulation study of Figures 4a-c in D. Brzyski, A.Gossmann,
# W. Su, M. Bogdan (2016), "Group SLOPE - adaptive selection
# of groups of predictors" (preprint available at https://arxiv.org/abs/1610.04960)
#---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# palette to be used in all plots
library(RColorBrewer)
pal = brewer.pal(3, "Set1")

# legend labels to be used in all plots
math.labs <- list(bquote(list("gSLOPE", q == 0.05)),
                  bquote(list("gSLOPE", q == 0.1)),
                  bquote("gLASSO"["CV"]))

# load the simulation results
load("../RData/figure_4_simulation_results.RData")

#-- Figure analogous to 4a

# some global simulation parameters
n.group <- length(group.id)
max.relevant.groups <- max(n.relevant)

# get means and error bars for FDP and power
results.summary <- results %>%
  select(-replication) %>%
  group_by(n.relevant) %>%
  summarize_all(funs(mean, lwr = (mean(.) - 2*sd(.)/sqrt(n.replications)),
                     upr = (mean(.) + 2*sd(.)/sqrt(n.replications))))

# pre-process for plotting
FDR.results <- results.summary %>% select(n.relevant, ends_with("FDP_mean")) %>%
  gather(scenario, FDR, -n.relevant) %>%
  mutate(scenario = gsub("\\.FDP_mean", "", scenario))
lwr <- results.summary %>% select(n.relevant, ends_with("FDP_lwr")) %>%
  gather(scenario, lwr, -n.relevant) %>%
  mutate(scenario = gsub("\\.FDP_lwr", "", scenario))
upr <- results.summary %>% select(n.relevant, ends_with("FDP_upr")) %>%
  gather(scenario, upr, -n.relevant) %>%
  mutate(scenario = gsub("\\.FDP_upr", "", scenario))
FDR.results <- FDR.results %>% left_join(lwr) %>% left_join(upr)

# draw a figure analogous to 4a
ggplot(FDR.results) +
  geom_segment(mapping = aes(x = 0, xend = max.relevant.groups, y = 0.1, yend = 0.1),
               linetype = 3, color = "black") +
  geom_segment(mapping = aes(x = 0, xend = max.relevant.groups, y = 0.05, yend = 0.05),
               linetype = 3, color = "black") +
  geom_line(mapping = aes(x = n.relevant, y = FDR, color = scenario)) +
  geom_point(mapping = aes(x = n.relevant, y = FDR, shape = scenario, color = scenario)) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr,
                            fill = scenario), alpha=0.25) +
  scale_shape_manual(breaks = c("lambda.05", "lambda.1", "grpLASSO"),
                     labels = c(math.labs[1], math.labs[2], math.labs[3]),
                     values = c(16, 17, 18)) +
  scale_color_manual(breaks = c("lambda.05", "lambda.1", "grpLASSO"),
                     labels = c(math.labs[1], math.labs[2], math.labs[3]),
                     values = pal) +
  scale_fill_manual(breaks = c("lambda.05", "lambda.1", "grpLASSO"),
                    labels = c(math.labs[1], math.labs[2], math.labs[3]),
                    values = pal) +
  xlab("Number of significant groups") + ylab("Estimated gFDR +/- 2SE") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.6, 0.5),
        legend.background = element_rect(color = "black"))

ggsave(file = "Figure_4a_like.png", width = 12, height = 12, units = "cm")

#--- Figure analogous to 4b

# pre-process for plotting
power.results <- results.summary %>% select(n.relevant, ends_with("power_mean")) %>%
  gather(scenario, power, -n.relevant) %>%
  mutate(scenario = gsub("\\.power_mean", "", scenario))
lwr <- results.summary %>% select(n.relevant, ends_with("power_lwr")) %>%
  gather(scenario, lwr, -n.relevant) %>%
  mutate(scenario = gsub("\\.power_lwr", "", scenario))
upr <- results.summary %>% select(n.relevant, ends_with("power_upr")) %>%
  gather(scenario, upr, -n.relevant) %>%
  mutate(scenario = gsub("\\.power_upr", "", scenario))
power.results <- power.results %>% left_join(lwr) %>% left_join(upr)

# plot estimated power with error bars
ggplot(power.results) +
  geom_line(mapping = aes(x = n.relevant, y = power, color = scenario)) +
  geom_point(mapping = aes(x = n.relevant, y = power, color = scenario, shape = scenario)) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr,
                            fill = scenario), alpha=0.25) +
  scale_shape_manual(breaks = c("lambda.05", "lambda.1", "grpLASSO"),
                     labels = c(math.labs[1], math.labs[2], math.labs[3]),
                     values = c(16, 17, 18)) +
  scale_color_manual(breaks = c("lambda.05", "lambda.1", "grpLASSO"),
                     labels = c(math.labs[1], math.labs[2], math.labs[3]),
                     values = pal) +
  scale_fill_manual(breaks = c("lambda.05", "lambda.1", "grpLASSO"),
                    labels = c(math.labs[1], math.labs[2], math.labs[3]),
                    values = pal) +
  xlab("Number of significant groups") + ylab("Estimated power +/- 2SE") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.6, 0.7),
        legend.background = element_rect(color = "black"))

ggsave(file = "Figure_4b_like.png", width = 12, height = 12, units = "cm")

#--- Figure 4c

data_frame(len = sapply(group.id, FUN = length)) %>%
  ggplot(aes(len)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey") +
  xlab("Group size") + ylab("Frequency") +
  theme_bw()

ggsave(file = "Figure_4c.png", width = 12, height = 12, units = "cm")
