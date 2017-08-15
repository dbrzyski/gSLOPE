#---------------------------------------------------------------------------
# This code reproduces the figures from simulation study of Figures 2a-e in
# D. Brzyski, A.Gossmann, W. Su, M. Bogdan (2016), "Group SLOPE -
# - adaptive selection of groups of predictors"
# Preprint available at https://arxiv.org/abs/1610.04960
#---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# palette to be used in all plots
library(RColorBrewer)
pal = brewer.pal(3, "Set1")

# legend labels to be used in all plots
math.labs <- list(bquote(list(q == 0.05, lambda^.("max"))),
                  bquote(list(q == 0.1, lambda^.("max"))),
                  bquote(list(q == 0.05, lambda^.("mean"))),
                  bquote(list(q == 0.1, lambda^.("mean"))))

#--- Figure 2a

# Load simulation results for Figure 2a
load("../RData/figure_2a_simulation_results.RData")

# some global simulation parameters
n.group <- 1000
n.replications <- max(results$replication)
max.relevant.groups <- max(results$n.relevant)

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

# draw Figure 2a
FDR.results %>% ggplot() +
  geom_segment(mapping = aes(x = 0, xend = max.relevant.groups, y = 0.1,
                             yend = 0.1*(n.group - max.relevant.groups)/n.group),
               linetype = 3, color = "black") +
  geom_segment(mapping = aes(x = 0, xend = max.relevant.groups, y = 0.05,
                             yend = 0.05*(n.group - max.relevant.groups)/n.group),
               linetype = 3, color = "black") +
  geom_line(mapping = aes(x = n.relevant, y = FDR, linetype = scenario), color = pal[3]) +
  geom_point(mapping = aes(x = n.relevant, y = FDR, shape = scenario), color = pal[3]) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr, group = scenario),
              fill = pal[3], alpha=0.25) +
  scale_shape_manual(breaks = c("lambda.max.1", "lambda.max.05"),
                     labels = c(math.labs[2], math.labs[1]),
                     values = c(16, 17)) +
  scale_linetype_manual(breaks = c("lambda.max.1", "lambda.max.05"),
                        labels = c(math.labs[2], math.labs[1]),
                        values = c(1, 2)) +
  xlab("Number of significant groups") + ylab("Estimated gFDR +/- 2SE") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.85),
        legend.background = element_rect(color = "black"))

ggsave(file = "Figure_2a.png", width = 12, height = 12, units = "cm")


#--- Figures 2bc

# load the simulation results for Figures 2b-e
load("../RData/figure_2bcde_simulation_results.RData")

# some global simulation parameters
n.group <- 1000
n.replications <- max(results$replication)
max.relevant.groups <- max(results$n.relevant)

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

# draw Figure 2b
ggplot(FDR.results[grep("max", FDR.results$scenario), ]) +
  geom_segment(mapping = aes(x = 0, xend = max.relevant.groups, y = 0.1,
                             yend = 0.1*(n.group - max.relevant.groups)/n.group),
               linetype = 3, color = "black") +
  geom_segment(mapping = aes(x = 0, xend = max.relevant.groups, y = 0.05,
                             yend = 0.05*(n.group - max.relevant.groups)/n.group),
               linetype = 3, color = "black") +
  geom_line(mapping = aes(x = n.relevant, y = FDR, linetype = scenario), color = pal[1]) +
  geom_point(mapping = aes(x = n.relevant, y = FDR, shape = scenario), color = pal[1]) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr, group = scenario),
              fill = pal[1], alpha=0.25) +
  scale_shape_manual(breaks = c("lambda.max.1", "lambda.max.05"),
                     labels = c(math.labs[2], math.labs[1]),
                     values = c(16, 17)) +
  scale_linetype_manual(breaks = c("lambda.max.1", "lambda.max.05"),
                        labels = c(math.labs[2], math.labs[1]),
                        values = c(1, 2)) +
  xlab("Number of significant groups") + ylab("Estimated gFDR +/- 2SE") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.85),
        legend.background = element_rect(color = "black"))

ggsave(file = "Figure_2b.png", width = 12, height = 12, units = "cm")

# draw Figure 2c
ggplot(FDR.results[grep("mean", FDR.results$scenario), ]) +
  geom_segment(mapping = aes(x = 0, xend = max.relevant.groups, y = 0.1,
                             yend = 0.1*(n.group - max.relevant.groups)/n.group),
               linetype = 3, color = "black") +
  geom_segment(mapping = aes(x = 0, xend = max.relevant.groups, y = 0.05,
                             yend = 0.05*(n.group - max.relevant.groups)/n.group),
               linetype = 3, color = "black") +
  geom_line(mapping = aes(x = n.relevant, y = FDR, linetype = scenario), color = pal[2]) +
  geom_point(mapping = aes(x = n.relevant, y = FDR, shape = scenario), color = pal[2]) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr, group = scenario),
              fill = pal[2], alpha=0.25) +
  scale_shape_manual(breaks = c("lambda.mean.1", "lambda.mean.05"),
                     labels = c(math.labs[4], math.labs[3]),
                     values = c(16, 17)) +
  scale_linetype_manual(breaks = c("lambda.mean.1", "lambda.mean.05"),
                        labels = c(math.labs[4], math.labs[3]),
                        values = c(1, 2)) +
  xlab("Number of significant groups") + ylab("Estimated gFDR +/- 2SE") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.85),
        legend.background = element_rect(color = "black"))

ggsave(file = "Figure_2c.png", width = 12, height = 12, units = "cm")

#--- Figure 2d

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
  scale_shape_manual(breaks = c("lambda.max.05", "lambda.max.1",
                                "lambda.mean.05", "lambda.mean.1"),
                     labels = math.labs,
                     values = c(17, 16, 17, 16)) +
  scale_color_manual(breaks = c("lambda.max.05", "lambda.max.1",
                                "lambda.mean.05", "lambda.mean.1"),
                     labels = math.labs,
                     values = c(pal[1], pal[1], pal[2], pal[2])) +
  scale_fill_manual(breaks = c("lambda.max.05", "lambda.max.1",
                               "lambda.mean.05", "lambda.mean.1"),
                    labels = math.labs,
                    values = c(pal[1], pal[1], pal[2], pal[2])) +
  xlab("Number of significant groups") + ylab("Estimated power +/- 2SE") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.3),
        legend.background = element_rect(color = "black"))

ggsave(file = "Figure_2d.png", width = 12, height = 12, units = "cm")

#--- Figure 2e

library(grpSLOPE)

# grouping structure used in the simulations
group <- c(rep(1:200, each=3),
           rep(201:400, each=4),
           rep(401:600, each=5),
           rep(601:800, each=6),
           rep(801:1000, each=7))
group.id <- getGroupID(group)
group.length <- sapply(group.id, FUN=length)
n.group <- length(group.id)
p <- sum(group.length)

# compute the lambda sequences
lambda.mean.1 <- lambdaGroupSLOPE(fdr = 0.1, group = group,
                                  wt = sqrt(group.length), method="mean")
lambda.mean.05 <- lambdaGroupSLOPE(fdr = 0.05, group = group,
                                   wt = sqrt(group.length), method="mean")
lambda.max.1 <- lambdaGroupSLOPE(method = "max", fdr = 0.1,
                                 group = group, wt = sqrt(group.length))
lambda.max.05 <- lambdaGroupSLOPE(method = "max", fdr = 0.05,
                                  group = group, wt = sqrt(group.length))

lambdas_df <- data_frame(lambdas = c(lambda.max.05, lambda.max.1,
                                     lambda.mean.05, lambda.mean.1),
                         setting = factor(rep(1:4, each = n.group)),
                         index = rep(1:n.group, 4))

# plot figure 2e
lambdas_df %>% filter(index <= 500) %>%
  ggplot(aes(index, lambdas, color = setting, linetype = setting)) +
    geom_line() +
    scale_color_manual(breaks = 1:4, labels = math.labs,
                       values = c(pal[1], pal[1], pal[2], pal[2])) +
    scale_linetype_manual(breaks =1:4, labels = math.labs,
                          values = c(2, 1, 2, 1)) +
    xlab(expression(i)) + ylab(expression(lambda[i])) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.8, 0.8),
          legend.background = element_rect(color = "black"))

ggsave(file = "Figure_2e.png", width = 12, height = 12, units = "cm")
