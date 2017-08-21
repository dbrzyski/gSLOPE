#---------------------------------------------------------------------------
# This code reproduces the figures from simulation study of Figures 3a-c in
# D. Brzyski, A.Gossmann, W. Su, M. Bogdan (2016), "Group SLOPE -
# - adaptive selection of groups of predictors"
# Preprint available at https://arxiv.org/abs/1610.04960
#---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

#--- Figures 3 a,b,c

# load the simulation results
load("../RData/figure_3_simulation_results.RData")
results <- rename(results, n_relevant = n.relevant)

# prepare data for plotting
fractions_df <- results %>%
  group_by(n_relevant) %>%
  summarize_all(sum) %>% tbl_df() %>%
  gather(grp_size, count, -n_relevant) %>%
  mutate(weight = gsub("^(.+)_wt_\\d$", "\\1", grp_size)) %>%
  mutate(grp_size = gsub("^.+_wt_(\\d)$", "\\1", grp_size)) %>%
  group_by(n_relevant, weight) %>%
  mutate(total_count_by_n_relevant = sum(count)) %>%
  mutate(fraction = count / total_count_by_n_relevant) %>%
  tbl_df()

fractions_df %>%
  mutate(weight = factor(weight, levels = c("sqrt", "one", "len"),
                         labels = c("w[i] == sqrt(l[i])",
                                    "w[i] == 1",
                                    "w[i] == l[i]"),
                         ordered = TRUE)) %>%
  mutate(grp_size = paste0(grp_size, "-el. groups")) %>%
  ggplot(aes(n_relevant, fraction, color = grp_size)) +
    geom_line() + geom_point(aes(shape = grp_size)) +
    facet_wrap(~weight, labeller = label_parsed) +
    labs(y = "Fraction of groups in STRG",
         x = "Number of relevant groups") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank())

ggsave("Figure_3abc.png", width = 8, height = 4)
