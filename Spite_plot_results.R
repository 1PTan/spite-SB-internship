# This code has been written by Pedro Batista Tan, for his major internship
# of the master program in Systems biology and Bioinformatics at Vrije
# Universiteit Amsterdam. The internship was under supervision from Rutger Hermsen,
# from the Theoretical biology group of the University of Utrecht.

# The script uses simulation data from the model parsed in csv files from python scripts
# It has been used to generate the final figures in the report

library(tidyverse)
library(latex2exp)
library(patchwork)

#--------------------------------------------------------------
# Fig 1 Stepwise model parameter sweep

# Import csv results
results_uu <- read.csv("Results_csv/Spite_results_df_uu.csv")

# Convert values below -1 to -1, since they can freely drift in negative values
# this was done just to have a better dynamic range for the color scale
results_uu$average_mean_spite <- pmax(-1, results_uu$average_mean_spite)

# Separate numeric diff_sigma and well-mixed ("i -> infty")
results_uu_n <- results_uu %>% filter(diff_sigma != "i")
results_uu_i <- results_uu %>% filter(diff_sigma == "i")

# Convert sigma columns to numeric
results_uu_n$spite_sigma <- results_uu_n$spite_sigma %>% as.character() %>% as.numeric()
results_uu_n$comp_sigma <- results_uu_n$comp_sigma %>% as.character() %>% as.numeric()
results_uu_n$diff_sigma <- results_uu_n$diff_sigma %>% as.character() %>% as.numeric()

results_uu_i$spite_sigma <- results_uu_i$spite_sigma %>% as.character() %>% as.numeric()
results_uu_i$comp_sigma <- results_uu_i$comp_sigma %>% as.character() %>% as.numeric()

# Create columns for the sigma ratios in the numeric dataframe
results_uu_n <- results_uu_n %>% mutate(
  spite_diff_ratio = spite_sigma/diff_sigma, 
  comp_diff_ratio = comp_sigma/diff_sigma)

# Add a column for the diff_sigma Latex labels
results_uu_n$diff_sigma_label <- factor(
  results_uu_n$diff_sigma,
  levels = c("4", "8", "16"),
  labels = c(TeX("$\\sigma_{d} = 4$"), 
             TeX("$\\sigma_{d} = 8$"), TeX("$\\sigma_{d} = 16$")))

results_uu_i$diff_sigma_label <- factor(
  results_uu_i$diff_sigma,
  levels = c("i"),
  labels = c(TeX("$\\sigma_{d} = \\infty$")))

# Fig 1A
results_uu_n  %>% ggplot(
  aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_mean_spite)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + 
  scale_fill_distiller(palette = 'RdBu', limit = c(-1,1) * max(abs(results_uu$average_mean_spite)))+
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       fill = TeX("Average $\\phi$")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 14), legend.position = "none", strip.text.x = element_text(
          size = 14))
png(filename="uu_sweep_cost0.005_f.png", width=600, height=600)
ggsave(file="uu_sweep_cost0.005_f.png", width=6, height=3, dpi=300)

# Fig 1B
results_uu_i %>% ggplot(
  aes(x=log2(as.numeric(as.character(spite_sigma))), y= log2(as.numeric(as.character(comp_sigma))), fill = average_mean_spite)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(-1,1) * max(abs(results_uu$average_mean_spite))) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       fill = TeX("Average $\\phi$")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 14), strip.text.x = element_text(
          size = 14))
png(filename="uu_sweep_i_cost0.005_f.png", width=350, height=500)
ggsave(file="uu_sweep_i_cost0.005_f.png", width=3.5, height=2.5, dpi=300)

# Fig 1C
# Plot the number of individuals depending on sigma_c
results_uu_n %>% ggplot(aes(x = log2(comp_sigma), y = average_n_individuals)) + geom_point() +
  scale_x_continuous(labels = c("0", "1", "2", "3", "4", "5", "6"), 
                     breaks = c(0,1,2,3,4,5,6)) + 
scale_y_continuous(
                   breaks = c(0,1,2,3,4,5)*10000, limits = c(0, 50000)) + 
  labs(x = TeX('$log_2(\\sigma_c)$'), y = "Population Size N") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=14), legend.title = element_text(size=14), 
        legend.text = element_text(size = 14))
png(filename="uu_sweep_nind_f.png", width=250, height=250)
ggsave(file="uu_sweep_nind_f.png", width=2.5, height=2.5, dpi=300)


#--------------------------------------------------------------
# Fig 2 Effects of diffusion, selfspite and selfspite + resistance

# The following steps are done for three sets of results: 
# Import csv results
# Convert diff sigma column to numeric
# Create a column containing the combination of spite and competition sigmas

results_uu_diff <- read.csv("Results_csv/Spite_results_df_uu_diff.csv")
results_uu_diff$diff_sigma <- results_uu_diff$diff_sigma %>% as.character() %>% as.numeric()
results_uu_diff <- results_uu_diff %>% mutate(spite_comp =  paste("spite:", spite_sigma, "comp:", comp_sigma))

results_selfspite <- read.csv("Results_csv/Spite_results_uu_selfspite.csv")
results_selfspite$diff_sigma <- results_selfspite$diff_sigma %>% as.character() %>% as.numeric()
results_selfspite <- results_selfspite %>% mutate(spite_comp =  paste("spite:", spite_sigma, "comp:", comp_sigma))

results_selfspite_res <- read.csv("Results_csv/Spite_results_selfspite_res.csv")
results_selfspite_res$diff_sigma <- results_selfspite_res$diff_sigma %>% as.character() %>% as.numeric()
results_selfspite_res <- results_selfspite_res %>% mutate(spite_comp =  paste("spite:", spite_sigma, "comp:", comp_sigma))

# Fig 2A
results_uu_diff %>% ggplot(aes(x = diff_sigma, y = average_mean_spite, color = spite_comp)) +
  geom_point() + geom_hline(yintercept = 0, color = "black") + geom_line() +
  labs(x = TeX('$\\sigma_d$'), y =  TeX('Average $\\phi$'),
       color = "") +
  scale_x_continuous(trans='log2', labels = c("1", "2", "4", "8", "16"), 
                     breaks = c(1,2,4,8,16), limits = c(1, 16)) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14)) + ylim(c(-1.25, 1))
#png(filename="uu_diff_sweep_cost0.005.png", width=600, height=600)
#ggsave(file="uu_diff_sweep_cost0.005.png", width=6, height=3, dpi=300)

# Fig 2B
results_selfspite %>% ggplot(aes(x = diff_sigma, y = average_mean_spite, color = spite_comp)) +
  geom_point() + geom_hline(yintercept = 0, color = "black") + geom_line() +
  labs(x = TeX('$\\sigma_d$'), y =  TeX('Average $\\phi$'),
       color = "") +
  scale_x_continuous(trans='log2', labels = c("1", "2", "4", "8", "16"), 
                     breaks = c(1,2,4,8,16), limits = c(1, 16)) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14)) + ylim(c(-2.25, 1))
#png(filename="uu_diff_sweep_cost0.005_selfspite.png", width=600, height=600)
#ggsave(file="uu_diff_sweep_cost0.005_selfspite.png", width=6, height=3, dpi=300)

# Fig 2C
results_selfspite_res %>% ggplot(aes(x = diff_sigma, 
                  y = average_mean_spite, color = spite_comp, 
                  group = interaction(cost, spite_comp))) +
  geom_point(aes(shape = as.factor(cost))) + geom_line() +
  geom_hline(yintercept = 0, color = "black") + 
  labs(x = TeX('$\\sigma_d$'), y =  TeX('Average $\\phi$'), shape = "Cost",
       color = "") +
  scale_x_continuous(trans='log2', labels = c("1", "2", "4", "8", "16"), 
                     breaks = c(1,2,4,8,16), limits = c(1, 16)) + 
  scale_y_continuous(labels = c("-2", "0", "2", "4", "6", "8"), 
                     breaks = c(-2,0,2,4,6, 8), limits = c(-2, 8.5)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))
#png(filename="uu_diff_sweep_selfspite_res.png", width=600, height=600)
#ggsave(file="uu_diff_sweep_selfspite_res.png", width=6, height=3, dpi=300)

# Fig 2 A,B,C Combined plot
results_uu_diff %>% ggplot(aes(x = diff_sigma, y = average_mean_spite, color = spite_comp)) +
  geom_point() + geom_hline(yintercept = 0, color = "black") + geom_line() +
  labs(x = TeX('$\\sigma_d$'), y =  TeX('Average $\\phi$'),
       color = "", title = "No Self-spite") + 
  scale_x_continuous(trans='log2', labels = c("1", "2", "4", "8", "16"), 
                     breaks = c(1,2,4,8,16), limits = c(1, 16)) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), legend.position = "none", strip.text.x = element_text(
          size = 14)) + ylim(c(-2.25, 1)) +
  
results_selfspite %>% ggplot(aes(x = diff_sigma, y = average_mean_spite, color = spite_comp)) +
  geom_point() + geom_hline(yintercept = 0, color = "black") + geom_line() +
  labs(x = TeX('$\\sigma_d$'), y =  TeX('Average $\\phi$'),
       color = "", title = "Self-spite") +  
  scale_x_continuous(trans='log2',        labels = c("1", "2", "4", "8", "16"), 
                     breaks = c(1,2,4,8,16), limits = c(1, 16)) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), legend.position = "none", strip.text.x = element_text(
          size = 14)) + ylim(c(-2.25, 1)) +
  
  
results_selfspite_res %>% ggplot(aes(x = diff_sigma, 
                                     y = average_mean_spite, color = spite_comp, 
                                     group = interaction(cost, spite_comp))) +
  geom_point(aes(shape = as.factor(cost))) + geom_line() + 
  geom_hline(yintercept = 0, color = "black") + 
  labs(x = TeX('$\\sigma_d$'), y =  TeX('Average $\\phi$'), shape = "Cost",
       color = "", title = "Resistance") +
  scale_x_continuous(trans='log2', labels = c("1", "2", "4", "8", "16"), 
                     breaks = c(1,2,4,8,16), limits = c(1, 16)) + 
  scale_y_continuous(labels = c("-2", "0", "2", "4", "6", "8"), 
                     breaks = c(-2,0,2,4,6, 8), limits = c(-2, 9)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14)) +  guides(color = guide_legend(override.aes = list(shape = NA, size = 1)))
png(filename="uu_diff_combined_f.png", width=800, height=250)
ggsave(file="uu_diff_combined_f.png", width=8, height=2.5, dpi=300)

#-------------------------------------------------------------
# Fig 3 Dichotomous Model Sweep

# Import csv results
results_dic_sweep_l <- read.csv("Results_csv/Results_dichotomous_l.csv")

# Convert columns to factors
results_dic_sweep_l$spite_expression_prob <- results_dic_sweep_l$spite_expression_prob %>% as.factor()
results_dic_sweep_l$kd <- results_dic_sweep_l$kd %>% as.factor()
results_dic_sweep_l$cost <- results_dic_sweep_l$cost %>% as.factor()

# Separate numeric diff_sigma and well-mixed ("i -> infty")
results_dic_sweep_n <- results_dic_sweep_l %>% filter(diff_sigma != "i")
results_dic_sweep_i <- results_dic_sweep_l %>% filter(diff_sigma == "i")

# Convert sigma columns to numeric
results_dic_sweep_n$spite_sigma <- results_dic_sweep_n$spite_sigma %>% as.character() %>% as.numeric()
results_dic_sweep_n$comp_sigma <- results_dic_sweep_n$comp_sigma %>% as.character() %>% as.numeric()
results_dic_sweep_n$diff_sigma <- results_dic_sweep_n$diff_sigma %>% as.character() %>% as.numeric()

results_dic_sweep_i$spite_sigma <- results_dic_sweep_i$spite_sigma %>% as.character() %>% as.numeric()
results_dic_sweep_i$comp_sigma <- results_dic_sweep_i$comp_sigma %>% as.character() %>% as.numeric()

# Create columns for the sigma ratios in the numeric dataframe
results_dic_sweep_n <- results_dic_sweep_n %>% mutate(
  spite_diff_ratio = spite_sigma/diff_sigma, 
  comp_diff_ratio = comp_sigma/diff_sigma)

# Add a column for the diff_sigma Latex labels
results_dic_sweep_n$diff_sigma_label <- factor(
  results_dic_sweep_n$diff_sigma,
  levels = c("4", "8", "16"),
  labels = c(TeX("$\\sigma_{d} = 4$"), 
             TeX("$\\sigma_{d} = 8$"), TeX("$\\sigma_{d} = 16$")))

results_dic_sweep_i$diff_sigma_label <- factor(
  results_dic_sweep_i$diff_sigma,
  levels = c("i"),
  labels = c(TeX("$\\sigma_{d} = \\infty$")))

# Plot the number of individuals depending on sigma_c
results_dic_sweep_n %>% ggplot(aes(x = log2(comp_sigma), y = average_n_individuals, color = average_spite_freq)) + 
  geom_point(aes(shape = as.factor(spite_expression_prob)), size = 3) +
  scale_x_continuous(labels = c("0", "1", "2", "3", "4", "5", "6"), breaks = c(0,1,2,3,4,5,6)) + 
  labs(x = TeX('$log_2(\\sigma_c)$'), y = "N individuals",
       colour = "Carrier Frequency", shape = "Spite expression") + ylim(0, 50000) + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))
png(filename="dichotomous_diff_sweep_nind_l.png", width=500, height=350)
ggsave(file="dichotomous_diff_sweep_nind_l.png", width=5, height=3.5, dpi=300)

# Figure 3A
#  Plot the dichotomous diffusion sweep with Cost 0.1 
# Individual plots are combined with patchwork "+"
results_dic_sweep_n  %>% filter(cost==0.1) %>% filter(spite_expression_prob==0.001) %>% 
  ggplot(aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_dic_sweep_n  %>% filter(cost==0.1) %>% filter(spite_expression_prob==0.01) %>% 
  ggplot(aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_dic_sweep_n  %>% filter(cost==0.1) %>% filter(spite_expression_prob==0.05) %>% 
  ggplot(aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'), 
       title = "Spite expression 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_dic_sweep_n  %>% filter(cost==0.1) %>% filter(spite_expression_prob==0.2) %>%
  ggplot(aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))
png(filename="dichotomous_diff_sweep_kd0.001_cost0.1_l.png", width=800, height=400)
ggsave(file="dichotomous_diff_sweep_kd0.001_cost0.1_l.png", width=8, height=4, dpi=300)

# Figure 3B
#  Plot the well-mixed dichotomous sweep with Cost 0.1
# Individual plots are combined with patchwork "+"

# Dichotomous diffusion sweep Cost 0.1 Well mixed diff sigma = "i" 
results_dic_sweep_i  %>% filter(cost==0.1) %>% filter(spite_expression_prob==0.001) %>% ggplot(
  aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_dic_sweep_i  %>% filter(cost==0.1) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_dic_sweep_i  %>% filter(cost==0.1) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_dic_sweep_i  %>% filter(cost==0.1) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))
png(filename="dichotomous_diff_i_kd0.001_cost0.1_l.png", width=450, height=400)
ggsave(file="dichotomous_diff_i_kd0.001_cost0.1_l.png", width=4.5, height=4, dpi=300)


#-------------------------------------------------------------
# Fig 4  Different initial frequencies f_0

# Import csv results. The csv also contains results for cost = 0.01
results_freq_l <- read.csv("Results_csv/Results_dichotomous_l_freq.csv")

# Separate numeric diff_sigma and well-mixed ("i -> infty")
results_freq_n <- results_freq_l %>% filter(diff_sigma != "i")
results_freq_i <- results_freq_l %>% filter(diff_sigma == "i")

# Convert sigma columns to numeric
results_freq_n$spite_sigma <- results_freq_n$spite_sigma %>% as.character() %>% as.numeric()
results_freq_n$comp_sigma <- results_freq_n$comp_sigma %>% as.character() %>% as.numeric()
results_freq_n$diff_sigma <- results_freq_n$diff_sigma %>% as.character() %>% as.numeric()

results_freq_i$spite_sigma <- results_freq_i$spite_sigma %>% as.character() %>% as.numeric()
results_freq_i$comp_sigma <- results_freq_i$comp_sigma %>% as.character() %>% as.numeric()

# Create columns for the sigma ratios in the numeric dataframe
results_freq_n <- results_freq_n %>% mutate(
  spite_diff_ratio = spite_sigma/diff_sigma, 
  comp_diff_ratio = comp_sigma/diff_sigma)

# Latex labels for diff_sigma
results_freq_n$diff_sigma_label <- factor(
  results_freq_n$diff_sigma,
  levels = c("4", "8", "16"),
  labels = c(TeX("$\\sigma_{d} = 4$"), 
             TeX("$\\sigma_{d} = 8$"), TeX("$\\sigma_{d} = 16$")))

results_freq_i$diff_sigma_label <- factor(
  results_freq_i$diff_sigma,
  levels = c("i"),
  labels = c(TeX("$\\sigma_{d} = \\infty$")))

# Dichotomous diffusion sweep n individuals freq
results_freq_n %>% ggplot(aes(x = log2(comp_sigma), y = average_n_individuals, color = average_spite_freq)) + 
  geom_point(aes(shape = as.factor(spite_expression_prob)), size = 2) +
  scale_x_continuous(labels = c("0", "1", "2", "3", "4", "5", "6"), breaks = c(0,1,2,3,4,5,6)) + 
  labs(x = TeX('$log_2(\\sigma_c)$'), y = "N individuals",
       colour = "Carrier Frequency", shape = "Spite expression") + ylim(0, 50000)
#png(filename="dichotomous_diff_sweep_nind_freq_l.png", width=600, height=400)
#ggsave(file="dichotomous_diff_sweep_nind_freq_l.png", width=6, height=4, dpi=300)

#----------------------------------
# Fig 4A
# Starting freq 0.25 cost 0.1
results_freq_n  %>% filter(initial_spite_frequency==0.25) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.001) %>% ggplot(
  aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none")  +
  
  results_freq_n  %>% filter(initial_spite_frequency==0.25) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none")  +
  
  results_freq_n  %>% filter(initial_spite_frequency==0.25) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'), 
       title = "Spite expression 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none")  +
  
  results_freq_n  %>% filter(initial_spite_frequency==0.25) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.2", fill = "Carrier Frequency")+ 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14)) 
png(filename="dichotomous_diff_sweep_f0_0.25_cost0.1_l.png", width=800, height=400)
ggsave(file="dichotomous_diff_sweep_f0_0.25_cost0.1_l.png", width=8, height=4, dpi=300)

# Fig 4B
# Starting freq 0.25 cost 0.1
# Well mixed diff sigma = "i"
results_freq_i  %>% filter(initial_spite_frequency==0.25) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.001) %>% ggplot(
  aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.001", fill = "Carrier Frequency")  + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none")  +
  
  results_freq_i  %>% filter(initial_spite_frequency==0.25) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none")  +
  
  results_freq_i  %>% filter(initial_spite_frequency==0.25) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none")  +
  
  results_freq_i  %>% filter(initial_spite_frequency==0.25) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.2", fill = "Carrier Frequency")  + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14)) 
png(filename="dichotomous_diff_i_f0_0.25_cost0.1_l.png", width=4.5, height=400)
ggsave(file="dichotomous_diff_i_f0_0.25_cost0.1_l.png", width=4.5, height=4, dpi=300)

#----------------------------------
# Fig 4C
# Starting freq 0.02 cost 0.1
results_freq_n  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.001) %>% ggplot(
  aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_freq_n  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_freq_n  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'), 
       title = "Spite expression 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_freq_n  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))
png(filename="dichotomous_diff_sweep_f0_0.02_cost0.1_l.png", width=800, height=400)
ggsave(file="dichotomous_diff_sweep_f0_0.02_cost0.1_l.png", width=8, height=4, dpi=300)

# Fig 4D
# Starting freq 0.02 cost 0.1
# Well mixed diff sigma = "i"
results_freq_i  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.001) %>% ggplot(
  aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_freq_i  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_freq_i  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.05", fill = "Carrier Frequency")  + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_freq_i  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))
png(filename="dichotomous_diff_i_f0_0.02_cost0.1_l.png", width=450, height=400)
ggsave(file="dichotomous_diff_i_f0_0.02_cost0.1_l.png", width=4.5, height=4, dpi=300)

#-------------------------------
# Fig 5 higher kd 0.002

# Import csv results
results_kd <- read.csv("Results_csv/Results_dichotomous_l_kd.csv")

# Separate numeric diff_sigma and well-mixed ("i -> infty")
results_kd_n <- results_kd %>% filter(diff_sigma != "i")
results_kd_i <- results_kd %>% filter(diff_sigma == "i")

# Convert sigma columns to numeric
results_kd_n$spite_sigma <- results_kd_n$spite_sigma %>% as.character() %>% as.numeric()
results_kd_n$comp_sigma <- results_kd_n$comp_sigma %>% as.character() %>% as.numeric()
results_kd_n$diff_sigma <- results_kd_n$diff_sigma %>% as.character() %>% as.numeric()

results_kd_i$spite_sigma <- results_kd_i$spite_sigma %>% as.character() %>% as.numeric()
results_kd_i$comp_sigma <- results_kd_i$comp_sigma %>% as.character() %>% as.numeric()

# Create columns for the sigma ratios in the numeric dataframe
results_kd_n <- results_kd_n %>% mutate(
  spite_diff_ratio = spite_sigma/diff_sigma, 
  comp_diff_ratio = comp_sigma/diff_sigma)


# Create columns for the sigma ratios in the numeric dataframe
results_kd_n$diff_sigma_label <- factor(
  results_kd_n$diff_sigma,
  levels = c("4", "8", "16"),
  labels = c(TeX("$\\sigma_{d} = 4$"), 
             TeX("$\\sigma_{d} = 8$"), TeX("$\\sigma_{d} = 16$")))

results_kd_i$diff_sigma_label <- factor(
  results_kd_i$diff_sigma,
  levels = c("i"),
  labels = c(TeX("$\\sigma_{d} = \\infty$")))

# Plot the number of individuals depending on sigma_c
results_kd_n %>% ggplot(aes(x = log2(comp_sigma), y = average_n_individuals, color = average_spite_freq)) + geom_point() +
  scale_x_continuous(labels = c("0", "1", "2", "3", "4", "5", "6"), breaks = c(0,1,2,3,4,5,6)) + 
  labs(x = TeX('$log_2(\\sigma_c)$'), y = "N individuals",
       colour = "Carrier Frequency") + ylim(0, 50000)

#png(filename="dichotomous_diff_sweep_nind_kd.png", width=600, height=400)
#ggsave(file="dichotomous_diff_sweep_nind_kd.png", width=6, height=4, dpi=300)

# Fig 5A
# kd 0.002 cost 0.1 Combined plots patchwork
results_kd_n  %>% filter(cost==0.1) %>% filter(kd==0.002) %>% filter(spite_expression_prob==0.001) %>% ggplot(
  aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kd_n  %>% filter(cost==0.1) %>% filter(kd==0.002) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kd_n  %>% filter(cost==0.1) %>% filter(kd==0.002) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'), 
       title = "Spite expression 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kd_n  %>% filter(cost==0.1) %>% filter(kd==0.002) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))
png(filename="dichotomous_diff_sweep_kd0.002_cost0.1_l.png", width=800, height=400)
ggsave(file="dichotomous_diff_sweep_kd0.002_cost0.1_l.png", width=8, height=4, dpi=300)

# Fig 5B
# Well mixed diff sigma = "i"
# kd 0.002 cost 0.1 Combined plots patchwork
results_kd_i  %>% filter(cost==0.1) %>% filter(kd==0.002) %>% filter(spite_expression_prob==0.001) %>% ggplot(
  aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title =element_text(size=14), legend.title = element_text(size=14), 
        legend.text = element_text(size = 12), legend.position = "none", 
        strip.text.x = element_text(size = 14)) +
  
  results_kd_i %>% filter(cost==0.1) %>% filter(kd==0.002) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title =element_text(size=14), legend.title = element_text(size=14), 
        legend.text = element_text(size = 12), legend.position = "none", 
        strip.text.x = element_text(size = 14)) +
  
  results_kd_i  %>% filter(cost==0.1) %>% filter(kd==0.002) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title =element_text(size=14), legend.title = element_text(size=14), 
        legend.text = element_text(size = 12), legend.position = "none", 
        strip.text.x = element_text(size = 14)) +
  
  results_kd_i  %>% filter(cost==0.1) %>% filter(kd==0.002) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title =element_text(size=14), legend.title = element_text(size=14), 
        legend.text = element_text(size = 12), 
        strip.text.x = element_text(size = 14))
png(filename="dichotomous_diff_i_kd0.002_cost0.1_l.png", width=450, height=400)
ggsave(file="dichotomous_diff_i_kd0.002_cost0.1_l.png", width=4.5, height=4, dpi=300)

#--------------------------------------------
# Fig 6

# import csv results higher density kc 0.25 and higher kd
results_kc <- read.csv("Results_csv/Results_dichotomous_l_kc_kd.csv")

# Separate numeric diff_sigma and well-mixed ("i -> infty")
results_kc_n <- results_kc %>% filter(diff_sigma != "i")
results_kc_i <- results_kc %>% filter(diff_sigma == "i")

# Convert sigma columns to numeric
results_kc_n$spite_sigma <- results_kc_n$spite_sigma %>% as.character() %>% as.numeric()
results_kc_n$comp_sigma <- results_kc_n$comp_sigma %>% as.character() %>% as.numeric()
results_kc_n$diff_sigma <- results_kc_n$diff_sigma %>% as.character() %>% as.numeric()

results_kc_i$spite_sigma <- results_kc_i$spite_sigma %>% as.character() %>% as.numeric()
results_kc_i$comp_sigma <- results_kc_i$comp_sigma %>% as.character() %>% as.numeric()

# Create columns for the sigma ratios in the numeric dataframe
results_kc_n <- results_kc_n %>% mutate(
  spite_diff_ratio = spite_sigma/diff_sigma, 
  comp_diff_ratio = comp_sigma/diff_sigma)

# Create columns for the sigma ratios in the numeric dataframe
results_kc_n$diff_sigma_label <- factor(
  results_kc_n$diff_sigma,
  levels = c("4", "8", "16"),
  labels = c(TeX("$\\sigma_{d} = 4$"), 
             TeX("$\\sigma_{d} = 8$"), TeX("$\\sigma_{d} = 16$")))

results_kc_i$diff_sigma_label <- factor(
  results_kc_i$diff_sigma,
  levels = c("i"),
  labels = c(TeX("$\\sigma_{d} = \\infty$")))

# Plot the number of individuals depending on sigma_c
# Dichotomous diffusion sweep n individuals kc
results_kc_n %>% filter(kd == 0.001) %>% ggplot(aes(x = log2(comp_sigma), y = average_n_individuals,
                                                    color = as.factor(spite_expression_prob))) + 
  geom_point() +
  scale_x_continuous(labels = c("0", "1", "2", "3", "4", "5", "6"), breaks = c(0,1,2,3,4,5,6)) + 
  labs(x = TeX('$log_2(\\sigma_c)$'), y = "Population Size N",
       colour = "Spite expression") + ylim(0, 100000) + 
  theme(legend.title = element_text(size=10), legend.text = element_text(size=10))
#png(filename="dichotomous_diff_sweep_nind_kc.png", width=350, height=200)
#ggsave(file="dichotomous_diff_sweep_nind_kc.png", width=3.5, height=2, dpi=300)

#-----------------------------------
# Fig 6A
# higher density kc 0.25 kd 0.001 cost 0.1
results_kc_n  %>% filter(initial_spite_frequency==0.5) %>% filter(cost == 0.1) %>% 
  filter(kd == 0.001) %>% filter(spite_expression_prob==0.001) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_n  %>% filter(initial_spite_frequency==0.5) %>% filter(cost == 0.1) %>%
  filter(kd == 0.001) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_n  %>% filter(initial_spite_frequency==0.5) %>% filter(cost == 0.1) %>%
  filter(kd == 0.001) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'), 
       title = "Spite expression 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none")+
  
  results_kc_n  %>% filter(initial_spite_frequency==0.5) %>% filter(cost == 0.1) %>%
  filter(kd == 0.001) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))
png(filename="dichotomous_diff_sweep_cost0.1_kc0.25_kd0.001.png", width=800, height=400)
ggsave(file="dichotomous_diff_sweep_cost0.1_kc0.25_kd0.001.png", width=8, height=4, dpi=300)

# Fig 6B
# higher density kc 0.25 kd 0.001 cost 0.1
# Well mixed diff sigma = "i"
results_kc_i  %>% filter(cost==0.1) %>% filter(kd==0.001) %>% filter(spite_expression_prob==0.001) %>% ggplot(
  aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_i  %>% filter(cost==0.1) %>% filter(kd==0.001) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_i  %>% filter(cost==0.1) %>% filter(kd==0.001) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_i  %>% filter(cost==0.1) %>% filter(kd==0.001) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))

png(filename="dichotomous_diff_i_cost0.1_kc0.25_kd0.001.png", width=450, height=400)
ggsave(file="dichotomous_diff_i_cost0.1_kc0.25_kd0.001.png", width=4.5, height=4, dpi=300)

#--------------------------------
# Fig 6C
# higher density kc 0.25 kd 0.005 cost 0.1
results_kc_n  %>% filter(initial_spite_frequency==0.5) %>% filter(cost == 0.1) %>% 
  filter(kd == 0.005) %>% filter(spite_expression_prob==0.001) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_n  %>% filter(initial_spite_frequency==0.5) %>% filter(cost == 0.1) %>%
  filter(kd == 0.005) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_n  %>% filter(initial_spite_frequency==0.5) %>% filter(cost == 0.1) %>%
  filter(kd == 0.005) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'), 
       title = "Spite expression 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_n  %>% filter(initial_spite_frequency==0.5) %>% filter(cost == 0.1) %>%
  filter(kd == 0.005) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))
png(filename="dichotomous_diff_sweep_cost0.1_kc0.25_kd0.005.png", width=800, height=400)
ggsave(file="dichotomous_diff_sweep_cost0.1_kc0.25_kd0.005.png", width=8, height=4, dpi=300)

# Fig 6D
# higher density kc 0.25 kd 0.005 cost 0.1
# Well mixed diff sigma = "i"

results_kc_i  %>% filter(cost==0.1) %>% filter(kd==0.005) %>% filter(spite_expression_prob==0.001) %>% ggplot(
  aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.001", fill = "Carrier Frequency")  + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_i  %>% filter(cost==0.1) %>% filter(kd==0.005) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.01", fill = "Carrier Frequency")   + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_i  %>% filter(cost==0.1) %>% filter(kd==0.005) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp.  0.05", fill = "Carrier Frequency")  + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_i  %>% filter(cost==0.1) %>% filter(kd==0.005) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.2", fill = "Carrier Frequency")  + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))

png(filename="dichotomous_diff_i_cost0.1_kc0.25_kd0.005.png", width=450, height=400)
ggsave(file="dichotomous_diff_i_cost0.1_kc0.25_kd0.005.png", width=4.5, height=4, dpi=300)

#----------------------------------------
# Fig 7

# import csv results higher density kc 0.25   kd 0.001  lower f_0 freq 0.02
results_kc_freq <- read.csv("Results_csv/Results_dichotomous_l_kc_freq.csv")

# Separate numeric diff_sigma and well-mixed ("i -> infty")
results_kc_freq_n <- results_kc_freq %>% filter(diff_sigma != "i")
results_kc_freq_i <- results_kc_freq %>% filter(diff_sigma == "i")

# Convert sigma columns to numeric
results_kc_freq_n$spite_sigma <- results_kc_freq_n$spite_sigma %>% as.character() %>% as.numeric()
results_kc_freq_n$comp_sigma <- results_kc_freq_n$comp_sigma %>% as.character() %>% as.numeric()
results_kc_freq_n$diff_sigma <- results_kc_freq_n$diff_sigma %>% as.character() %>% as.numeric()

results_kc_freq_i$spite_sigma <- results_kc_freq_i$spite_sigma %>% as.character() %>% as.numeric()
results_kc_freq_i$comp_sigma <- results_kc_freq_i$comp_sigma %>% as.character() %>% as.numeric()

# Create columns for the sigma ratios in the numeric dataframe
results_kc_freq_n <- results_kc_freq_n %>% mutate(
  spite_diff_ratio = spite_sigma/diff_sigma, 
  comp_diff_ratio = comp_sigma/diff_sigma)

# Add a column for the diff_sigma Latex labels
results_kc_freq_n$diff_sigma_label <- factor(
  results_kc_freq_n$diff_sigma,
  levels = c("4", "8", "16"),
  labels = c(TeX("$\\sigma_{d} = 4$"), 
             TeX("$\\sigma_{d} = 8$"), TeX("$\\sigma_{d} = 16$")))

results_kc_freq_i$diff_sigma_label <- factor(
  results_kc_freq_i$diff_sigma,
  levels = c("i"),
  labels = c(TeX("$\\sigma_{d} = \\infty$")))

# Plot the number of individuals depending on sigma_c
# Dichotomous diffusion sweep n individuals kc  freq 0.02
results_kc_freq_n %>% ggplot(aes(x = log2(comp_sigma), y = average_n_individuals, color = average_spite_freq)) + 
  geom_point(aes(shape = as.factor(spite_expression_prob)), size = 3) +
  scale_x_continuous(labels = c("0", "1", "2", "3", "4", "5", "6"), breaks = c(0,1,2,3,4,5,6)) + 
  labs(x = TeX('$log_2(\\sigma_c)$'), y = "N individuals",
       colour = "Carrier Frequency", shape = "Spite expression") + ylim(0, 100000)
#png(filename="dichotomous_diff_sweep_nind_kc_freq.png", width=600, height=400)
#ggsave(file="dichotomous_diff_sweep_nind_kc_freq.png", width=6, height=4, dpi=300)

#--------------------------------
# Fig 7A
# higher density kc 0.25 freq 0.02 kd 0.001 cost 0.1
results_kc_freq_n  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>% 
  filter(kd == 0.001) %>% filter(spite_expression_prob==0.001) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_freq_n  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>%
  filter(kd == 0.001) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_freq_n  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>%
  filter(kd == 0.001) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'), 
       title = "Spite expression 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_freq_n  %>% filter(initial_spite_frequency==0.02) %>% filter(cost == 0.1) %>%
  filter(kd == 0.001) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_diff_ratio), y= log2(comp_diff_ratio), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  labs(x = TeX('$log_2(\\sigma_s/\\sigma_d)$'), y =  TeX('$log_2(\\sigma_c/\\sigma_d)$'),
       title = "Spite expression 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))
png(filename="dichotomous_diff_sweep_cost0.1_kc0.25_kd0.001_f0_0.02.png", width=800, height=400)
ggsave(file="dichotomous_diff_sweep_cost0.1_kc0.25_kd0.001_f0_0.02.png", width=8, height=4, dpi=300)

# Fig 7B
# higher density kc 0.25
# kd 0.001 cost 0.1 Well mixed diff sigma = "i"
results_kc_freq_i  %>% filter(cost==0.1) %>% filter(kd==0.001) %>% filter(spite_expression_prob==0.001) %>% ggplot(
  aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.001", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_freq_i  %>% filter(cost==0.1) %>% filter(kd==0.001) %>% filter(spite_expression_prob==0.01) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.01", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_freq_i  %>% filter(cost==0.1) %>% filter(kd==0.001) %>% filter(spite_expression_prob==0.05) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.05", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14), legend.position = "none") +
  
  results_kc_freq_i  %>% filter(cost==0.1) %>% filter(kd==0.001) %>% filter(spite_expression_prob==0.2) %>% ggplot(
    aes(x=log2(spite_sigma), y= log2(comp_sigma), fill = average_spite_freq)) + 
  geom_tile() + facet_wrap(~diff_sigma_label, labeller = label_parsed) + scale_fill_distiller(palette = 'RdBu', limit = c(0,1)) +
  scale_x_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) + 
  scale_y_continuous(labels = c("1", "2", "3", "4", "5", "6", "7"), breaks = c(1,2,3,4,5,6,7)) +
  labs(x = TeX('$log_2(\\sigma_s)$'), y =  TeX('$log_2(\\sigma_c)$'),
       title = "Exp. 0.2", fill = "Carrier Frequency") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title =element_text(size=15), legend.title = element_text(size=14), 
        legend.text = element_text(size = 15), strip.text.x = element_text(
          size = 14))

png(filename="dichotomous_diff_i_cost0.1_kc0.25_kd0.001_f0_0.02.png", width=450, height=400)
ggsave(file="dichotomous_diff_i_cost0.1_kc0.25_kd0.001_f0_0.02.png", width=4.5, height=4, dpi=300)
