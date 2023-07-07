rm(list=ls())
set.seed(821)

# environment ====
library(functions)
library(ggpubr)
library(EpiViz)
library(ggforestplot)
library(ggforestplot)
library(wesanderson)
library(tidyverse)
library(rlang)
library(plyr)
library(data.table)
library(cowplot)
source("scripts/my_forestplot.R")
source("scripts/colour_palette.R")

# data ====
data_main <- read.table("analysis/002_analysis/proteomics/results_formatted.txt", header = T, sep = "\t")
data_sensitivity <- read.table("analysis/002_analysis/proteomics/results_2y-followup_formatted.txt", header = T, sep = "\t")
data <- bind_rows(data_main,data_sensitivity)
data$ID <- paste0(data$exposure, "_", data$cancer, "_", data$sex, "_", data$model, "_", data$followup)
data$followup <- as.character(data$followup)
data$followup[is.na(data$followup)] <- "all"
data$followup <- as.factor(data$followup)
psignif <- 0.05/61
ci <- 0.95

# directional consistency across analyses ====
data <- select(data, ID, exposure, cancer, sex, followup, model, analysis, b)
data <- subset(data, model == "model 5")

## make data frmaes ====
### sex-combined ====
#### overall
data1 <- subset(data, sex == "sex-combined")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, followup == "above")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_combined_overall_above <- directional_consistency(data1)
data_combined_overall_above$exposure <- "combined_overall_above"

data1 <- subset(data, sex == "sex-combined")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, followup == "below")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_combined_overall_below <- directional_consistency(data1)
data_combined_overall_below$exposure <- "combined_overall_below"

data1 <- subset(data, sex == "sex-combined")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, followup == "all")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_combined_overall_all <- directional_consistency(data1)
data_combined_overall_all$exposure <- "combined_overall_all"

#### distal
data1 <- subset(data, sex == "sex-combined")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, followup == "above")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_combined_distal_above <- directional_consistency(data1)
data_combined_distal_above$exposure <- "combined_distal_above"

data1 <- subset(data, sex == "sex-combined")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, followup == "below")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_combined_distal_below <- directional_consistency(data1)
data_combined_distal_below$exposure <- "combined_distal_below"

data1 <- subset(data, sex == "sex-combined")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, followup == "all")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_combined_distal_all <- directional_consistency(data1)
data_combined_distal_all$exposure <- "combined_distal_all"

#### proximal
data1 <- subset(data, sex == "sex-combined")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, followup == "above")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_combined_proximal_above <- directional_consistency(data1)
data_combined_proximal_above$exposure <- "combined_proximal_above"

data1 <- subset(data, sex == "sex-combined")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, followup == "below")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_combined_proximal_below <- directional_consistency(data1)
data_combined_proximal_below$exposure <- "combined_proximal_below"

data1 <- subset(data, sex == "sex-combined")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, followup == "all")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_combined_proximal_all <- directional_consistency(data1)
data_combined_proximal_all$exposure <- "combined_proximal_all"

### female ====
#### overall
data1 <- subset(data, sex == "female")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, followup == "above")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_female_overall_above <- directional_consistency(data1)
data_female_overall_above$exposure <- "female_overall_above"

data1 <- subset(data, sex == "female")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, followup == "below")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_female_overall_below <- directional_consistency(data1)
data_female_overall_below$exposure <- "female_overall_below"

data1 <- subset(data, sex == "female")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, followup == "all")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_female_overall_all <- directional_consistency(data1)
data_female_overall_all$exposure <- "female_overall_all"

#### distal
data1 <- subset(data, sex == "female")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, followup == "above")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_female_distal_above <- directional_consistency(data1)
data_female_distal_above$exposure <- "female_distal_above"

data1 <- subset(data, sex == "female")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, followup == "below")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_female_distal_below <- directional_consistency(data1)
data_female_distal_below$exposure <- "female_distal_below"

data1 <- subset(data, sex == "female")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, followup == "all")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_female_distal_all <- directional_consistency(data1)
data_female_distal_all$exposure <- "female_distal_all"

#### proximal
data1 <- subset(data, sex == "female")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, followup == "above")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_female_proximal_above <- directional_consistency(data1)
data_female_proximal_above$exposure <- "female_proximal_above"

data1 <- subset(data, sex == "female")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, followup == "below")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_female_proximal_below <- directional_consistency(data1)
data_female_proximal_below$exposure <- "female_proximal_below"

data1 <- subset(data, sex == "female")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, followup == "all")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_female_proximal_all <- directional_consistency(data1)
data_female_proximal_all$exposure <- "female_proximal_all"


### male ====
#### overall
data1 <- subset(data, sex == "male")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, followup == "above")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_male_overall_above <- directional_consistency(data1)
data_male_overall_above$exposure <- "male_overall_above"

data1 <- subset(data, sex == "male")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, followup == "below")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_male_overall_below <- directional_consistency(data1)
data_male_overall_below$exposure <- "male_overall_below"

data1 <- subset(data, sex == "male")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, followup == "all")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_male_overall_all <- directional_consistency(data1)
data_male_overall_all$exposure <- "male_overall_all"

#### distal
data1 <- subset(data, sex == "male")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, followup == "above")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_male_distal_above <- directional_consistency(data1)
data_male_distal_above$exposure <- "male_distal_above"

data1 <- subset(data, sex == "male")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, followup == "below")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_male_distal_below <- directional_consistency(data1)
data_male_distal_below$exposure <- "male_distal_below"

data1 <- subset(data, sex == "male")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, followup == "all")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_male_distal_all <- directional_consistency(data1)
data_male_distal_all$exposure <- "male_distal_all"

#### proximal
data1 <- subset(data, sex == "male")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, followup == "above")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_male_proximal_above <- directional_consistency(data1)
data_male_proximal_above$exposure <- "male_proximal_above"

data1 <- subset(data, sex == "male")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, followup == "below")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_male_proximal_below <- directional_consistency(data1)
data_male_proximal_below$exposure <- "male_proximal_below"

data1 <- subset(data, sex == "male")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, followup == "all")
data1 <- as.data.frame(pivot_wider(data1, names_from = analysis, values_from = b))
data1 <- data1[,c(1,7,8)]
data_male_proximal_all <- directional_consistency(data1)
data_male_proximal_all$exposure <- "male_proximal_all"

## plot ====
data_direction <- rbind(data_combined_overall_above, data_combined_overall_below, data_combined_overall_all,
                        data_combined_distal_above, data_combined_distal_below, data_combined_distal_all,
                        data_combined_proximal_above, data_combined_proximal_below, data_combined_proximal_all,
                        data_female_overall_above, data_female_overall_below, data_female_overall_all,
                        data_female_distal_above, data_female_distal_below, data_female_distal_all,
                        data_female_proximal_above, data_female_proximal_below, data_female_proximal_all,
                        data_male_overall_above, data_male_overall_below, data_male_overall_all,
                        data_male_distal_above, data_male_distal_below, data_male_distal_all,
                        data_male_proximal_above, data_male_proximal_below, data_male_proximal_all)

data_direction$exposure <- factor(data_direction$exposure, levels = c("combined_overall_above", "combined_overall_below", "combined_overall_all",
                                                                      "combined_distal_above", "combined_distal_below", "combined_distal_all",
                                                                      "combined_proximal_above", "combined_proximal_below", "combined_proximal_all",
                                                                      "female_overall_above", "female_overall_below", "female_overall_all",
                                                                      "female_distal_above", "female_distal_below", "female_distal_all",
                                                                      "female_proximal_above", "female_proximal_below", "female_proximal_all",
                                                                      "male_overall_above", "male_overall_below", "male_overall_all",
                                                                      "male_distal_above", "male_distal_below", "male_distal_all",
                                                                      "male_proximal_above", "male_proximal_below", "male_proximal_all"))

pdf("analysis/002_analysis/proteomics/figures/006_main-sensitivity_comparison_m5.pdf",
    width = 8, height = 6, pointsize = 10)
ggplot(data = data_direction,
       aes(x = exposure)) +
  geom_bar(aes(fill = direction_group))  +
  scale_fill_manual(values = discrete_wes_pal) +
  guides(fill = guide_legend(override.aes = list(size = 5),
                             title = "",
                             label.hjust = 0,
                             label.vjust = 0.5)) +
  xlab("") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")
dev.off()
table(data_direction$direction_group, data_direction$exposure)

# main results in manuscript comparison with sensitivity ====
data <- read.table("analysis/002_analysis/proteomics/results_formatted.txt", header = T, sep = "\t")
data$cancer <- factor(data$cancer, levels = c("overall", "distal", "proximal"))
data$sex <- factor(data$sex, levels = c("male", "female", "sex-combined"))
data$model <- factor(data$model, levels = c("model 5", "model 4", "model 3", "model 2", "model 1"))
data$followup <- factor(data$followup, levels = c("below", "above"))
data$ID <- paste0(data$exposure, "_", data$cancer, "_", data$sex, "_", data$model)
plot_data <- subset(data, model == "model 5")
plot_data <- plot_data[is.na(plot_data$followup),]
plot_data_above <- subset(plot_data, lower_ci > 1)
plot_data_below <- subset(plot_data, upper_ci < 1)
plot_data <- rbind(plot_data_below, plot_data_above)
plot_data <- plot_data[order(plot_data$exposure),]
ID_main <- plot_data$ID

data_main <- read.table("analysis/002_analysis/proteomics/results_formatted.txt", header = T, sep = "\t")
data_sensitivity <- read.table("analysis/002_analysis/proteomics/results_2y-followup_formatted.txt", header = T, sep = "\t")
data <- bind_rows(data_main,data_sensitivity)
data$ID <- paste0(data$exposure, "_", data$cancer, "_", data$sex, "_", data$model)
data$followup <- as.character(data$followup)
data$followup[is.na(data$followup)] <- "all"
data$followup <- as.factor(data$followup)

plot_data <- data[data$ID %in% ID_main, ]
plot_data <- subset(plot_data, followup == "all")

psignif <- 0.05/60
ci <- 0.95

pdf("analysis/002_analysis/proteomics/figures/006_main-sensitivity_comparison_main-highlighted_m5.pdf",
    width = 8, height = 12, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = sex,
              shape = analysis,
              logodds = T, 
              space = 0.9) +
  ggforce::facet_col(
    facets = ~cancer,
    scales = "free_y",
    space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank()) 
dev.off()

# median followup results in manuscript comparison with sensitivity ====
data <- read.table("analysis/002_analysis/proteomics/results_formatted.txt", header = T, sep = "\t")
data$cancer <- factor(data$cancer, levels = c("overall", "distal", "proximal"))
data$sex <- factor(data$sex, levels = c("male", "female", "sex-combined"))
data$model <- factor(data$model, levels = c("model 5", "model 4", "model 3", "model 2", "model 1"))
data$followup <- factor(data$followup, levels = c("below", "above"))
data$ID <- paste0(data$exposure, "_", data$cancer, "_", data$sex, "_", data$model)
plot_data <- subset(data, model == "model 5")
plot_data <- plot_data[!is.na(plot_data$followup),]
plot_data_above <- subset(plot_data, lower_ci > 1)
plot_data_below <- subset(plot_data, upper_ci < 1)
plot_data <- rbind(plot_data_below, plot_data_above)
ID_main <- plot_data$ID

data_main <- read.table("analysis/002_analysis/proteomics/results_formatted.txt", header = T, sep = "\t")
data_sensitivity <- read.table("analysis/002_analysis/proteomics/results_2y-followup_formatted.txt", header = T, sep = "\t")
data <- bind_rows(data_main,data_sensitivity)
data$ID <- paste0(data$exposure, "_", data$cancer, "_", data$sex, "_", data$model)
data$followup <- as.character(data$followup)
data$followup[is.na(data$followup)] <- "all"
data$followup <- as.factor(data$followup)

plot_data <- data[data$ID %in% ID_main, ]
plot_data <- subset(plot_data, followup != "all")

psignif <- 0.05/60
ci <- 0.95

pdf("analysis/002_analysis/proteomics/figures/006_main-sensitivity_comparison_followup-highlighted_m5.pdf",
    width = 10, height = 16, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = followup,
              shape = analysis,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             rows = vars(sex),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

# interesting proteins ====
data_main <- read.table("analysis/002_analysis/proteomics/results_formatted.txt", header = T, sep = "\t")
data_sensitivity <- read.table("analysis/002_analysis/proteomics/results_2y-followup_formatted.txt", header = T, sep = "\t")
data <- bind_rows(data_main,data_sensitivity)
data$ID <- paste0(data$exposure, "_", data$cancer, "_", data$sex, "_", data$model)
data$followup <- as.character(data$followup)
data$followup[is.na(data$followup)] <- "all"
data$followup <- as.factor(data$followup)

plot_data <- subset(data, model == "model 5")
proteins <- c("Olk_Il7", "Olk_Pdcd1", "Olk_Tnfrsf4", "Olk_Tnfrsf9")
plot_data <- plot_data[plot_data$exposure %in% proteins, ]
plot_data$followup <- factor(plot_data$followup, levels = c("below", "all", "above"))
plot_data <- plot_data %>%
  mutate(ci_solid = ifelse(lower_ci > 1 | upper_ci < 1, 0, 1))

pdf("analysis/002_analysis/proteomics/figures/006_main-sensitivity_comparison_followup-highlighted_interesting_m5.pdf",
    width = 12, height = 6, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = ci_solid,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = followup,
              shape = analysis,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             rows = vars(sex),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(panel.spacing = unit(2, "lines"))
dev.off()
