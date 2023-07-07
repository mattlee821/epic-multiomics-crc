rm(list=ls())
set.seed(821)

# environment ====
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
data <- read.table("analysis/002_analysis/metabolomics/positive/results_formatted.txt", header = T, sep = "\t")
data$cancer <- factor(data$cancer, levels = c("overall", "distal", "proximal"))
data$sex <- factor(data$sex, levels = c("male", "female", "sex-combined"))
data$model <- factor(data$model, levels = c("model 5", "model 4", "model 3", "model 2", "model 1"))
data$followup <- factor(data$followup, levels = c("below", "above"))
data$ID <- paste0(data$exposure, "_", data$cancer, "_", data$sex, "_", data$model, "_", data$followup)

# multiple testing ====
data_threshold <- subset(data, model == "model 5")
data_threshold <- subset(data_threshold, p < 0.05)

psignif <- 0.05/233
ci <- 0.95

# model comparison ====
## sex-combined ====
plot_data <- subset(data, sex == "sex-combined")
plot_data <- plot_data[is.na(plot_data$followup),]

pdf("analysis/002_analysis/metabolomics/positive/figures/001_model_comparison_sex-combined.pdf",
    width = 16, height = 130, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = model,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
                     scales = "free_y",
                     space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()

## female ====
plot_data <- subset(data, sex == "female")
plot_data <- plot_data[is.na(plot_data$followup),]

pdf("analysis/002_analysis/metabolomics/positive/figures/001_model_comparison_female.pdf",
    width = 16, height = 130, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = model,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()

## male ====
plot_data <- subset(data, sex == "male")
plot_data <- plot_data[is.na(plot_data$followup),]

pdf("analysis/002_analysis/metabolomics/positive/figures/001_model_comparison_male.pdf",
    width = 16, height = 130, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = model,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()


# combined, male, female - all cancers ====
plot_data <- subset(data, model == "model 5")
plot_data <- plot_data[is.na(plot_data$followup),]

pdf("analysis/002_analysis/metabolomics/positive/figures/002_sex_site-specific_m5.pdf",
    width = 16, height = 22, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = sex,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()

# combined, male, female - all cancers - highlighted ====
plot_data <- subset(data, model == "model 5")
plot_data <- plot_data[is.na(plot_data$followup),]
plot_data_above <- subset(plot_data, lower_ci > 1)
plot_data_below <- subset(plot_data, upper_ci < 1)
plot_data_null <- rbind(plot_data_below, plot_data_above)
plot_data_null$shape <- "non-null" 
ID <- plot_data_null$ID
plot_data <- plot_data[!plot_data$ID %in% ID, ]
plot_data$shape <- "null"
plot_data <- rbind(plot_data, plot_data_null)

pdf("analysis/002_analysis/metabolomics/positive/figures/002_sex_site-specific_highlight_m5.pdf",
    width = 16, height = 22, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = sex,
              shape = shape,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()

# combined, male, female - all cancers - highlighted only ====
plot_data <- subset(data, model == "model 5")
plot_data <- plot_data[is.na(plot_data$followup),]
plot_data_above <- subset(plot_data, lower_ci > 1)
plot_data_below <- subset(plot_data, upper_ci < 1)
plot_data <- rbind(plot_data_below, plot_data_above)
plot_data <- plot_data[order(plot_data$exposure),]

pdf("analysis/002_analysis/metabolomics/positive/figures/002_sex_site-specific_highlight-only_m5.pdf",
    width = 12, height = 8, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = sex,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank()) 
dev.off()

pdf("analysis/002_analysis/metabolomics/positive/figures/002_sex_site-specific_highlight-only_m5_2.pdf",
    width = 6, height = 8, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = sex,
              logodds = T, 
              space = 0.9) +
  ggforce::facet_col(
    facets = ~cancer,
    scales = "free_y",
    space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank()) 
dev.off()

# directional consitency ====
## combined - all ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "sex-combined")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(4,5)]
data1$exposure <- "combined_all"
data_direction_combined_all <- data1

## combined - distal proximal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "sex-combined")
data1 <- subset(data1, cancer != "overall")

data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "combined_distal-proximal"
data_direction_combined_dp <- data1

## combined - overall distal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "sex-combined")
data1 <- subset(data1, cancer != "proximal")

data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "combined_overall-distal"
data_direction_combined_od <- data1

## combined - overall proximal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "sex-combined")
data1 <- subset(data1, cancer != "distal")

data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "combined_overall-proximal"
data_direction_combined_op <- data1

## male - all ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "male")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(4,5)]
data1$exposure <- "male_all"
data_direction_male_all <- data1

## male - distal proximal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "male")
data1 <- subset(data1, cancer != "overall")

data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "male_distal-proximal"
data_direction_male_dp <- data1

## male - overall distal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "male")
data1 <- subset(data1, cancer != "proximal")

data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "male_overall-distal"
data_direction_male_od <- data1

## male - overall proximal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "male")
data1 <- subset(data1, cancer != "distal")

data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "male_overall-proximal"
data_direction_male_op <- data1

## female - all ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "female")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(4,5)]
data1$exposure <- "female_all"
data_direction_female_all <- data1

## female - distal proximal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "female")
data1 <- subset(data1, cancer != "overall")

data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "female_distal-proximal"
data_direction_female_dp <- data1

## male - overall distal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "female")
data1 <- subset(data1, cancer != "proximal")

data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "female_overall-distal"
data_direction_female_od <- data1

## male - overall proximal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, sex == "female")
data1 <- subset(data1, cancer != "distal")

data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, cancer, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = cancer, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "female_overall-proximal"
data_direction_female_op <- data1

## join and plot ====
data_direction <- rbind(data_direction_combined_all, data_direction_combined_dp,
                        data_direction_combined_od, data_direction_combined_op,
                        data_direction_male_all, data_direction_male_dp,
                        data_direction_male_od, data_direction_male_op,
                        data_direction_female_all, data_direction_female_dp,
                        data_direction_female_od, data_direction_female_op)

data_direction$exposure <- factor(data_direction$exposure, levels = c("combined_all", "combined_distal-proximal",
                                                                      "combined_overall-distal", "combined_overall-proximal",
                                                                      "male_all", "male_distal-proximal",
                                                                      "male_overall-distal", "male_overall-proximal",
                                                                      "female_all", "female_distal-proximal",
                                                                      "female_overall-distal", "female_overall-proximal"))
                                  
pdf("analysis/002_analysis/metabolomics/positive/figures/003_site-specific_sex_m5_directional-consistency.pdf",
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

# directional consisteny - across sites within sex ====
## directional consistency combined, male, female - overall ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "overall")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(4,5)]
data1$exposure <- "overall_all"
data_direction_overall_all <- data1

## directional consistency male, female - overall ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, sex != "sex-combined")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "overall_male-female"
data_direction_overall_mf <- data1

## directional consistency combined, male - overall ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, sex != "female")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "overall_all-male"
data_direction_overall_am <- data1

## directional consistency combined, female - overall ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "overall")
data1 <- subset(data1, sex != "male")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "overall_all-female"
data_direction_overall_af <- data1

## directional consistency combined, male, female - distal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "distal")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(4,5)]
data1$exposure <- "distal_all"
data_direction_distal_all <- data1

## directional consistency male, female - distal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, sex != "sex-combined")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "distal_male-female"
data_direction_distal_mf <- data1

## directional consistency combined, male - overall ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, sex != "female")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "distal_all-male"
data_direction_distal_am <- data1

## directional consistency combined, female - overall ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, sex != "male")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "distal_all-female"
data_direction_distal_af <- data1

## directional consistency combined, male, female - proximal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "distal")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(4,5)]
data1$exposure <- "proximal_all"
data_direction_proximal_all <- data1

## directional consistency male, female - proximal ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "distal")
data1 <- subset(data1, sex != "sex-combined")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "proximal_male-female"
data_direction_proximal_mf <- data1

## directional consistency combined, male - overall ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, sex != "female")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "proximal_all-male"
data_direction_proximal_am <- data1

## directional consistency combined, female - overall ====
data1 <- subset(data, model == "model 5")
data1 <- subset(data1, cancer == "proximal")
data1 <- subset(data1, sex != "male")
data1 <- data1[is.na(data1$followup),]
data1 <- select(data1, exposure, sex, b)
data1 <- as.data.frame(pivot_wider(data1, names_from = sex, values_from = b))

rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                                     direction == 2 ~ "negative",
                                     direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- data1[,c(3,4)]
data1$exposure <- "proximal_all-female"
data_direction_proximal_af <- data1

## join and plot ====
data_direction <- rbind(data_direction_overall_all, data_direction_overall_mf,
                        data_direction_overall_am, data_direction_overall_af,
                        data_direction_distal_all, data_direction_distal_mf,
                        data_direction_distal_am, data_direction_distal_af,
                        data_direction_proximal_all, data_direction_proximal_mf,
                        data_direction_proximal_am, data_direction_proximal_af)

data_direction$exposure <- factor(data_direction$exposure, levels = c("overall_all", "overall_male-female",
                                                                      "overall_all-male", "overall_all-female",
                                                                      "distal_all", "distal_male-female",
                                                                      "distal_all-male", "distal_all-female",
                                                                      "proximal_all", "proximal_male-female",
                                                                      "proximal_all-male", "proximal_all-female"))

pdf("analysis/002_analysis/metabolomics/positive/figures/003_sex_site-specific_m5_directional-consistency.pdf",
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










# median followup time ====
## data ====
data_analysis <- fread("data/processed_data.txt")
data_analysis$Length_Bld <- data_analysis$Length_Bld/365
data_analysis_male <- subset(data_analysis, Sex == 1)
data_analysis_female <- subset(data_analysis, Sex == 2)
data_analysis_combined <- data_analysis

### site specific
proximal <- c("C180", "C181", "C182", "C183", "C184", "C185")
distal <- c("C186", "C187")
overlapping <- "C188"
unspecified <- "C189"
data_analysis_proximal <- data_analysis[data_analysis$Siteclrt %in% proximal, ]
data_analysis_proximal <- data_analysis_proximal$Match_Caseset
data_analysis_proximal_combined <- data_analysis[data_analysis$Match_Caseset %in% data_analysis_proximal, ]
data_analysis_proximal_male <- subset(data_analysis_proximal_combined, Sex == 1)
data_analysis_proximal_female <- subset(data_analysis_proximal_combined, Sex == 2)

data_analysis_distal <- data_analysis[data_analysis$Siteclrt %in% distal, ]
data_analysis_distal <- data_analysis_distal$Match_Caseset
data_analysis_distal_combined <- data_analysis[data_analysis$Match_Caseset %in% data_analysis_distal, ]
data_analysis_distal_male <- subset(data_analysis_distal_combined, Sex == 1)
data_analysis_distal_female <- subset(data_analysis_distal_combined, Sex == 2)

data_analysis_overlapping <- data_analysis[data_analysis$Siteclrt %in% overlapping, ]
data_analysis_overlapping <- data_analysis_overlapping$Match_Caseset
data_analysis_overlapping_combined <- data_analysis[data_analysis$Match_Caseset %in% data_analysis_overlapping, ]
data_analysis_overlapping_male <- subset(data_analysis_overlapping_combined, Sex == 1)
data_analysis_overlapping_female <- subset(data_analysis_overlapping_combined, Sex == 2)

data_analysis_unspecified <- data_analysis[data_analysis$Siteclrt %in% unspecified, ]
data_analysis_unspecified <- data_analysis_unspecified$Match_Caseset
data_analysis_unspecified_combined <- data_analysis[data_analysis$Match_Caseset %in% data_analysis_unspecified, ]
data_analysis_unspecified_male <- subset(data_analysis_unspecified_combined, Sex == 1)
data_analysis_unspecified_female <- subset(data_analysis_unspecified_combined, Sex == 2)

### follow-up time data 
data_analysis_combined <- median(data_analysis_combined$Length_Bld)
data_analysis_proximal_combined <- median(data_analysis_proximal_combined$Length_Bld)
data_analysis_distal_combined <- median(data_analysis_distal_combined$Length_Bld)

data_analysis_male <- median(data_analysis_male$Length_Bld)
data_analysis_proximal_male <- median(data_analysis_proximal_male$Length_Bld)
data_analysis_distal_male <- median(data_analysis_distal_male$Length_Bld)

data_analysis_female <- median(data_analysis_female$Length_Bld)
data_analysis_proximal_female <- median(data_analysis_proximal_female$Length_Bld)
data_analysis_distal_female <- median(data_analysis_distal_female$Length_Bld)

follow_up <- data.frame(
  sex = c("overall", "distal", "proximal"),
  combined = c(round(data_analysis_combined, 2), round(data_analysis_distal_combined, 2), round(data_analysis_proximal_combined, 2)),
  female = c(round(data_analysis_female, 2), round(data_analysis_distal_female, 2), round(data_analysis_proximal_female, 2)),
  male = c(round(data_analysis_male, 2), round(data_analysis_distal_male, 2), round(data_analysis_proximal_male, 2))
)
write.table(follow_up, "analysis/tables/median-followup.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


### combined, male, female - all cancers - highlighted ====
plot_data <- subset(data, model == "model 5")
plot_data <- plot_data[!is.na(plot_data$followup),]
table(plot_data$followup)
plot_data_above <- subset(plot_data, lower_ci > 1)
plot_data_below <- subset(plot_data, upper_ci < 1)
plot_data_null <- rbind(plot_data_below, plot_data_above)
plot_data_null$shape <- "non-null" 
ID <- plot_data_null$ID
plot_data <- plot_data[!plot_data$ID %in% ID, ]
plot_data$shape <- "null"
plot_data <- rbind(plot_data, plot_data_null)

#### combined ====
plot_data_combined <- subset(plot_data, sex == "sex-combined")
pdf("analysis/002_analysis/metabolomics/positive/figures/004_combined_site-specific_followup_m5_highlight.pdf",
    width = 16, height = 22, pointsize = 10)
my_forestplot(df = plot_data_combined,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = followup,
              shape = shape,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()

#### female ====
plot_data_female <- subset(plot_data, sex == "female")
pdf("analysis/002_analysis/metabolomics/positive/figures/004_female_site-specific_followup_m5_highlight.pdf",
    width = 16, height = 22, pointsize = 10)
my_forestplot(df = plot_data_female,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = followup,
              shape = shape,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()

#### male ====
plot_data_male <- subset(plot_data, sex == "male")
pdf("analysis/002_analysis/metabolomics/positive/figures/004_male_site-specific_followup_m5_highlight.pdf",
    width = 16, height = 22, pointsize = 10)
my_forestplot(df = plot_data_male,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = followup,
              shape = shape,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()


### combined, male, female - all cancers - highlighted only ====
plot_data <- subset(data, model == "model 5")
plot_data <- plot_data[!is.na(plot_data$followup),]
plot_data_above <- subset(plot_data, lower_ci > 1)
plot_data_below <- subset(plot_data, upper_ci < 1)
plot_data <- rbind(plot_data_below, plot_data_above)

pdf("analysis/002_analysis/metabolomics/positive/figures/004_sex_site-specific_followup_highlight_m5.pdf",
    width = 12, height = 10, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = sex,
              shape = followup,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()

pdf("analysis/002_analysis/metabolomics/positive/figures/004_sex_site-specific_followup_highlight_m5_2.pdf",
    width = 8, height = 12, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = p,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = sex,
              shape = followup,
              logodds = T, 
              space = 0.9) +
  ggforce::facet_col(
    facets = ~cancer,
    scales = "free_y",
    space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank()) 
dev.off()

# interesting proteins ====
plot_data <- subset(data, model == "model 5")
proteins <- c("Olk_Il7", "Olk_Pdcd1", "Olk_Tnfrsf4", "Olk_Tnfrsf9")
plot_data <- plot_data[plot_data$exposure %in% proteins, ]
plot_data$followup <- as.character(plot_data$followup)
plot_data$followup[is.na(plot_data$followup)] <- "all"
plot_data$followup <- as.factor(plot_data$followup)
table(plot_data$followup)
plot_data$followup <- factor(plot_data$followup, levels = c("below", "all", "above"))
plot_data <- plot_data %>%
  mutate(ci_solid = ifelse(lower_ci > 1 | upper_ci < 1, 0, 1))

pdf("analysis/002_analysis/metabolomics/positive/figures/005_sex_site-specific_followup_highlight_interesting_m5.pdf",
    width = 12, height = 6, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = ci_solid,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = sex,
              shape = followup,
              logodds = T, 
              space = 0.9) +
  facet_grid(cols = vars(cancer),
             scales = "free_y",
             space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(panel.spacing = unit(2, "lines"))
dev.off()
