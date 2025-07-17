library(tidyverse)
library(ggplot2)
path <- "/home/callum/rotation2/plate_reader/CFU_data/"
image_path <- "/home/callum/rotation2/plate_reader/images/"

#file_name <- "2025-05-16-Silwet02"


# First the approach for the CFU survival assay
### Needs updating if more strains are added ###

file_name <- "2025-06-19-chelators01"
data <- read.csv(paste(path,file_name,".csv",sep=""))
data$experiment <- "M9_01"
data$pretreatment <- "M9"
# drop notes column
data <- subset(data, select = -X)

file_name <- "2025-06-26-chelators02"
data2 <- read.csv(paste(path,file_name,".csv",sep=""))
data2$experiment <- "M9_02"
data2$pretreatment <- "M9"
# drop notes columns
data2 <- subset(data2, select = -c(X,X.1))

file_name <- "2025-06-27-chelators03"
data3 <- read.csv(paste(path,file_name,".csv",sep=""))
data3$experiment <- "M9ctr_01"
data3$pretreatment <- "M9-citr"


all_data <- rbind(data, data2, data3)
# Generate count column (used 3 uL and want per mL so *1/0.003,
# 10edilution factor cells present)
all_data$count <- all_data$Number * 1/0.003 * 10 ** all_data$Dilution_factor
# log count column
all_data$logged <- log10(all_data$count)

head(all_data)
# Get the blank value for each experiment
all_data_blanks <- all_data %>%
  filter(Conc == 0) %>%
  group_by(experiment) %>%
  summarise(
    blank_value = mean(subset(logged, Conc == 0)),
  )
# Recombine
all_data_with_blanks <- all_data %>%
  left_join(all_data_blanks, by = "experiment")
# Get mean and standard deviation form each expeirment
summary_data <- all_data_with_blanks %>%
  group_by(Strain, Conc, Condition, pretreatment, experiment, blank_value) %>%
  summarise(
    mean_value = mean(logged, na.rm = TRUE),
    sd_value = sd(logged, na.rm = TRUE),
    .groups = "drop"
  )
# Get the log fold change
summary_data$logFC <- summary_data$mean_value - summary_data$blank_value
head(summary_data, n = 20)
head(all_data_with_blanks)

# COmplete the data so consistent bar spacing used
## Create a complete version of summary_data with all combinations of:
## Conc, pretreatment, and experiment
summary_data <- summary_data %>%
  complete(experiment, nesting(Condition, Conc),
           fill = list(logFC = 0, sd_value = NA))
## The nesting means it only fills where there is an existing combination
## of chelator and concentration, so we don't create lots of unused citrate
## concentrations

# Same thing applied to the individual replicate points
all_data_with_blanks <- all_data_with_blanks %>%
  complete(experiment, nesting(Condition, Conc),
           fill = list(logged = NA, replicate = 1))
head(summary_data)
# Plotting time
position_dodge_val <- position_dodge(width = 0.7)

ggplot(data = subset(summary_data, 
    Condition %in% c("Citrate", "EDTA") 
    #& pretreatment %in% c("M9","M9-citr")
    ))+
    # Plot mean
    geom_bar(stat = "identity",
             aes(x = as.factor(Conc), 
                 y = logFC,
                 fill = experiment,
                 group = experiment),
              position = position_dodge_val,
              width = 0.6)+
    # Add baseline
    geom_hline(yintercept = 0)+
    geom_errorbar(
        aes(x = as.factor(Conc),
            ymin = (logFC - sd_value),
            ymax = (logFC + sd_value),
            group = experiment),
            width = 0.7,
            position = position_dodge_val)+
    geom_point(data = all_data_with_blanks,
        aes(x = as.factor(Conc),
            y = logged-blank_value,
            col = pretreatment,
            shape = as.factor(Rep),
            group = experiment
        ),
        col = "black",
        size = 4,
        position = position_dodge_val,
        show.legend = FALSE
    )+
    xlab("Chelator concentration (mM)")+
    ylab("log-fold change")+
    facet_wrap(vars(Condition), nrow = 2, scales = "free")+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.1,0.12),
        legend.background = element_rect(fill = NA)
    )

ggsave(
    "chelCFU_summary.png",
    plot = last_plot(),
    device = png,
    path = image_path,
    scale = 2,
    width = 14,
    height = 14,
    units = c("cm"),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
)

# Previous method for plotting individual figures for each experiment

# Trying logging the data before calculating mean/sd
data$logged <- log10(data$count)
logged_mean_CFU <- 7.809411
# Compute summary stats
summary_log_data <- data %>%
  group_by(Strain, Conc, Condition) %>%
  summarise(
    mean_value = mean(logged, na.rm = TRUE),
    sd_value = sd(logged, na.rm = TRUE),
    .groups = "drop"
  )
head(summary_log_data)
head(data)
blank <-  filter(data, Conc == 0)
head(blank)
logged_mean_CFU <- mean(blank$logged)
print(logged_mean_CFU)
ggplot(data = summary_log_data)+
    geom_bar(stat = "identity",
             aes(x = as.factor(Conc), 
                 y = mean_value-logged_mean_CFU,
                 fill = Condition),
              position = "dodge",)+
    geom_hline(yintercept = log10(mean_CFU/mean_CFU))+
    geom_errorbar(
        aes(x = as.factor(Conc),
            ymin = (mean_value - sd_value)-logged_mean_CFU,
            ymax = (mean_value + sd_value)-logged_mean_CFU,
            group = Strain),
            width = 0.7,
            position = "dodge")+
    geom_point(data = data,
        aes(x = as.factor(Conc),
            y = log10(count)-logged_mean_CFU,
            col = Strain,
            shape = as.factor(Rep)
        ),
        col = "black",
        size = 4,
        #position = position_dodge(0.8)
    )+
    xlab("Chelator concentration")+
    ylab("log-fold change")+
    facet_wrap(vars(Condition))+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.1,0.12),
        legend.background = element_rect(fill = NA)
    )

ggsave(
    "chelCFU3_logged.png",
    plot = last_plot(),
    device = png,
    path = image_path,
    scale = 2,
    width = 14,
    height = 8,
    units = c("cm"),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
)


# Silwet CFU figure 

file_name <- "2025-05-16-Silwet02"
slw_data <- read.csv(paste(path,file_name,".csv",sep=""))

head(slw_data)

# Do the logging and blank calculations

# Generate count column
slw_data$count <- slw_data$Number * 1/0.003 * 10 ** data$Dilution_factor
# log count column
slw_data$logged <- log10(slw_data$count)

# Get blank for each strain? probably useful
slw_data_blanks <- slw_data %>%
  filter(Conc == 0) %>%
  group_by(Strain) %>%
  summarise(
    blank_value = mean(subset(logged, Conc == 0)),
  )
head(slw_data_blanks)
# Recombine
slw_data_with_blanks <- slw_data %>%
  left_join(slw_data_blanks, by = "Strain")


# Get mean and standard deviation form each expeirment
slw_summary_data <- slw_data_with_blanks %>%
  group_by(Strain, Conc, blank_value) %>%
  summarise(
    mean_value = mean(logged, na.rm = TRUE),
    sd_value = sd(logged, na.rm = TRUE),
    .groups = "drop"
  )
# Get the log fold change
slw_summary_data$logFC <- slw_summary_data$mean_value - slw_summary_data$blank_value
head(slw_summary_data, n = 20)
head(all_data_with_blanks)

# Plotting time
position_dodge_val <- position_dodge(width = 0.7)

ggplot(data = slw_summary_data)+
    # Plot mean
    geom_bar(stat = "identity",
             aes(x = as.factor(Conc), 
                 y = logFC,
                 fill = Strain,
                 group = Strain),
              )+
    # Add baseline
    geom_hline(yintercept = 0)+
    geom_errorbar(
        aes(x = as.factor(Conc),
            ymin = (logFC - sd_value),
            ymax = (logFC + sd_value),
            ),
            width = 0.7,
            )+
    geom_point(data = slw_data_with_blanks,
        aes(x = as.factor(Conc),
            y = logged-blank_value,
            shape = as.factor(Rep),
        ),
        col = "black",
        size = 4,
        show.legend = FALSE
    )+
    xlab("Silwet concentration (v/v %)")+
    ylab("log-fold change")+
    facet_wrap(vars(Strain), nrow = 1)+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.9,0.12),
        legend.background = element_rect(fill = NA)
    )

ggsave(
    "slwtCFU_logFC.png",
    plot = last_plot(),
    device = png,
    path = image_path,
    scale = 2,
    width = 14,
    height = 8,
    units = c("cm"),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
)
