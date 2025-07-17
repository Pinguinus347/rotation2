library(tidyverse)
library(ggplot2)

path <- "/home/callum/rotation2/plate_reader/long_data/"
file_name <- "silwet01_r"
full_file <- "silwet01"
file_name2 <- "silwet02"
data <- read.csv(paste(path,file_name,".csv",sep=""))
data2 <- read.csv(paste(path,file_name2,".csv",sep=""))
head(data2)

head(data)



data_full <- read.csv(paste(path,full_file,".csv",sep=""))
head(data)
tail(data)
blank_subset <- subset(data, species %in% "blank")
head(blank_subset)



# 3. normalise for blank by subtraction or division
# Get blank values as a baseline
blank_baseline <- data %>%
  filter(species == "blank") %>%
  select(Time, silwet_prc, blank_value = value)
head(blank_baseline)
# Step 2: Join with the original data
adjusted_data <- data %>%
  left_join(blank_baseline, by = c("Time", "silwet_prc")) %>%
  filter(species != "blank") %>%
  mutate(sub_value = value - blank_value)
head(adjusted_data)

# 2. Compute summary stats
summary_data <- adjusted_data %>%
  group_by(Time, species, silwet_prc, experiment) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    mean_sub = mean(value-blank_value, na.rm = TRUE),
    sd_sub = sd(value-blank_value, na.rm = TRUE),
    mean_div = mean(value/blank_value, na.rm = TRUE),
    sd_div = sd(value/blank_value, na.rm = TRUE),
    
    .groups = "drop"
  )
head(summary_data)


summary_data <- data %>%
  group_by(Time, species, silwet_prc) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),    
    .groups = "drop"
  )

times_hrs <- c(0,5,10,15)
times <- times_hrs * 60*60
times
summary_data$Time_rounded <- signif(summary_data$Time, digits = 4)
summary_data$Time_rounded
head(subset(summary_data, Time %in% times))
ggplot(data = subset(summary_data, Time_rounded %in% times),
    aes(x = as.factor(silwet_prc)))+
    geom_point(aes(
            y = mean_value,
            col = as.factor(species)
        ),
        size = 2)+
    facet_wrap(vars(Time_rounded/(60*60)))+

    geom_errorbar(
        aes(
            ymin = mean_value - sd_value,
            ymax = mean_value + sd_value,
            col = as.factor(species),
        ),
        width = 0.5
    )+
    xlab("Silwet percentage")+
    ylab("Mean OD value")
    geom_ribbon(aes(
            ymin = mean_sub - sd_sub, ymax = mean_sub + sd_sub,
            fill = as.factor(silwet_prc),
        ),
        alpha = 0.2
    )+
    ylim(0,NA)

summary_data <- data2 %>%
  group_by(Time, species, silwet_prc) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),    
    .groups = "drop"
  )
times_hrs <- c(0,5,10,15)
times <- times_hrs * 60*60
times
summary_data$Time_rounded <- signif(summary_data$Time, digits = 4)
summary_data$Time_rounded
head(subset(summary_data, Time %in% times))
ggplot(data = subset(summary_data, Time_rounded %in% times),
    aes(x = as.factor(silwet_prc)))+
    geom_point(aes(
            y = mean_value,
            col = as.factor(species)
        ),
        size = 2)+
    facet_wrap(vars(Time_rounded/(60*60)))+

    geom_errorbar(
        aes(
            ymin = mean_value - sd_value,
            ymax = mean_value + sd_value,
            col = as.factor(species),
        ),
        width = 0.5
    )+
    xlab("Silwet percentage")+
    ylab("Mean OD value")


# Plotting the subtraction adjusted data

ggplot(data = summary_data %>% filter(species != "blank"),
    aes(x = Time/(60*60)))+
    geom_line(aes(
            y = mean_sub,
            col = as.factor(silwet_prc)
        ),
        size = 2)+
    geom_ribbon(aes(
            ymin = mean_sub - sd_sub, ymax = mean_sub + sd_sub,
            fill = as.factor(silwet_prc),
        ),
        alpha = 0.2
    )+
    ylim(0,NA)+
    facet_wrap(vars(species), nrow = 2)
summary_data


# Getting the WT growth curve
wt_baseline <- summary_data %>%
  filter(silwet_prc == 0) %>%
  select(Time, species, wt_value = mean_div, wt_sub_value = mean_sub)
# Step 2: Join with the original data
adjusted_summary <- summary_data %>%
  left_join(wt_baseline, by = c("Time", "species")) #%>%
  #filter(silwet_prc != 0)
head(adjusted_summary)
adjusted_summary %>% select(wt_value,wt_sub_value)
# Plotting the division adjusted data

ggplot(data = summary_data %>% filter(species != "blank"),
    aes(x = Time/(60*60)))+
    geom_line(aes(
            y = mean_div,
            col = as.factor(silwet_prc)
        ),
        linewith = 2)+
    geom_ribbon(aes(
            ymin = mean_div - sd_div, ymax = mean_div + sd_div,
            fill = as.factor(silwet_prc),
        ),
        alpha = 0.3
    )+
    facet_wrap(vars(species), nrow = 2)

# Plotting the division and WT adjusted data

ggplot(data = adjusted_summary %>% filter(species != "blank"),
    aes(x = Time/(60*60)))+
    geom_line(aes(
            y = mean_div/wt_value,
            col = as.factor(silwet_prc)
        ),
        size = 2)+
    geom_ribbon(aes(
            ymin = (mean_div - sd_div)/wt_value, ymax = (mean_div + sd_div)/wt_value,
            fill = as.factor(silwet_prc),
        ),
        alpha = 0.3
    )+
    facet_wrap(vars(species), nrow = 2)+
    xlab("Time (hrs)")+
    ylab("(OD600/blank)/WT")+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.92,0.92),
        legend.background = element_rect(fill = NA)
    )

# Plotting the subtraction and WT adjusted data

ggplot(data = adjusted_summary %>% filter(species != "blank"),
    aes(x = Time/(60*60)))+
    geom_line(aes(
            y = mean_sub/wt_sub_value,
            col = as.factor(silwet_prc)
        ),
        size = 2)+
    geom_ribbon(aes(
            ymin = (mean_sub - sd_sub)/wt_sub_value, ymax = (mean_sub + sd_sub)/wt_sub_value,
            fill = as.factor(silwet_prc),
        ),
        alpha = 0.3
    )+
    facet_wrap(vars(species), nrow = 2)+
    ylim(0,NA)+
    xlab("Time (hrs)")+
    ylab("(OD600 - blank)/(WT - blank)")+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.92,0.92),
        legend.background = element_rect(fill = NA)
    )

# Plotting all data to identify anomalies
tail(data_full)
ggplot(data = data_full, aes(
        x = Time/(60*60),
        
    ))+
    geom_line(aes(y = (value),
        #col = as.factor(silwet_prc),
        #group = interaction(silwet_prc,column))
    ))+
    facet_wrap(vars(interaction(column, row)), ncol = 9, nrow = 8)+
    xlab("Time (hrs)")+
    ylab("OD600")+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.92,0.92),
        legend.background = element_rect(fill = NA)
    )
ggplot(data = data2, aes(
        x = Time/(60*60),
        
    ))+
    geom_line(aes(y = (value),
        #col = as.factor(silwet_prc),
        #group = interaction(silwet_prc,column))
    ))+
    facet_wrap(vars(interaction(column, row)), ncol = 12, nrow = 8)+
    xlab("Time (hrs)")+
    ylab("OD600")+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.92,0.92),
        legend.background = element_rect(fill = NA)
    )
# Just plotting the blank
ggplot(data = subset(data_full, species %in% "blank"), aes(
        x = Time/(60*60),
        
    ))+
    geom_line(aes(y = (value),
        col = as.factor(silwet_prc),
        #group = interaction(silwet_prc,column))
    ))+
    geom_line(aes(y=Temp/(mean(Temp))))+
    #facet_wrap(vars(interaction(column, row)), ncol = 9, nrow = 8)+
    xlab("Time (hrs)")+
    ylab("OD600")+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.92,0.92),
        legend.background = element_rect(fill = NA)
    )


library(dplyr)

# Compute rate of change (growth rate) by group
growth_data <- data %>%
  arrange(species, silwet_prc, column, Time) %>%
  group_by(species, silwet_prc, column) %>%
  mutate(
    delta_time = Time - lag(Time),
    delta_value = value - lag(value),
    growth_rate = delta_value / delta_time
  ) %>%
  ungroup()
head(growth_data)
ggplot(growth_data, aes(x = Time / (60 * 60), y = growth_rate)) +
  geom_line(aes(
    col = as.factor(silwet_prc),
    group = interaction(silwet_prc, column)
  )) +
  ylim(-0.000005,NA)+
  facet_wrap(vars(species)) +
  labs(y = "Growth rate", x = "Time (hours)")
head(summary_data)
summary_growth <- summary_data %>%
  arrange(species, silwet_prc, Time) %>%
  group_by(species, silwet_prc) %>%
  mutate(
    delta_time = Time - lag(Time),
    delta_value = mean_sub - lag(mean_sub),
    growth_rate = delta_value / delta_time
  ) %>%
  ungroup()
head(growth_data)
ggplot(summary_growth, aes(x = Time / (60 * 60), y = growth_rate)) +
  geom_line(aes(
    col = as.factor(silwet_prc),
    group = interaction(silwet_prc)
  )) +
  facet_wrap(vars(species)) +
  labs(y = "Growth rate", x = "Time (hours)")
