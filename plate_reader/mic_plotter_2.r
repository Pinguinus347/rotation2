library(tidyverse)
library(ggplot2)
library(dplyr)
path <- "/home/callum/rotation2/plate_reader/long_data/"
image_path <- "/home/callum/rotation2/plate_reader/images/"
file_name <- "silwet01"
file_name2 <- "silwet02"
data <- read.csv(paste(path,file_name,".csv",sep=""))
data2 <- read.csv(paste(path,file_name2,".csv",sep=""))

# Plotting all data to identify anomalies
tail(data)
ggplot(data = data, aes(
        x = Time/(60*60),
        
    ))+
    geom_line(aes(y = (value),
        #col = as.factor(silwet_prc),
        #group = interaction(silwet_prc,column))
    ))+
    facet_wrap(vars(interaction(column, row)), nrow = 8)+
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

# Remove anomalies
data <- data %>% filter(variable != "E8")
head(data)
data$variable


total <- rbind(data, data2)
tail(total)

# 2. Compute summary stats
summary_data <- total %>%
  group_by(Time, species, silwet_prc, experiment) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    .groups = "drop"
  )
head(summary_data)
# 3. Plot the average OD at various timepoints
# times_hrs <- c(0,10,15,20)
times_hrs <- c(15)
times <- times_hrs * 60*60
# Must round times down as measurements are not taken directly on the hour
summary_data$Time_rounded <- signif(summary_data$Time, digits = 4)
summary_data$Time_rounded
head(subset(summary_data, Time %in% times))
ggplot(data = subset(summary_data, Time_rounded %in% times),
    aes(x = as.factor(silwet_prc)))+
    geom_point(aes(
            y = mean_value,
            col = as.factor(species),
            group = interaction(experiment),
            #shape = experiment,
        ),
        position = position_dodge(width = 0.5),
        size = 2)+
    facet_wrap(vars(Time_rounded/(60*60)))+

    geom_errorbar(
        aes(
            ymin = mean_value - sd_value,
            ymax = mean_value + sd_value,
            col = as.factor(species),
            group = interaction(experiment),
            #linetype = experiment
        ),
        position = position_dodge(width = 0.5),
        width = 0.5
    )+
    ylim(NA,0.6)+
    xlab("Silwet percentage")+
    ylab("Mean OD value")+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.12,0.9),
        legend.background = element_rect(fill = NA)
    )

ggsave(
    "slwt15h.png",
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


file_name <- "divion01"
data1 <- read.csv(paste(path,file_name,".csv",sep=""))
file_name <- "divion02"
data2 <- read.csv(paste(path,file_name,".csv",sep=""))
file_name <- "divion03"
data3 <- read.csv(paste(path,file_name,".csv",sep=""))
file_name <- "divion04"
data4 <- read.csv(paste(path,file_name,".csv",sep=""))
# Plotting all data to identify anomalies
tail(data)
ggplot(data = data4, aes(
        x = Time/(60*60),
        
    ))+
    geom_line(aes(y = (value),
        #col = as.factor(silwet_prc),
        #group = interaction(silwet_prc,column))
    ),
    col = "blue")+
    geom_line(data = data3, aes(x = Time/(60*60), y = (value),
        #col = as.factor(silwet_prc),
        #group = interaction(silwet_prc,column))
    ),
    col = "black")+
    geom_line(data = data2, aes(x = Time/(60*60), y = (value),
        #col = as.factor(silwet_prc),
        #group = interaction(silwet_prc,column))
    ),
    col = "red")+
    # geom_line(data = data1, aes(x = Time/(60*60), y = (value),
    #     #col = as.factor(silwet_prc),
    #     #group = interaction(silwet_prc,column))
    # ),
    # col = "green")+
    #xlim(0,24)+
    geom_vline(xintercept = 24, linetype = 2)+
    facet_wrap(vars(interaction(column, row)), nrow = 8)+
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


# Filter data3
data3 <- data3 %>% filter(column != "7")
data3 <- data3 %>% filter(column != "8")
data3 <- data3 %>% filter(variable != "D1")
data3 <- data3 %>% filter(variable != "D2")
data3 <- data3 %>% filter(variable != "D3")
data3 <- data3 %>% filter(variable != "D4")
data3 <- data3 %>% filter(variable != "D5")
data3 <- data3 %>% filter(variable != "D6")

head(data)
# data <- data %>% filter(variable != "B8")
# data <- data %>% filter(variable != "B9")
# data <- data %>% filter(variable != "B11")
# data <- data %>% filter(variable != "B12")
# data <- data %>% filter(variable != "E11")
# data <- data %>% filter(variable != "E12")


head(data)
data$variable
# 2. Compute summary stats
total <- rbind(data2, data3, data4)
summary_data <- total %>%
  group_by(Time, species, ion_conc, ion_type, experiment) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    .groups = "drop"
  )
head(summary_data)
# 3. Plot the average OD at various timepoints
times_hrs <- c(0,10,15,20)
times <- times_hrs * 60*60
# Must round the time to hours to plot on the hour as measurments are not made exactly on the hour
summary_data$Time_hrs <- summary_data$Time/(60*60)
summary_data$Time_hrs <- signif(summary_data$Time_hrs, digits = 3)
print(summary_data$Time_hrs)
# Must round times down as measurements are not taken directly on the hour
summary_data$Time_rounded <- signif(summary_data$Time, digits = 4)
summary_data$Time_rounded
head(subset(summary_data, Time %in% times))
ggplot(data = subset(summary_data, Time_hrs %in% 15),# & ion_type %in% "Mg"),
    aes(x = as.factor(ion_conc)))+
    geom_point(aes(
            y = mean_value,
            col = as.factor(species),
            group = interaction(experiment),
            #shape = as.factor(experiment),
        ),
        position = position_dodge(width = 0.5),
        size = 2)+
    facet_wrap(vars(ion_type))+

    geom_errorbar(
        aes(
            ymin = mean_value - sd_value,
            ymax = mean_value + sd_value,
            col = as.factor(species),
            group = interaction(experiment),
            #linetype = as.factor(ion_type)
        ),
        position = position_dodge(width = 0.5),
        width = 0.5
    )+
    #ylim(NA,0.6)+
    xlab("Major divalent ion concentration (mM)")+
    ylab("Mean OD value")+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.12,0.9),
        legend.background = element_rect(fill = NA)
    )

ggsave(
    "divion15h.png",
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


head(summary_data)
ggplot(data = subset(summary_data, ion_type %in% "Mg"),
    aes(x = Time/(60*60)))+
    geom_line(aes(
            y = mean_value,
            col = as.factor(ion_conc),
            linetype = interaction(experiment),
            #linetype = as.factor(chelator_type),
        ),
        size = 2)+
    geom_ribbon(aes(
            ymin = mean_value - sd_value, ymax = mean_value + sd_value,
            fill = as.factor(ion_conc),
            group = interaction(experiment, ion_conc)
            #group = interaction(ion_conc),
            #linetype = as.factor(ion_conc),
        ),
        alpha = 0.2
    )+
    ylim(0,NA)+
    xlim(0,24)+
    facet_wrap(vars(species), nrow = 2)


file_name <- "ctrEDTA01"
data1 <- read.csv(paste(path,file_name,".csv",sep=""))
file_name <- "ctrEDTA02"
data2 <- read.csv(paste(path,file_name,".csv",sep=""))
file_name <- "ctrEDTA03"
data3 <- read.csv(paste(path,file_name,".csv",sep=""))

# Plotting all data to identify anomalies
head(data)
ggplot(data = data3, aes(
        x = Time/(60*60),
    ))+
    geom_line(aes(y = (value),
        #col = as.factor(silwet_prc),
        #group = interaction(silwet_prc,column))
    ),
    col = "black")+
    #xlim(0,24)+
    geom_vline(xintercept = 24, linetype = 2)+
    facet_wrap(vars(interaction(column, row)), nrow = 8)+
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

#data <- data %>% filter(column != "12")
head(data)

total <- rbind(data1, data2, data3)
tail(total)

# 2. Compute summary stats
summary_data <- data3 %>%
  group_by(Time, species, chelator_conc, chelator_type,ion_conc, ion_type, experiment) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

# data$variable
# # 2. Compute summary stats
# summary_data <- data %>%
#   group_by(Time, species, chelator_conc, chelator_type, ion_conc, ion_type, experiment) %>%
#   summarise(
#     mean_value = mean(value, na.rm = TRUE),
#     sd_value = sd(value, na.rm = TRUE),
#     .groups = "drop"
#   )
# head(summary_data)
# 3. Plot the average OD at various timepoints
times_hrs <- c(15)
times <- times_hrs * 60*60
# Must round the time to hours to plot on the hour as measurments are not made exactly on the hour
summary_data$Time_hrs <- summary_data$Time/(60*60)
summary_data$Time_hrs <- signif(summary_data$Time_hrs, digits = 3)
print(summary_data$Time_hrs)
# Must round times down as measurements are not taken directly on the hour
summary_data$Time_rounded <- signif(summary_data$Time, digits = 4)
summary_data$Time_rounded
head(subset(summary_data, Time %in% times))

# summary_data$species <- factor(summary_data$species,
#                                levels = c("DC3000", "blank"),
#                                ordered = TRUE)

ggplot(data = subset(summary_data, Time_hrs %in% times_hrs),
    aes(x = as.factor(chelator_conc)))+
    geom_point(aes(
            y = mean_value,
            col = interaction(chelator_type,species),
            group = interaction(ion_conc),
            #shape = as.factor(ion_conc),
        ),
        position = position_dodge(width = 0.5),
        size = 2,
        show.legend = FALSE)+
    facet_wrap(vars(ion_conc), nrow=3)+

    geom_errorbar(
        aes(
            ymin = mean_value - sd_value,
            ymax = mean_value + sd_value,
            col = interaction(chelator_type,species),
            group = interaction(ion_conc),
            #linetype = as.factor(ion_conc)
        ),
        position = position_dodge(width = 0.5),
        width = 0.5,
        show.legend = FALSE
    )+
    #ylim(NA,0.6)+
    xlab("Chelator concentration (mM)")+
    ylab("Mean OD value")+
    theme(
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA),
        text = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = c(0.12,0.9),
        legend.background = element_rect(fill = NA),
        
    )

ggsave(
    "1448achltrs.png",
    plot = last_plot(),
    device = png,
    path = image_path,
    scale = 2,
    width = 8,
    height = 10,
    units = c("cm"),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
)
head(summary_data)
ggplot(data = summary_data,
    aes(x = Time/(60*60)))+
    geom_line(aes(
            y = mean_value,
            col = interaction(chelator_conc,species),
            #group = interaction(ion_conc),
            #linetype = as.factor(chelator_type),
        ),
        size = 2)+
    geom_ribbon(aes(
            ymin = mean_value - sd_value, ymax = mean_value + sd_value,
            fill = interaction(chelator_conc,species),
            #group = interaction(ion_conc),
            #linetype = as.factor(ion_conc),
        ),
        alpha = 0.2
    )+
    ylim(0,NA)+
    facet_wrap(vars(ion_conc,chelator_type), nrow = 2)

head(data1)













ggplot(data = subset(data1, row %in% "F" & column %in% c(4,5,6)))+
    geom_line(aes(y = (value), x = Time/(60*60),
        col = as.factor(column),
        #group = interaction(silwet_prc,column))
    ), show.legend = FALSE, linewidth = 1.5
    )+
    xlim(0,24)+
    #geom_vline(xintercept = 24, linetype = 2)+
    #facet_wrap(vars(interaction(column, row)), nrow = 8)+
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
ggsave(
    "ExampleOD600.png",
    plot = last_plot(),
    device = png,
    path = image_path,
    scale = 2,
    width = 5,
    height = 6,
    units = c("cm"),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
)
