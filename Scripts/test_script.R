install.packages("tidyverse")
library(tidyverse)

setwd("C:/Users/krinem/OneDrive - UKCEH/fastStructure")
fam_data <- read.table("newLEAFgenotypes_Mothers+Native_NoWCRegenOrDups.FAM", header=F, sep = "\t")
class(fam_data) #its useful to keep checking how your data is being treated - here its a DF

# Extract sample and family information
samplelist <- fam_data[, 2]
familylist <- fam_data[, 1]

# Create a data frame for storing the results
all_data <- tibble(sample = character(),
                   family = character(),
                   k = numeric(),
                   Q = character(),
                   value = numeric())
my_k <- 2:4

# Loop over values of k
for (k in c(2:4)) {
  # Read the data from the corresponding .Q file
  data <- read_delim(paste0("newLEAFgenotypes_Mothers+Native_NoWCRegenOrDups.", k, ".meanQ"),
                     col_names = paste0("Q", seq(1:k)),
                     delim = "  ")
  
  # Add sample, family, and k information to the data
  data$sample <- samplelist
  data$family <- familylist
  data$k <- k
  
  # Convert data from wide to long format
  data %>% gather(Q, value, -sample, -family, -k) -> data
  
  # Append the current iteration to the full dataset
  all_data <- rbind(all_data, data)
}


all_data %>%
  filter(k == 3) %>%
 ggplot(.,aes(x=sample,y=value,fill=factor(Q))) +
  geom_bar(stat="identity",position="stack",  width = 1) +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c("1","2", "3")) +
  facet_grid(~family, scales = 'free', space = 'free')

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

all_data %>%
  mutate(family = factor(family))  %>%
  ggplot(.,aes(x=factor(sample),y=value,fill=factor(Q))) +
  geom_bar(stat="identity",position="stack", width = 1) +
  xlab("Population") + ylab("Ancestry") +
  facet_grid(k ~ family, scales = "free_x", space = "free", labeller = labeller(k = grp.labs)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none") +
  scale_fill_brewer(palette="Spectral") +
  scale_x_discrete(labels = all_data$family) +  
  #facet_wrap(~k,ncol=1) + 
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    panel.grid = element_blank()
  ) 

# Assuming 'family' is the variable representing populations, replace it with your actual variable
all_data %>%
  ggplot(aes(x = sample, y = value, fill = factor(Q))) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  xlab("family") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none") +
  scale_fill_brewer(palette = "Spectral") +
  scale_x_discrete(labels = function(x) gsub(".*_", "", x)) +
  facet_wrap(~k, ncol = 1)

library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)

all_data %>%
  filter(k == 3) %>%
  ggplot(., aes(factor(sample), value, fill = factor(Q))) +
  geom_bar(color = "white", linewidth = 0.1, stat = "identity", position = "stack", width = 1) +
  facet_grid(~fct_inorder(family), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "", title = "K=3", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_text(color = "black"),
        panel.spacing.x = unit(0.1, "lines"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(angle = 45, vjust = 1, hjust = 1, size = 9, margin = margin(0, 20, 0, 15)),
  strip.background = element_rect(fill = "transparent", color = NA)) +
  scale_fill_gdocs(guide = FALSE)




cols = c("#b98d3e","#9970c1","#64a860","#cc545e")

Q_levels <- levels(factor(all_data$Q))  # Get the levels of Q

# Ensure the order of levels in Q is consistent
Q <- factor(all_data$Q, levels = Q_levels)



all_data <- mutate(all_data, color = cols[as.numeric(Q_levels) == as.numeric(Q)])





all_data %>%
  ggplot(., aes(sample, value, fill = factor(Q))) + 
  geom_bar(color = "black", stat = "identity", width = 1) +
  scale_fill_manual(values = cols) +
  labs(y = "Ancestry") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_grid(k ~ family, switch = "x", 
             scales = "free_x",
             space = "free", labeller = labeller(k = grp.labs)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),  # Hide x-axis labels
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(0.1, "lines"),
    panel.grid = element_blank(),
    strip.text.y.right = element_text(angle = 360, vjust = 0.5, hjust = 1),
    strip.text.x.bottom = element_text(angle = 80, vjust = 1, hjust = 1, margin = margin(0, 30, 0, 0)),
    strip.background.x = element_rect(fill = "transparent", color = NA)) +
    coord_cartesian(clip = "off")

#Problems
#I cant get the x axis labels to centre
#I cnt get the labels to stop cropping each other
#colours not consistent between plots?