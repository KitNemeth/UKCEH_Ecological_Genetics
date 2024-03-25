setwd("C:/Users/krinem/OneDrive - UKCEH")

library(knitr)
library(dplyr)
library(tidyverse)
library(readxl)
library(magrittr)
library(kableExtra)
library(genetics)
library(xtable)
library(reshape2)
library(janitor)
library(diffdf)
library(parallel)
library(here)
library(AGHmatrix)
library(svglite)

# we load the GRM
Gmat_VanRaden <- readRDS("Gmat_VanRaden.rds")

# Function to extract the genetic relatedness coefficients in a vector
extract_relatedness_coeffs_from_GRM <- function(GRM){
  # Extract upper triangle of the matrix
  upper_triangle <- upper.tri(GRM, diag=TRUE)
  
  # Set lower triangle elements to NA
  GRM[!upper_triangle] <- NA
  
  # Convert the resulting matrix to a dataframe
  df <- as.data.frame(as.table(GRM))
  
  # Remove rows with NA values (optional if you want to remove NA values)
  df <- na.omit(df)
  
  # Rename the columns
  names(df) <- c("Ind1", "Ind2", "r")
  
  df
}

df <- readRDS("FilteredGenomicData.rds")

# Get unique population codes
population_codes <- unique(df$PopulationCode)

# Initialize an empty list to store FieldCode for each population
FieldCode_list <- list()

# Loop through each population code
for (population_code in population_codes) {
  # Subset the data for the current population code
  subset_data <- subset(df, PopulationCode == population_code)
  # Get the FieldCode column for the subset
  FieldCode <- subset_data[["FieldCode"]]
  # Store the FieldCode in the list with the population code as the list name
  FieldCode_list[[population_code]] <- FieldCode
}

# Loop through each population
for (population_code in population_codes) {
  # Check if population code exists in FieldCode_list
  if (population_code %in% names(FieldCode_list)) {
    # Get the FieldCode list for the current population code
    FieldCode_list_current <- FieldCode_list[[population_code]]
    
    # Subset the df_r dataframe based on the current population's FieldCode list
    df_rsubset <- subset(df_r, Ind1 %in% FieldCode_list_current & Ind2 %in% FieldCode_list_current)
    
    # Define the thresholds to separate the degrees of relatedness
    KING_thresholds <- c(0.088, 0.177, 0.354, 0.707)
    
    # We count the number of individual pairs within each degree of relatedness 
    r_counts <- df_rsubset %>%
      filter(!Ind1 == Ind2) %>% 
      summarize(
        identical = sum(r > KING_thresholds[[4]]),
        first_degree = sum(r <= KING_thresholds[[4]] & r > KING_thresholds[[3]]),
        second_degree = sum(r <= KING_thresholds[[3]] & r > KING_thresholds[[2]]),
        third_degree = sum(r <= KING_thresholds[[2]] & r > KING_thresholds[[1]]),
        unrelated = sum(r <= KING_thresholds[[1]])
      )
    
    # Prepare the dataset for the plot 
    dfplot <- df_rsubset %>% 
      filter(!Ind1 == Ind2) %>%
      mutate(Index = row_number())
    
    # Define the breaks and labels for the y-scale
    y_breaks <- c(0, 0.25, 0.5, 0.75, 1, KING_thresholds)
    y_labels_colors <- ifelse(y_breaks %in% KING_thresholds, "orange", "black")
    
    # Plot
    plot <- ggplot(dfplot, aes(x = Index, 
                               y = r, 
                               color = cut(r, breaks = c(-Inf, KING_thresholds, Inf),
                                           labels = c(paste0("Unrelated (n = ", r_counts$unrelated, ")"), 
                                                      paste0("3rd degree (n = ", r_counts$third_degree, ")"), 
                                                      paste0("2nd degree (n = ", r_counts$second_degree, ")"), 
                                                      paste0("1st degree (n = ", r_counts$first_degree, ")"), 
                                                      paste0("Identical (n = ", r_counts$identical, ")"))))) +
      geom_point() +
      geom_hline(yintercept = KING_thresholds, linetype = "dashed", color = "orange") +
      labs(x = "Pairs of individuals", y = "Genetic relatedness coefficient") +
      theme_bw() +
      scale_color_manual(values = c("blue", "green", "deeppink", "purple", "yellow")) +
      scale_y_continuous(breaks = y_breaks) +
      theme(axis.text.y = element_text(angle = 0, hjust = 1, color = y_labels_colors),
            legend.position = c(0.2, 0.85),
            legend.title = element_blank()) +
      ggtitle(paste("Population Code:", population_code))
    
    # Define filename with population code
    out <- paste("plot_", gsub(" ", "_", population_code), ".svg", sep="")
    
    # Save the plot as an SVG file
    ggsave(filename = out, plot, device = "svg")
  } else {
    # Print a message if the population code does not exist in the FieldCode_list
    message("Population code '", population_code, "' does not exist in the FieldCode_list.")
  }
}
