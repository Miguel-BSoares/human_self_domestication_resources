###############FROM HERE WE COMPARE ALLELIC FREQUENCIES##################
#########################################################################
library(ggplot2)
library(BSDA)
library(qqman)
library(dplyr)
library(gtools)
library(tidyr)

#requires a file 
aggro_loci_position<-read.table("loci_position_aggro",header = TRUE,sep = " ")
aggro_loci_position_df<-as.data.frame(aggro_loci_position)

amh_freq_txt<-read.table("amh_freq.afreq",header = FALSE,sep = "\t")
nea_freq_txt<-read.table("neander_freq.afreq",header = FALSE,sep = "\t")

species<-cbind(aggro_loci_position[,1],aggro_loci_position[,2],amh_freq_txt[,6], nea_freq_txt[,6])
species_df<-as.data.frame(species)
species_df
colnames(species_df) <- c("loci", "position","amh", "neander")
print(species_df)

sorted_higher_maf <- species_df[order(as.numeric(as.character(species_df_filtered$neander))), ]
sorted<-species_df[order(species_df$neander),]

#stats# binomial test
# Step 0: Add the depth information  
# Merge the coverage data based on gene and position
species_df <- merge(species_df, aggro_sum_loci, by = c("loci","position"))
species_df <- species_df[ -c(7) ] #remove unecessary columns - example
str(species_df)

# Ensure columns are numeric
species_df$amh <- as.numeric(species_df$amh)
species_df$neander <- as.numeric(species_df$neander)
species_df$alt <- as.numeric(species_df$alt)

# Filter rows where coverage is greater than xx
species_df_filtered <- species_df[species_df$alt > 2, ]

# Log-transform the coverage (to stabilize variance)
species_df_filtered$log_alt <- log(species_df_filtered$alt)

# Now normalize the frequencies using the log-transformed coverage
species_df_filtered$amh_normalized_log <- species_df_filtered$amh / species_df_filtered$log_alt
species_df_filtered$neander_normalized_log <- species_df_filtered$neander / species_df_filtered$log_alt

# Perform the comparison (neander > amh) only on rows where coverage > 5
comparison_result <- species_df_filtered$neander_normalized_log > species_df_filtered$amh_normalized_log

# Count the number of successes (neander > amh)
num_successes <- sum(comparison_result)
num_trials <- length(comparison_result)

# Perform the binomial test
binom_test_10 <- binom.test(num_successes, num_trials, p = 0.5, alternative = "two.sided") # change according to different depths 
binom_test_2 <- binom.test(num_successes, num_trials, p = 0.5, alternative = "two.sided") # change according to different depths
binom_test_5 <- binom.test(num_successes, num_trials, p = 0.5, alternative = "two.sided") # change according to different depths
binom_test_10 <- binom.test(num_successes, num_trials, p = 0.5, alternative = "two.sided") # change according to different depths

### Plotting
species_df_long <- pivot_longer(species_df, cols = c("amh", "neander"),
                        names_to = "species", values_to = "value")
species_df_long$value <- as.numeric(species_df_long$value)
ggplot(species_df_long, aes(x = value, fill = species)) +
  geom_histogram(binwidth = 0.01, color = "black", alpha = 0.6, position = "identity") +
  labs(
    title = "MAF distribution aggro_db",
    x = "Value",
    y = "Frequency"
  ) +
  theme_minimal() +
  scale_fill_manual(
  values = c("amh" = "blue", "neander" = "red"),
  labels = c("amh" = "AMH", "neander" = "NEA")) +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    legend.position = "right",  # Position legend at the top
    legend.title = element_text(size = 12, face = "bold"),  # Customize legend title font
    legend.text = element_text(size = 10),  # Customize legend text size
    legend.key.size = unit(1, "cm"),  # Customize size of the legend key
    legend.key.width = unit(1, "cm"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)# Customize width of the legend key
  )

###Identify loci with alternative allelic frequencies above 0.5
neander_above_0_5 <- species_df_long[species_df_long$species == "neander" & species_df_long$value > 0.5, ]
print(neander_above_0_5)

amh_above_0_5 <- species_df_long[species_df_long$species == "amh" & species_df_long$value > 0.5, ]
print(amh_above_0_5)

################################################################
############## AGE DATA ###################3

##frequencies must have been obtained from vcf files 

bin_a<-read.table("bin_a_freq.afreq",header = FALSE,sep = "\t")
bin_b<-read.table("bin_b_freq.afreq",header = FALSE,sep = "\t")
bin_c<-read.table("bin_c_freq.afreq",header = FALSE,sep = "\t")
bin_d<-read.table("bin_d_freq.afreq",header = FALSE,sep = "\t")
bin_e<-read.table("bin_e_freq.afreq",header = FALSE,sep = "\t")
bin_f<-read.table("bin_f_freq.afreq",header = FALSE,sep = "\t")
bin_g<-read.table("bin_g_freq.afreq",header = FALSE,sep = "\t")
aggro_loci_position<-read.table("loci_position_aggro",header = TRUE,sep = " ")
aggro_loci_position_df<-as.data.frame(aggro_loci_position)
age_bins<-cbind(aggro_loci_position[,1],aggro_loci_position[,2],bin_a[,6],bin_b[,6],bin_c[,6],bin_d[,6],bin_e[,6],bin_f[,6],bin_g[,6])
age_bins_df<-as.data.frame(age_bins)
age_bins_df

colnames(age_bins_df) <- c("gene","position", "200-2000ybp", "3000-4800ybp","5000ybp-5800ybp","6000ybp-6800ybp","7000ybp-7900ybp","8000ybp-8500ybp","9000ybp-10800ybp")
str(age_bins_df)

age_bins_df[, 2:9]<-sapply(age_bins_df[, 2:9], as.numeric)

mean_values <- rowMeans(age_bins_df[, 3:9],na.rm = TRUE)  # Mean of columns 3 to 9 for each row
std_dev_values <- apply(age_bins_df[, 3:9], 1, sd,na.rm = TRUE)  # Standard deviation of columns 3 to 9 for each row

if (any(is.na(mean_values)) | any(is.na(std_dev_values))) {
  stop("There are missing values in the mean or standard deviation calculations.")
}

z_scores <- sweep(age_bins_df[, 3:9], 1, mean_values, FUN = "-")  # Subtract the mean
z_scores <- sweep(z_scores, 1, std_dev_values, FUN = "/")  # Divide by standard deviation

# Combine SNP names (column 1) with the Z-scores
z_score_df <- cbind(age_bins_df[, 1:2, drop = FALSE], z_scores)

# View the resulting data frame with SNP names and Z-scores
print(z_score_df)

####
mean_z_scores <- rowMeans(z_score_df[, 3:9], na.rm = TRUE)
std_dev_z_scores <- apply(z_score_df[, 3:9], 1, sd, na.rm = TRUE)


# View the updated data frame
print(z_score_df)
z_score_df<-z_score_df[with(z_score_df, order(z_score_df$gene, z_score_df$position)), ] #order

# Round all numeric values in the data frame to 2 decimal places
z_score_df[, columns_to_check] <- round(z_score_df[, columns_to_check], 2)

# Function to count values exceeding upper limit and return count
count_exceeding_limit <- function(row_values, upper_limit) {
  # Round values to 2 decimal places for consistency
  row_values <- round(row_values, 2)
  
  # Compare values to the upper limit
  comparisons <- row_values > upper_limit
  
  # Return the count of TRUE values (those exceeding the upper limit)
  return(sum(comparisons))
}

# Apply the comparison function row-wise for the entire data frame
z_score_df$CountExceeding1.25 <- apply(z_score_df, 1, function(row) {
  # Convert the row to a vector
  row_values <- as.numeric(row[columns_to_check])
  
  # Apply the function for upper limit 1.25
  count_exceeding_limit(row_values, 1.25)
})

z_score_df$CountExceeding2 <- apply(z_score_df, 1, function(row) {
  # Convert the row to a vector
  row_values <- as.numeric(row[columns_to_check])
  
  # Apply the function for upper limit 2
  count_exceeding_limit(row_values, 2)
})

z_score_df$CountExceeding2.5 <- apply(z_score_df, 1, function(row) {
  # Convert the row to a vector
  row_values <- as.numeric(row[columns_to_check])
  
  # Apply the function for upper limit 2.5
  count_exceeding_limit(row_values, 2.5)
})

# View the updated dataframe with the new columns
tail(z_score_df)  # View the first few rows to ensure the columns are added correctly
z_score_df
