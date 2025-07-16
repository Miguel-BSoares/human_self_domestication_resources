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

species_df <- merge(species_df, aggro_sum_loci, by = c("loci","position"))
species_df <- species_df[ -c(7) ] #remove unecessary columns - example
str(species_df)

species_df$amh <- as.numeric(species_df$amh)
species_df$neander <- as.numeric(species_df$neander)
species_df$alt <- as.numeric(species_df$alt)

species_df_filtered <- species_df[species_df$alt > 2, ] ## choose depth at this step
species_df_filtered$log_alt <- log(species_df_filtered$alt)

species_df_filtered$amh_normalized_log <- species_df_filtered$amh / species_df_filtered$log_alt
species_df_filtered$neander_normalized_log <- species_df_filtered$neander / species_df_filtered$log_alt

comparison_result <- species_df_filtered$neander_normalized_log > species_df_filtered$amh_normalized_log ## check depth
num_successes <- sum(comparison_result)
num_trials <- length(comparison_result)

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

# z scores and new table
z_scores <- sweep(age_bins_df[, 3:9], 1, mean_values, FUN = "-") 
z_scores <- sweep(z_scores, 1, std_dev_values, FUN = "/")  
z_score_df <- cbind(age_bins_df[, 1:2, drop = FALSE], z_scores) # Combine SNP names (column 1) with the Z-scores
print(z_score_df)

####
mean_z_scores <- rowMeans(z_score_df[, 3:9], na.rm = TRUE)
std_dev_z_scores <- apply(z_score_df[, 3:9], 1, sd, na.rm = TRUE)

print(z_score_df)
z_score_df<-z_score_df[with(z_score_df, order(z_score_df$gene, z_score_df$position)), ] #order
z_score_df[, columns_to_check] <- round(z_score_df[, columns_to_check], 2) # Round all numeric values in the data frame to 2 decimal places

# Function to count values exceeding upper limit and return count
count_exceeding_limit <- function(row_values, upper_limit) {
  row_values <- round(row_values, 2)
  comparisons <- row_values > upper_limit
  return(sum(comparisons))
}

# Apply the comparison function row-wise for the entire data frame - can be looped for different upper limits
z_score_df$CountExceeding1.25 <- apply(z_score_df, 1, function(row) {
  row_values <- as.numeric(row[columns_to_check])
  count_exceeding_limit(row_values, 1.25) ## upper limit
})

z_score_df$CountExceeding2 <- apply(z_score_df, 1, function(row) {
  row_values <- as.numeric(row[columns_to_check])
  count_exceeding_limit(row_values, 2) ## upper limit
})

z_score_df$CountExceeding2.5 <- apply(z_score_df, 1, function(row) {
  row_values <- as.numeric(row[columns_to_check])
    count_exceeding_limit(row_values, 2.5) ## upper limit
})
