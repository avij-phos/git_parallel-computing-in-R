
# Benchmark the performance of parallel methods (car, furrr, mirai) 


tryCatch({
  
  library(tictoc)
  library(vegan)
  library(car)
  library(future)
  library(furrr)
  library(ggplot2)
  library(dplyr)
  library(mirai)
  library(forcats)
  library(boot)
  cat("All libraries loaded successfully.\n")
}, error = function(e) {
  cat("ERROR: Failed to load one or more required libraries.\n")
  cat("Please ensure all packages (tictoc, vegan, car, future, furrr, ggplot2, dplyr, mirai, forcats) are installed.\n")
  stop(e)
})

# Checking current directory and file
cat("Current working directory:", getwd(), "\n")
cat("File exists:", file.exists("feature_table.txt"), "\n")

# Load the data
file_path <- "feature_table.txt"

if(file.exists(file_path)) {
  feature_table_raw <- read.table(
    file = file_path,
    sep = "\t",
    header = FALSE,
    row.names = 1,
    check.names = FALSE,
    fill = TRUE,
    quote = "",
    comment.char = "",
    skip = 2
  )
  cat("SUCCESS: Data loaded! Dimensions:", dim(feature_table_raw), "\n")
  cat("First few column names:", head(colnames(feature_table_raw)), "\n")
} else {
  cat("ERROR: File not found. Available files:\n")
  print(list.files())
}

cat("SUCCESS: Data loaded! Dimensions:", dim(feature_table_raw), "\n")

# 1. Check the actual structure
cat("\n--- Data Structure ---\n")
cat("Number of rows (features):", nrow(feature_table_raw), "\n")
cat("Number of columns:", ncol(feature_table_raw), "\n")
cat("Column names:", colnames(feature_table_raw), "\n")

# 2. Process the data -3 samples (V2, V3, V4) and taxonomy (V5)
taxonomy_data <- feature_table_raw[, 4, drop = FALSE]
feature_counts <- feature_table_raw[, 1:3]

cat("\nTaxonomy column dimensions:", dim(taxonomy_data), "\n")
cat("Feature counts dimensions:", dim(feature_counts), "\n")

# 3. Transpose so samples are rows and features are columns
count_matrix <- as.data.frame(t(feature_counts))

# 4. Convert to numeric and set sample names
count_matrix[] <- lapply(count_matrix, as.numeric)
#Naming the samples
rownames(count_matrix) <- c("Sample_V2", "Sample_V3", "Sample_V4")

cat("\nFinal count matrix dimensions:", dim(count_matrix), "\n")
cat("Samples (rows):", nrow(count_matrix), "\n") 
cat("Features (columns):", ncol(count_matrix), "\n")

# 5. Check sample totals
cat("\nTotal reads per sample:\n")
sample_totals <- rowSums(count_matrix)
print(sample_totals)

# 6. Check the data
cat("\nFirst few rows of count matrix:\n")
print(head(count_matrix[, 1:5]))

# Calculate diversity metrics
cat("\n--- Calculating Diversity Metrics ---\n")

#shannon diversity for each sample
shannon_diversity <- vegan::diversity(count_matrix, index = 'shannon')
cat("Shannon diversity for each sample:\n")
print(shannon_diversity)

# Simpson diversity
simpson_diversity <- diversity(count_matrix, index = 'simpson')
cat("\nSimpson diversity:\n")
print(simpson_diversity)

# Species richness
species_richness <- apply(count_matrix > 0, 1, sum)
cat("\nSpecies richness (number of features present):\n")
print(species_richness)

# Evenness (Pielou's evenness)
evenness <- ifelse(species_richness > 1, shannon_diversity / log(species_richness), 0)
cat("\nPielou's evenness:\n")
print(evenness)

# Create summary table
diversity_summary <- data.frame(
  Sample = rownames(count_matrix),
  Total_Reads = sample_totals,
  Shannon = shannon_diversity,
  Simpson = simpson_diversity,
  Richness = species_richness,
  Evenness = evenness
)
cat("\n=== DIVERSITY SUMMARY ===\n")
print(diversity_summary)

# Defining the Shannon diversity function for bootstrap
calculate_shannon <- function(data, indices) {
  resampled_data <- data[indices, ]
  shannon_values <- vegan::diversity(resampled_data, index = 'shannon')
  return(mean(shannon_values))
}

# Bootstrap analysis 
cat("\n--- Running Bootstrap Analysis ---\n")
set.seed(123)
# Running bootstrap 3 samples
boot_results <- boot::boot(count_matrix, statistic = calculate_shannon, R = 1000)
cat("\nBootstrap Results:\n")
print(boot_results)

#parallel benchmarking

# Define constants for parallel benchmarking
N_BOOTSTRAPS <- 1000 
N_CORES <- 4

cat("\n=== PARALLEL BENCHMARKING ===\n")
#standard bootstrap
cat("\n1.standard sequentail bootstrap:\n")
tic("Sequential")
boot_seq <- boot::boot(count_matrix, statistic = calculate_shannon, R = 1000)
time_seq <- toc()

# 2. car package Bootstrap:
cat("\n2. car package Bootstrap:\n")
tic("CAR")
boot_car <- boot(
  data = count_matrix, 
  statistic = calculate_shannon, 
  R = N_BOOTSTRAPS, 
  parallel = "multicore", 
  ncpus = N_CORES 
)
time_car <- toc()

# Method 3: Using future + furrr for parallel processing
cat("\n3. Future/Furrr Parallel Bootstrap:\n")
# Set up the parallel plan using N_CORES
plan(multisession, workers = N_CORES) 

tic("Future/Furrr")

# Using future_map_dbl from furrr package
boot_future <- future_map_dbl(1:N_BOOTSTRAPS, function(i) {
  indices <- sample(1:nrow(count_matrix), replace = TRUE) 
  calculate_shannon(count_matrix, indices)
}, .options = furrr_options(seed = TRUE))

time_future <- toc()

# Method 4: Using mirai for parallel processing

cat("\n4. Mirai Parallel Bootstrap (Simplified):\n") 

# Pre-generate indices
set.seed(123)
indices_list <- lapply(1:N_BOOTSTRAPS, function(i) {
  sample(1:nrow(count_matrix), size = nrow(count_matrix), replace = TRUE)
})

tic("Mirai")

# Start mirai daemons
daemons(n = N_CORES)

# Converting count_matrix to matrix for better compatibility
count_matrix_as_matrix <- as.matrix(count_matrix)

# Create multiple mirai tasks
mirai_tasks <- list()
for (i in 1:N_BOOTSTRAPS){
  mirai_tasks[[i]] <- mirai(
    {
      library(vegan)
      # Ensure we have valid data
      if (nrow(data) == 0 || ncol(data) == 0) {
        return(NA)
      }
      # Calculate diversity
      resampled_data <- data[indices, , drop = FALSE]
      shannon_values <- diversity(resampled_data, index = 'shannon')
      mean(shannon_values)
    },
    data = count_matrix_as_matrix,
    indices = indices_list[[i]]
  )
}

# Retrieve results
mirai_results <- lapply(1:N_BOOTSTRAPS, function(i) {
  result_obj <- call_mirai(mirai_tasks[[i]])
  if (is.numeric(result_obj$data)) result_obj$data else NA
})

boot_mirai <- unlist(mirai_results)
success_count <- sum(!is.na(boot_mirai))

cat("Successful tasks:", success_count, "/", N_BOOTSTRAPS, "\n")
time_mirai <- toc()

daemons(0)
#comparing results
cat("\n=== BENCHMARK RESULTS ===\n")
cat("sequentail bootstrap mean:", mean(boot_seq$t), "\n")
cat("CAR bootstrap mean:", mean(boot_car$t), "\n")
cat("Future/Furr bootstrap mean:", mean(boot_future), "\n")
cat("Mirai bootstrap mean:", round(mean(boot_mirai, na.rm = TRUE), 5), "\n")


# Creating benchmark results dataframe BEFORE calculating speedup
benchmark_results <- data.frame(
  Method = c("Sequential", "CAR", "Future/Furrr", "Mirai"),
  Time_seconds = c(
    time_seq$toc - time_seq$tic,
    time_car$toc - time_car$tic, 
    time_future$toc - time_future$tic,
    time_mirai$toc - time_mirai$tic
  )
)

#calculating speedup
sequential_time <- time_seq$toc - time_seq$tic
benchmark_results$Speedup <- sequential_time / benchmark_results$Time_seconds
# Calculate speedup compared to sequential
sequential_time <- time_seq$toc - time_seq$tic
benchmark_results$Speedup <- sequential_time / benchmark_results$Time_seconds

cat("\n=== PERFORMANCE COMPARISON ===\n")
print(benchmark_results)

# Create speedup visualization
p_speedup <- ggplot(benchmark_results, aes(x = fct_reorder(Method, Speedup), y = Speedup, fill = Method)) +
  geom_col() +
  geom_text(aes(label = round(Speedup, 2)), vjust = -0.5) +
  labs(
    title = "Parallel Speedup Comparison",
    x = "Method",
    y = "Speedup (vs Sequential)",
    subtitle = "Higher is better"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_speedup)
ggsave("benchmark_speedup.png", p_speedup, width = 10, height = 6, dpi = 300)
s



