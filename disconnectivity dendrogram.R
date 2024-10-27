
#---------------------------------------------------------

# Load necessary libraries
library(jsonlite)
library(dendextend)  # For more control over dendrograms

# Define the directory path
directory_path <- '~/Desktop/IPM/Phasetransitionproject/Dendrogram data for R /'

# List all JSON files in the directory
json_files <- list.files(directory_path, pattern = "\\.json$", full.names = TRUE)

# Initialize a list to store the content of each JSON file
json_data_list <- list()

# Loop through each file and read the JSON content
for (file in json_files) {
  json_data <- fromJSON(file)
  
  # Extract file name for plot title
  name <- sub(".*DENDROGRAM(.*)\\.json$", "\\1", file)
  print(name)
  
  # Extract the Distance Matrix, LM, and LM Energy
  LM_energy <- json_data$`LM energy`
  LM <- json_data$`LM`
  
  # Convert LM to a character vector (if needed)
  LM <- as.character(LM)
  
  # Define your PATTERNS that should correspond to LM
  # Replace this with your actual PATTERNS vector
  PATTERNS <- json_data$`PATTERNS`
  
  # Check if LM and PATTERNS have the same length
  if (length(PATTERNS) != length(LM)) {
    stop("The length of the PATTERNS list must match the length of the LM list.")
  }
  
  # Number of features
  num_features <- length(LM)
  
  # Initialize an empty distance matrix
  distance_matrix <- matrix(0, nrow = num_features, ncol = num_features)
  
  # Create a map from feature names to indices
  feature_index <- setNames(seq_along(LM), LM)
  distance_dictionary <- json_data$`distance dictionary`
  
  # Fill the distance matrix using the converted_dict
  for (key in names(distance_dictionary)) {
    # Extract the feature pairs from the string key
    pair <- strsplit(gsub("[()']", "", key), ", ")[[1]]
    f1 <- pair[1]
    f2 <- pair[2]
    
    # Get the corresponding distance
    distance <- as.numeric(distance_dictionary[[key]])
    
    # Get the indices for the features in the distance matrix
    idx1 <- feature_index[[f1]]
    idx2 <- feature_index[[f2]]
    
    # Fill the distance matrix symmetrically
    distance_matrix[idx1, idx2] <- distance
    distance_matrix[idx2, idx1] <- distance
  }
  
  # Set row and column names for the distance matrix
  rownames(distance_matrix) <- LM
  colnames(distance_matrix) <- LM
  
  # Set custom heights using the LM labels
  custom_heights <- setNames(LM_energy, LM)
  
  # Perform hierarchical clustering
  hc <- hclust(as.dist(distance_matrix), method = "average")
  
  # Convert the hclust object to a dendrogram
  dend <- as.dendrogram(hc)
  
  # Function to set custom leaf heights using LM labels
  set_custom_leaf_height <- function(dend, custom_heights) {
    dend <- dendrapply(dend, function(node) {
      if (is.leaf(node)) {
        label <- attr(node, "label")
        if (!is.null(custom_heights[label])) {
          attr(node, "height") <- custom_heights[label]
        }
      }
      return(node)
    })
    return(dend)
  }
  
  # Apply custom hanging heights
  dend <- set_custom_leaf_height(dend, custom_heights)
  
  # Substitute LM labels with PATTERNS
  matching_indices <- match(labels(dend), LM)
  
  # Check for any NA values before substitution
  if (any(is.na(matching_indices))) {
    stop("Some labels in the dendrogram do not have a corresponding entry in LM.")
  }
  combined_labels <- paste0(PATTERNS[matching_indices], " (", LM[matching_indices], ")")
  
  # Update labels to the combined PATTERNS and LM
  labels(dend) <- combined_labels
  
  # Adjust the font size of the labels using dendextend
  dend <- dend %>% set("labels_cex", 0.5)  # 0.5 scales the label size down to 50%
  
  # Plot the customized dendrogram with combined labels and smaller font size
  plot(dend, main = name, ylab = 'Energy', xlab = 'Local Minima', ylim = c(-0.4, 1))
  
  # Get the order of the leaves in the dendrogram
  leaf_order <- order.dendrogram(dend)
  
  # Get the corresponding LM labels for the leaves in the dendrogram order
  ordered_LM <- LM[leaf_order]
  
  # Retrieve the x-coordinates for each leaf using the leaf order
  x_coords <- 1:length(ordered_LM)
  
  # Plot the orange circles at the terminal points of the leaves
  # The y-coordinates are taken from the custom heights corresponding to the ordered LM labels
  y_coords <- custom_heights[ordered_LM]
  
  # Plot the points on the dendrogram
  points(x_coords, y_coords, pch = 18, col = "orange", cex = 1)
  # Store each file's content
  json_data_list[[file]] <- json_data
  
  # Debug information
  print(name)
  print(LM)
  print(PATTERNS)
}

# End of the script

