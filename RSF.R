#
#
# Problem : distribute receives 0 length selections...
#
#

#
# Utils
#

character_to_numeric <- function(x) as.numeric(strsplit(x, ",")[[1]])

normalize <- function(sequence){
  sequence <- sapply(sequence, function(x, m) x - m, m = mean(sequence))
  sequence <- sapply(sequence, function(x, s) x / s, s = sd(sequence))
  return(sequence)
}

#
# Normalize all series in a df
#

normalize_df <- function(df){
  for (k in 1:nrow(df)){
    df[k,] <- normalize(df[k,])
  }
  return(df)
}

install.packages("data.tree")
library(data.tree)

#
# Calculate accuracy
#

accuracy <- function(results, ts){
  return(sum(results[!is.na(results[])] == ts[!is.na(results[])]) / length(ts[!is.na(results[])]))
}

load_df <- function(path){
  df <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
  
  # Use lapply to convert list of character to list of numeric
  # Use as.data.frame.vector to build data.frame containing vector classes instead of many columns
  
  df <- as.data.frame.vector(lapply(df[,1], character_to_numeric))
  
  # Use apply to get labels ([1]) from list (x) of vectors (x[[1]])
  
  y <- apply(df, 1, function(x) x[[1]][1])
  
  # Use lapply to get all elements other than labels from list of vector
  
  df <- as.data.frame.vector(lapply(df[,1], function(x) normalize(x[-1])))
  
  return(list(df, y))
}

#
# Implement the random shapelet forest algorithm
#

make_leaf <- function(labels, filter = NULL){
  leaf <- Node$new(sprintf("L%d", leaf_id))
  leaf_id <<- leaf_id + 1
  if (is.null(filter)) {
    leaf$label <- unique(labels)
  } else {
    leaf$label <- unique(labels[filter])
  }
  return(leaf)
}

#
# Generate the tree nodes by splitting according to best split (bs)
#

distribute <- function(data, labels, selection, distance_df, bs, min_ssize, max_ssize, nb_shapelets){
  # Filter to split distances < t and >= t
  best_shapelet_ds <- distance_df[,bs["max_ig_shapelet"]]
  left <- best_shapelet_ds < bs["threshold"]
  right <- best_shapelet_ds >= bs["threshold"]
  
  # Creating the root or new node
  node_is_pure <- length(unique(labels[selection])) == 1
  if (node_is_pure) {
    node <- make_leaf(labels[selection])
  } else {
    node <- Node$new(sprintf("B%d", branch_id))
    branch_id <<- branch_id + 1
    node$bs <- bs
    
    # Attach the child nodes to the current node
    if (sum(left) > 0){
      left_child <- random_shapelet_tree(data, labels, selection[left], min_ssize, max_ssize, nb_shapelets)
      left_child$id <- 1
      node$AddChildNode(left_child)
    }

    if (sum(right) > 0){
      right_child <- random_shapelet_tree(data, labels, selection[right], min_ssize, max_ssize, nb_shapelets)
      right_child$id <- 2
      node$AddChildNode(right_child)
    }
  }
  # Return the current node
  return(node)
}

#
# Compute the entropy of sub_ds
#

entropy <- function(labels, selection){
  label_sample <- labels[selection]
  label_list <- unique(label_sample)
  sample_card <- length(label_sample)
  p_label <- c()
  for (k in 1:length(label_list)){
    p_label[k] <- sum(label_sample == label_list[k]) / sample_card
  }
  return(-sum(p_label * log2(p_label)))
}

#
# Compute the information gain for a split in a and b ensembles
#

information_gain <- function(labels, selection, split_point){
  l_ds <- length(selection)
  l_a <- length(selection[1:split_point])
  l_b <- length(selection[(split_point + 1):l_ds])
  
  h_ds <- entropy(labels, selection)
  h_a <- entropy(labels, selection[1:split_point])
  h_b <- entropy(labels, selection[(split_point + 1):l_ds])
  
  return(h_ds - ((l_a/l_ds) * h_a + (l_b/l_ds) * h_b))
}

#
# Compute the distance between two subsequences
#

shapelet_distance <- function(s1, s2){
  return(sum((s1 - s2)^2))
}

#
# Calculate the distance between a shapelet and a series
#
# NOTE: doesn't handle cases when series_size % size != 0, which causes small losses of informations
#

ts_shapelet_distance <- function(candidate, series){
  series_size <- length(series)
  size <- length(candidate)
  # distances <- data.frame()
  distances <- c()
  k <- 1
  start <- 1
  while (start + size <= series_size + 1){
    # portion_distance <- data.frame("portion" = k, "distance" = shapelet_distance(candidate, series[start:(start - 1 + size)]))
    # distances <- rbind(distances, portion_distance)
    distances[[k]] <- shapelet_distance(candidate, series[start:(start - 1 + size)])
    k <- k + 1
    start <- start + size
  }
  
  # distances <- distances[order(distances[,"distance"]),]
  return(min(distances))
}

#
# Calculate information gain of each possible split for each shapelet to determine best threshold and split index
#

best_split <- function(labels, selection, distance_df){
  bs <- c(NA, NA, NA, NA, NA, NA)
  names(bs) <- c("max_ig_confidence", "series", "max_ig", "max_ig_index", "max_ig_shapelet", "threshold")
  for (k in 1:ncol(distance_df)){
    # order the selection by distance
    selection <- selection[order(distance_df[,k], decreasing = TRUE)]
    ds <- as.data.frame(distance_df[order(distance_df[,k], decreasing = TRUE)])
    l <- 1
    while (l <= nrow(ds)){
      is_best_split <- FALSE

      # Find split and calculate IG
      if (l < nrow(ds)){
        ig <- information_gain(labels, selection, split_point = l)
      } else {
        ig <- 0
      }

      if (is.na(bs["max_ig"])){
        is_best_split <- TRUE
      } else if (ig > bs["max_ig"]){
        is_best_split <- TRUE
      }

      if (is_best_split){
        # Store the best split informations
        bs["max_ig"] <- ig
        bs["max_ig_confidence"] <- ig * l
        bs["max_ig_index"] <- l
        bs["max_ig_shapelet"] <- k
        
        # Save the id of the actual series (in original series)
        bs["series"] <- selection[l]

        if (l < nrow(ds)){
          bs["threshold"] <- mean(ds[l:(l + 1), ])
        } else {
          bs["threshold"] <- ds[l, ]
        }
      }
      l <- l + 1
    }
  }
  return(bs)
}

#
# Extract shapelet
#

extract_shapelet <- function(ts, min_ssize, max_ssize){
  ssize <- sample(min_ssize:max_ssize, 1)
  sstart <- sample(1:(length(ts) - ssize), 1)
  return(ts[sstart:(sstart + ssize - 1)])
}

#
# Generate a candidate shapelet from the samples series
#

sample_shapelet <- function(data, selection, min_ssize, max_ssize){
  ts_id <- sample(1:nrow(data), 1)
  ts <- data[ts_id,]
  
  # Objective :
  # modify sample shapelet to sample until :
  # X% of the best similar portions of new_sample are not common with any of 
  # X% of the best similar portions of EACH shapelet from the same TS that generates a BS

  # The problem :
  # The distance between the sample shapelets and each TS are currently calculated AFTER the shapelets are sampled
  # The calculation of the distance between S and TS doesn't save positions of the portions
  # The information of the BS found for this sample aren't accessible -> SOLVED : global variable bs.df

  # calculate the distance between new_sample and "mother TS" to avoid twin shapelets
  
  # find brother shapelets among best split shapelets
  brothers <- bs.df["series" == ts_id,]
  shapelet <- extract_shapelet(ts[[1]], min_ssize, max_ssize)
  return(shapelet)
}

#
# Randomly generate nb_shapelets from the sample of series
#

generate_candidates <- function(data, selection, min_ssize, max_ssize, nb_shapelets){
  candidate_list <- list()
  for (k in 1:nb_shapelets){
    candidate_list[[k]] <- sample_shapelet(data, selection, min_ssize, max_ssize)
  }
  return(candidate_list)
}

#
# Compute the distances between all ts and candidate shapelets and store in a data.frame
#

candidate_distances <- function(data, selection, candidate_list){
  selection_size <- length(selection)
  distance_df <- matrix(nrow = selection_size, ncol = length(candidate_list))
  for (k in 1:selection_size){
    # The distance set for TS no. k is saved at index k
    distance_df[k,] <- sapply(candidate_list, ts_shapelet_distance, series = data[selection[k],][[1]])
  }
  return(distance_df)
}

#
# Generate a random shapelet tree
#

random_shapelet_tree <- function(data, labels, selection, min_ssize, max_ssize, nb_shapelets){
  candidate_list <- generate_candidates(data, selection, min_ssize, max_ssize, nb_shapelets)
  distance_df <- candidate_distances(data, selection, candidate_list)
  bs <- best_split(labels, selection, distance_df)
  node <- distribute(data, labels, selection, distance_df, bs, min_ssize, max_ssize, nb_shapelets)
  
  # Store bs for comparison with sampled shapelets
  bs.df <<- rbind(bs.df, bs)
  names(bs.df) <<- names(bs)
  
  # Save shapelet in node for classification
  node$shapelet <- candidate_list[[bs["max_ig_shapelet"]]]
  return(node)
}

#
# Generate the random shapelet forest
#

random_shapelet_forest <- function(data, labels, trees, min_ssize, max_ssize, nb_shapelets){
  rsf <- list()
  avg_b <- c()
  avg_l <- c()
  for (k in 1:trees){
    # initialize node counters
    branch_id <<- 1
    leaf_id <<- 1

    # initialize bs.df
    bs.df <<- data.frame()
    
    rsf[[k]] <- random_shapelet_tree(data, labels, sample(nrow(data), nrow(data)), min_ssize, max_ssize, nb_shapelets)

    avg_b[k] <- branch_id
    avg_l[k] <- leaf_id
    
    # Display info. about the RSF
    if (display_rst_info && k %% 5 == 0){
      print(sprintf("RSF size = %d, avg. branch = %.1f, avg. leaf = %.1f", k, mean(avg_b), mean(avg_l)))
    }
  }
  return(rsf)
}

#
# Return left or right child from id 1 or 2
#

get_child_node_by_id <- function(rst, id){
  for (k in 1:length(rst$children)){
    if (rst$children[[k]]$id == id){
      return(rst$children[[k]])
    }
  }
}

#
# Return the class label predicted by rsf for example
#

classify_example <- function(rst, example){
  if (!is.null(rst$label)){
    return(rst$label)
  }

  dist <- ts_shapelet_distance(rst$shapelet, example[[1]])
  if (dist < rst$bs["threshold"]){
    left_child <- get_child_node_by_id(rst, 1)
    return(classify_example(left_child, example))
  } else {
    right_child <- get_child_node_by_id(rst, 2)
    return(classify_example(right_child, example))
  }
}

#
# Define the hyper parameters according to Hills and al. (2014)
#

search_hyper_parameters <- function(data, sample_size, max_iteration, nb_shapelets){
  data_sample <- data[sample(nrow(data), sample_size),]
  shapelet_metrics <- data.frame()

  for (k in 1:max_iteration){
    # for each series of the sample
    for (l in 1:sample_size){
      # for each required shapelet
      for (m in 1:nb_shapelets){
        max_ssize <- sample(1:length(data_sample[[l]]), 1)
        min_ssize <- sample(1:max_ssize, 1)
        shapelet <- extract_shapelet(data_sample[[l]], min_ssize, max_ssize)
        distance <- ts_shapelet_distance(shapelet, data_sample[[l]])

        # Insert shapelet informations in df
        shapelet_data <- data.frame("series" = l, "max_ssize" = max_ssize, "min_ssize" = min_ssize, "distance" = distance)
        shapelet_metrics <- rbind(shapelet_metrics, shapelet_data)
      }
    }
  }

  # Only take the best candidate of each series
  bs <- data.frame()
  for (l in 1:sample_size){
    series_bs_list <- shapelet_metrics[shapelet_metrics[,"series"] == l,]
    series_bs <- series_bs_list[which.min(as.matrix(series_bs_list["distance"])),]
    bs <- rbind(bs, series_bs)
  }
  
  bs <- data.frame(bs["max_ssize"] - bs["min_ssize"])
  names(bs) <- c("length")
  bs <- as.data.frame(bs[order(bs["length"], decreasing = TRUE),])[c(floor(sample_size/4), 3 * floor(sample_size/4)),]
  
  # min_ssize = longest shapelet of the first quartile
  min_ssize <- bs[[2]]
  if (min_ssize == 0){
    min_ssize <- 1
  }

  # max_ssize = shortest shapelet of the fourth quartile
  max_ssize <- bs[[1]]
  if (max_ssize == 0){
    max_ssize <- 1
  }
  
  return(data.frame(min_ssize = c(min_ssize), max_ssize = c(max_ssize)))
}

#
# Generate a random shapelet forest from df and return classification of test_df
#

classify_series <- function(rsf, test_df) {
  results <- c()
  row_results <- c()
  for (k in 1:nrow(test_df)){
    for (l in 1:length(rsf)){
      rst <- rsf[[l]]
      row_results[l] <- classify_example(rst, test_df[k,])
    }
    results[k] <- names(which.max(table(row_results)))
    
    # Display percentage of examples classified
    if (display_classified_percent && (100 * k/nrow(test_df)) %% 10 == 0){
      print(sprintf("%.2f%% examples classified", 100 * k/nrow(test_df)))
    }
  }
  return(results)
}

#
# Verbosity parameters
#

display_hp <<- TRUE
display_rst_info <<- TRUE
display_classified_percent <<- FALSE

#
# Determine min./max. shapelet size
# Generate random shapelet forest
# Perform classification with rsf
# Calculate accuracy
#

perform_test <- function(){
  hp <- search_hyper_parameters(df, sample_size = 20, max_iteration = 20, nb_shapelets = 5)
  
  # Display hyperparameters
  if (display_hp){
    print(hp)
  }
  
  rsf <- random_shapelet_forest(df, y, trees = 5, hp$min_ssize, hp$max_ssize, nb_shapelets = 5)
  results <- classify_series(rsf, test_df)
  return(accuracy(results, test_y))
}

#
# Load training set
#

df_list <- load_df("C:/Users/Louis/Documents/Datascience/data mining/UCR_TS_Archive_2015/FISH/FISH_TRAIN")
df <- df_list[[1]]
y <- df_list[[2]]

#
# Load testing set
#

test_df_list <- load_df("C:/Users/Louis/Documents/Datascience/data mining/UCR_TS_Archive_2015/FISH/FISH_TEST")
test_df <- test_df_list[[1]]
test_y <- test_df_list[[2]]

accuracy <- perform_test()
print(accuracy)