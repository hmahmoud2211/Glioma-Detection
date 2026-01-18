# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required libraries
required_packages <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "GSVA", 
                       "caret", "dplyr", "tidyr", "xgboost", "randomForest",
                       "e1071", "glmnet", "nnet", "class", "MASS",
                       "Biobase", "GSEABase", "ggplot2", "pheatmap")

# Install Bioconductor packages if needed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 1. Load and prepare count data
prepare_data <- function(counts_path, meta_path) {
  tryCatch({
    # Load count data
    counts <- read.csv(counts_path, row.names = 1)
    
    # Load metadata
    meta_data <- read.csv(meta_path)
    
    # Clean up cancer type names and convert to factor
    meta_data$cancer_type <- as.factor(meta_data$Sample_characteristics_ch1)
    
    # Remove any potential duplicates in sample names
    meta_data <- meta_data[!duplicated(meta_data$ID_REF), ]
    
    # Ensure count matrix and metadata match
    common_samples <- intersect(colnames(counts), meta_data$ID_REF)
    if(length(common_samples) == 0) {
      stop("No matching samples between count data and metadata")
    }
    
    counts <- counts[, common_samples]
    meta_data <- meta_data[meta_data$ID_REF %in% common_samples, ]
    
    # Print summary
    cat("\nData summary:\n")
    cat("Number of genes:", nrow(counts), "\n")
    cat("Number of samples:", ncol(counts), "\n")
    cat("\nCancer types distribution:\n")
    print(table(meta_data$cancer_type))
    
    return(list(counts = counts, 
                meta_data = meta_data))
  }, error = function(e) {
    stop("Error in data preparation: ", conditionMessage(e))
  })
}

# 2. Perform DESeq2 Analysis
run_deseq2 <- function(counts, meta_data) {
  tryCatch({
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
      countData = round(counts), # ensure counts are integers
      colData = meta_data,
      design = ~ cancer_type
    )
    
    # Filter low count genes
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    # Set reference level as the most common cancer type
    most_common <- names(sort(table(meta_data$cancer_type), decreasing = TRUE)[1])
    dds$cancer_type <- relevel(dds$cancer_type, ref = most_common)
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Get results for each cancer type comparison
    cancer_types <- levels(meta_data$cancer_type)
    all_results <- list()
    
    for(cancer in cancer_types) {
      if(cancer != most_common) {
        res <- results(dds, contrast = c("cancer_type", cancer, most_common))
        res <- as.data.frame(res)
        res$GeneID <- rownames(res)
        all_results[[cancer]] <- res
      }
    }
    
    return(list(dds = dds, results = all_results))
  }, error = function(e) {
    stop("Error in DESeq2 analysis: ", conditionMessage(e))
  })
}

# 3. Perform GO Analysis
perform_go_analysis <- function(deseq_results) {
  tryCatch({
    go_results <- lapply(names(deseq_results), function(cancer) {
      res <- deseq_results[[cancer]]
      # Get significant genes (padj < 0.05 and |log2FoldChange| > 1)
      sig_genes <- res$GeneID[!is.na(res$padj) & 
                                res$padj < 0.05 & 
                                abs(res$log2FoldChange) > 1]
      
      if(length(sig_genes) == 0) {
        return(NULL)
      }
      
      # Run GO enrichment
      go_res <- enrichGO(gene = sig_genes,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
      
      return(go_res)
    })
    
    names(go_results) <- names(deseq_results)
    go_results <- go_results[!sapply(go_results, is.null)]
    
    if(length(go_results) == 0) {
      stop("No significant GO terms found")
    }
    
    return(go_results)
  }, error = function(e) {
    stop("Error in GO analysis: ", conditionMessage(e))
  })
}

# 4. Perform GSVA Analysis
perform_gsva_analysis <- function(dds, go_results) {
  tryCatch({
    # Get normalized counts
    norm_counts <- counts(dds, normalized = TRUE)
    
    # Create ExpressionSet
    expr_set <- ExpressionSet(assayData = norm_counts)
    
    # Process each GO result set
    gsva_results <- lapply(names(go_results), function(cancer_type) {
      go_res <- go_results[[cancer_type]]
      
      # Skip if no enriched terms
      if (nrow(go_res@result) == 0) {
        return(NULL)
      }
      
      # Extract gene sets and create unique names
      go_terms <- go_res@result$ID[1:min(50, nrow(go_res@result))]  # Limit to top 50 terms
      go_genes <- go_res@geneSets[go_terms]
      
      # Create unique names
      unique_names <- paste0(cancer_type, "_", seq_along(go_terms))
      
      # Create gene sets
      go_sets <- mapply(
        function(genes, name) {
          if (length(genes) >= 5) {  # Only keep sets with at least 5 genes
            GeneSet(as.character(genes), 
                    setName = name,
                    setIdentifier = name)
          } else {
            NULL
          }
        },
        go_genes,
        unique_names,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )
      
      # Remove NULL entries
      go_sets <- Filter(Negate(is.null), go_sets)
      
      if (length(go_sets) == 0) {
        return(NULL)
      }
      
      # Convert to GeneSetCollection
      go_sets <- GeneSetCollection(go_sets)
      
      # Run GSVA
      gsva_params <- GSVA::gsvaParam(
        expr = expr_set,
        geneSets = go_sets,
        minSize = 5,
        maxSize = 500
      )
      
      GSVA::gsva(gsva_params)
    })
    
    names(gsva_results) <- names(go_results)
    gsva_results <- gsva_results[!sapply(gsva_results, is.null)]
    
    if (length(gsva_results) == 0) {
      stop("No valid GSVA results obtained")
    }
    
    return(gsva_results)
    
  }, error = function(e) {
    stop("Error in GSVA analysis: ", conditionMessage(e))
  })
}

# ============================================================
# HELPER FUNCTIONS FOR HIERARCHICAL CLASSIFICATION
# ============================================================

# Extract glioma TYPE from cancer_type label
extract_glioma_type <- function(cancer_type) {
  type <- tolower(as.character(cancer_type))
  
  # Check for glioblastoma first (most specific)
  if (grepl("glioblastoma", type)) {
    return("Glioblastoma")
  }
  # Check for oligodendroastrocytoma (before checking for astrocytoma or oligodendroglioma)
  if (grepl("oligodendroastrocytoma|oligoastrocytoma", type)) {
    return("Oligodendroastrocytoma")
  }
  # Check for oligodendroglioma
  if (grepl("oligodendroglioma", type)) {
    return("Oligodendroglioma")
  }
  # Check for astrocytoma
  if (grepl("astrocytoma", type)) {
    return("Astrocytoma")
  }
  
  return("Unknown")
}

# Extract glioma GRADE from cancer_type label  
# Simplified to 3 classes for better accuracy: Low (II), High (III/IV), Recurrent
extract_glioma_grade <- function(cancer_type) {
  type <- tolower(as.character(cancer_type))
  
  # Check for recurrent first
  if (grepl("^recurrent", type)) {
    return("Recurrent")
  }
  # Check for glioblastoma (Grade IV) or anaplastic (Grade III) = High grade
  if (grepl("glioblastoma", type) || grepl("anaplastic", type)) {
    return("High_Grade")
  }
  # Default is Low Grade (Grade II)
  return("Low_Grade")
}

# ============================================================
# SIMPLIFIED HIGH-ACCURACY CLASSIFICATION APPROACH
# ============================================================
# Key principles:
# 1. NO synthetic oversampling (use class weights instead)
# 2. Simple but effective feature selection
# 3. Focus on robust models (RF, XGBoost, SVM)
# 4. Simple majority voting ensemble

# Feature selection using Random Forest importance
select_rf_features <- function(X, y, top_n = 200) {
  cat("    Feature selection using RF importance...\n")
  
  n_features <- ncol(X)
  if (n_features <= top_n) {
    cat("    Keeping all", n_features, "features\n")
    return(1:n_features)
  }
  
  # Quick RF to get feature importance
  rf_quick <- randomForest(X, y, ntree = 100, importance = TRUE)
  importance_scores <- importance(rf_quick)[, "MeanDecreaseGini"]
  
  # Select top features by importance
  top_n <- min(top_n, n_features)
  selected_idx <- order(importance_scores, decreasing = TRUE)[1:top_n]
  
  cat("    Selected top", length(selected_idx), "features by RF importance\n")
  return(selected_idx)
}

# Calculate balanced class weights
calculate_class_weights <- function(y) {
  class_counts <- table(y)
  total <- length(y)
  n_classes <- length(class_counts)
  
  # Balanced weights
  weights <- total / (n_classes * class_counts)
  return(weights)
}

# Train a single classifier - OPTIMIZED FOR HIGH ACCURACY
train_single_classifier <- function(X_train, y_train, X_test, y_test, task_name = "Classification") {
  
  model_results_list <- list()
  accuracy_scores <- c()
  all_predictions <- list()
  
  cat("\n========== Training", task_name, "Models ==========\n")
  cat("  Training samples:", nrow(X_train), "\n")
  cat("  Test samples:", nrow(X_test), "\n")
  cat("  Features:", ncol(X_train), "\n")
  cat("  Classes:", length(levels(y_train)), "\n")
  cat("  Class distribution:\n")
  print(table(y_train))
  
  # ----------------------
  # PREPROCESSING
  # ----------------------
  
  # Fix class labels (replace spaces with underscores for caret compatibility)
  original_levels <- levels(y_train)
  clean_levels <- make.names(original_levels)
  level_mapping <- setNames(original_levels, clean_levels)
  
  y_train_clean <- factor(make.names(as.character(y_train)), levels = clean_levels)
  y_test_clean <- factor(make.names(as.character(y_test)), levels = clean_levels)
  
  # Scale features
  X_train_scaled <- scale(X_train)
  scale_center <- attr(X_train_scaled, "scaled:center")
  scale_scale <- attr(X_train_scaled, "scaled:scale")
  scale_scale[scale_scale == 0 | is.na(scale_scale)] <- 1
  X_test_scaled <- scale(X_test, center = scale_center, scale = scale_scale)
  
  # Replace NA with 0
  X_train_scaled[is.na(X_train_scaled)] <- 0
  X_test_scaled[is.na(X_test_scaled)] <- 0
  
  # Class weights for imbalanced data
  class_weights <- calculate_class_weights(y_train_clean)
  cat("  Class weights:", paste(round(class_weights, 2), collapse = ", "), "\n")
  
  # ----------------------
  # Use all features (no selection for now - it doesn't help)
  # ----------------------
  X_train_selected <- X_train_scaled
  X_test_selected <- X_test_scaled
  
  # ----------------------
  # Use caret for proper cross-validated training
  # ----------------------
  set.seed(42)  # For reproducible CV
  
  # Create seeds for caret (needed for reproducibility)
  seeds <- vector(mode = "list", length = 26)  # 5 folds * 5 repeats + 1 final
  for(i in 1:25) seeds[[i]] <- sample.int(1000, 50)  # 50 = max tuning params
  seeds[[26]] <- sample.int(1000, 1)  # For final model
  
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 5,  # More repeats for stability
    classProbs = TRUE,
    verboseIter = FALSE,
    sampling = "down",  # Downsample majority class
    seeds = seeds
  )
  
  # ----------------------
  # 1. Random Forest with caret tuning
  # ----------------------
  cat("\n  [1/5] Random Forest (CV-tuned)...\n")
  tryCatch({
    rf_grid <- expand.grid(mtry = c(floor(sqrt(ncol(X_train_selected))), 
                                     floor(ncol(X_train_selected)/5),
                                     floor(ncol(X_train_selected)/10),
                                     floor(ncol(X_train_selected)/3)))
    
    rf_fit <- train(
      x = X_train_selected, y = y_train_clean,
      method = "rf",
      trControl = ctrl,
      tuneGrid = rf_grid,
      ntree = 2000,  # Increased for better performance
      importance = TRUE
    )
    
    rf_pred_clean <- predict(rf_fit, X_test_selected)
    rf_pred <- factor(level_mapping[as.character(rf_pred_clean)], levels = original_levels)
    rf_conf <- confusionMatrix(rf_pred, y_test)
    rf_acc <- rf_conf$overall["Accuracy"]
    accuracy_scores["RandomForest"] <- rf_acc
    all_predictions[["RandomForest"]] <- rf_pred
    model_results_list[["RandomForest"]] <- list(
      model = rf_fit$finalModel,
      predictions = rf_pred,
      confusion_matrix = rf_conf,
      importance = varImp(rf_fit)
    )
    cat("    Accuracy:", round(rf_acc, 4), "| Best mtry:", rf_fit$bestTune$mtry, "\n")
  }, error = function(e) {
    cat("    Failed:", conditionMessage(e), "\n")
  })
  
  # ----------------------
  # 2. XGBoost with caret tuning
  # ----------------------
  cat("\n  [2/4] XGBoost (direct training)...\n")
  tryCatch({
    # Convert class labels to numeric for xgboost
    y_numeric <- as.numeric(y_train_clean) - 1
    num_classes <- length(unique(y_train_clean))
    
    # Create DMatrix
    dtrain <- xgb.DMatrix(data = as.matrix(X_train_selected), label = y_numeric)
    dtest <- xgb.DMatrix(data = as.matrix(X_test_selected))
    
    # Set parameters
    params <- list(
      objective = "multi:softmax",
      num_class = num_classes,
      max_depth = 4,
      eta = 0.1,
      subsample = 0.8,
      colsample_bytree = 0.8,
      min_child_weight = 3
    )
    
    # Train with early stopping using internal CV
    xgb_model <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = 150,
      verbose = 0,
      early_stopping_rounds = 20,
      watchlist = list(train = dtrain)
    )
    
    # Predict
    xgb_pred_numeric <- predict(xgb_model, dtest)
    xgb_pred_clean <- factor(levels(y_train_clean)[xgb_pred_numeric + 1], levels = levels(y_train_clean))
    xgb_pred <- factor(level_mapping[as.character(xgb_pred_clean)], levels = original_levels)
    
    xgb_conf <- confusionMatrix(xgb_pred, y_test)
    xgb_acc <- xgb_conf$overall["Accuracy"]
    accuracy_scores["XGBoost"] <- xgb_acc
    all_predictions[["XGBoost"]] <- xgb_pred
    model_results_list[["XGBoost"]] <- list(
      model = xgb_model,
      predictions = xgb_pred,
      confusion_matrix = xgb_conf
    )
    cat("    Accuracy:", round(xgb_acc, 4), "\n")
  }, error = function(e) {
    cat("    Failed:", conditionMessage(e), "\n")
  })
  
  # ----------------------
  # 3. SVM with caret tuning (try both Radial and Linear)
  # ----------------------
  cat("\n  [3/5] SVM (CV-tuned)...\n")
  tryCatch({
    # Try Radial SVM
    svm_radial_grid <- expand.grid(
      C = c(0.1, 1, 10, 100),
      sigma = c(0.001, 0.01, 0.1)
    )
    
    svm_radial_fit <- train(
      x = X_train_selected, y = y_train_clean,
      method = "svmRadial",
      trControl = ctrl,
      tuneGrid = svm_radial_grid
    )
    
    svm_radial_pred_clean <- predict(svm_radial_fit, X_test_selected)
    svm_radial_pred <- factor(level_mapping[as.character(svm_radial_pred_clean)], levels = original_levels)
    svm_radial_acc <- confusionMatrix(svm_radial_pred, y_test)$overall["Accuracy"]
    
    # Try Linear SVM (often better for high-dimensional data)
    svm_linear_grid <- expand.grid(C = c(0.01, 0.1, 1, 10))
    
    svm_linear_fit <- train(
      x = X_train_selected, y = y_train_clean,
      method = "svmLinear",
      trControl = ctrl,
      tuneGrid = svm_linear_grid
    )
    
    svm_linear_pred_clean <- predict(svm_linear_fit, X_test_selected)
    svm_linear_pred <- factor(level_mapping[as.character(svm_linear_pred_clean)], levels = original_levels)
    svm_linear_acc <- confusionMatrix(svm_linear_pred, y_test)$overall["Accuracy"]
    
    # Use the better SVM
    if (svm_radial_acc >= svm_linear_acc) {
      svm_pred <- svm_radial_pred
      svm_fit <- svm_radial_fit
      svm_acc <- svm_radial_acc
      kernel_used <- "Radial"
    } else {
      svm_pred <- svm_linear_pred
      svm_fit <- svm_linear_fit
      svm_acc <- svm_linear_acc
      kernel_used <- "Linear"
    }
    
    svm_conf <- confusionMatrix(svm_pred, y_test)
    accuracy_scores["SVM"] <- svm_acc
    all_predictions[["SVM"]] <- svm_pred
    model_results_list[["SVM"]] <- list(
      model = svm_fit$finalModel,
      predictions = svm_pred,
      confusion_matrix = svm_conf
    )
    cat("    Accuracy:", round(svm_acc, 4), "| Best kernel:", kernel_used, 
        "| C:", svm_fit$bestTune$C, "\n")
  }, error = function(e) {
    cat("    Failed:", conditionMessage(e), "\n")
  })
  
  # ----------------------
  # 4. Elastic Net (glmnet)
  # ----------------------
  cat("\n  [4/5] Elastic Net (CV-tuned)...\n")
  tryCatch({
    glmnet_grid <- expand.grid(
      alpha = c(0, 0.5, 1),
      lambda = 10^seq(-4, 0, length = 20)
    )
    
    glmnet_fit <- train(
      x = X_train_selected, y = y_train_clean,
      method = "glmnet",
      trControl = ctrl,
      tuneGrid = glmnet_grid,
      family = "multinomial"
    )
    
    glmnet_pred_clean <- predict(glmnet_fit, X_test_selected)
    glmnet_pred <- factor(level_mapping[as.character(glmnet_pred_clean)], levels = original_levels)
    glmnet_conf <- confusionMatrix(glmnet_pred, y_test)
    glmnet_acc <- glmnet_conf$overall["Accuracy"]
    accuracy_scores["ElasticNet"] <- glmnet_acc
    all_predictions[["ElasticNet"]] <- glmnet_pred
    model_results_list[["ElasticNet"]] <- list(
      model = glmnet_fit$finalModel,
      predictions = glmnet_pred,
      confusion_matrix = glmnet_conf
    )
    cat("    Accuracy:", round(glmnet_acc, 4), "\n")
  }, error = function(e) {
    cat("    Failed:", conditionMessage(e), "\n")
  })
  
  # ----------------------
  # 5. K-Nearest Neighbors
  # ----------------------
  cat("\n  [6/6] KNN (CV-tuned)...\n")
  tryCatch({
    knn_grid <- expand.grid(k = c(3, 5, 7, 9))
    
    knn_fit <- train(
      x = X_train_selected, y = y_train_clean,
      method = "knn",
      trControl = ctrl,
      tuneGrid = knn_grid
    )
    
    knn_pred_clean <- predict(knn_fit, X_test_selected)
    knn_pred <- factor(level_mapping[as.character(knn_pred_clean)], levels = original_levels)
    knn_conf <- confusionMatrix(knn_pred, y_test)
    knn_acc <- knn_conf$overall["Accuracy"]
    accuracy_scores["KNN"] <- knn_acc
    all_predictions[["KNN"]] <- knn_pred
    model_results_list[["KNN"]] <- list(
      model = knn_fit$finalModel,
      predictions = knn_pred,
      confusion_matrix = knn_conf
    )
    cat("    Accuracy:", round(knn_acc, 4), "| Best k:", knn_fit$bestTune$k, "\n")
  }, error = function(e) {
    cat("    Failed:", conditionMessage(e), "\n")
  })
  
  # ----------------------
  # ENSEMBLE: Weighted Majority Voting (weights based on accuracy)
  # ----------------------
  cat("\n  [ENSEMBLE] Top-3 Weighted Majority Voting...\n")
  tryCatch({
    if (length(all_predictions) >= 2) {
      # Get weights from accuracy scores (exclude any failed models)
      valid_models <- names(all_predictions)
      valid_acc <- accuracy_scores[valid_models]
      
      # Only use top 3 models for ensemble
      top_n <- min(3, length(valid_acc))
      top_models <- names(sort(valid_acc, decreasing = TRUE))[1:top_n]
      valid_acc <- valid_acc[top_models]
      
      cat("    Using top", top_n, "models:", paste(top_models, collapse = ", "), "\n")
      
      # Normalize weights (higher accuracy = higher weight)
      weights <- valid_acc / sum(valid_acc)
      
      # Weighted voting using only top models
      pred_matrix <- sapply(all_predictions[top_models], as.character)
      
      ensemble_pred <- apply(pred_matrix, 1, function(row) {
        # Count weighted votes
        vote_counts <- sapply(original_levels, function(lev) {
          sum(weights[row == lev])
        })
        original_levels[which.max(vote_counts)]
      })
      ensemble_pred <- factor(ensemble_pred, levels = original_levels)
      
      ensemble_conf <- confusionMatrix(ensemble_pred, y_test)
      ensemble_acc <- ensemble_conf$overall["Accuracy"]
      accuracy_scores["Ensemble"] <- ensemble_acc
      all_predictions[["Ensemble"]] <- ensemble_pred
      model_results_list[["Ensemble"]] <- list(
        model = list(models = top_models),
        predictions = ensemble_pred,
        confusion_matrix = ensemble_conf
      )
      cat("    Ensemble Accuracy:", round(ensemble_acc, 4), "\n")
    }
  }, error = function(e) {
    cat("    Ensemble failed:", conditionMessage(e), "\n")
  })
  
  # Select best model
  if (length(accuracy_scores) == 0) {
    stop("No models were successfully trained for ", task_name)
  }
  
  accuracy_df <- data.frame(
    Model = names(accuracy_scores),
    Accuracy = as.numeric(accuracy_scores)
  )
  accuracy_df <- accuracy_df[order(accuracy_df$Accuracy, decreasing = TRUE), ]
  
  best_model_name <- accuracy_df$Model[1]
  best_accuracy <- accuracy_df$Accuracy[1]
  
  cat("\n  *** Best", task_name, "Model:", best_model_name, "with Accuracy:", round(best_accuracy, 4), "***\n")
  
  return(list(
    best_model_name = best_model_name,
    best_model = model_results_list[[best_model_name]],
    all_models = model_results_list,
    accuracy_comparison = accuracy_df,
    preprocessing = list(
      scale_center = scale_center,
      scale_scale = scale_scale
    )
  ))
}

# ============================================================
# 5. HIERARCHICAL TWO-STAGE CLASSIFICATION
# ============================================================
# Stage 1: Predict glioma TYPE (Astrocytoma, Oligodendroglioma, Oligodendroastrocytoma, Glioblastoma)
# Stage 2: Predict glioma GRADE (Grade II, Grade III, Grade IV, Recurrent)
# Final: Combine predictions for full cancer_type classification

train_model <- function(gsva_results, meta_data) {
  tryCatch({
    cat("\n##########################################################\n")
    cat("# HIERARCHICAL TWO-STAGE GLIOMA CLASSIFICATION PIPELINE  #\n")
    cat("##########################################################\n")
    
    # Combine GSVA scores from all comparisons
    all_scores <- do.call(rbind, lapply(gsva_results, function(x) {
      if (!is.null(x)) exprs(x) else NULL
    }))
    
    # Remove any duplicate row names
    all_scores <- all_scores[!duplicated(rownames(all_scores)), ]
    
    # Prepare the feature matrix
    X <- t(all_scores)
    
    # Prepare the original target variable
    y_original <- factor(meta_data$cancer_type)
    
    # Extract TYPE and GRADE labels
    y_type <- factor(sapply(meta_data$cancer_type, extract_glioma_type))
    y_grade <- factor(sapply(meta_data$cancer_type, extract_glioma_grade))
    
    # Ensure X and y have matching samples
    common_samples <- intersect(rownames(X), meta_data$ID_REF)
    X <- X[common_samples, ]
    y_original <- y_original[match(common_samples, meta_data$ID_REF)]
    y_type <- y_type[match(common_samples, meta_data$ID_REF)]
    y_grade <- y_grade[match(common_samples, meta_data$ID_REF)]
    
    cat("\n========== Data Summary ==========\n")
    cat("Total samples:", length(common_samples), "\n")
    cat("Total features:", ncol(X), "\n")
    cat("\nGlioma TYPE distribution:\n")
    print(table(y_type))
    cat("\nGlioma GRADE distribution:\n")
    print(table(y_grade))
    cat("\nOriginal cancer_type distribution:\n")
    print(table(y_original))
    
    # Split data - use 85/15 split for more training data (small dataset)
    set.seed(123)
    train_index <- createDataPartition(y_original, p = 0.85, list = FALSE)
    
    X_train <- X[train_index, ]
    X_test <- X[-train_index, ]
    
    y_original_train <- y_original[train_index]
    y_original_test <- y_original[-train_index]
    
    y_type_train <- y_type[train_index]
    y_type_test <- y_type[-train_index]
    
    y_grade_train <- y_grade[train_index]
    y_grade_test <- y_grade[-train_index]
    
    # ============================================================
    # STAGE 1: GLIOMA TYPE CLASSIFICATION
    # ============================================================
    cat("\n##########################################################\n")
    cat("# STAGE 1: GLIOMA TYPE CLASSIFICATION                     #\n")
    cat("# (Astrocytoma, Oligodendroglioma, Oligodendroastrocytoma, Glioblastoma)\n")
    cat("##########################################################\n")
    
    type_results <- train_single_classifier(
      X_train, y_type_train, X_test, y_type_test, 
      task_name = "TYPE"
    )
    
    # ============================================================
    # STAGE 2: GLIOMA GRADE CLASSIFICATION
    # ============================================================
    cat("\n##########################################################\n")
    cat("# STAGE 2: GLIOMA GRADE CLASSIFICATION                    #\n")
    cat("# (Grade II, Grade III, Grade IV, Recurrent)              #\n")
    cat("##########################################################\n")
    
    grade_results <- train_single_classifier(
      X_train, y_grade_train, X_test, y_grade_test, 
      task_name = "GRADE"
    )
    
    # ============================================================
    # COMBINED PREDICTION: Reconstruct full cancer_type
    # ============================================================
    cat("\n##########################################################\n")
    cat("# COMBINED HIERARCHICAL PREDICTION                        #\n")
    cat("##########################################################\n")
    
    # Get best predictions from each stage
    type_pred <- type_results$best_model$predictions
    grade_pred <- grade_results$best_model$predictions
    
    # Combine predictions to reconstruct the full cancer_type
    combined_pred <- mapply(function(type, grade) {
      type <- as.character(type)
      grade <- as.character(grade)
      
      # Handle Glioblastoma specially
      if (type == "Glioblastoma") {
        if (grade == "Recurrent") {
          return("recurrent Glioblastomas")
        } else if (grade == "Grade_IV_Primary") {
          return("primary Glioblastomas")
        } else if (grade == "Grade_IV_Secondary") {
          return("secondary Glioblastomas")
        } else {
          return("primary Glioblastomas")  # Default for GBM
        }
      }
      
      # Build the cancer type string for other types
      type_lower <- tolower(type)
      if (type_lower == "oligodendroastrocytoma") {
        type_plural <- "oligodendroastrocytomas"
      } else {
        type_plural <- paste0(type_lower, "s")
      }
      
      if (grade == "Recurrent") {
        # Check if it was anaplastic recurrent
        if (grepl("anaplastic", as.character(y_original_test[1]), ignore.case = TRUE)) {
          return(paste("recurrent anaplastic", type_plural))
        }
        return(paste("recurrent", type_plural))
      } else if (grade == "Grade_III") {
        return(paste("anaplastic", type_plural))
      } else {
        return(type_plural)
      }
    }, type_pred, grade_pred)
    
    combined_pred <- factor(combined_pred, levels = levels(y_original))
    
    # Calculate combined accuracy
    # Handle cases where combined prediction doesn't match exact labels
    combined_correct <- sum(combined_pred == y_original_test, na.rm = TRUE)
    combined_acc <- combined_correct / length(y_original_test)
    
    cat("\nCombined Hierarchical Accuracy:", round(combined_acc, 4), "\n")
    
    # ============================================================
    # ALSO TRAIN DIRECT FULL CLASSIFICATION FOR COMPARISON
    # ============================================================
    cat("\n##########################################################\n")
    cat("# DIRECT FULL CLASSIFICATION (for comparison)            #\n")
    cat("##########################################################\n")
    
    direct_results <- train_single_classifier(
      X_train, y_original_train, X_test, y_original_test, 
      task_name = "DIRECT (Full cancer_type)"
    )
    
    # ============================================================
    # FINAL COMPARISON AND SUMMARY
    # ============================================================
    cat("\n##########################################################\n")
    cat("# FINAL MODEL COMPARISON SUMMARY                         #\n")
    cat("##########################################################\n")
    
    final_comparison <- data.frame(
      Approach = c(
        paste("Stage1_TYPE:", type_results$best_model_name),
        paste("Stage2_GRADE:", grade_results$best_model_name),
        "Combined_Hierarchical",
        paste("Direct:", direct_results$best_model_name)
      ),
      Accuracy = c(
        type_results$accuracy_comparison$Accuracy[1],
        grade_results$accuracy_comparison$Accuracy[1],
        combined_acc,
        direct_results$accuracy_comparison$Accuracy[1]
      )
    )
    
    cat("\n")
    print(final_comparison)
    
    # Determine overall best approach
    hierarchical_effective_acc <- min(
      type_results$accuracy_comparison$Accuracy[1],
      grade_results$accuracy_comparison$Accuracy[1]
    )
    direct_acc <- direct_results$accuracy_comparison$Accuracy[1]
    
    cat("\n========== RECOMMENDATION ==========\n")
    cat("TYPE Classification Accuracy:   ", round(type_results$accuracy_comparison$Accuracy[1], 4), "\n")
    cat("GRADE Classification Accuracy:  ", round(grade_results$accuracy_comparison$Accuracy[1], 4), "\n")
    cat("Direct Full Classification:     ", round(direct_acc, 4), "\n")
    
    # Create combined accuracy table for all models
    all_accuracy <- rbind(
      data.frame(
        Model = paste0("TYPE_", type_results$accuracy_comparison$Model),
        Accuracy = type_results$accuracy_comparison$Accuracy,
        Stage = "Type"
      ),
      data.frame(
        Model = paste0("GRADE_", grade_results$accuracy_comparison$Model),
        Accuracy = grade_results$accuracy_comparison$Accuracy,
        Stage = "Grade"
      ),
      data.frame(
        Model = paste0("DIRECT_", direct_results$accuracy_comparison$Model),
        Accuracy = direct_results$accuracy_comparison$Accuracy,
        Stage = "Direct"
      )
    )
    all_accuracy <- all_accuracy[order(all_accuracy$Accuracy, decreasing = TRUE), ]
    
    cat("\n========== ALL MODELS RANKED ==========\n")
    print(all_accuracy)
    
    # Return comprehensive results
    return(list(
      # Hierarchical results
      type_classifier = type_results,
      grade_classifier = grade_results,
      combined_predictions = combined_pred,
      combined_accuracy = combined_acc,
      
      # Direct classification results
      direct_classifier = direct_results,
      
      # Best model info (use direct classifier as the main one for compatibility)
      best_model_name = paste0("Hierarchical_", type_results$best_model_name, "_", grade_results$best_model_name),
      model = list(
        type_model = type_results$best_model$model,
        grade_model = grade_results$best_model$model
      ),
      predictions = combined_pred,
      actual = y_original_test,
      confusion_matrix = direct_results$best_model$confusion_matrix,
      importance = type_results$best_model$importance,
      feature_matrix = X,
      
      # All models for reference
      all_models = list(
        type = type_results$all_models,
        grade = grade_results$all_models,
        direct = direct_results$all_models
      ),
      accuracy_comparison = all_accuracy,
      
      # Stage-specific results
      type_accuracy = type_results$accuracy_comparison,
      grade_accuracy = grade_results$accuracy_comparison,
      final_comparison = final_comparison
    ))
    
  }, error = function(e) {
    stop("Error in model training: ", conditionMessage(e))
  })
}

# Legacy function for backwards compatibility - trains only direct classification
train_model_direct <- function(gsva_results, meta_data) {
  tryCatch({
    # Combine GSVA scores from all comparisons
    all_scores <- do.call(rbind, lapply(gsva_results, function(x) {
      if (!is.null(x)) exprs(x) else NULL
    }))
    
    # Remove any duplicate row names
    all_scores <- all_scores[!duplicated(rownames(all_scores)), ]
    
    # Prepare the feature matrix
    X <- t(all_scores)
    
    # Prepare the target variable
    y <- factor(meta_data$cancer_type)
    
    # Ensure X and y have matching samples
    common_samples <- intersect(rownames(X), meta_data$ID_REF)
    X <- X[common_samples, ]
    y <- y[match(common_samples, meta_data$ID_REF)]
    
    # Split data
    set.seed(123)
    train_index <- createDataPartition(y, p = 0.7, list = FALSE)
    X_train <- X[train_index, ]
    X_test <- X[-train_index, ]
    y_train <- y[train_index]
    y_test <- y[-train_index]
    
    # Store all model results
    model_results_list <- list()
    accuracy_scores <- c()
    
    cat("\n========== Training Multiple Models ==========\n")
    
    # ----------------------
    # 1. Random Forest
    # ----------------------
    cat("\n[1/5] Training Random Forest...\n")
    tryCatch({
      rf_model <- randomForest(X_train, y_train, ntree = 500, importance = TRUE)
      rf_pred <- predict(rf_model, X_test)
      rf_conf <- confusionMatrix(rf_pred, y_test)
      rf_acc <- rf_conf$overall["Accuracy"]
      accuracy_scores["RandomForest"] <- rf_acc
      model_results_list[["RandomForest"]] <- list(
        model = rf_model,
        predictions = rf_pred,
        confusion_matrix = rf_conf,
        importance = importance(rf_model)
      )
      cat("  Random Forest Accuracy:", round(rf_acc, 4), "\n")
    }, error = function(e) {
      cat("  Random Forest failed:", conditionMessage(e), "\n")
    })
    
    # ----------------------
    # 2. XGBoost
    # ----------------------
    cat("\n[2/5] Training XGBoost...\n")
    tryCatch({
      y_numeric_train <- as.numeric(y_train) - 1
      y_numeric_test <- as.numeric(y_test) - 1
      num_classes <- length(levels(y))
      class_labels <- levels(y)
      
      dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_numeric_train)
      dtest <- xgb.DMatrix(data = as.matrix(X_test), label = y_numeric_test)
      
      params <- list(
        objective = "multi:softmax",
        num_class = num_classes,
        eval_metric = "mlogloss",
        eta = 0.1,
        max_depth = 6,
        subsample = 0.8,
        colsample_bytree = 0.8
      )
      
      xgb_model <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = 100,
        watchlist = list(train = dtrain, test = dtest),
        early_stopping_rounds = 10,
        verbose = 0
      )
      
      xgb_pred_numeric <- predict(xgb_model, dtest)
      xgb_pred <- factor(class_labels[xgb_pred_numeric + 1], levels = class_labels)
      xgb_conf <- confusionMatrix(xgb_pred, y_test)
      xgb_acc <- xgb_conf$overall["Accuracy"]
      accuracy_scores["XGBoost"] <- xgb_acc
      model_results_list[["XGBoost"]] <- list(
        model = xgb_model,
        predictions = xgb_pred,
        confusion_matrix = xgb_conf,
        importance = xgb.importance(feature_names = colnames(X_train), model = xgb_model),
        class_labels = class_labels
      )
      cat("  XGBoost Accuracy:", round(xgb_acc, 4), "\n")
    }, error = function(e) {
      cat("  XGBoost failed:", conditionMessage(e), "\n")
    })
    
    # ----------------------
    # 3. Support Vector Machine (SVM)
    # ----------------------
    cat("\n[3/5] Training SVM...\n")
    tryCatch({
      svm_model <- svm(X_train, y_train, kernel = "radial", probability = TRUE)
      svm_pred <- predict(svm_model, X_test)
      svm_conf <- confusionMatrix(svm_pred, y_test)
      svm_acc <- svm_conf$overall["Accuracy"]
      accuracy_scores["SVM"] <- svm_acc
      model_results_list[["SVM"]] <- list(
        model = svm_model,
        predictions = svm_pred,
        confusion_matrix = svm_conf,
        importance = NULL
      )
      cat("  SVM Accuracy:", round(svm_acc, 4), "\n")
    }, error = function(e) {
      cat("  SVM failed:", conditionMessage(e), "\n")
    })
    
    # ----------------------
    # 4. Neural Network (nnet)
    # ----------------------
    cat("\n[4/5] Training Neural Network...\n")
    tryCatch({
      # Scale data for neural network
      X_train_scaled <- scale(X_train)
      X_test_scaled <- scale(X_test, center = attr(X_train_scaled, "scaled:center"),
                              scale = attr(X_train_scaled, "scaled:scale"))
      
      # Use PCA to reduce dimensionality and avoid too many weights
      n_features <- ncol(X_train_scaled)
      max_components <- min(30, n_features, nrow(X_train_scaled) - 1)  # Limit to 30 components
      
      pca_model <- prcomp(X_train_scaled, center = FALSE, scale. = FALSE)
      X_train_pca <- pca_model$x[, 1:max_components, drop = FALSE]
      X_test_pca <- predict(pca_model, X_test_scaled)[, 1:max_components, drop = FALSE]
      
      # Calculate appropriate hidden layer size (avoid too many weights)
      n_classes <- length(levels(y_train))
      hidden_size <- min(10, max(3, floor(max_components / 2)))
      
      nnet_model <- nnet(X_train_pca, class.ind(y_train), 
                         size = hidden_size, softmax = TRUE, maxit = 500, 
                         trace = FALSE, MaxNWts = 5000, decay = 0.01)
      nnet_pred_prob <- predict(nnet_model, X_test_pca)
      nnet_pred <- factor(levels(y)[apply(nnet_pred_prob, 1, which.max)], levels = levels(y))
      nnet_conf <- confusionMatrix(nnet_pred, y_test)
      nnet_acc <- nnet_conf$overall["Accuracy"]
      accuracy_scores["NeuralNetwork"] <- nnet_acc
      model_results_list[["NeuralNetwork"]] <- list(
        model = nnet_model,
        predictions = nnet_pred,
        confusion_matrix = nnet_conf,
        importance = NULL,
        scaling = list(center = attr(X_train_scaled, "scaled:center"),
                       scale = attr(X_train_scaled, "scaled:scale")),
        pca_model = pca_model,
        n_components = max_components
      )
      cat("  Neural Network Accuracy:", round(nnet_acc, 4), "(using", max_components, "PCA components)\n")
    }, error = function(e) {
      cat("  Neural Network failed:", conditionMessage(e), "\n")
    })
    
    # ----------------------
    # 5. Elastic Net (glmnet)
    # ----------------------
    cat("\n[5/5] Training Elastic Net...\n")
    tryCatch({
      # Filter classes with sufficient observations (at least 2 per class in training)
      class_counts_train <- table(y_train)
      valid_classes <- names(class_counts_train[class_counts_train >= 2])
      
      if (length(valid_classes) < 2) {
        stop("Not enough classes with sufficient observations for Elastic Net")
      }
      
      # Subset to valid classes only
      train_idx_valid <- y_train %in% valid_classes
      test_idx_valid <- y_test %in% valid_classes
      
      X_train_glm <- X_train[train_idx_valid, , drop = FALSE]
      y_train_glm <- droplevels(y_train[train_idx_valid])
      X_test_glm <- X_test[test_idx_valid, , drop = FALSE]
      y_test_glm <- droplevels(y_test[test_idx_valid])
      
      if (length(levels(y_train_glm)) < 2) {
        stop("Need at least 2 classes for classification")
      }
      
      cat("  Using", length(levels(y_train_glm)), "classes with sufficient observations\n")
      
      # Use cross-validation to find best lambda
      cv_glmnet <- cv.glmnet(as.matrix(X_train_glm), y_train_glm, family = "multinomial", 
                              alpha = 0.5, type.measure = "class", nfolds = min(5, min(class_counts_train[valid_classes])))
      glmnet_pred <- predict(cv_glmnet, as.matrix(X_test_glm), s = "lambda.min", type = "class")
      glmnet_pred <- factor(glmnet_pred[,1], levels = levels(y_test_glm))
      glmnet_conf <- confusionMatrix(glmnet_pred, y_test_glm)
      glmnet_acc <- glmnet_conf$overall["Accuracy"]
      accuracy_scores["ElasticNet"] <- glmnet_acc
      
      # Get coefficients for importance
      coef_list <- coef(cv_glmnet, s = "lambda.min")
      model_results_list[["ElasticNet"]] <- list(
        model = cv_glmnet,
        predictions = glmnet_pred,
        confusion_matrix = glmnet_conf,
        importance = coef_list
      )
      cat("  Elastic Net Accuracy:", round(glmnet_acc, 4), "\n")
    }, error = function(e) {
      cat("  Elastic Net failed:", conditionMessage(e), "\n")
    })
    
    # ----------------------
    # Select Best Model
    # ----------------------
    cat("\n========== Model Comparison Summary ==========\n")
    accuracy_df <- data.frame(
      Model = names(accuracy_scores),
      Accuracy = as.numeric(accuracy_scores)
    )
    accuracy_df <- accuracy_df[order(accuracy_df$Accuracy, decreasing = TRUE), ]
    print(accuracy_df)
    
    best_model_name <- accuracy_df$Model[1]
    best_accuracy <- accuracy_df$Accuracy[1]
    
    cat("\n*** Best Model:", best_model_name, "with Accuracy:", round(best_accuracy, 4), "***\n")
    
    best_model_result <- model_results_list[[best_model_name]]
    
    return(list(
      best_model_name = best_model_name,
      model = best_model_result$model,
      predictions = best_model_result$predictions,
      actual = y_test,
      confusion_matrix = best_model_result$confusion_matrix,
      importance = best_model_result$importance,
      feature_matrix = X,
      all_models = model_results_list,
      accuracy_comparison = accuracy_df
    ))
    
  }, error = function(e) {
    stop("Error in model training: ", conditionMessage(e))
  })
}

# Main execution function
main <- function(counts_path, meta_path,
                 output_dir = dirname(counts_path),
                 viz_dir = file.path(output_dir, "visualizations"),
                 run_viz = FALSE) {
  results_dir <- file.path(output_dir, "analysis_results")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. Load data
  cat("Loading and preparing data...\n")
  data <- prepare_data(counts_path, meta_path)
  
  # 2. Run DESeq2
  cat("Performing DESeq2 analysis...\n")
  deseq_results <- run_deseq2(data$counts, data$meta_data)
  
  # Save DESeq2 results
  for(cancer in names(deseq_results$results)) {
    write.csv(deseq_results$results[[cancer]], 
              file.path(output_dir, paste0("deseq2_", make.names(cancer), ".csv")))
  }
  
  # 3. GO Analysis
  cat("Performing GO analysis...\n")
  go_results <- perform_go_analysis(deseq_results$results)
  
  # 4. GSVA Analysis
  cat("Performing GSVA analysis...\n")
  gsva_results <- perform_gsva_analysis(deseq_results$dds, go_results)
  
  # 5. Train Model (Hierarchical: Type + Grade)
  cat("Training prediction models (Hierarchical: Type + Grade)...\n")
  model_results <- train_model(gsva_results, data$meta_data)
  
  # Save results and create visualizations
  saveRDS(list(
    deseq_results = deseq_results,
    go_results = go_results,
    gsva_results = gsva_results,
    model_results = model_results
  ), file.path(results_dir, "complete_analysis.rds"))

  # Save all accuracy comparisons
  write.csv(model_results$accuracy_comparison,
            file.path(results_dir, "model_accuracy.csv"),
            row.names = FALSE)
  
  # Save type-specific accuracy
  write.csv(model_results$type_accuracy,
            file.path(results_dir, "model_accuracy_TYPE.csv"),
            row.names = FALSE)
  
  # Save grade-specific accuracy
  write.csv(model_results$grade_accuracy,
            file.path(results_dir, "model_accuracy_GRADE.csv"),
            row.names = FALSE)
  
  # Save final comparison
  write.csv(model_results$final_comparison,
            file.path(results_dir, "hierarchical_comparison.csv"),
            row.names = FALSE)
  
  writeLines(paste0("BestModel=", model_results$best_model_name),
             file.path(results_dir, "best_model.txt"))
  
  
  # Feature Selection
  all_scores <- do.call(rbind, lapply(gsva_results, function(x) {
    if (!is.null(x)) exprs(x) else NULL
  }))

  if (!is.null(all_scores) && nrow(all_scores) > 0) {
    all_scores <- all_scores[!duplicated(rownames(all_scores)), , drop = FALSE]
    X_all <- t(all_scores)

    vars_before <- apply(X_all, 2, var)
    top_n <- min(200, length(vars_before))
    top_features <- names(sort(vars_before, decreasing = TRUE))[1:top_n]
    X_top <- X_all[, top_features, drop = FALSE]

    pdf(file.path(viz_dir, "FeatureSelection_variance_distribution.pdf"))
    hist(vars_before, breaks = 50, col = "steelblue",
         main = "Variance distribution of all GSVA pathways",
         xlab = "Variance")
    if (length(vars_before) >= top_n) {
      abline(v = sort(vars_before, decreasing = TRUE)[top_n],
             col = "red", lwd = 2)
    }
    dev.off()

    pdf(file.path(viz_dir, "FeatureSelection_top200_pathways.pdf"))
    barplot(
      sort(vars_before, decreasing = TRUE)[1:top_n],
      main = paste0("Top ", top_n, " most variable pathways"),
      ylab = "Variance",
      border = NA,
      col = colorRampPalette(c("navy", "white", "firebrick3"))(top_n)
    )
    dev.off()

    pdf(file.path(viz_dir, "FeatureSelection_heatmap_selected_features.pdf"))
    pheatmap(
      X_top,
      show_rownames = FALSE,
      show_colnames = FALSE,
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      main = paste0("Heatmap of selected GSVA pathways (top ", top_n, ")"),
      border_color = NA
    )
    dev.off()
  } else {
    warning("No GSVA scores available for feature selection plots.")
  }
  
  
  
  # Create summary plots
  pdf(file.path(results_dir, "analysis_summary.pdf"))
  
  # Plot confusion matrix
  plot(model_results$confusion_matrix$table, 
       main = paste("Confusion Matrix -", model_results$best_model_name),
       xlab = "Predicted",
       ylab = "Actual")
  
  # Plot model comparison
  par(mar = c(8, 4, 4, 2))
  barplot(model_results$accuracy_comparison$Accuracy,
          names.arg = model_results$accuracy_comparison$Model,
          main = "Model Accuracy Comparison",
          ylab = "Accuracy",
          col = ifelse(model_results$accuracy_comparison$Model == model_results$best_model_name, 
                       "darkgreen", "steelblue"),
          las = 2,
          ylim = c(0, 1))
  abline(h = max(model_results$accuracy_comparison$Accuracy), col = "red", lty = 2)

  # Plot top features importance (if available)
  if (!is.null(model_results$importance)) {
    par(mar = c(10, 4, 4, 2))
    
    if (model_results$best_model_name == "RandomForest") {
      imp_df <- as.data.frame(model_results$importance)
      imp_df$feature <- rownames(imp_df)
      imp_df <- imp_df[order(imp_df$MeanDecreaseAccuracy, decreasing = TRUE), ]
      barplot(head(imp_df$MeanDecreaseAccuracy, 20),
              names.arg = head(imp_df$feature, 20),
              main = "Top 20 Important Features (Random Forest)",
              las = 2)
    } else if (model_results$best_model_name == "XGBoost") {
      imp_df <- as.data.frame(model_results$importance)
      imp_df <- imp_df[order(imp_df$Gain, decreasing = TRUE), ]
      barplot(head(imp_df$Gain, 20),
              names.arg = head(imp_df$Feature, 20),
              main = "Top 20 Important Features (XGBoost Gain)",
              las = 2)
    }
  }
  
  dev.off()

  if (isTRUE(run_viz)) {
    run_visualizations(file.path(results_dir, "complete_analysis.rds"), viz_dir)
  }

  return(model_results)
}

############################################################
## NEW: Feature selection helper for machine learning
############################################################

feature_select_gsva <- function(gsva_results, meta_data,
                                top_n = 200, remove_nzv = TRUE) {
  # 1) Combine GSVA scores as in train_model()
  all_scores <- do.call(rbind, lapply(gsva_results, function(x) {
    if (!is.null(x)) exprs(x) else NULL
  }))
  
  if (is.null(all_scores) || nrow(all_scores) == 0) {
    stop("No GSVA scores found for feature selection.")
  }
  
  # 2) Remove duplicated pathway names
  all_scores <- all_scores[!duplicated(rownames(all_scores)), ]
  
  # 3) Transpose to samples × pathways
  X <- t(all_scores)
  
  # 4) Align samples with metadata
  common_samples <- intersect(rownames(X), meta_data$ID_REF)
  if (length(common_samples) == 0) {
    stop("No matching samples between GSVA matrix and metadata in feature_select_gsva().")
  }
  
  X <- X[common_samples, , drop = FALSE]
  
  # 5) Optional: remove near-zero variance features
  if (remove_nzv) {
    nzv_idx <- caret::nearZeroVar(X)
    if (length(nzv_idx) > 0) {
      X <- X[, -nzv_idx, drop = FALSE]
    }
  }
  
  # 6) Keep top_n most variable features
  if (ncol(X) > top_n) {
    vars <- apply(X, 2, var)
    top_features <- names(sort(vars, decreasing = TRUE))[1:top_n]
    X <- X[, top_features, drop = FALSE]
  }
  
  # 7) Return feature-selected matrix and info
  return(list(
    X = X,
    samples = rownames(X),
    features = colnames(X)
  ))
}


parse_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--counts") {
      parsed$counts <- args[i + 1]
      i <- i + 2
      next
    }
    if (key == "--meta") {
      parsed$meta <- args[i + 1]
      i <- i + 2
      next
    }
    if (key == "--out") {
      parsed$out <- args[i + 1]
      i <- i + 2
      next
    }
    if (key == "--viz") {
      parsed$viz <- TRUE
      i <- i + 1
      next
    }
    stop("Unknown argument: ", key)
  }

  if (is.null(parsed$counts) || is.null(parsed$meta)) {
    stop("Usage: Rscript glioma_analysis.R --counts <counts.csv> --meta <meta.csv> [--out <output_dir>] [--viz]")
  }
  if (is.null(parsed$out)) {
    parsed$out <- dirname(parsed$counts)
  }
  if (is.null(parsed$viz)) {
    parsed$viz <- FALSE
  }

  return(parsed)
}

run_cli <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  tryCatch({
    results <- main(
      counts_path = args$counts,
      meta_path = args$meta,
      output_dir = args$out,
      run_viz = args$viz
    )

    cat("\n##########################################################\n")
    cat("# ANALYSIS COMPLETE - HIERARCHICAL CLASSIFICATION        #\n")
    cat("##########################################################\n")
    
    cat("\n========== HIERARCHICAL MODEL SUMMARY ==========\n")
    cat("\n--- Stage 1: TYPE Classification ---\n")
    cat("Best TYPE Model:", results$type_classifier$best_model_name, "\n")
    cat("TYPE Accuracy:", round(results$type_classifier$accuracy_comparison$Accuracy[1], 4), "\n")
    
    cat("\n--- Stage 2: GRADE Classification ---\n")
    cat("Best GRADE Model:", results$grade_classifier$best_model_name, "\n")
    cat("GRADE Accuracy:", round(results$grade_classifier$accuracy_comparison$Accuracy[1], 4), "\n")
    
    cat("\n--- Combined Hierarchical Approach ---\n")
    cat("Combined Accuracy:", round(results$combined_accuracy, 4), "\n")
    
    cat("\n--- Direct Full Classification (for comparison) ---\n")
    cat("Best Direct Model:", results$direct_classifier$best_model_name, "\n")
    cat("Direct Accuracy:", round(results$direct_classifier$accuracy_comparison$Accuracy[1], 4), "\n")
    
    cat("\n========== FINAL COMPARISON ==========\n")
    print(results$final_comparison)
    
    cat("\n========== ALL MODELS RANKED BY ACCURACY ==========\n")
    print(results$accuracy_comparison)
    
    cat("\n========== TYPE CLASSIFICATION DETAILS ==========\n")
    print(results$type_classifier$best_model$confusion_matrix)
    
    cat("\n========== GRADE CLASSIFICATION DETAILS ==========\n")
    print(results$grade_classifier$best_model$confusion_matrix)
    
  }, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
    traceback()
  })
}

if (!interactive()) {
  run_cli()
}


run_visualizations <- function(analysis_rds_path, viz_dir) {
  ############################################################
  ## VISUALIZATIONS FOR GLIOMA PROJECT
  ## Run this AFTER your main pipeline has finished
  ############################################################

## 1) Setup ------------------------------------------------

  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

# Libraries (most are already installed from your main script)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(Biobase)

  # Load the analysis results from the provided RDS file
  analysis <- readRDS(analysis_rds_path)

# Unpack
deseq_results <- analysis$deseq_results   # list: $dds, $results (DE tables)
go_results    <- analysis$go_results      # named list of enrichResult
gsva_results  <- analysis$gsva_results    # named list of ExpressionSet (GSVA)
model_results <- analysis$model_results   # best model + all_models etc.

dds       <- deseq_results$dds
de_tables <- deseq_results$results        # list of DESeq2 result data.frames per cancer_type
meta_data <- as.data.frame(colData(dds))


############################################################
## AUTOMATIC GRADE 2 GLIOMA VISUALIZATIONS
############################################################

# Your dataset contains 3 Grade-II glioma groups:
grade2_levels <- c(
  "astrocytomas",
  "oligodendrogliomas",
  "oligodendroastrocytomas"
)

# Keep only the Grade-II labels that actually exist
grade2_levels <- grade2_levels[grade2_levels %in% levels(meta_data$cancer_type)]

cat("Detected Grade 2 subtypes:\n")
print(grade2_levels)

for (g_label in grade2_levels) {
  
  if (!g_label %in% names(de_tables)) {
    cat("Skipping", g_label, "(no DESeq2 table found)\n")
    next
  }
  
  cat("\n=== Making Grade 2 plots for:", g_label, "===\n")
  
  res_g <- de_tables[[g_label]]
  
  # Mark significant genes
  res_g$significant <- with(
    res_g,
    !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1
  )
  
  # -------------------------------------------
  # 1) Histogram of log2FC for this Grade-II subtype
  # -------------------------------------------
  pdf(file.path(viz_dir, paste0("Grade2_", g_label, "_log2FC_histogram.pdf")))
  print(
    ggplot(res_g, aes(x = log2FoldChange, fill = significant)) +
      geom_histogram(bins = 60, alpha = 0.7) +
      scale_fill_manual(values = c("grey70", "firebrick3")) +
      labs(
        title = paste("log2FC distribution for Grade-II subtype:", g_label),
        x = "log2 fold change",
        y = "Number of genes",
        fill = "Significant"
      ) +
      theme_minimal()
  )
  dev.off()
  
  # -------------------------------------------
  # 2) Density plots of top DE genes for this Grade-II subtype
  # -------------------------------------------
  
  norm_counts <- counts(dds, normalized = TRUE)
  
  # take top genes sorted by adjusted p-value
  top_genes_g <- res_g[order(res_g$padj), "GeneID"]
  top_genes_g <- head(top_genes_g[!is.na(top_genes_g)], 5)
  
  expr_df <- as.data.frame(t(norm_counts[top_genes_g, , drop = FALSE]))
  expr_df$sample <- rownames(expr_df)
  
  expr_df <- merge(
    expr_df,
    meta_data[, c("ID_REF", "cancer_type")],
    by.x = "sample",
    by.y = "ID_REF"
  )
  
  expr_long <- tidyr::pivot_longer(
    expr_df,
    cols = all_of(top_genes_g),
    names_to = "gene",
    values_to = "expr"
  )
  
  expr_long$group <- ifelse(
    expr_long$cancer_type == g_label,
    g_label,
    "Other"
  )
  
  pdf(file.path(viz_dir, paste0("Grade2_", g_label, "_top_genes_density.pdf")))
  print(
    ggplot(expr_long, aes(x = log2(expr + 1), color = group)) +
      geom_density() +
      facet_wrap(~ gene, scales = "free") +
      labs(
        title = paste("Top DE genes for Grade-II subtype:", g_label),
        x = "log2(normalized counts + 1)",
        y = "Density",
        color = "Group"
      ) +
      theme_minimal()
  )
  dev.off()
}


############################################################
## 2) DESeq2: Volcano plot (genes) ------------------------
############################################################

plot_volcano <- function(deseq_res_df, cancer_type,
                         lfc_thresh = 1, padj_thresh = 0.05, top_n = 10) {
  df <- deseq_res_df
  df$negLog10Padj <- -log10(df$padj)
  
  df$significant <- with(
    df,
    !is.na(padj) & padj < padj_thresh & abs(log2FoldChange) > lfc_thresh
  )
  
  # Top significant genes for labeling
  top_genes <- df[df$significant, ]
  top_genes <- top_genes[order(top_genes$padj), ]
  top_genes <- head(top_genes, top_n)
  
  ggplot(df, aes(x = log2FoldChange, y = negLog10Padj)) +
    geom_point(aes(color = significant), alpha = 0.6) +
    scale_color_manual(values = c("grey70", "red")) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
    geom_text(
      data = top_genes,
      aes(label = GeneID),
      vjust = 1.2, size = 3
    ) +
    labs(
      title = paste("Volcano plot -", cancer_type),
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = "Significant"
    ) +
    theme_minimal()
}

# Save ALL DESeq2 volcano plots to one PDF
pdf(file.path(viz_dir, "DESeq2_volcano_plots.pdf"))
for (cancer in names(de_tables)) {
  p <- plot_volcano(de_tables[[cancer]], cancer_type = cancer)
  print(p)
}
dev.off()



############################################################
## PCA PLOT WITH CUSTOM CANCER-TYPE COLORS + PDF OUTPUT
############################################################

# 0) Build color map for the cancer types that actually exist
ct_levels <- sort(unique(meta_data$cancer_type))

preferred_colors <- c(
  # Grade IV (GBM)
  "primary Glioblastomas"              = "#E74C3C",
  "secondary Glioblastomas"            = "#943126",
  
  # Grade III (anaplastic)
  "anaplastic astrocytomas"            = "#F39C12",
  "anaplastic oligodendrogliomas"      = "#D68910",
  "anaplastic oligodendroastrocytomas" = "#B9770E",
  
  # Grade II
  "astrocytomas"                       = "#2ECC71",
  "oligodendrogliomas"                 = "#27AE60",
  "oligodendroastrocytomas"            = "#145A32",
  
  # Recurrent groups (blue palette)
  "recurrent astrocytomas"                     = "#5DADE2",
  "recurrent oligodendrogliomas"               = "#3498DB",
  "recurrent oligodendroastrocytomas"          = "#2E86C1",
  "recurrent anaplastic astrocytomas"          = "#1F618D",
  "recurrent anaplastic oligodendrogliomas"    = "#154360",
  "recurrent Glioblastomas"                    = "#1B4F72",
  "recurrent anaplastic oligodendroastrocytomas" = "#0B5345"
)

# keep only colors for levels that are actually present
base_colors <- preferred_colors[names(preferred_colors) %in% ct_levels]

# if any levels are missing colors, assign grey→black automatically
missing_levels <- setdiff(ct_levels, names(base_colors))
if (length(missing_levels) > 0) {
  extra_cols <- colorRampPalette(c("grey70", "black"))(length(missing_levels))
  names(extra_cols) <- missing_levels
  base_colors <- c(base_colors, extra_cols)
}

# 1) Get PCA data from DESeq2 (VST object 'vsd' already created earlier)
pca_df <- DESeq2::plotPCA(vsd, intgroup = "cancer_type", returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

# 2) Make PCA plot with manual colors
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cancer_type)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = base_colors, name = "cancer_type") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 8)
  )

# 3) Save PCA to PDF
pca_pdf <- file.path(viz_dir, "DESeq2_PCA_custom_colors.pdf")
pdf(pca_pdf, width = 8, height = 6)
print(p_pca)
dev.off()

############################################################
## END PCA BLOCK
############################################################

############################################################
## 3) DESeq2: Heatmap of top variable genes --------------
############################################################
############################################################
## CUSTOM HEATMAP WITH AUTOMATIC ANNOTATION COLORS
############################################################

# 0) Get all cancer_type levels actually present
ct_levels <- sort(unique(annotation_col$cancer_type))

# 1) Start with a "preferred" color map for known labels
preferred_colors <- c(
  # Grade IV (GBM)
  "primary Glioblastomas"              = "#E74C3C",
  "secondary Glioblastomas"            = "#943126",
  
  # Grade III (anaplastic)
  "anaplastic astrocytomas"            = "#F39C12",
  "anaplastic oligodendrogliomas"      = "#D68910",
  "anaplastic oligodendroastrocytomas" = "#B9770E",
  
  # Grade II
  "astrocytomas"                       = "#2ECC71",
  "oligodendrogliomas"                 = "#27AE60",
  "oligodendroastrocytomas"            = "#145A32",
  
  # Recurrent groups (blue palette)
  "recurrent astrocytomas"                     = "#5DADE2",
  "recurrent oligodendrogliomas"               = "#3498DB",
  "recurrent oligodendroastrocytomas"          = "#2E86C1",
  "recurrent anaplastic astrocytomas"          = "#1F618D",
  "recurrent anaplastic oligodendrogliomas"    = "#154360",
  "recurrent Glioblastomas"                    = "#1B4F72",
  "recurrent anaplastic oligodendroastrocytomas" = "#0B5345"
)

# 2) Keep ONLY colors that correspond to actually present levels
base_colors <- preferred_colors[names(preferred_colors) %in% ct_levels]

# 3) If there are any cancer_type levels without a color, assign them automatically
missing_levels <- setdiff(ct_levels, names(base_colors))

if (length(missing_levels) > 0) {
  extra_cols <- colorRampPalette(c("grey70", "black"))(length(missing_levels))
  names(extra_cols) <- missing_levels
  base_colors <- c(base_colors, extra_cols)
}

# 4) Build annotation_colors list for pheatmap
annotation_colors <- list(
  cancer_type = base_colors
)

# 5) Define BLUE → WHITE → RED heatmap colors
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# 6) Save to PDF
pdf(
  file.path(viz_dir, "DESeq2_top_genes_heatmap_BLUE_RED.pdf"),
  width = 11,
  height = 11
)

pheatmap(
  mat_top,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = heatmap_colors,
  main = paste("Top", n_genes, "Most Variable Genes (VST)"),
  fontsize = 10,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete"
)

dev.off()
############################################################
## END HEATMAP BLOCK
############################################################
#



############################################################
## 4) GO: Volcano-style plots & heatmap -------------------
############################################################

# For GO, we can treat enrichment like this:
#   x-axis  = GeneRatio (converted to numeric)
#   y-axis  = -log10(p.adjust)
# to get a "GO volcano" sense: strong + significant terms on the top-right.
############################################################
## Better GO visualization: dot plot
############################################################

go_dotplot_better <- function(go_res, cancer_type,
                              top_n = 15) {
  if (is.null(go_res) || nrow(go_res@result) == 0) {
    return(NULL)
  }
  
  df <- as.data.frame(go_res@result)
  
  # Convert "GeneRatio" like "10/100" to numeric
  df$GeneRatioNum <- sapply(df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
  
  # Significance as -log10(p.adjust)
  df$negLog10Padj <- -log10(df$p.adjust)
  
  # Keep top_n most significant terms
  df <- df[order(df$p.adjust), ]
  df_top <- head(df, top_n)
  
  # Order terms by significance (or GeneRatio, if you prefer)
  df_top$Description <- factor(
    df_top$Description,
    levels = rev(df_top$Description)   # so most significant on top
  )
  
  ggplot(df_top,
         aes(x = GeneRatioNum,
             y = Description)) +
    geom_point(aes(size = Count,
                   color = negLog10Padj)) +
    scale_color_gradient(
      low  = "steelblue",
      high = "firebrick3",
      name = "-log10 adj p"
    ) +
    scale_size_continuous(
      range = c(2, 8),
      name  = "Gene count"
    ) +
    labs(
      title = paste("GO dot plot -", cancer_type),
      x = "Gene ratio (DEGs in GO term / DEGs total)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title  = element_text(size = 12, face = "bold")
    )
}

############################################################
## GO dot plots (nicer visualization than volcano)
############################################################

go_dot_pdf <- file.path(viz_dir, "GO_dotplots_all_cancertypes.pdf")
pdf(go_dot_pdf, width = 9, height = 7)

for (ct in names(go_results)) {
  p <- go_dotplot_better(go_results[[ct]], cancer_type = ct, top_n = 15)
  if (!is.null(p)) print(p)
}

dev.off()



# ############################################################
## GSVA PATHWAY HEATMAP – MATCH DESeq2 STYLE
############################################################

# mat_paths and annotation_col_paths already created above:
#   mat_paths:   top pathways × samples
#   annotation_col_paths$cancer_type: subtype for each sample

# 1) Build consistent colors for cancer_type (same as DESeq2)
ct_levels <- sort(unique(annotation_col_paths$cancer_type))

preferred_colors <- c(
  # Grade IV (GBM)
  "primary Glioblastomas"              = "#E74C3C",
  "secondary Glioblastomas"            = "#943126",
  
  # Grade III (anaplastic)
  "anaplastic astrocytomas"            = "#F39C12",
  "anaplastic oligodendrogliomas"      = "#D68910",
  "anaplastic oligodendroastrocytomas" = "#B9770E",
  
  # Grade II
  "astrocytomas"                       = "#2ECC71",
  "oligodendrogliomas"                 = "#27AE60",
  "oligodendroastrocytomas"            = "#145A32",
  
  # Recurrent groups (blue palette)
  "recurrent astrocytomas"                     = "#5DADE2",
  "recurrent oligodendrogliomas"               = "#3498DB",
  "recurrent oligodendroastrocytomas"          = "#2E86C1",
  "recurrent anaplastic astrocytomas"          = "#1F618D",
  "recurrent anaplastic oligodendrogliomas"    = "#154360",
  "recurrent Glioblastomas"                    = "#1B4F72",
  "recurrent anaplastic oligodendroastrocytomas" = "#0B5345"
)

# keep only levels that exist in this matrix
base_colors <- preferred_colors[names(preferred_colors) %in% ct_levels]

# if any levels have no color yet, assign grey→black automatically
missing_levels <- setdiff(ct_levels, names(base_colors))
if (length(missing_levels) > 0) {
  extra_cols <- colorRampPalette(c("grey70", "black"))(length(missing_levels))
  names(extra_cols) <- missing_levels
  base_colors <- c(base_colors, extra_cols)
}

annotation_colors_gsva <- list(
  cancer_type = base_colors
)

# 2) Use SAME blue→white→red palette as DESeq2 heatmap
gsva_heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# 3) Save GSVA heatmap to PDF
gsva_pdf2 <- file.path(viz_dir, "GSVA_pathway_heatmap_BLUE_RED.pdf")
pdf(gsva_pdf2, width = 11, height = 11)

pheatmap(
  mat_paths,
  annotation_col   = annotation_col_paths,
  annotation_colors = annotation_colors_gsva,
  color            = gsva_heatmap_colors,
  main             = paste("Top", n_pathways, "most variable pathways (GSVA)"),
  show_rownames    = FALSE,
  show_colnames    = FALSE,
  border_color     = NA,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete"
)

dev.off()

############################################################
## END GSVA HEATMAP BLOCK
############################################################


############################################################
## PIPELINE OVERVIEW VISUALIZATIONS
############################################################

library(dplyr)
library(ggplot2)

pipeline_pdf <- file.path(viz_dir, "pipeline_overview.pdf")
pdf(pipeline_pdf, width = 10, height = 8)

############################################################
## 1) Number of significant DEGs per cancer_type
############################################################

library(ggplot2)

# Build summary of significant DEGs per cancer_type
deg_list <- lapply(names(de_tables), function(ct) {
  df <- de_tables[[ct]]
  
  # logical index of significant genes
  sig_idx <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1
  
  data.frame(
    cancer_type = ct,
    total_genes = nrow(df),
    sig_genes   = sum(sig_idx)
  )
})

deg_summary <- do.call(rbind, deg_list)

# Create a DEG summary PDF only
deg_pdf <- file.path(viz_dir, "pipeline_DEG_summary.pdf")
pdf(deg_pdf, width = 8, height = 6)

p_deg <- ggplot(deg_summary, aes(x = cancer_type, y = sig_genes)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Number of significant DE genes per cancer type",
    x = "Cancer type",
    y = "# of DE genes (padj < 0.05 & |log2FC| > 1)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_deg)

dev.off()


############################################################
## 2) Number of enriched GO terms per cancer_type
############################################################
# Make sure ggplot2 is loaded
library(ggplot2)

# Build GO summary without dplyr
go_list <- lapply(names(go_results), function(ct) {
  gr <- go_results[[ct]]
  
  if (is.null(gr)) {
    data.frame(cancer_type = ct, go_terms = 0)
  } else {
    df <- as.data.frame(gr@result)
    
    # significant GO terms: padj < 0.05
    sig_idx <- !is.na(df$p.adjust) & df$p.adjust < 0.05
    
    data.frame(
      cancer_type = ct,
      go_terms    = sum(sig_idx)
    )
  }
})

go_summary <- do.call(rbind, go_list)

# Create PDF in your Visualizations folder
go_pdf <- file.path(viz_dir, "pipeline_GO_summary.pdf")
pdf(go_pdf, width = 8, height = 6)

p_go <- ggplot(go_summary, aes(x = cancer_type, y = go_terms)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(
    title = "Number of enriched GO BP terms per cancer type",
    x = "Cancer type",
    y = "# of enriched GO terms (padj < 0.05)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_go)

dev.off()


############################################################
## 3) Number of GSVA pathways per cancer_type
############################################################

# Build GSVA summary without dplyr
gsva_list <- lapply(names(gsva_results), function(ct) {
  gs <- gsva_results[[ct]]
  
  if (is.null(gs)) {
    data.frame(cancer_type = ct, pathways = 0)
  } else {
    # Number of pathways = number of GSVA gene sets
    data.frame(
      cancer_type = ct,
      pathways = nrow(exprs(gs))
    )
  }
})

gsva_summary <- do.call(rbind, gsva_list)

# Create PDF
gsva_pdf <- file.path(viz_dir, "pipeline_GSVA_summary.pdf")
pdf(gsva_pdf, width = 8, height = 6)

p_gsva <- ggplot(gsva_summary, aes(x = cancer_type, y = pathways)) +
  geom_bar(stat = "identity", fill = "seagreen") +
  labs(
    title = "Number of GSVA pathways derived per cancer type",
    x = "Cancer type",
    y = "# of GSVA gene sets (pathways)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_gsva)

dev.off()



############################################################
## Faceted histograms: compare shapes across cancer_types
############################################################

# Combine all DE tables into one big data.frame
de_all_list <- lapply(names(de_tables), function(ct) {
  df <- de_tables[[ct]]
  df$cancer_type <- ct
  df
})
de_all <- do.call(rbind, de_all_list)
de_all <- de_all[!is.na(de_all$log2FoldChange), ]

logfc_facet_pdf <- file.path(viz_dir, "DESeq2_log2FC_histograms_faceted.pdf")
pdf(logfc_facet_pdf, width = 11, height = 8)

p_facet <- ggplot(de_all, aes(x = log2FoldChange)) +
  geom_histogram(aes(y = ..density..),
                 bins = 50,
                 fill = "steelblue",
                 alpha = 0.6) +
  geom_density(color = "darkred", size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ cancer_type, scales = "free_y") +
  labs(
    title = "log2FoldChange distributions across cancer types",
    x = "log2FoldChange (vs primary Glioblastomas)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8)
  )

print(p_facet)
dev.off()





############################################################
## 4) Model accuracy comparison (pipeline outcome)
############################################################

# accuracy_comparison is already in model_results from your pipeline
acc_df <- model_results$accuracy_comparison

# Create PDF for model accuracy plot
acc_pdf <- file.path(viz_dir, "pipeline_model_accuracy.pdf")
pdf(acc_pdf, width = 8, height = 6)

p_acc <- ggplot(acc_df, aes(x = Model, y = Accuracy)) +
  geom_bar(
    stat = "identity",
    aes(fill = Model == model_results$best_model_name)
  ) +
  scale_fill_manual(
    values = c("grey70", "darkgreen"),
    guide = FALSE
  ) +
  geom_hline(
    yintercept = max(acc_df$Accuracy),
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    title = "Model accuracies (final pipeline step)",
    x = "Model",
    y = "Accuracy on test set"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_acc)

dev.off()


############################################################
## NEW: nicer colour palette heatmaps
############################################################

# Colour palette: blue → white → red
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# 1) VST top genes heatmap (again, but with custom colors)
pdf(file.path(viz_dir, "DESeq2_top_genes_heatmap_nice_colors.pdf"))
pheatmap(
  mat_top,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  main = paste("Top", n_genes, "most variable genes (VST)"),
  color = heatmap_colors
)
dev.off()

# 2) GSVA pathways heatmap with same palette
pdf(file.path(viz_dir, "GSVA_pathway_heatmap_nice_colors.pdf"))
pheatmap(
  mat_paths,
  annotation_col = annotation_col_paths,
  show_rownames = FALSE,
  main = paste("Top", n_pathways, "most variable pathways (GSVA)"),
  color = heatmap_colors
)
dev.off()


############################################################
## NEW: QC plots for counts & samples
############################################################

qc_pdf <- file.path(viz_dir, "QC_counts_and_PCA.pdf")
pdf(qc_pdf, width = 10, height = 8)

# 1) Library size per sample
lib_sizes <- colSums(counts(dds))
lib_df <- data.frame(
  sample = names(lib_sizes),
  library_size = lib_sizes,
  cancer_type = meta_data$cancer_type[match(names(lib_sizes), rownames(meta_data))]
)

p_lib <- ggplot(lib_df, aes(x = sample, y = library_size, fill = cancer_type)) +
  geom_bar(stat = "identity") +
  labs(title = "Library size per sample", x = "Sample", y = "Total counts") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
print(p_lib)

# 2) Boxplot of log2 normalized counts
norm_counts <- counts(dds, normalized = TRUE)
log_norm <- log2(norm_counts + 1)
log_df <- as.data.frame(log_norm)
log_df$gene <- rownames(log_df)
log_long <- tidyr::pivot_longer(
  log_df,
  cols = -gene,
  names_to = "sample",
  values_to = "log2expr"
)

p_box <- ggplot(log_long, aes(x = sample, y = log2expr)) +
  geom_boxplot(outlier.size = 0.3) +
  labs(title = "Distribution of log2 normalized counts per sample",
       x = "Sample", y = "log2(normalized counts + 1)") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
print(p_box)

# 3) Sample distance heatmap
sample_dist <- dist(t(vsd_mat))
sample_dist_mat <- as.matrix(sample_dist)
rownames(sample_dist_mat) <- colnames(vsd_mat)
colnames(sample_dist_mat) <- colnames(vsd_mat)

pheatmap(
  sample_dist_mat,
  annotation_col = annotation_col,
  main = "Sample-to-sample distances (VST)",
  color = heatmap_colors
)

# 4) PCA (you already have one, but included for completeness)
pca_plot <- plotPCA(vsd, intgroup = "cancer_type")
print(pca_plot)

dev.off()

table(meta_data$cancer_type)
}
