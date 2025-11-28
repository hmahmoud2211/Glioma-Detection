# Load required libraries
required_packages <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "GSVA", 
                       "caret", "dplyr", "tidyr", "xgboost", "randomForest",
                       "e1071", "glmnet", "nnet",
                       "Biobase", "GSEABase", "ggplot2")

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

# 5. Build and Train Multiple Models - Select Best One
train_model <- function(gsva_results, meta_data) {
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
main <- function(counts_path, meta_path) {
  output_dir <- dirname(counts_path)
  
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
  
  # 5. Train Model
  cat("Training prediction model...\n")
  model_results <- train_model(gsva_results, data$meta_data)
  
  # Create results directory
  results_dir <- file.path(output_dir, "analysis_results")
  dir.create(results_dir, showWarnings = FALSE)
  
  # Save results and create visualizations
  saveRDS(list(
    deseq_results = deseq_results,
    go_results = go_results,
    gsva_results = gsva_results,
    model_results = model_results
  ), file.path(results_dir, "complete_analysis.rds"))
  
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
  
  return(model_results)
}

# Run the pipeline
tryCatch({
  counts_path <- "C:\\Users\\user\\Downloads\\Bioinfo\\Need_CSV\\GSE48865_raw_counts_GRCh38.p13_NCBI.csv"
  meta_path <- "C:/Users/user/Downloads/Bioinfo/Need_CSV/meta_data.csv"
  
  results <- main(counts_path, meta_path)
  
  cat("\nAnalysis complete!\n")
  cat("Best Model:", results$best_model_name, "\n")
  cat("Model Accuracy:", results$confusion_matrix$overall["Accuracy"], "\n")
  cat("\nAll Models Comparison:\n")
  print(results$accuracy_comparison)
  cat("\nDetailed Performance Metrics for Best Model:\n")
  print(results$confusion_matrix)
  
}, error = function(e) {
  cat("Error occurred:", conditionMessage(e), "\n")
  print(traceback())
})
