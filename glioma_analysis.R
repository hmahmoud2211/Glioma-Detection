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
  
  
  #Feature Selection
  # Compute variance before filtering
  vars_before <- apply(t(all_scores), 2, var)
  
  pdf(file.path(viz_dir, "FeatureSelection_variance_distribution.pdf"))
  hist(vars_before, breaks = 50, col = "steelblue",
       main = "Variance distribution of all GSVA pathways",
       xlab = "Variance")
  
  abline(v = sort(vars_before, decreasing = TRUE)[200],
         col = "red", lwd = 2)
  dev.off()
  pdf(file.path(viz_dir, "FeatureSelection_top200_pathways.pdf"))
  barplot(
    sort(vars_before, decreasing = TRUE)[1:200],
    main = "Top 200 most variable pathways",
    ylab = "Variance",
    border = NA,
    col = colorRampPalette(c("navy", "white", "firebrick3"))(200)
  )
  dev.off()
  pdf(file.path(viz_dir, "FeatureSelection_heatmap_selected_features.pdf"))
  pheatmap(
    X, 
    show_rownames = FALSE,
    show_colnames = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    main = "Heatmap of selected GSVA pathways (top 200)",
    border_color = NA
  )
  dev.off()
  
  
  
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

  return(model_results)  
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


# Run the pipeline
tryCatch({
  counts_path <- "G:/Glioma Project/Glioma-Detection/Need_CSV/GSE48865_raw_counts_GRCh38.p13_NCBI.csv"
  meta_path <- "G:/Glioma Project/Glioma-Detection/Need_CSV/meta_data.csv"
  
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


############################################################
## VISUALIZATIONS FOR GLIOMA PROJECT
## Run this AFTER your main pipeline has finished
############################################################

## 1) Setup ------------------------------------------------

# This is the ONLY explicit path we'll use (for saving plots)
viz_dir <- "G:\\Glioma Project\\Glioma-Detection\\Visualizations"
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

# Libraries (most are already installed from your main script)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(Biobase)

# Load the analysis results from the correct RDS file
analysis_rds_path <- "G:/Glioma Project/Glioma-Detection/Need_CSV/analysis_results/complete_analysis.rds"
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



