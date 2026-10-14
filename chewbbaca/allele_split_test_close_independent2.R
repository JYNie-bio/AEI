args <- commandArgs(trailingOnly = TRUE)

setwd("~/chewbbaca/")
# ==================== 配置 ====================
PROFILE_PATH <- "results_alleles.tsv"
METADATA_PATH <- "~/MIC_20230906_14880_addtraintest.txt"
#OUTPUT_PATH <- "test_near_vs_independent.txt"
THRESHOLD <- 3  # cgMLST 等位基因差异阈值，≤3 为近缘组

# ==================== 1. 读取等位基因谱 ====================
profiles <- read.delim(PROFILE_PATH, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
# 将非数值（LNF, PLOT5, INF- 等）转换为 NA，0 也视为缺失
profiles[] <- lapply(profiles, function(x) {
  x <- as.character(x)
  x <- sub("^INF-", "", x)              # INF-183 -> "183"
  x <- suppressWarnings(as.numeric(x))  # "183" -> 183; "LNF" -> NA
  x[x == 0] <- NA                       # 0 -> NA
  return(x)
})
#profiles[profiles == 0] <- NA
profile_mat <- as.matrix(profiles)

# ==================== 2. 读取元数据 ====================
meta <- read.delim(METADATA_PATH, check.names = FALSE, stringsAsFactors = FALSE)
id_col <- names(meta)[1]
antibiotic_cols <- grep("_traintest$", names(meta), value = TRUE)

# ==================== 3. 对每种抗生素分别分析 ====================
all_results <- list()

ab_col<-antibiotic_cols[as.numeric(args[1])]
  ab_name <- sub("_traintest$", "", ab_col)
OUTPUT_PATH <- paste0("test_near_vs_independent_",ab_name,".txt")

  train_ids <- meta$Genome_name[which(meta[[ab_col]] == "train")]
  test_ids  <- meta$Genome_name[which((meta[[ab_col]] == "test") | (meta[[ab_col]] == "uniquetest"))]
  
  # 过滤存在于等位基因谱中的样本
  train_ids <- intersect(train_ids, rownames(profile_mat))
  test_ids  <- intersect(test_ids, rownames(profile_mat))
  
  if (length(train_ids) == 0 || length(test_ids) == 0) next
  
  # 提取子矩阵
  train_mat <- profile_mat[train_ids, , drop = FALSE]
  test_mat  <- profile_mat[test_ids, , drop = FALSE]
  
  # 预计算训练矩阵的非 NA 掩码，加速后续循环
  train_not_na <- !is.na(train_mat)
  
  # 对每个测试样本计算到训练集的最小距离
  for (i in seq_along(test_ids)) {
    test_id  <- test_ids[i]
    test_vec <- test_mat[i, ]
    
    # 有效位点：训练样本非 NA 且测试样本非 NA
    #chewBBACA 的等位基因编号从 1 开始。0 通常代表缺失数据，而不是一个真实的等位基因。0 意味着：该位点没有被成功检出，或者被判定为不可靠（如 contig 末端、旁系同源、组装错误等）。
    # 广播 test_vec 为矩阵
    test_vec_mat <- matrix(test_vec, nrow = nrow(train_mat), ncol = ncol(train_mat), byrow = TRUE)
    
    # 有效位点：双方都非 NA
    valid <- train_not_na & !is.na(test_vec_mat)
    
    # 差异：仅在有效位点处比较
    diff <- (train_mat != test_vec_mat) & valid
    distances <- rowSums(diff, na.rm = TRUE)
    valid_counts <- rowSums(valid, na.rm = TRUE)
    
    # 只保留有效位点 > 0 的训练样本
    distances[valid_counts == 0] <- NA
    
    if (all(is.na(distances))) {
      min_dist <- NA
      closest_train <- NA
      group <- "no_valid_comparison"
    } else {
      min_dist <- min(distances, na.rm = TRUE)
      # closest_idx <- which.min(distances)
      # closest_train <- train_ids[closest_idx]
      closest_indices <- which(distances == min_dist)
      closest_train <- paste(train_ids[closest_indices], collapse = ";")
      group <- ifelse(min_dist <= THRESHOLD, "near", "independent")
    }
    
    all_results[[length(all_results) + 1]] <- data.frame(
      antibiotic = ab_name,
      test_id = test_id,
      min_distance_to_train = min_dist,
      closest_train_isolate = closest_train,
      group = group,
      stringsAsFactors = FALSE
    )
  }
  cat("  Done:", length(test_ids), "test samples processed\n")


# ==================== 4. 合并并保存结果 ====================
results_df <- do.call(rbind, all_results)
write.table(results_df, OUTPUT_PATH, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nResults saved to:", OUTPUT_PATH, "\n")

# 打印摘要
cat("\nSummary by antibiotic and group:\n")
print(table(results_df$antibiotic, results_df$group))
