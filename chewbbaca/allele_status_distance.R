# ==================== 配置 ====================
setwd("~/chewbbaca/")
PROFILE_PATH <- "results_alleles.tsv"
OUTPUT_PATH  <- "distance_matrix_all_samples_nie.txt"

# ==================== 1. 读取并清洗数据 ====================
cat("Reading allele profiles...\n")
profiles <- read.delim(PROFILE_PATH, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# 清洗：INF- 去前缀，非数值转 NA，0 转 NA
profiles[] <- lapply(profiles, function(x) {
  x <- as.character(x)
  x <- sub("^INF-", "", x)              # INF-183 -> "183"
  x <- suppressWarnings(as.numeric(x))  # "183" -> 183; "LNF" -> NA
  x[x == 0] <- NA                       # 0 -> NA
  return(x)
})

profile_mat <- as.matrix(profiles)
all_samples <- rownames(profile_mat)
n <- length(all_samples)
m <- ncol(profile_mat)

cat("Samples:", n, " Loci:", m, "\n")

# ==================== 2. 初始化距离矩阵 ====================
dist_mat <- matrix(0L, n, n, dimnames = list(all_samples, all_samples))

# ==================== 3. 分块计算两两距离 ====================
chunk_size <- 1000   # 每次处理 1000 个样本，可根据内存调整
start_time <- Sys.time()

for (i in 1:(n - 1)) {
  test_vec <- profile_mat[i, ]
  remaining <- (i + 1):n
  n_remaining <- length(remaining)
  
  for (chunk_start in seq(1, n_remaining, by = chunk_size)) {
    chunk_end <- min(chunk_start + chunk_size - 1, n_remaining)
    chunk_idx <- remaining[chunk_start:chunk_end]
    
    sub_mat <- profile_mat[chunk_idx, , drop = FALSE]
    
    # 有效位点：双方都非 NA
    # 显式广播 test_vec 为矩阵
    test_vec_mat <- matrix(test_vec, nrow = nrow(sub_mat), ncol = ncol(sub_mat), byrow = TRUE)
    
    valid <- !is.na(sub_mat) & !is.na(test_vec_mat)
    # 差异：在有效位点上等位基因不同
    diff <- (sub_mat != test_vec_mat) & valid
    dists <- rowSums(diff, na.rm = TRUE)
    
    # 填入对称位置
    dist_mat[i, chunk_idx] <- dists
    dist_mat[chunk_idx, i] <- dists
  }
  
  # 进度提示
  if (i %% 100 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat(sprintf("Processed %d / %d samples (%.1f min)\n", i, n - 1, elapsed))
  }
}

cat("Distance matrix computation completed.\n")

# ==================== 4. 保存结果 ====================
write.table(dist_mat, OUTPUT_PATH, sep = "\t", quote = FALSE, col.names = NA)
cat("Distance matrix saved to:", OUTPUT_PATH, "\n")
