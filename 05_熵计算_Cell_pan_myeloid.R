# 加载必要的库
library(dbscan)  # 用于 k-NN 计算
library(dplyr)   # 数据操作
library(ggplot2) # 数据可视化

# 创建示例数据
n_cells <- 5000  # 假设有 5000 个细胞
n_clusters <- 5  # 假设有 5 个簇
n_cancer_types <- 3  # 假设有 3 种癌症类型
set.seed(123)  # 设置随机种子以确保结果可复现

cell_data <- data.frame(
  UMAP_1 = rnorm(n_cells, mean = 0, sd = 1),  # 随机生成UMAP坐标
  UMAP_2 = rnorm(n_cells, mean = 0, sd = 1),
  cluster = sample(1:n_clusters, n_cells, replace = TRUE),  # 随机分配到5个簇
  cancer_type = sample(c("Cancer_A", "Cancer_B", "Cancer_C"), n_cells, replace = TRUE)  # 随机分配癌症类型
)

# 1. 随机采样 1000 个细胞（从 5000 个细胞中）
sampled_cells <- cell_data %>% sample_n(1000)

# 2. 构建 k-NN 图（k = 80）
# 确保 k 小于数据的行数（1000 行数据）
k <- 80  # 最近邻个数
knn_result <- kNN(as.matrix(sampled_cells[, c("UMAP_1", "UMAP_2")]), k = k)

# 检查 knn_result 的结果是否正确
# knn_result$id 是一个矩阵，行数等于采样的细胞数，列数等于 k
dim(knn_result$id)  # 应该是 1000 x 80

# 3. 计算每个细胞的邻居簇比例
# 获取邻居的簇信息（注意列数为 k，行数为细胞数）
neighbor_clusters <- apply(knn_result$id, 1, function(neighbors) {
  sampled_cells$cluster[neighbors]  # 提取邻居的 cluster 信息
})

# 将邻居信息转置为矩阵
neighbor_clusters <- t(neighbor_clusters)

# 检查 neighbor_clusters 的形状是否正确
dim(neighbor_clusters)  # 应该是 1000 x 80

# 计算每个细胞的香农熵
entropy_values <- apply(neighbor_clusters, 1, function(neighbor_cluster) {
  neighbor_cluster_counts <- table(neighbor_cluster)  # 统计邻居簇的数量
  p_t <- neighbor_cluster_counts / sum(neighbor_cluster_counts)  # 计算比例
  -sum(p_t * log2(p_t))  # 香农熵公式
})

# 4. 将熵值添加到 sampled_cells 数据框
sampled_cells$entropy <- entropy_values

# 5. 结果分析
# 5.1 查看熵值的总结信息
summary(sampled_cells$entropy)

# 5.2 可视化熵值的分布
ggplot(sampled_cells, aes(x = entropy)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Shannon Entropy",
       x = "Entropy",
       y = "Cell Count")

# 6. 比较不同癌症类型的熵值
# 6.1 按癌症类型分组，计算每组的熵值均值和标准差
grouped_entropy <- sampled_cells %>%
  group_by(cancer_type) %>%
  summarise(mean_entropy = mean(entropy, na.rm = TRUE),  # 计算均值
            sd_entropy = sd(entropy, na.rm = TRUE))  # 计算标准差

# 查看分组结果
print(grouped_entropy)

# 6.2 可视化分组熵值
ggplot(grouped_entropy, aes(x = cancer_type, y = mean_entropy)) +
  geom_bar(stat = "identity", fill = "coral") +
  geom_errorbar(aes(ymin = mean_entropy - sd_entropy, ymax = mean_entropy + sd_entropy), width = 0.2) +
  labs(title = "Mean Entropy by Cancer Type",
       x = "Cancer Type",
       y = "Mean Shannon Entropy") +
  theme_minimal()
