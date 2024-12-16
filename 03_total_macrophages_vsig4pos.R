scobj = schard::h5ad2seurat('Total_macrophages_harmony.h5ad')
Idents(scobj) = "leiden"
table(scobj@meta.data$leiden)
scobj <- NormalizeData(scobj, normalization.method = 'LogNormalize', scale.factor = 10000)
leiden_counts <- table(scobj@meta.data$leiden)
sorted_levels <- names(sort(leiden_counts, decreasing = TRUE))
scobj@meta.data$leiden <- factor(scobj@meta.data$leiden, levels = sorted_levels)
table(scobj@meta.data$leiden)
#DEG分析#####################################
deg <- FindAllMarkers(scobj,assay = "RNA",only.pos = TRUE,min.pct = 0.1, 
                      min.diff.pct = 0.05,logfc.threshold = 0.25,slot="data") %>%
  Add_Pct_Diff()%>%
  group_by(cluster)%>%
  arrange(desc(avg_log2FC))%>%
  arrange(cluster)
write_rds(deg,"Total_macrophages_27clusters_DEGs.rds")
filter_deg <- deg %>%
  group_by(cluster) %>%
  arrange(cluster, desc(pct_diff))

top100 = filter_deg %>% group_by(cluster) %>% 
  top_n(n = 100, wt = pct_diff)

tmp <- top100 %>%
  select(cluster, gene)
gene_matrix <- tmp %>%
  group_by(cluster) %>%
  mutate(row_id = row_number()) %>%  
  pivot_wider(names_from = cluster, values_from = gene) %>%  
  select(-row_id)
#作图#####################################
key <- c("VSIG4","TREM2","FOLR2","CD5L","CD168","CLEC2","CLEC4F","TIMD4")
key2 <- c("ABCA1", "ABCG1", "APOE", "AXL", "CETP", "CD14", "CD163", "CD207", "CD36", 
  "CD5L", "CD68", "CD9", "CETP", "CSTB", "DPP4", "ESAM", "FCGR1A", "FCN1", 
  "FOLR2", "GPX1", "GPNMB", "HBA1", "HBB", "HMOX1", "HP", "IGFBP7", "IL1B", 
  "LIPA", "MARCO", "MRC1", "MSR1", "NR1H3", "NPC1", "NPC2", "ORM1", "SAA1", 
  "SLC40A1", "SLC48A1", "SPP1", "THBD", "TIMD4", "TREM2", "VSIG4")

DotPlot(scobj,feature=key,dot.scale=6,cols="RdBu",
        scale.min=5,scale.max = 60)+
  scale_x_discrete("")+scale_y_discrete("")+RotatedAxis()+
  theme(axis.text.x=element_text(face="italic",size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(legend.title = element_text(size=6))+
  theme(legend.text=element_text(size=6))+
  guides(colour = guide_colourbar(title = "Average Expression",barwidth = .5, barheight = 3))
DimPlot(seuratObj3, reduction = "Xumap_", cols=c54,
        label = T, label.size = 4,group.by = 'leiden') #18个亚群
FeaturePlot(scobj,feature=c("VSIG4"),reduction="Xumap_",pt.size = 0.5)
#VSIG4表达的阳性亚群处理和保存h5ad#####################################
scobj@meta.data$leiden <- as.character(scobj@meta.data$leiden)
clusters <- c('2', '1', '4', '11', '16', '17', '23', '24')
table(scobj@meta.data$leiden)
vsig <- subset(scobj, subset = leiden %in% c('2', '1', '4', '11', '16', '17', '23', '24'))
# 假设原 Seurat 对象名称为 old_seurat_object
old_seurat_object<-vsig
old_seurat_object$X_index<-NULL
old_seurat_object$`_index`<-NULL
old_seurat_object$X_index.1<-NULL
old_seurat_object$a<-NULL

# Step 1: 提取计数矩阵和元数据
count_matrix <- GetAssayData(old_seurat_object, slot = "counts")
metadata <- old_seurat_object@meta.data
# Step 2: 创建新 Seurat 对象
options(Seurat.object.assay.version = "v3")
new_seurat_object <- CreateSeuratObject(counts = count_matrix)
# Step 3: 将元数据添加到新对象
# 确保 metadata 的行名与 count_matrix 的列名一致
rownames(metadata) <- colnames(count_matrix)
new_seurat_object <- AddMetaData(new_seurat_object, metadata = metadata)
sce<-new_seurat_object
obj_name<-"Total_macrophages_VSIG4_pos"
factor_cols <- sapply(sce@meta.data, is.factor)
sce@meta.data[factor_cols] <- lapply(sce@meta.data[factor_cols], as.character)

# step 2: 保存为 .h5seurat 临时文件
h5seurat_filename <- paste0(obj_name, ".h5seurat")
SaveH5Seurat(sce, filename = h5seurat_filename, overwrite = TRUE)

# step 3: 将 .h5seurat 文件转换为 .h5ad 文件
Convert(h5seurat_filename, dest = "h5ad", assay = "RNA", overwrite = TRUE)

#处理metadata###================================================================================
info <- read_xlsx("./良恶性整合/Nature_neutrophil_etiology.xlsx")
hcc <-subset(scobj,subset = condition == "HCC")
table(hcc@meta.data$orig.ident)
table(hcc@meta.data$etiology)
hcc@meta.data$etiology <- info$Virus[match(hcc@meta.data$orig.ident, info$Patient)]
table(scobj@meta.data$etiology)
hcc_cells <- rownames(hcc@meta.data)
scobj@meta.data[hcc_cells, "etiology"] <- hcc@meta.data$etiology
#age
hcc <-subset(scobj,subset = condition == "HCC")
table(scobj@meta.data$age)
table(info$Age)
table(hcc@meta.data$etiology)
hcc@meta.data$age <- info$Age[match(hcc@meta.data$orig.ident, info$Patient)]
table(scobj@meta.data$age)
hcc_cells <- rownames(hcc@meta.data)
scobj@meta.data[hcc_cells, "age"] <- hcc@meta.data$age
scobj@meta.data <- scobj@meta.data %>% mutate(age = recode(age, "NA" = "Unclear"))
table(scobj@meta.data$age)
tmp<-subset(scobj,subset = age=="HBV")
tmp@meta.data$orig.ident,tmp@meta.data$age
# 创建映射表，包含编号和年龄
age_map <- data.frame(
  orig.ident = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P13", "P14", "P15"),
  age = c(59, 71, 75, 68, 74, 72, 65, 48, 77, 77, 55, 56, 56, 76)
)
# 更新 tmp@meta.data$age 列
tmp@meta.data$age <- age_map$age[match(tmp@meta.data$batch, age_map$orig.ident)]
tmp_cells <- rownames(tmp@meta.data)
scobj@meta.data[tmp_cells, "age"] <- tmp@meta.data$age
table(scobj@meta.data$age)
scobj@meta.data$age2<-scobj@meta.data$age
scobj@meta.data$age <- as.character(scobj@meta.data$age)
scobj@meta.data <- scobj@meta.data %>%
  mutate(
    age2 = case_when(
      age %in% c("<30", "<30 years", "30-60", "30-60 years") ~ "<60",  # 归为 <60
      age %in% c(">60", ">60 years") ~ ">=60",  # 归为 >=60
      age %in% c(">0 years", "Unclear") ~ "Unclear",  # 归为 Unclear
      grepl("^\\d+$", age) & as.numeric(age) < 60 ~ "<60",  # 只转换纯数字，且小于60
      grepl("^\\d+$", age) & as.numeric(age) >= 60 ~ ">=60",  # 只转换纯数字，且大于等于60
      TRUE ~  "Unclear"
    )
  )

# 检查新列 age2 的分布
table(scobj@meta.data$age2)
#检查NA#####################################
na <-subset(scobj,subset = etiology == "NA")
na@meta.data$new_column <- sub("_[^_]+$", "", na@meta.data$X_index) #A004 ICC癌旁
na_cells <- rownames(na@meta.data)
scobj@meta.data[na_cells, "etiology"] <- "Unclear"
table(scobj@meta.data$etiology)
cir <-subset(scobj,subset = etiology == "Cirrhosis")
unique_sample_names <- unique(sub("_[^_]+$", "", rownames(cir@meta.data)))
unique_sample_names
scobj@meta.data$etiology2<-scobj@meta.data$etiology
scobj@meta.data <- scobj@meta.data %>%
  mutate(
    etiology2 = case_when(
      etiology %in% c("Cirrhosis", "NBNC") ~ "Unclear",
      etiology %in% c("HBV", "HCV") ~ "Virus",
      etiology %in% c("PBC", "PSC") ~ "Autoimmune",
      TRUE ~ etiology 
    )
  )
table(scobj@meta.data$etiology2)
table(scobj@meta.data$condition)
table(scobj@meta.data$group)
scobj@meta.data <- scobj@meta.data %>% mutate(group = recode(group, "NA" = "HCC"))

#VSIG4阳性亚群的特征#####################################

c1marker <- FindMarkers(object = scobj,ident.1 = 1,ident.2 = 2,
                        logfc.threshold = 0.25,test.use = "wilcox",
                        only.pos = TRUE)
c2marker <- FindMarkers(object = scobj,ident.1 = 2,ident.2 = 1,
                        logfc.threshold = 0.25,test.use = "wilcox",
                        only.pos = TRUE)
vsig <- NormalizeData(vsig, normalization.method = 'LogNormalize', scale.factor = 10000)
key <- c("VSIG4","TREM2","FOLR2","CD5L","CD168","CLEC2","CLEC4F","TIMD4")
DotPlot(scobj,feature=key,dot.scale=6,cols="RdBu",
       scale.max = 60)+
  scale_x_discrete("")+scale_y_discrete("")+RotatedAxis()+
  theme(axis.text.x=element_text(face="italic",size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(legend.title = element_text(size=6))+
  theme(legend.text=element_text(size=6))+
  guides(colour = guide_colourbar(title = "Average Expression",barwidth = .5, barheight = 3))
DimPlot(vsig, reduction = "Xumap_", cols=c54,
        label = T, label.size = 4,group.by = 'leiden') #18个亚群
FeaturePlot(scobj,feature=key,reduction="Xumap_",pt.size = 0.5,raster=FALSE)
