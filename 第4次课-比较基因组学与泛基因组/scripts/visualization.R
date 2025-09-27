#!/usr/bin/env Rscript
# 泛基因组数据可视化脚本
# 用于生成泛基因组分析的各种图表
# 
# 作者: 微生物基因组学课程组
# 版本: 1.0.0
# 日期: 2025

# 加载必需的R包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(reshape2)
  library(RColorBrewer)
  library(pheatmap)
  library(ape)
  library(ggtree)
  library(gridExtra)
})

# 设置参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  input_dir <- "."
  output_dir <- "figures"
} else {
  input_dir <- args[1]
  output_dir <- ifelse(length(args) > 1, args[2], "figures")
}

# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("泛基因组可视化脚本 v1.0.0\n")
cat("输入目录:", input_dir, "\n")
cat("输出目录:", output_dir, "\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

# 函数：读取基因存在/缺失矩阵
read_gene_presence_matrix <- function(file_path) {
  cat("读取基因存在/缺失矩阵...\n")
  
  if (!file.exists(file_path)) {
    stop("文件不存在: ", file_path)
  }
  
  # 读取数据
  data <- read_csv(file_path, show_col_types = FALSE)
  
  # 获取菌株列（跳过前14列的注释信息）
  strain_cols <- colnames(data)[15:ncol(data)]
  
  cat("✅ 数据读取完成\n")
  cat("   基因家族数量:", nrow(data), "\n")
  cat("   菌株数量:", length(strain_cols), "\n")
  
  return(list(data = data, strain_cols = strain_cols))
}

# 函数：创建泛基因组增长曲线
create_pangenome_curve <- function(gene_data, strain_cols, output_dir) {
  cat("生成泛基因组增长曲线...\n")
  
  n_strains <- length(strain_cols)
  
  # 计算累积基因数
  pangenome_sizes <- c()
  core_sizes <- c()
  
  for (i in 1:n_strains) {
    current_strains <- strain_cols[1:i]
    
    # 计算泛基因组大小（至少在一个菌株中存在）
    pangenome_count <- 0
    core_count <- 0
    
    for (row_idx in 1:nrow(gene_data$data)) {
      presence_count <- sum(!is.na(gene_data$data[row_idx, current_strains]) & 
                           gene_data$data[row_idx, current_strains] != "")
      
      if (presence_count > 0) {
        pangenome_count <- pangenome_count + 1
      }
      
      if (presence_count == i) {
        core_count <- core_count + 1
      }
    }
    
    pangenome_sizes <- c(pangenome_sizes, pangenome_count)
    core_sizes <- c(core_sizes, core_count)
  }
  
  # 创建数据框
  curve_data <- data.frame(
    n_genomes = 1:n_strains,
    pangenome = pangenome_sizes,
    core = core_sizes
  )
  
  # 重塑数据用于绘图
  curve_long <- melt(curve_data, id.vars = "n_genomes", 
                     variable.name = "gene_type", value.name = "gene_count")
  
  # 创建图表
  p <- ggplot(curve_long, aes(x = n_genomes, y = gene_count, color = gene_type)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    scale_color_manual(values = c("pangenome" = "#e74c3c", "core" = "#2ecc71"),
                       labels = c("泛基因组", "核心基因组")) +
    labs(title = "泛基因组增长曲线",
         x = "基因组数量",
         y = "基因数量",
         color = "基因类型") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "bottom")
  
  # 保存图表
  curve_file <- file.path(output_dir, "pangenome_curve.png")
  ggsave(curve_file, p, width = 10, height = 6, dpi = 300)
  cat("泛基因组增长曲线已保存到:", curve_file, "\n")
  
  return(curve_data)
}

# 函数：创建基因组成饼图
create_genome_composition_pie <- function(gene_data, strain_cols, output_dir) {
  cat("生成基因组成饼图...\n")
  
  n_strains <- length(strain_cols)
  
  # 分类基因
  core_count <- 0
  accessory_count <- 0
  unique_count <- 0
  
  for (row_idx in 1:nrow(gene_data$data)) {
    presence_count <- sum(!is.na(gene_data$data[row_idx, strain_cols]) & 
                         gene_data$data[row_idx, strain_cols] != "")
    
    if (presence_count == n_strains) {
      core_count <- core_count + 1
    } else if (presence_count == 1) {
      unique_count <- unique_count + 1
    } else {
      accessory_count <- accessory_count + 1
    }
  }
  
  # 创建饼图数据
  pie_data <- data.frame(
    category = c("核心基因", "辅助基因", "独特基因"),
    count = c(core_count, accessory_count, unique_count),
    percentage = c(core_count, accessory_count, unique_count) / 
                 (core_count + accessory_count + unique_count) * 100
  )
  
  # 添加标签
  pie_data$label <- paste0(pie_data$category, "\n", 
                          pie_data$count, " (", 
                          round(pie_data$percentage, 1), "%)")
  
  # 创建饼图
  p <- ggplot(pie_data, aes(x = "", y = count, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("核心基因" = "#2ecc71", 
                                "辅助基因" = "#f39c12", 
                                "独特基因" = "#e74c3c")) +
    labs(title = "泛基因组组成",
         fill = "基因类型") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5))
  
  # 保存图表
  pie_file <- file.path(output_dir, "core_accessory_pie.png")
  ggsave(pie_file, p, width = 8, height = 6, dpi = 300)
  cat("基因组成饼图已保存到:", pie_file, "\n")
  
  return(pie_data)
}

# 函数：创建基因存在/缺失热图
create_presence_absence_heatmap <- function(gene_data, strain_cols, output_dir) {
  cat("生成基因存在/缺失热图...\n")
  
  # 选择前100个基因用于展示（避免图片过大）
  n_genes_to_show <- min(100, nrow(gene_data$data))
  
  # 创建二进制矩阵
  binary_matrix <- matrix(0, nrow = n_genes_to_show, ncol = length(strain_cols))
  rownames(binary_matrix) <- gene_data$data$Gene[1:n_genes_to_show]
  colnames(binary_matrix) <- strain_cols
  
  for (i in 1:n_genes_to_show) {
    for (j in 1:length(strain_cols)) {
      if (!is.na(gene_data$data[i, strain_cols[j]]) && 
          gene_data$data[i, strain_cols[j]] != "") {
        binary_matrix[i, j] <- 1
      }
    }
  }
  
  # 创建热图
  heatmap_file <- file.path(output_dir, "presence_absence_heatmap.png")
  
  png(heatmap_file, width = 12, height = 10, units = "in", res = 300)
  
  pheatmap(binary_matrix,
           color = c("white", "#2c3e50"),
           border_color = "grey90",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = FALSE,
           show_colnames = TRUE,
           main = "基因存在/缺失热图 (前100个基因)",
           fontsize = 10,
           legend_labels = c("缺失", "存在"))
  
  dev.off()
  
  cat("基因存在/缺失热图已保存到:", heatmap_file, "\n")
}

# 函数：可视化系统发育树
visualize_phylogenetic_tree <- function(input_dir, output_dir) {
  cat("可视化系统发育树...\n")
  
  # 查找系统发育树文件
  tree_files <- list.files(input_dir, pattern = "*.newick", full.names = TRUE)
  
  if (length(tree_files) == 0) {
    cat("⚠️  未找到系统发育树文件 (.newick)\n")
    return(NULL)
  }
  
  for (tree_file in tree_files) {
    tryCatch({
      # 读取系统发育树
      tree <- read.tree(tree_file)
      
      # 创建树的可视化
      p <- ggtree(tree) + 
        geom_tiplab(size = 3) +
        geom_nodepoint(color = "#2c3e50", size = 2) +
        labs(title = paste("系统发育树:", basename(tree_file))) +
        theme_tree2() +
        theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
      
      # 保存图表
      tree_plot_file <- file.path(output_dir, 
                                  paste0("phylogenetic_tree_", 
                                        tools::file_path_sans_ext(basename(tree_file)), 
                                        ".png"))
      ggsave(tree_plot_file, p, width = 10, height = 8, dpi = 300)
      cat("系统发育树图已保存到:", tree_plot_file, "\n")
      
    }, error = function(e) {
      cat("⚠️  无法处理树文件", basename(tree_file), ":", e$message, "\n")
    })
  }
}

# 函数：创建综合报告图
create_summary_plot <- function(curve_data, pie_data, output_dir) {
  cat("生成综合报告图...\n")
  
  # 重新创建曲线图（用于组合）
  curve_long <- melt(curve_data, id.vars = "n_genomes", 
                     variable.name = "gene_type", value.name = "gene_count")
  
  p1 <- ggplot(curve_long, aes(x = n_genomes, y = gene_count, color = gene_type)) +
    geom_line(size = 1) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c("pangenome" = "#e74c3c", "core" = "#2ecc71"),
                       labels = c("泛基因组", "核心基因组")) +
    labs(title = "泛基因组增长曲线",
         x = "基因组数量",
         y = "基因数量",
         color = "基因类型") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "bottom")
  
  # 重新创建饼图（用于组合）
  p2 <- ggplot(pie_data, aes(x = "", y = count, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("核心基因" = "#2ecc71", 
                                "辅助基因" = "#f39c12", 
                                "独特基因" = "#e74c3c")) +
    labs(title = "基因组成分布",
         fill = "基因类型") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  
  # 组合图表
  combined_plot <- grid.arrange(p1, p2, ncol = 2, 
                               top = "泛基因组分析综合报告")
  
  # 保存组合图
  summary_file <- file.path(output_dir, "pangenome_summary.png")
  ggsave(summary_file, combined_plot, width = 14, height = 6, dpi = 300)
  cat("综合报告图已保存到:", summary_file, "\n")
}

# 主函数
main <- function() {
  tryCatch({
    # 查找基因存在/缺失矩阵文件
    gene_presence_file <- file.path(input_dir, "gene_presence_absence.csv")
    
    if (!file.exists(gene_presence_file)) {
      stop("未找到基因存在/缺失矩阵文件: gene_presence_absence.csv")
    }
    
    # 读取数据
    gene_data <- read_gene_presence_matrix(gene_presence_file)
    
    # 生成各种图表
    curve_data <- create_pangenome_curve(gene_data, gene_data$strain_cols, output_dir)
    pie_data <- create_genome_composition_pie(gene_data, gene_data$strain_cols, output_dir)
    create_presence_absence_heatmap(gene_data, gene_data$strain_cols, output_dir)
    
    # 可视化系统发育树
    visualize_phylogenetic_tree(input_dir, output_dir)
    
    # 创建综合报告
    create_summary_plot(curve_data, pie_data, output_dir)
    
    cat("\n✅ 所有图表生成完成！\n")
    cat("图表文件保存在:", output_dir, "\n")
    
    # 列出生成的文件
    cat("\n生成的图表文件:\n")
    plot_files <- list.files(output_dir, pattern = "*.png", full.names = FALSE)
    for (file in plot_files) {
      cat("  -", file, "\n")
    }
    
  }, error = function(e) {
    cat("❌ 错误:", e$message, "\n")
    quit(status = 1)
  })
}

# 运行主函数
main()