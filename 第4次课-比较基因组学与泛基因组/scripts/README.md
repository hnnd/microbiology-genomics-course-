# 泛基因组分析脚本说明

本目录包含第4次课实践操作所需的所有分析脚本。

## 脚本文件列表

### 1. pangenome_analysis.py
**功能**: 泛基因组统计分析
**用途**: 
- 质量控制检查
- 基本统计分析
- 生成统计报告

**使用方法**:
```bash
# 质量控制
python3 scripts/pangenome_analysis.py --mode qc --input gff_files/

# 统计分析
python3 scripts/pangenome_analysis.py --mode stats --input ./
```

**输出文件**:
- `results/quality_report.txt` - 质量控制报告
- `results/pangenome_statistics.txt` - 统计分析报告
- `results/pangenome_stats.json` - 详细统计数据

### 2. core_phylogeny.py
**功能**: 核心基因组系统发育分析
**用途**:
- 构建系统发育树
- 计算进化距离
- 生成系统发育报告

**使用方法**:
```bash
# 使用邻接法构建系统发育树
python3 scripts/core_phylogeny.py --input core_gene_alignment.aln --method nj

# 使用UPGMA方法
python3 scripts/core_phylogeny.py --input core_gene_alignment.aln --method upgma
```

**输出文件**:
- `results/phylogenetic_tree_nj.newick` - 系统发育树文件
- `results/phylogenetic_tree_nj.png` - 系统发育树图片
- `results/distance_matrix.txt` - 距离矩阵
- `results/phylogeny_report.txt` - 系统发育分析报告

### 3. functional_analysis.py
**功能**: 功能分类分析
**用途**:
- 基因功能分类
- 菌株特异性基因分析
- 辅助基因组深入分析

**使用方法**:
```bash
# 基础功能分析
python3 scripts/functional_analysis.py --input gene_presence_absence.csv --mode basic

# 辅助基因组分析
python3 scripts/functional_analysis.py --input gene_presence_absence.csv --mode accessory
```

**输出文件**:
- `results/functional_classification.txt` - 功能分类报告
- `results/strain_specific_genes.txt` - 菌株特异性基因
- `results/accessory_functions.txt` - 辅助基因组分析
- `results/functional_pie_charts.png` - 功能分类饼图
- `results/functional_comparison.png` - 功能分类比较图

### 4. visualization.R
**功能**: 数据可视化
**用途**:
- 生成泛基因组增长曲线
- 创建基因组成饼图
- 绘制基因存在/缺失热图
- 可视化系统发育树

**使用方法**:
```bash
# 使用默认参数
Rscript scripts/visualization.R

# 指定输入和输出目录
Rscript scripts/visualization.R . figures/
```

**输出文件**:
- `figures/pangenome_curve.png` - 泛基因组增长曲线
- `figures/core_accessory_pie.png` - 基因组成饼图
- `figures/presence_absence_heatmap.png` - 基因存在/缺失热图
- `figures/phylogenetic_tree_*.png` - 系统发育树图片
- `figures/pangenome_summary.png` - 综合报告图

## 依赖软件和包

### Python依赖
```bash
pip install pandas numpy matplotlib seaborn biopython
```

### R依赖
```r
install.packages(c("ggplot2", "dplyr", "readr", "reshape2", 
                   "RColorBrewer", "pheatmap", "ape", "ggtree", "gridExtra"))
```

### 外部软件
- Roary (>= 3.13.0)
- Prokka (>= 1.14.6)
- FastTree (>= 2.1.11)

## 完整分析流程

1. **数据预处理**
   ```bash
   # 使用Prokka注释基因组
   prokka --outdir strain1_annotation --prefix strain1 strain1.gbff
   
   # 收集GFF文件
   mkdir gff_files
   cp */*.gff gff_files/
   ```

2. **质量控制**
   ```bash
   python3 scripts/pangenome_analysis.py --mode qc --input gff_files/
   ```

3. **泛基因组分析**
   ```bash
   # 运行Roary
   roary -p 4 -e -n -v gff_files/*.gff
   
   # 统计分析
   python3 scripts/pangenome_analysis.py --mode stats --input ./
   ```

4. **系统发育分析**
   ```bash
   python3 scripts/core_phylogeny.py --input core_gene_alignment.aln --method nj
   ```

5. **功能分析**
   ```bash
   python3 scripts/functional_analysis.py --input gene_presence_absence.csv --mode basic
   ```

6. **数据可视化**
   ```bash
   Rscript scripts/visualization.R
   ```

## 故障排除

### 常见问题

1. **Python模块导入错误**
   - 确保已安装所有必需的Python包
   - 检查Python版本 (需要Python 3.7+)

2. **R包加载失败**
   - 安装缺失的R包
   - 检查R版本 (需要R 4.0+)

3. **文件路径错误**
   - 确保在正确的工作目录中运行脚本
   - 检查输入文件是否存在

4. **内存不足**
   - 减少并行线程数
   - 处理较小的数据集

### 获取帮助

如果遇到问题，可以：
1. 查看脚本的帮助信息：`python3 script.py --help`
2. 检查错误日志和输出信息
3. 联系课程老师或同学

## 版本信息

- 创建日期: 2025年
- 版本: 1.0.0
- 作者: 微生物基因组学课程组

## 许可证

本脚本集合遵循教育用途开源协议，仅供学习和研究使用。