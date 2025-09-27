# 第5次课：微生物系统发育基因组学

## 课程概览

本次课程专注于微生物系统发育基因组学的理论基础和实践应用，涵盖全基因组系统发育方法、分子钟分析和发散时间估算等核心内容。

## 文件结构

```
第5次课-微生物系统发育基因组学/
├── README.md                                    # 本文件
├── 理论课件-微生物系统发育基因组学.md              # Marp格式理论课件
├── 实践操作-全基因组系统发育分析.md               # 实践操作指南
├── scripts/                                    # 分析脚本目录
│   ├── ani_calculator.py                       # ANI计算脚本
│   ├── phylogeny_builder.py                    # 系统发育树构建脚本
│   ├── divergence_estimator.py                 # 发散时间估算脚本
│   └── visualization_generator.py              # 可视化生成脚本
├── data/                                       # 数据文件目录
├── results/                                    # 分析结果目录
└── figures/                                    # 图表输出目录
```

## 理论课程内容

### 主要知识点

1. **全基因组系统发育方法**
   - 16S rRNA vs 全基因组系统发育
   - ANI、AAI、dDDH等相似性指标
   - 核心基因组系统发育分析

2. **分子钟与发散时间**
   - 分子钟理论基础
   - 校准点选择策略
   - 贝叶斯和最大似然方法

3. **系统发育信号评估**
   - 基因树与物种树冲突
   - 水平基因转移检测
   - 系统发育网络分析

### 课件特色

- 使用Marp格式，支持现代化演示
- 包含丰富的图表和可视化内容
- 设计了互动讨论环节
- 提供知识检查和测验

## 实践操作内容

### 操作流程

1. **基因组相似性分析** (40分钟)
   - FastANI计算ANI值
   - MASH距离计算
   - 相似性聚类分析

2. **全基因组系统发育分析** (60分钟)
   - 核心基因识别
   - 多序列比对
   - 系统发育树构建

3. **发散时间估算** (20分钟)
   - 分子钟分析
   - 结果可视化

### 脚本功能

#### ani_calculator.py
- 批量计算基因组间ANI值
- 生成相似性矩阵和热图
- 执行层次聚类分析
- 生成分析报告

#### phylogeny_builder.py
- 识别核心基因
- 执行多序列比对
- 连接多基因序列
- 构建系统发育树

#### divergence_estimator.py
- 估算分子进化速率
- 计算节点发散时间
- 创建时间校准树
- 生成时间轴可视化

#### visualization_generator.py
- 生成综合ANI分析图
- 可视化系统发育树
- 创建发散时间时间轴
- 生成分析仪表板

## 软件依赖

### 必需软件
- FastANI (>= 1.3)
- MASH (>= 2.3)
- IQ-TREE (>= 2.0)
- MUSCLE (可选)

### Python依赖包
```bash
pip install biopython pandas numpy matplotlib seaborn scipy
```

## 使用说明

### 环境准备
```bash
# 创建工作目录
mkdir -p ~/genomics_course/lesson5
cd ~/genomics_course/lesson5

# 复制脚本文件
cp -r scripts/ .
mkdir -p data results figures
```

### 快速开始
```bash
# 1. 准备基因组文件列表
ls *.fna > genome_list.txt

# 2. 运行ANI分析
python3 scripts/ani_calculator.py --cluster

# 3. 构建系统发育树
python3 scripts/phylogeny_builder.py --step all

# 4. 估算发散时间
python3 scripts/divergence_estimator.py

# 5. 生成可视化
python3 scripts/visualization_generator.py
```

### 分步执行
```bash
# 仅计算ANI
python3 scripts/ani_calculator.py

# 仅识别核心基因
python3 scripts/phylogeny_builder.py --step core_genes

# 仅执行序列比对
python3 scripts/phylogeny_builder.py --step alignment

# 仅构建系统发育树
python3 scripts/phylogeny_builder.py --step tree
```

## 预期输出

### 主要结果文件
- `results/ani_matrix.txt` - ANI相似性矩阵
- `results/concatenated_alignment.fasta` - 连接的多基因比对
- `results/concatenated_alignment.fasta.treefile` - 系统发育树
- `results/divergence_times.txt` - 发散时间估算
- `results/timetree.nwk` - 时间校准树

### 可视化图表
- `figures/ani_heatmap.png` - ANI热图
- `figures/phylogenetic_tree_visualization.png` - 系统发育树
- `figures/divergence_timeline.png` - 发散时间轴
- `figures/phylogenomics_dashboard.png` - 综合分析仪表板

## 教学目标

完成本次课程后，学生应能够：

1. **理解**全基因组系统发育分析的原理和优势
2. **掌握**ANI、AAI等基因组相似性指标的计算和解释
3. **应用**分子钟方法估算物种发散时间
4. **分析**系统发育信号与水平基因转移的影响
5. **评估**不同系统发育方法的适用性和局限性

## 扩展学习

### 推荐阅读
1. Konstantinidis & Tiedje. Genomic insights that advance the species definition for prokaryotes. PNAS. 2005
2. Parks et al. A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. Nature Biotechnology. 2018
3. Richter & Rosselló-Móra. Shifting the genomic gold standard for the prokaryotic species definition. PNAS. 2009

### 在线资源
- [GTDB数据库](https://gtdb.ecogenomic.org/)
- [FastANI文档](https://github.com/ParBLiSS/FastANI)
- [IQ-TREE教程](http://www.iqtree.org/doc/)

### 进阶练习
1. 使用更大的基因组数据集进行分析
2. 比较不同系统发育方法的结果
3. 整合宏基因组数据进行群落系统发育分析

## 技术支持

如遇到技术问题，请：
1. 检查软件版本和依赖关系
2. 查看脚本中的错误处理信息
3. 参考官方文档和FAQ
4. 联系课程助教或老师

## 版本信息

- 创建日期：2025年
- 版本：v1.0.0
- 最后更新：2025年
- 维护者：课程开发团队

---

*本课程材料遵循开源协议，欢迎改进和分享*