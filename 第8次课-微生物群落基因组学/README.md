# 第8次课：微生物群落基因组学

## 课程概述

本次课程是微生物基因组学课程的最后一次课，重点介绍微生物群落基因组学的理论基础和实践方法。课程涵盖宏基因组vs单菌基因组、群落功能分析和环境基因组学应用等内容。

## 文件结构

```
第8次课-微生物群落基因组学/
├── README.md                                    # 本文件
├── 理论课件-微生物群落基因组学.md                  # Marp格式理论课件
├── 实践操作-宏基因组与群落功能分析.md              # 实践操作指南
├── images/                                      # SVG图像文件
│   ├── mag_analysis_workflow.svg               # MAG分析流程图
│   ├── community_function_analysis.svg         # 群落功能分析流程图
│   ├── nitrogen_cycle_genes.svg                # 氮循环基因网络图
│   └── multi_omics_integration.svg             # 多组学整合分析图
├── scripts/                                    # Python分析脚本
│   ├── mag_quality_assessor.py                 # MAG质量评估脚本
│   ├── community_function_analyzer.py          # 群落功能分析脚本
│   ├── pathway_reconstructor.py                # 代谢通路重建脚本
│   └── visualization_generator.py              # 可视化生成脚本
├── data/                                       # 数据文件目录
├── results/                                    # 分析结果目录
└── figures/                                    # 图表输出目录
```

## 理论课件内容

### 主要章节
1. **宏基因组vs单菌基因组**
   - 培养vs非培养微生物
   - MAG技术原理和应用
   - 单细胞基因组技术

2. **微生物群落功能**
   - 群落代谢网络
   - 功能冗余与互补
   - 关键物种识别

3. **环境基因组学应用**
   - 生物地球化学循环
   - 污染修复基因
   - 极端环境适应
   - 气候变化响应

4. **前沿技术与应用前景**
   - 新兴技术发展
   - 多组学整合分析
   - 实际应用前景

### 课件特色
- 使用Marp格式，支持现代化演示
- 包含丰富的SVG矢量图表
- 互动环节和思考问题
- 前沿技术和应用案例

## 实践操作内容

### 操作目标
- 掌握MAG基因组质量评估方法
- 理解群落功能分析流程
- 应用宏基因组数据进行功能预测
- 分析环境微生物群落的生态功能

### 主要步骤
1. **MAG基因组质量评估** (40分钟)
   - 基本信息统计
   - CheckM质量评估
   - 物种分类注释

2. **群落功能分析** (50分钟)
   - 物种组成分析
   - 功能基因定量
   - 代谢通路重建

3. **结果分析与可视化** (30分钟)
   - 代谢网络重建
   - 群落功能可视化
   - 统计分析

### 分析工具
- **CheckM**: MAG质量评估
- **GTDB-Tk**: 物种分类注释
- **MetaPhlAn**: 物种组成分析
- **HUMAnN**: 功能基因定量
- **Python脚本**: 自定义分析和可视化

## 脚本说明

### mag_quality_assessor.py
- **功能**: 批量评估MAG基因组质量
- **输入**: MAG基因组FASTA文件
- **输出**: 质量评估报告和统计结果
- **主要指标**: 完整性、污染率、应变异质性

### community_function_analyzer.py
- **功能**: 分析微生物群落功能特征
- **输入**: 宏基因组测序数据
- **输出**: 功能谱、多样性指数、基因家族表
- **分析内容**: KEGG通路、COG分类、Pfam结构域

### pathway_reconstructor.py
- **功能**: 重建代谢通路和网络
- **输入**: 功能注释结果
- **输出**: 代谢网络、通路完整性、关键代谢物
- **分析方法**: 网络拓扑分析、中心性计算

### visualization_generator.py
- **功能**: 生成分析结果可视化图表
- **输入**: 各种分析结果文件
- **输出**: PNG格式图表文件
- **图表类型**: 热图、柱状图、网络图、仪表板

## 使用说明

### 环境要求
- Python 3.8+
- 必需软件包: pandas, numpy, matplotlib, seaborn, networkx
- 生物信息学软件: CheckM, GTDB-Tk, MetaPhlAn, HUMAnN

### 运行步骤
1. **准备环境**
   ```bash
   # 创建工作目录
   mkdir -p ~/genomics_course/lesson8
   cd ~/genomics_course/lesson8
   
   # 复制脚本文件
   cp scripts/* .
   ```

2. **运行MAG质量评估**
   ```bash
   python3 mag_quality_assessor.py --mag_dir data/mag_genomes --output_dir results
   ```

3. **运行群落功能分析**
   ```bash
   python3 community_function_analyzer.py --input metagenome_sample.fastq.gz --output_dir results
   ```

4. **运行代谢通路重建**
   ```bash
   python3 pathway_reconstructor.py --input_dir results/functional_analysis --output_dir results
   ```

5. **生成可视化图表**
   ```bash
   python3 visualization_generator.py --results_dir results --output_dir figures
   ```

### 预期结果
- MAG质量评估报告
- 群落功能谱分析
- 代谢网络重建结果
- 多种可视化图表
- 综合分析报告

## 学习目标

完成本次课程后，学生应能够：

1. **理解核心概念**
   - 宏基因组学与传统微生物学的区别
   - MAG技术的原理和应用
   - 群落功能分析的基本方法

2. **掌握分析技能**
   - MAG基因组质量评估
   - 群落功能注释和定量
   - 代谢网络重建和分析
   - 结果可视化和解释

3. **应用实践能力**
   - 独立完成宏基因组分析项目
   - 评估微生物群落的生态功能
   - 解释分析结果的生物学意义

## 扩展学习

### 进阶主题
- 时间序列宏基因组分析
- 空间宏基因组学
- 多组学数据整合
- 机器学习在宏基因组学中的应用

### 相关资源
- **数据库**: IMG/M, MGnify, GTDB
- **软件工具**: Anvi'o, MEGAN, QIIME2
- **在线平台**: MG-RAST, Galaxy
- **文献资源**: Nature Microbiology, Microbiome

## 课程总结

第8次课作为微生物基因组学课程的总结，展示了从单个基因组到群落基因组的研究进展。通过理论学习和实践操作，学生将全面掌握现代微生物基因组学的核心技术和应用方法，为未来的科研工作奠定坚实基础。

---

**课程信息**
- 课程名称：微生物基因组学
- 第8次课：微生物群落基因组学
- 时长：4小时（理论2小时 + 实践2小时）
- 难度：⭐⭐⭐⭐☆