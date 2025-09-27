# 第7次课：功能基因挖掘与应用

## 课程概述

本次课程专注于微生物功能基因的挖掘、分析和应用，涵盖次级代谢产物基因簇、抗生素抗性基因和毒力因子的系统性分析方法。

## 文件结构

```
第7次课-功能基因挖掘与应用/
├── README.md                                    # 课程说明文档
├── 理论课件-功能基因挖掘与应用.md                  # Marp格式理论课件
├── 实践操作-次级代谢与抗性基因分析.md              # 实践操作指南
├── images/                                      # SVG图像文件
│   ├── enzyme_mining_workflow.svg               # 酶基因挖掘工作流程
│   ├── biosynthetic_gene_cluster.svg            # 生物合成基因簇结构
│   ├── resistance_mechanisms.svg                # 抗性机制网络图
│   └── gene_application_pipeline.svg            # 基因应用开发流程
├── scripts/                                     # Python分析脚本
│   ├── bgc_analyzer.py                          # BGC分析脚本
│   ├── resistance_annotator.py                  # 抗性基因注释脚本
│   ├── virulence_scanner.py                     # 毒力因子扫描脚本
│   └── network_builder.py                       # 功能网络构建脚本
├── data/                                        # 数据文件目录
├── results/                                     # 分析结果目录
└── figures/                                     # 可视化图表目录
```

## 主要内容

### 理论部分 (2小时)

1. **生物技术基因挖掘策略**
   - 酶基因功能预测方法
   - 新酶活性发现策略
   - 酶工程改造靶点识别
   - 工业酶基因挖掘案例

2. **次级代谢产物基因簇分析**
   - BGC结构特征和组织模式
   - PKS/NRPS基因识别方法
   - 产物结构预测技术
   - 新化合物发现策略

3. **抗性与毒力机制解析**
   - 抗生素抗性基因分类
   - 抗性机制多样性分析
   - 毒力因子进化特征
   - 病原性岛识别方法

4. **功能基因发现与验证**
   - 高通量筛选策略
   - 功能验证实验设计
   - 应用开发流程
   - 前沿技术发展趋势

### 实践部分 (2小时)

1. **次级代谢基因簇分析** (45分钟)
   - antiSMASH在线分析
   - 本地BGC预测脚本
   - 基因簇结果解读

2. **抗性基因注释分析** (35分钟)
   - CARD数据库搜索
   - 抗性机制分类
   - 抗性谱预测

3. **毒力因子分析** (30分钟)
   - VFDB数据库搜索
   - 毒力岛预测
   - 致病性评估

4. **功能网络构建** (30分钟)
   - 数据整合
   - 网络构建
   - 结果可视化

## 使用说明

### 环境要求

- Python 3.8+
- BLAST+ 2.10.0+
- HMMER 3.3+
- Prodigal 2.6.3+
- Diamond 2.0.0+

### 快速开始

1. **查看理论课件**
   ```bash
   # 使用Marp渲染课件
   marp 理论课件-功能基因挖掘与应用.md --pdf
   ```

2. **运行实践操作**
   ```bash
   # 创建工作目录
   mkdir -p ~/genomics_course/lesson7
   cd ~/genomics_course/lesson7
   
   # 运行BGC分析
   python3 scripts/bgc_analyzer.py --input genome.fna --output results/bgc_analysis
   
   # 运行抗性基因注释
   python3 scripts/resistance_annotator.py --proteins proteins.faa --genome genome.fna --output results/resistance_analysis
   
   # 运行毒力因子扫描
   python3 scripts/virulence_scanner.py --genome genome.fna --proteins proteins.faa --output results/virulence_analysis
   
   # 构建功能网络
   python3 scripts/network_builder.py --mode integrate --bgc results/bgc_analysis/cluster_summary.txt --resistance results/resistance_analysis/resistance_genes.txt --virulence results/virulence_analysis/virulence_factors.txt --output results/
   ```

## 学习目标

完成本次课程后，学生应能够：

- 掌握生物技术相关功能基因的挖掘策略和预测方法
- 理解次级代谢产物生物合成基因簇的结构特征和分析方法
- 分析微生物抗性和毒力机制的分子基础
- 应用功能基因挖掘技术解决实际生物技术问题

## 重要工具和数据库

### 在线工具
- [antiSMASH](https://antismash.secondarymetabolites.org) - 次级代谢基因簇预测
- [CARD](https://card.mcmaster.ca) - 抗生素抗性基因数据库
- [VFDB](http://www.mgc.ac.cn/VFs) - 毒力因子数据库

### 本地工具
- BLAST+ - 序列相似性搜索
- HMMER - 蛋白质结构域搜索
- Prodigal - 基因预测
- Diamond - 快速序列比对

## 扩展资源

### 推荐阅读
1. Medema MH, et al. antiSMASH: rapid identification, annotation and analysis of secondary metabolite biosynthesis gene clusters. *Nucleic Acids Res*. 2011
2. McArthur AG, et al. The comprehensive antibiotic resistance database. *Antimicrob Agents Chemother*. 2013
3. Chen L, et al. VFDB: a reference database for bacterial virulence factors. *Nucleic Acids Res*. 2005

### 相关课程
- 第3次课：基因预测与功能注释系统
- 第4次课：比较基因组学与泛基因组
- 第6次课：微生物代谢网络重建

## 技术支持

如遇到技术问题，请：
1. 查看实践操作指南中的常见问题部分
2. 检查软件环境和依赖包安装
3. 联系课程助教或同学讨论
4. 参考工具官方文档

---

*课程材料版本：v1.0.0*  
*最后更新：2025年*