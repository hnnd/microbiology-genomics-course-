# 第6次课：微生物代谢网络重建

## 课程内容概述

本次课程包含理论课件和实践操作两部分：

### 理论课件
- **文件**: `理论课件-微生物代谢网络重建.md`
- **格式**: Marp演示文稿
- **内容**: 代谢重建理论、营养需求预测、代谢模型构建、网络分析方法

### 实践操作
- **文件**: `实践操作-代谢通路重建与分析.md`
- **时长**: 2小时
- **内容**: KEGG数据库使用、代谢网络重建、营养需求预测、网络可视化

## 目录结构

```
第6次课-微生物代谢网络重建/
├── 理论课件-微生物代谢网络重建.md     # Marp格式课件
├── 实践操作-代谢通路重建与分析.md     # 实践操作指南
├── images/                           # 课件SVG图片
│   ├── metabolic_network_topology.svg
│   ├── amino_acid_synthesis.svg
│   ├── growth_factors.svg
│   ├── fba_workflow.svg
│   └── metabolic_modules.svg
├── scripts/                          # Python分析脚本
│   ├── kegg_pathway_analyzer.py      # KEGG通路分析器
│   ├── metabolic_reconstructor.py    # 代谢网络重建器
│   ├── nutrition_predictor.py        # 营养需求预测器
│   └── network_visualizer.py         # 网络可视化器
├── data/                            # 示例数据文件
│   ├── additional_reactions.txt      # 额外反应数据
│   └── reference_metabolism.json     # 参考代谢数据
├── results/                         # 分析结果输出目录
└── figures/                         # 图表输出目录
```

## 使用说明

### 1. 环境准备

确保安装了以下Python包：
```bash
pip install biopython pandas matplotlib seaborn networkx requests
```

### 2. 数据准备

下载示例数据（大肠杆菌基因组注释）：
```bash
# 在课程目录下运行
mkdir -p ~/genomics_course/lesson6
cd ~/genomics_course/lesson6

# 下载数据文件（实践操作指南中有详细说明）
```

### 3. 分析流程

按照以下顺序运行分析脚本：

#### 步骤1：提取酶基因信息
```bash
python3 scripts/kegg_pathway_analyzer.py --extract-enzymes
```

#### 步骤2：KEGG通路映射
```bash
python3 scripts/kegg_pathway_analyzer.py --pathway-mapping
```

#### 步骤3：代谢网络重建
```bash
python3 scripts/metabolic_reconstructor.py --auto-reconstruct
```

#### 步骤4：网络质量验证
```bash
python3 scripts/metabolic_reconstructor.py --validate-network
```

#### 步骤5：营养需求预测
```bash
python3 scripts/nutrition_predictor.py --predict-requirements
```

#### 步骤6：氨基酸合成分析
```bash
python3 scripts/nutrition_predictor.py --amino-acid-synthesis
python3 scripts/nutrition_predictor.py --generate-aa-report
```

#### 步骤7：网络可视化
```bash
python3 scripts/network_visualizer.py --visualize-network
python3 scripts/network_visualizer.py --pathway-heatmap
python3 scripts/network_visualizer.py --nutrition-summary
```

### 4. 结果文件

分析完成后，将在以下目录生成结果：

- `results/`: 分析数据和报告
  - `enzyme_genes.txt`: 提取的酶基因列表
  - `pathway_analysis/`: KEGG通路分析结果
  - `metabolic_network.json`: 重建的代谢网络
  - `nutrition_prediction.txt`: 营养需求预测报告
  
- `figures/`: 可视化图表
  - `metabolic_network.png`: 代谢网络图
  - `pathway_heatmap.png`: 通路完整性热图
  - `nutrition_summary.png`: 营养需求总结图

## 脚本功能说明

### kegg_pathway_analyzer.py
- 从基因组注释中提取酶基因
- 映射到KEGG代谢通路
- 评估通路完整性

### metabolic_reconstructor.py  
- 基于酶基因重建代谢网络
- 验证网络质量
- 识别代谢间隙
- 导出SBML格式模型

### nutrition_predictor.py
- 预测碳源利用能力
- 分析氨基酸合成能力
- 评估维生素需求
- 生成营养需求报告

### network_visualizer.py
- 生成代谢网络可视化
- 创建通路完整性热图
- 制作营养需求总结图
- 绘制特定通路详图

## 注意事项

1. **KEGG API限制**: 脚本包含延时机制避免API限制，如遇问题可使用 `--slow-mode` 参数

2. **数据质量**: 分析结果依赖于基因组注释质量，建议使用高质量的注释数据

3. **网络连接**: 部分功能需要访问KEGG数据库，确保网络连接正常

4. **计算资源**: 大型基因组的网络重建可能需要较多内存和时间

## 扩展练习

1. **参数优化**: 尝试不同的分析参数，观察对结果的影响
2. **多菌株比较**: 分析不同微生物的代谢网络差异
3. **实验验证**: 设计实验验证营养需求预测结果
4. **方法改进**: 基于文献改进预测算法

## 技术支持

如遇到技术问题，请：
1. 检查Python环境和依赖包
2. 查看脚本运行日志
3. 参考实践操作指南中的故障排除部分
4. 联系课程老师或同学

---

**课程目标**: 通过本次实践，学生将掌握微生物代谢网络重建的完整流程，理解从基因组到代谢型的预测方法，并能够应用这些方法分析微生物的营养需求和代谢能力。