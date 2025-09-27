# 第3次课：基因预测与功能注释系统

## 课程概述

本次课程深入探讨微生物基因组的基因预测和功能注释方法，包括原核基因预测优化、特殊基因预测、非编码RNA预测和自动注释流程设计。

## 文件结构

```
第3次课-基因预测与功能注释系统/
├── README.md                                    # 本文件
├── 理论课件-基因预测与功能注释系统.md              # Marp格式理论课件
├── 实践操作-基因注释与非编码RNA分析.md             # 实践操作指南
├── images/                                      # 课件图片资源
│   ├── start_codon_usage.svg                   # 起始密码子使用频率图
│   ├── signal_peptide_structure.svg            # 信号肽结构图
│   ├── trna_structure.svg                      # tRNA二级结构图
│   ├── riboswitch_mechanism.svg                # 核糖开关机制图
│   └── annotation_pipeline.svg                 # 注释流水线图
└── scripts/                                     # Python分析脚本
    ├── annotation_comparison.py                 # 注释结果比较脚本
    ├── ncrna_analyzer.py                       # 非编码RNA分析脚本
    ├── quality_assessor.py                     # 注释质量评估脚本
    └── visualization_generator.py               # 可视化生成脚本
```

## 理论课件内容

### 主要知识点
1. **原核基因预测优化**
   - 起始密码子识别
   - 开放阅读框重叠处理
   - 基因模型训练与优化
   - 假基因识别方法

2. **特殊基因预测**
   - 信号肽预测算法
   - 跨膜蛋白识别
   - 分泌蛋白预测
   - 抗菌肽基因识别

3. **非编码RNA预测**
   - tRNA基因结构与预测
   - rRNA基因识别挑战
   - 小RNA（sRNA）预测
   - 核糖开关识别

4. **自动注释流程设计**
   - 注释流水线构建
   - 数据库选择策略
   - 结果整合方法
   - 质量控制标准

## 实践操作内容

### 操作目标
- 掌握多种基因注释工具的使用和结果比较
- 理解基因预测和功能注释的质量评估方法
- 应用非编码RNA预测工具进行全面基因组注释
- 分析注释结果的一致性和差异性，进行手工修正

### 主要步骤
1. **多工具基因注释比较** (45分钟)
   - Prokka快速注释
   - Prodigal基因预测
   - 注释结果比较

2. **非编码RNA预测** (45分钟)
   - tRNA基因预测
   - rRNA基因预测
   - 小RNA预测
   - 非编码RNA分析整合

3. **注释质量评估与优化** (30分钟)
   - 注释完整性评估
   - 与参考注释比较
   - 手工注释修正示例

## 脚本功能说明

### annotation_comparison.py
- **功能**：比较不同注释工具的基因预测结果
- **输入**：Prokka和Prodigal的GFF文件
- **输出**：比较统计报告和可视化图表
- **主要分析**：基因数量统计、重叠分析、差异基因识别

### ncrna_analyzer.py
- **功能**：整合和分析非编码RNA预测结果
- **输入**：tRNAscan-SE、RNAmmer、Infernal结果文件
- **输出**：ncRNA分类统计和分布图
- **主要分析**：tRNA类型统计、rRNA识别、小RNA功能分类

### quality_assessor.py
- **功能**：评估基因注释的质量指标
- **输入**：基因组序列和注释文件
- **输出**：质量评估报告和综合评分
- **主要指标**：基因密度、编码密度、注释覆盖率、链偏好性

### visualization_generator.py
- **功能**：生成各种可视化图表
- **输出**：基因组注释地图、功能分类图、质量仪表板
- **图表类型**：柱状图、饼图、热图、雷达图

## 使用说明

### 环境要求
- Python 3.7+
- Prokka 1.14.6+
- Prodigal 2.6.3+
- tRNAscan-SE 2.0.9+
- RNAmmer 1.2+
- Infernal 1.1.4+

### 快速开始
1. 准备工作环境和数据
2. 运行基因注释工具
3. 执行Python分析脚本
4. 查看结果报告和图表

### 数据文件
- **示例基因组**：大肠杆菌K-12 MG1655
- **参考注释**：NCBI官方注释结果
- **数据库**：Rfam、UniProt、KEGG等

## 学习成果

完成本次课程后，学生将能够：
- 熟练使用主流基因注释工具
- 理解不同预测方法的优缺点
- 评估注释质量并进行优化
- 设计完整的注释流水线

## 扩展练习

1. **基础扩展**
   - 尝试不同参数设置
   - 比较多个物种基因组
   - 测试其他注释工具

2. **进阶挑战**
   - 开发自定义注释流程
   - 整合多种证据来源
   - 实现批量处理功能

3. **创新探索**
   - 机器学习方法应用
   - 实验数据整合
   - 新工具性能评估

## 参考资料

### 核心文献
1. Seemann T. Prokka: rapid prokaryotic genome annotation. *Bioinformatics*. 2014
2. Hyatt D, et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics*. 2010
3. Lowe TM, Eddy SR. tRNAscan-SE: a program for improved detection of transfer RNA genes. *Nucleic Acids Res*. 2007

### 在线资源
- Prokka GitHub: https://github.com/tseemann/prokka
- Rfam数据库: http://rfam.xfam.org/
- NCBI基因组: https://www.ncbi.nlm.nih.gov/genome/

## 技术支持

如遇到问题，请：
1. 检查软件版本和依赖
2. 查看错误日志
3. 参考官方文档
4. 联系课程助教

---

*本课程材料遵循开源协议，欢迎改进和分享*