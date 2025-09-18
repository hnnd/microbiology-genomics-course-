# 软件环境配置

## 核心软件列表

### 基因组组装
- SPAdes
- Unicycler
- Flye

### 基因注释
- Prokka
- RAST
- antiSMASH

### 比较基因组学
- Roary
- Panaroo
- OrthoFinder

### 系统发育分析
- FastANI
- BEAST
- IQ-TREE

### 代谢分析
- KEGG Mapper
- ModelSEED
- Cytoscape

### 质量评估
- CheckM
- BUSCO
- QUAST

## 安装指南

### 方式1：Conda环境
```bash
# 创建课程专用环境
conda create -n microbial-genomics python=3.8
conda activate microbial-genomics

# 安装核心软件包
conda install -c bioconda spades prokka roary checkm-genome
```

### 方式2：Docker容器
```bash
# 拉取预配置镜像
docker pull microbial-genomics:latest

# 运行容器
docker run -it -v $(pwd):/workspace microbial-genomics:latest
```

### 方式3：云平台
- Galaxy平台使用指南
- RAST在线平台
- KBase平台教程

## 环境验证

提供软件安装验证脚本和测试数据。