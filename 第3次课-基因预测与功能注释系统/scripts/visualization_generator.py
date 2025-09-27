#!/usr/bin/env python3
"""
可视化生成脚本
为基因注释结果生成各种可视化图表
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict, Counter
import matplotlib.patches as patches

def create_genome_annotation_map():
    """
    创建基因组注释地图
    """
    print("生成基因组注释地图...")
    
    # 这里创建一个示例的基因组注释可视化
    # 实际使用时需要解析真实的GFF文件
    
    fig, ax = plt.subplots(figsize=(15, 8))
    
    # 模拟基因组长度（大肠杆菌约4.6Mb）
    genome_length = 4600000
    
    # 模拟一些基因位置和类型
    np.random.seed(42)
    n_genes = 4300
    
    gene_starts = sorted(np.random.randint(0, genome_length, n_genes))
    gene_lengths = np.random.normal(1000, 300, n_genes)
    gene_lengths = np.clip(gene_lengths, 200, 5000)  # 限制基因长度范围
    
    # 基因类型分布
    gene_types = np.random.choice(['CDS', 'tRNA', 'rRNA', 'pseudogene'], 
                                 n_genes, p=[0.95, 0.03, 0.005, 0.015])
    
    # 颜色映射
    type_colors = {
        'CDS': '#3498db',
        'tRNA': '#e74c3c', 
        'rRNA': '#f39c12',
        'pseudogene': '#95a5a6'
    }
    
    # 绘制基因组主干
    ax.plot([0, genome_length], [0, 0], 'k-', linewidth=3, alpha=0.7)
    
    # 绘制基因
    y_positions = {'CDS': 0.1, 'tRNA': 0.2, 'rRNA': 0.3, 'pseudogene': -0.1}
    
    for start, length, gene_type in zip(gene_starts, gene_lengths, gene_types):
        end = start + length
        y_pos = y_positions[gene_type]
        
        # 绘制基因矩形
        rect = patches.Rectangle((start, y_pos-0.02), length, 0.04, 
                               facecolor=type_colors[gene_type], 
                               alpha=0.7, edgecolor='none')
        ax.add_patch(rect)
    
    # 设置坐标轴
    ax.set_xlim(0, genome_length)
    ax.set_ylim(-0.3, 0.4)
    ax.set_xlabel('基因组位置 (bp)')
    ax.set_ylabel('基因类型')
    
    # 设置y轴标签
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()))
    
    # 添加标题和图例
    ax.set_title('基因组注释概览图', fontsize=16, fontweight='bold')
    
    # 创建图例
    legend_elements = [patches.Patch(facecolor=color, label=gene_type) 
                      for gene_type, color in type_colors.items()]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # 添加刻度标记
    major_ticks = np.arange(0, genome_length + 1, 1000000)
    ax.set_xticks(major_ticks)
    ax.set_xticklabels([f'{int(x/1000000)}M' for x in major_ticks])
    
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/genome_annotation_map.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_functional_classification_chart():
    """
    创建功能分类图表
    """
    print("生成功能分类图表...")
    
    # 模拟COG功能分类数据
    cog_categories = {
        'J': ('翻译、核糖体结构和生物发生', 180),
        'A': ('RNA处理和修饰', 5),
        'K': ('转录', 310),
        'L': ('复制、重组和修复', 190),
        'B': ('染色质结构和动力学', 2),
        'D': ('细胞周期控制、细胞分裂、染色体分配', 45),
        'Y': ('核结构', 0),
        'V': ('防御机制', 85),
        'T': ('信号转导机制', 240),
        'M': ('细胞壁/膜/包膜生物发生', 230),
        'N': ('细胞运动', 120),
        'Z': ('细胞骨架', 8),
        'W': ('胞外结构', 15),
        'U': ('胞内运输、分泌和囊泡运输', 95),
        'O': ('翻译后修饰、蛋白质周转、分子伴侣', 180),
        'C': ('能量产生和转换', 250),
        'G': ('碳水化合物运输和代谢', 280),
        'E': ('氨基酸运输和代谢', 290),
        'F': ('核苷酸运输和代谢', 95),
        'H': ('辅酶运输和代谢', 140),
        'I': ('脂质运输和代谢', 120),
        'P': ('无机离子运输和代谢', 220),
        'Q': ('次级代谢物生物合成、运输和分解代谢', 60),
        'R': ('一般功能预测', 450),
        'S': ('功能未知', 320)
    }
    
    # 创建图表
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 1. COG功能分类柱状图
    categories = list(cog_categories.keys())
    descriptions = [cog_categories[cat][0] for cat in categories]
    counts = [cog_categories[cat][1] for cat in categories]
    
    # 按数量排序
    sorted_data = sorted(zip(categories, descriptions, counts), key=lambda x: x[2], reverse=True)
    top_15 = sorted_data[:15]  # 只显示前15个
    
    top_cats = [item[0] for item in top_15]
    top_descs = [item[1][:20] + '...' if len(item[1]) > 20 else item[1] for item in top_15]
    top_counts = [item[2] for item in top_15]
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(top_cats)))
    bars = ax1.barh(range(len(top_cats)), top_counts, color=colors, alpha=0.7)
    
    ax1.set_yticks(range(len(top_cats)))
    ax1.set_yticklabels([f'{cat}: {desc}' for cat, desc in zip(top_cats, top_descs)])
    ax1.set_xlabel('基因数量')
    ax1.set_title('COG功能分类统计（前15类）')
    ax1.grid(axis='x', alpha=0.3)
    
    # 添加数值标签
    for i, (bar, count) in enumerate(zip(bars, top_counts)):
        ax1.text(bar.get_width() + 5, bar.get_y() + bar.get_height()/2,
                str(count), ha='left', va='center')
    
    # 2. 功能分类饼图（按大类合并）
    metabolism_cats = ['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q']
    information_cats = ['J', 'A', 'K', 'L', 'B']
    cellular_cats = ['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O']
    unknown_cats = ['R', 'S']
    
    metabolism_count = sum([cog_categories[cat][1] for cat in metabolism_cats])
    information_count = sum([cog_categories[cat][1] for cat in information_cats])
    cellular_count = sum([cog_categories[cat][1] for cat in cellular_cats])
    unknown_count = sum([cog_categories[cat][1] for cat in unknown_cats])
    
    pie_data = [metabolism_count, information_count, cellular_count, unknown_count]
    pie_labels = ['代谢', '信息存储与处理', '细胞过程与信号', '功能预测/未知']
    pie_colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12']
    
    wedges, texts, autotexts = ax2.pie(pie_data, labels=pie_labels, colors=pie_colors,
                                      autopct='%1.1f%%', startangle=90)
    ax2.set_title('功能大类分布')
    
    plt.tight_layout()
    plt.savefig('figures/functional_classification.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_annotation_quality_dashboard():
    """
    创建注释质量仪表板
    """
    print("生成注释质量仪表板...")
    
    # 模拟质量指标数据
    quality_metrics = {
        '基因预测准确率': 95.2,
        '功能注释覆盖率': 78.5,
        'tRNA预测完整性': 98.8,
        'rRNA预测完整性': 100.0,
        '基因密度合理性': 92.1,
        '编码密度合理性': 89.7
    }
    
    fig = plt.figure(figsize=(16, 12))
    
    # 创建网格布局
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # 为每个指标创建仪表盘
    positions = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    
    for i, (metric, value) in enumerate(quality_metrics.items()):
        if i < len(positions):
            row, col = positions[i]
            ax = fig.add_subplot(gs[row, col], projection='polar')
            
            # 创建仪表盘
            theta = np.linspace(0, np.pi, 100)
            r = np.ones_like(theta)
            
            # 背景弧
            ax.plot(theta, r, color='lightgray', linewidth=10, alpha=0.3)
            
            # 数值弧
            value_theta = np.linspace(0, np.pi * value / 100, int(value))
            value_r = np.ones_like(value_theta)
            
            # 根据数值选择颜色
            if value >= 90:
                color = '#2ecc71'  # 绿色
            elif value >= 70:
                color = '#f39c12'  # 橙色
            else:
                color = '#e74c3c'  # 红色
            
            ax.plot(value_theta, value_r, color=color, linewidth=10)
            
            # 设置样式
            ax.set_ylim(0, 1.2)
            ax.set_theta_zero_location('W')
            ax.set_theta_direction(1)
            ax.set_thetagrids([0, 45, 90, 135, 180], ['0%', '25%', '50%', '75%', '100%'])
            ax.set_rgrids([])
            ax.set_title(f'{metric}\n{value:.1f}%', pad=20, fontsize=10, fontweight='bold')
            
            # 添加数值标签
            ax.text(np.pi/2, 0.5, f'{value:.1f}%', ha='center', va='center', 
                   fontsize=12, fontweight='bold', color=color)
    
    # 添加总体评分
    overall_score = np.mean(list(quality_metrics.values()))
    ax_overall = fig.add_subplot(gs[2, :])
    
    # 创建水平进度条
    bar_width = 0.3
    ax_overall.barh(0, overall_score, height=bar_width, color='#3498db', alpha=0.7)
    ax_overall.barh(0, 100, height=bar_width, color='lightgray', alpha=0.3)
    
    ax_overall.set_xlim(0, 100)
    ax_overall.set_ylim(-bar_width, bar_width)
    ax_overall.set_xlabel('质量得分')
    ax_overall.set_title(f'总体注释质量评分: {overall_score:.1f}/100', 
                        fontsize=14, fontweight='bold')
    
    # 添加评分等级
    if overall_score >= 90:
        grade = "优秀"
        grade_color = '#2ecc71'
    elif overall_score >= 80:
        grade = "良好"
        grade_color = '#3498db'
    elif overall_score >= 70:
        grade = "中等"
        grade_color = '#f39c12'
    else:
        grade = "需要改进"
        grade_color = '#e74c3c'
    
    ax_overall.text(overall_score + 2, 0, f'{grade}', ha='left', va='center',
                   fontsize=12, fontweight='bold', color=grade_color)
    
    # 移除y轴
    ax_overall.set_yticks([])
    ax_overall.spines['top'].set_visible(False)
    ax_overall.spines['right'].set_visible(False)
    ax_overall.spines['left'].set_visible(False)
    
    plt.suptitle('基因注释质量评估仪表板', fontsize=16, fontweight='bold')
    plt.savefig('figures/quality_dashboard.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_comparison_heatmap():
    """
    创建工具比较热图
    """
    print("生成工具比较热图...")
    
    # 模拟不同工具的性能数据
    tools = ['Prokka', 'Prodigal', 'GeneMark', 'Glimmer', 'RAST']
    metrics = ['速度', '准确率', '易用性', '功能注释', '特殊基因', '输出格式']
    
    # 评分矩阵（1-5分）
    scores = np.array([
        [5, 4, 5, 5, 4, 5],  # Prokka
        [5, 5, 4, 2, 3, 3],  # Prodigal
        [3, 5, 3, 3, 4, 4],  # GeneMark
        [4, 4, 3, 2, 3, 3],  # Glimmer
        [2, 4, 5, 5, 4, 5]   # RAST
    ])
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 创建热图
    im = ax.imshow(scores, cmap='RdYlGn', aspect='auto', vmin=1, vmax=5)
    
    # 设置坐标轴
    ax.set_xticks(np.arange(len(metrics)))
    ax.set_yticks(np.arange(len(tools)))
    ax.set_xticklabels(metrics)
    ax.set_yticklabels(tools)
    
    # 旋转x轴标签
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    # 添加数值标签
    for i in range(len(tools)):
        for j in range(len(metrics)):
            text = ax.text(j, i, scores[i, j], ha="center", va="center",
                          color="black", fontweight='bold')
    
    # 添加颜色条
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('评分 (1-5)', rotation=-90, va="bottom")
    
    ax.set_title("基因注释工具性能比较", fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('figures/tools_comparison_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()

def generate_all_visualizations():
    """
    生成所有可视化图表
    """
    print("开始生成可视化图表...")
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 创建输出目录
    os.makedirs('figures', exist_ok=True)
    
    # 生成各种图表
    create_genome_annotation_map()
    create_functional_classification_chart()
    create_annotation_quality_dashboard()
    create_comparison_heatmap()
    
    print("所有可视化图表生成完成！")
    print("图表保存位置: figures/目录")
    print("生成的图表:")
    print("  - genome_annotation_map.png: 基因组注释概览图")
    print("  - functional_classification.png: 功能分类统计图")
    print("  - quality_dashboard.png: 质量评估仪表板")
    print("  - tools_comparison_heatmap.png: 工具比较热图")

if __name__ == "__main__":
    generate_all_visualizations()