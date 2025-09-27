#!/usr/bin/env python3
"""
基因注释结果比较脚本
比较不同工具的基因预测结果，分析差异和重叠情况
"""

import os
import sys
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

def parse_gff(gff_file, feature_type="CDS"):
    """
    解析GFF文件，提取指定类型的特征
    """
    features = []
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 9 and parts[2] == feature_type:
                    features.append({
                        'seqid': parts[0],
                        'source': parts[1],
                        'type': parts[2],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'score': parts[5],
                        'strand': parts[6],
                        'phase': parts[7],
                        'attributes': parts[8]
                    })
    except FileNotFoundError:
        print(f"警告: 文件 {gff_file} 未找到")
        return []
    
    return features

def calculate_overlap(gene1, gene2, min_overlap=0.5):
    """
    计算两个基因的重叠程度
    """
    start1, end1 = gene1['start'], gene1['end']
    start2, end2 = gene2['start'], gene2['end']
    
    # 计算重叠区域
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    if overlap_start >= overlap_end:
        return 0.0
    
    overlap_length = overlap_end - overlap_start
    gene1_length = end1 - start1
    gene2_length = end2 - start2
    
    # 计算重叠比例（相对于较短基因）
    min_length = min(gene1_length, gene2_length)
    overlap_ratio = overlap_length / min_length if min_length > 0 else 0.0
    
    return overlap_ratio

def compare_annotations():
    """
    主要比较函数
    """
    print("开始基因注释结果比较分析...")
    
    # 检查输入文件
    prokka_gff = "prokka_results/ecoli_prokka.gff"
    prodigal_gff = "prodigal_results.gff"
    reference_gff = "ecoli_reference.gff"
    
    # 解析GFF文件
    print("解析GFF文件...")
    prokka_genes = parse_gff(prokka_gff)
    prodigal_genes = parse_gff(prodigal_gff)
    reference_genes = parse_gff(reference_gff)
    
    # 基本统计
    stats = {
        'Prokka': len(prokka_genes),
        'Prodigal': len(prodigal_genes),
        'Reference': len(reference_genes)
    }
    
    print(f"基因预测数量统计:")
    for tool, count in stats.items():
        print(f"  {tool}: {count} 个基因")
    
    # 计算重叠情况
    print("\n计算基因重叠情况...")
    
    # Prokka vs Prodigal
    prokka_prodigal_overlap = 0
    for pg in prokka_genes:
        for pdg in prodigal_genes:
            if pg['seqid'] == pdg['seqid'] and calculate_overlap(pg, pdg) >= 0.5:
                prokka_prodigal_overlap += 1
                break
    
    # Prokka vs Reference
    prokka_ref_overlap = 0
    for pg in prokka_genes:
        for rg in reference_genes:
            if pg['seqid'] == rg['seqid'] and calculate_overlap(pg, rg) >= 0.5:
                prokka_ref_overlap += 1
                break
    
    # Prodigal vs Reference
    prodigal_ref_overlap = 0
    for pdg in prodigal_genes:
        for rg in reference_genes:
            if pdg['seqid'] == rg['seqid'] and calculate_overlap(pdg, rg) >= 0.5:
                prodigal_ref_overlap += 1
                break
    
    # 生成比较报告
    os.makedirs('results', exist_ok=True)
    
    with open('results/annotation_comparison.txt', 'w') as f:
        f.write("基因注释工具比较报告\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("1. 基因预测数量统计\n")
        f.write("-" * 30 + "\n")
        for tool, count in stats.items():
            f.write(f"{tool:12}: {count:6} 个基因\n")
        
        f.write(f"\n2. 基因重叠分析（重叠阈值 >= 50%）\n")
        f.write("-" * 40 + "\n")
        f.write(f"Prokka vs Prodigal  : {prokka_prodigal_overlap:6} 个重叠基因 ({prokka_prodigal_overlap/min(stats['Prokka'], stats['Prodigal'])*100:.1f}%)\n")
        f.write(f"Prokka vs Reference : {prokka_ref_overlap:6} 个重叠基因 ({prokka_ref_overlap/min(stats['Prokka'], stats['Reference'])*100:.1f}%)\n")
        f.write(f"Prodigal vs Reference: {prodigal_ref_overlap:6} 个重叠基因 ({prodigal_ref_overlap/min(stats['Prodigal'], stats['Reference'])*100:.1f}%)\n")
        
        f.write(f"\n3. 差异分析\n")
        f.write("-" * 20 + "\n")
        f.write(f"Prokka 特有基因    : {stats['Prokka'] - prokka_prodigal_overlap:6} 个\n")
        f.write(f"Prodigal 特有基因  : {stats['Prodigal'] - prokka_prodigal_overlap:6} 个\n")
        f.write(f"参考注释缺失基因   : {stats['Reference'] - max(prokka_ref_overlap, prodigal_ref_overlap):6} 个\n")
    
    # 生成可视化图表
    create_comparison_plots(stats, prokka_prodigal_overlap, prokka_ref_overlap, prodigal_ref_overlap)
    
    print(f"\n比较分析完成！")
    print(f"详细报告已保存到: results/annotation_comparison.txt")
    print(f"可视化图表已保存到: figures/目录")

def create_comparison_plots(stats, pp_overlap, pr_overlap, pdr_overlap):
    """
    创建比较可视化图表
    """
    os.makedirs('figures', exist_ok=True)
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 1. 基因数量比较柱状图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    tools = list(stats.keys())
    counts = list(stats.values())
    colors = ['#3498db', '#e74c3c', '#2ecc71']
    
    bars = ax1.bar(tools, counts, color=colors, alpha=0.7)
    ax1.set_ylabel('基因数量')
    ax1.set_title('不同工具基因预测数量比较')
    ax1.grid(axis='y', alpha=0.3)
    
    # 添加数值标签
    for bar, count in zip(bars, counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10,
                str(count), ha='center', va='bottom')
    
    # 2. 重叠情况热图
    overlap_matrix = [
        [stats['Prokka'], pp_overlap, pr_overlap],
        [pp_overlap, stats['Prodigal'], pdr_overlap],
        [pr_overlap, pdr_overlap, stats['Reference']]
    ]
    
    overlap_df = pd.DataFrame(overlap_matrix, 
                             index=['Prokka', 'Prodigal', 'Reference'],
                             columns=['Prokka', 'Prodigal', 'Reference'])
    
    sns.heatmap(overlap_df, annot=True, fmt='d', cmap='Blues', ax=ax2)
    ax2.set_title('基因重叠矩阵')
    
    plt.tight_layout()
    plt.savefig('figures/annotation_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. 韦恩图样式的重叠分析
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # 计算重叠比例
    total_genes = max(stats.values())
    prokka_ratio = stats['Prokka'] / total_genes
    prodigal_ratio = stats['Prodigal'] / total_genes
    overlap_ratio = pp_overlap / total_genes
    
    categories = ['Prokka特有', '共同预测', 'Prodigal特有']
    sizes = [stats['Prokka'] - pp_overlap, pp_overlap, stats['Prodigal'] - pp_overlap]
    colors = ['#3498db', '#9b59b6', '#e74c3c']
    
    wedges, texts, autotexts = ax.pie(sizes, labels=categories, colors=colors, 
                                     autopct='%1.1f%%', startangle=90)
    ax.set_title('Prokka vs Prodigal 基因预测重叠分析')
    
    plt.savefig('figures/overlap_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    compare_annotations()