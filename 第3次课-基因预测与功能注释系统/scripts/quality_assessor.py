#!/usr/bin/env python3
"""
基因注释质量评估脚本
评估基因预测和功能注释的质量指标
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import re
from Bio import SeqIO

def parse_genome_fasta(fasta_file):
    """
    解析基因组FASTA文件，获取基本信息
    """
    genome_info = {
        'total_length': 0,
        'gc_content': 0,
        'n_count': 0,
        'contigs': 0
    }
    
    try:
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq))
            genome_info['contigs'] += 1
        
        # 合并所有序列
        full_sequence = ''.join(sequences)
        genome_info['total_length'] = len(full_sequence)
        
        # 计算GC含量
        gc_count = full_sequence.count('G') + full_sequence.count('C')
        genome_info['gc_content'] = gc_count / len(full_sequence) * 100
        
        # 计算N的数量
        genome_info['n_count'] = full_sequence.count('N')
        
    except FileNotFoundError:
        print(f"警告: 基因组文件 {fasta_file} 未找到")
        return None
    
    return genome_info

def parse_prokka_gff(gff_file):
    """
    解析Prokka GFF文件，提取注释信息
    """
    features = defaultdict(list)
    gene_info = []
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 9:
                    feature_type = parts[2]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    attributes = parts[8]
                    
                    features[feature_type].append({
                        'start': start,
                        'end': end,
                        'length': end - start + 1,
                        'strand': strand,
                        'attributes': attributes
                    })
                    
                    # 提取基因功能信息
                    if feature_type == 'CDS':
                        product = "hypothetical protein"
                        if "product=" in attributes:
                            product_match = re.search(r'product=([^;]+)', attributes)
                            if product_match:
                                product = product_match.group(1)
                        
                        gene_info.append({
                            'start': start,
                            'end': end,
                            'length': end - start + 1,
                            'strand': strand,
                            'product': product
                        })
    
    except FileNotFoundError:
        print(f"警告: GFF文件 {gff_file} 未找到")
        return {}, []
    
    return features, gene_info

def calculate_quality_metrics(genome_info, features, gene_info):
    """
    计算注释质量指标
    """
    metrics = {}
    
    if not genome_info:
        return metrics
    
    # 基本统计
    total_genes = len(gene_info)
    total_cds = len(features.get('CDS', []))
    total_trna = len(features.get('tRNA', []))
    total_rrna = len(features.get('rRNA', []))
    
    metrics['total_genes'] = total_genes
    metrics['total_cds'] = total_cds
    metrics['total_trna'] = total_trna
    metrics['total_rrna'] = total_rrna
    
    # 基因密度
    genome_length = genome_info['total_length']
    metrics['gene_density'] = total_genes / (genome_length / 1000)  # 每kb基因数
    
    # 编码密度
    total_coding_length = sum([gene['length'] for gene in gene_info])
    metrics['coding_density'] = total_coding_length / genome_length * 100  # 编码序列占比
    
    # 平均基因长度
    if gene_info:
        metrics['avg_gene_length'] = sum([gene['length'] for gene in gene_info]) / len(gene_info)
        metrics['median_gene_length'] = sorted([gene['length'] for gene in gene_info])[len(gene_info)//2]
    
    # 功能注释统计
    hypothetical_count = 0
    annotated_count = 0
    
    for gene in gene_info:
        product = gene['product'].lower()
        if 'hypothetical' in product or 'unknown' in product or 'putative' in product:
            hypothetical_count += 1
        else:
            annotated_count += 1
    
    metrics['hypothetical_proteins'] = hypothetical_count
    metrics['annotated_proteins'] = annotated_count
    metrics['annotation_coverage'] = annotated_count / total_genes * 100 if total_genes > 0 else 0
    
    # 基因长度分布
    gene_lengths = [gene['length'] for gene in gene_info]
    metrics['short_genes'] = len([l for l in gene_lengths if l < 300])  # <100 aa
    metrics['medium_genes'] = len([l for l in gene_lengths if 300 <= l <= 1500])  # 100-500 aa
    metrics['long_genes'] = len([l for l in gene_lengths if l > 1500])  # >500 aa
    
    # 链偏好性
    plus_strand = len([gene for gene in gene_info if gene['strand'] == '+'])
    minus_strand = len([gene for gene in gene_info if gene['strand'] == '-'])
    metrics['strand_bias'] = abs(plus_strand - minus_strand) / total_genes * 100 if total_genes > 0 else 0
    
    return metrics

def assess_annotation_quality():
    """
    主要质量评估函数
    """
    print("开始基因注释质量评估...")
    
    # 解析输入文件
    print("解析基因组文件...")
    genome_info = parse_genome_fasta("ecoli_genome.fna")
    
    print("解析注释文件...")
    features, gene_info = parse_prokka_gff("prokka_results/ecoli_prokka.gff")
    
    if not genome_info or not gene_info:
        print("错误: 无法读取必要的输入文件")
        return
    
    # 计算质量指标
    print("计算质量指标...")
    metrics = calculate_quality_metrics(genome_info, features, gene_info)
    
    # 生成质量评估报告
    os.makedirs('results', exist_ok=True)
    
    with open('results/quality_assessment.txt', 'w') as f:
        f.write("基因注释质量评估报告\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("1. 基因组基本信息\n")
        f.write("-" * 30 + "\n")
        f.write(f"基因组长度: {genome_info['total_length']:,} bp\n")
        f.write(f"GC含量: {genome_info['gc_content']:.2f}%\n")
        f.write(f"序列数量: {genome_info['contigs']}\n")
        f.write(f"N碱基数量: {genome_info['n_count']}\n")
        
        f.write(f"\n2. 基因预测统计\n")
        f.write("-" * 30 + "\n")
        f.write(f"总基因数: {metrics['total_genes']}\n")
        f.write(f"编码基因(CDS): {metrics['total_cds']}\n")
        f.write(f"tRNA基因: {metrics['total_trna']}\n")
        f.write(f"rRNA基因: {metrics['total_rrna']}\n")
        f.write(f"基因密度: {metrics['gene_density']:.2f} 基因/kb\n")
        f.write(f"编码密度: {metrics['coding_density']:.2f}%\n")
        
        f.write(f"\n3. 基因长度分析\n")
        f.write("-" * 30 + "\n")
        f.write(f"平均基因长度: {metrics.get('avg_gene_length', 0):.0f} bp\n")
        f.write(f"中位基因长度: {metrics.get('median_gene_length', 0):.0f} bp\n")
        f.write(f"短基因(<300bp): {metrics['short_genes']} ({metrics['short_genes']/metrics['total_genes']*100:.1f}%)\n")
        f.write(f"中等基因(300-1500bp): {metrics['medium_genes']} ({metrics['medium_genes']/metrics['total_genes']*100:.1f}%)\n")
        f.write(f"长基因(>1500bp): {metrics['long_genes']} ({metrics['long_genes']/metrics['total_genes']*100:.1f}%)\n")
        
        f.write(f"\n4. 功能注释质量\n")
        f.write("-" * 30 + "\n")
        f.write(f"有功能注释: {metrics['annotated_proteins']} ({metrics['annotation_coverage']:.1f}%)\n")
        f.write(f"假设蛋白: {metrics['hypothetical_proteins']} ({100-metrics['annotation_coverage']:.1f}%)\n")
        
        f.write(f"\n5. 其他质量指标\n")
        f.write("-" * 30 + "\n")
        f.write(f"链偏好性: {metrics['strand_bias']:.2f}%\n")
        
        f.write(f"\n6. 质量评价\n")
        f.write("-" * 20 + "\n")
        
        # 质量评价标准
        quality_score = 0
        comments = []
        
        # 基因密度评价（大肠杆菌正常范围：0.9-1.1 基因/kb）
        if 0.8 <= metrics['gene_density'] <= 1.2:
            quality_score += 20
            comments.append("✓ 基因密度正常")
        else:
            comments.append("⚠ 基因密度异常")
        
        # 编码密度评价（大肠杆菌正常范围：85-95%）
        if 80 <= metrics['coding_density'] <= 95:
            quality_score += 20
            comments.append("✓ 编码密度正常")
        else:
            comments.append("⚠ 编码密度异常")
        
        # 功能注释覆盖率评价
        if metrics['annotation_coverage'] >= 70:
            quality_score += 25
            comments.append("✓ 功能注释覆盖率良好")
        elif metrics['annotation_coverage'] >= 50:
            quality_score += 15
            comments.append("△ 功能注释覆盖率中等")
        else:
            comments.append("✗ 功能注释覆盖率较低")
        
        # tRNA数量评价（大肠杆菌约86个）
        if 80 <= metrics['total_trna'] <= 95:
            quality_score += 15
            comments.append("✓ tRNA数量正常")
        else:
            comments.append("⚠ tRNA数量异常")
        
        # rRNA数量评价（大肠杆菌约7个）
        if 5 <= metrics['total_rrna'] <= 10:
            quality_score += 10
            comments.append("✓ rRNA数量正常")
        else:
            comments.append("⚠ rRNA数量异常")
        
        # 链偏好性评价
        if metrics['strand_bias'] <= 10:
            quality_score += 10
            comments.append("✓ 链分布均衡")
        else:
            comments.append("⚠ 链分布不均衡")
        
        f.write(f"总体质量得分: {quality_score}/100\n")
        f.write("质量评价详情:\n")
        for comment in comments:
            f.write(f"  {comment}\n")
        
        # 质量等级
        if quality_score >= 80:
            quality_level = "优秀"
        elif quality_score >= 60:
            quality_level = "良好"
        elif quality_score >= 40:
            quality_level = "中等"
        else:
            quality_level = "需要改进"
        
        f.write(f"\n质量等级: {quality_level}\n")
    
    # 生成可视化图表
    create_quality_plots(genome_info, metrics, gene_info)
    
    print(f"\n质量评估完成！")
    print(f"详细报告: results/quality_assessment.txt")
    print(f"可视化图表: figures/目录")

def create_quality_plots(genome_info, metrics, gene_info):
    """
    创建质量评估可视化图表
    """
    os.makedirs('figures', exist_ok=True)
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 1. 基因长度分布直方图
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    gene_lengths = [gene['length'] for gene in gene_info]
    ax1.hist(gene_lengths, bins=50, color='#3498db', alpha=0.7, edgecolor='black')
    ax1.set_xlabel('基因长度 (bp)')
    ax1.set_ylabel('频数')
    ax1.set_title('基因长度分布')
    ax1.axvline(metrics.get('avg_gene_length', 0), color='red', linestyle='--', 
               label=f'平均长度: {metrics.get("avg_gene_length", 0):.0f} bp')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # 2. 功能注释饼图
    annotation_data = [metrics['annotated_proteins'], metrics['hypothetical_proteins']]
    annotation_labels = ['有功能注释', '假设蛋白']
    colors = ['#2ecc71', '#e74c3c']
    
    wedges, texts, autotexts = ax2.pie(annotation_data, labels=annotation_labels, 
                                      colors=colors, autopct='%1.1f%%', 
                                      startangle=90)
    ax2.set_title('功能注释覆盖率')
    
    # 3. 基因类型统计
    gene_types = ['CDS', 'tRNA', 'rRNA']
    gene_counts = [metrics['total_cds'], metrics['total_trna'], metrics['total_rrna']]
    colors3 = ['#3498db', '#f39c12', '#9b59b6']
    
    bars = ax3.bar(gene_types, gene_counts, color=colors3, alpha=0.7)
    ax3.set_ylabel('数量')
    ax3.set_title('基因类型分布')
    ax3.grid(axis='y', alpha=0.3)
    
    # 添加数值标签
    for bar, count in zip(bars, gene_counts):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(gene_counts)*0.01,
                str(count), ha='center', va='bottom')
    
    # 4. 质量指标雷达图
    categories = ['基因密度', '编码密度', '注释覆盖率', 'tRNA完整性', 'rRNA完整性']
    
    # 标准化指标到0-100范围
    values = [
        min(metrics['gene_density'] / 1.0 * 100, 100),  # 基因密度标准化
        metrics['coding_density'],  # 编码密度已经是百分比
        metrics['annotation_coverage'],  # 注释覆盖率已经是百分比
        min(metrics['total_trna'] / 86 * 100, 100),  # tRNA完整性（以86为标准）
        min(metrics['total_rrna'] / 7 * 100, 100)   # rRNA完整性（以7为标准）
    ]
    
    # 创建雷达图
    angles = [n / float(len(categories)) * 2 * 3.14159 for n in range(len(categories))]
    angles += angles[:1]  # 闭合图形
    values += values[:1]
    
    ax4.plot(angles, values, 'o-', linewidth=2, color='#3498db')
    ax4.fill(angles, values, alpha=0.25, color='#3498db')
    ax4.set_xticks(angles[:-1])
    ax4.set_xticklabels(categories)
    ax4.set_ylim(0, 100)
    ax4.set_title('注释质量综合评估')
    ax4.grid(True)
    
    plt.tight_layout()
    plt.savefig('figures/quality_assessment.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. 基因组特征概览
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 创建基因组特征汇总表
    summary_data = {
        '指标': ['基因组长度', 'GC含量', '总基因数', '基因密度', '编码密度', '注释覆盖率'],
        '数值': [f"{genome_info['total_length']:,} bp", 
                f"{genome_info['gc_content']:.2f}%",
                f"{metrics['total_genes']:,}",
                f"{metrics['gene_density']:.2f} 基因/kb",
                f"{metrics['coding_density']:.2f}%",
                f"{metrics['annotation_coverage']:.1f}%"],
        '评价': ['正常', '正常', '正常', '正常', '正常', '良好']
    }
    
    # 创建表格
    table_data = list(zip(summary_data['指标'], summary_data['数值'], summary_data['评价']))
    table = ax.table(cellText=table_data, 
                    colLabels=['指标', '数值', '评价'],
                    cellLoc='center',
                    loc='center',
                    colWidths=[0.3, 0.4, 0.3])
    
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 2)
    
    # 设置表格样式
    for i in range(len(summary_data['指标']) + 1):
        for j in range(3):
            cell = table[(i, j)]
            if i == 0:  # 表头
                cell.set_facecolor('#3498db')
                cell.set_text_props(weight='bold', color='white')
            else:
                if j == 2:  # 评价列
                    if '正常' in summary_data['评价'][i-1] or '良好' in summary_data['评价'][i-1]:
                        cell.set_facecolor('#d5f4e6')
                    else:
                        cell.set_facecolor('#ffeaa7')
    
    ax.axis('off')
    ax.set_title('基因组注释质量汇总', fontsize=16, fontweight='bold', pad=20)
    
    plt.savefig('figures/quality_summary.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    assess_annotation_quality()