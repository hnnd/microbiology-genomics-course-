#!/usr/bin/env python3
"""
非编码RNA分析脚本
整合tRNA、rRNA和小RNA预测结果，进行分类统计和可视化
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import re

def parse_trna_results(trna_file):
    """
    解析tRNAscan-SE结果文件
    """
    trnas = []
    
    try:
        with open(trna_file, 'r') as f:
            lines = f.readlines()
            
        # 跳过头部信息，找到结果开始行
        start_idx = 0
        for i, line in enumerate(lines):
            if line.startswith('Sequence') and 'tRNA' in line:
                start_idx = i + 2  # 跳过表头
                break
        
        for line in lines[start_idx:]:
            if line.strip() and not line.startswith('-'):
                parts = line.strip().split()
                if len(parts) >= 9:
                    trnas.append({
                        'sequence': parts[0],
                        'trna_num': parts[1],
                        'start': int(parts[2]),
                        'end': int(parts[3]),
                        'aa_type': parts[4],
                        'anticodon': parts[5],
                        'intron_start': parts[6] if parts[6] != '0' else None,
                        'intron_end': parts[7] if parts[7] != '0' else None,
                        'score': float(parts[8])
                    })
    
    except FileNotFoundError:
        print(f"警告: tRNA结果文件 {trna_file} 未找到")
        return []
    
    return trnas

def parse_rrna_results(rrna_file):
    """
    解析RNAmmer结果文件（GFF格式）
    """
    rrnas = []
    
    try:
        with open(rrna_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 9:
                    # 从attributes中提取rRNA类型
                    attributes = parts[8]
                    rrna_type = "unknown"
                    if "16s_rRNA" in attributes:
                        rrna_type = "16S"
                    elif "23s_rRNA" in attributes:
                        rrna_type = "23S"
                    elif "5s_rRNA" in attributes:
                        rrna_type = "5S"
                    
                    rrnas.append({
                        'seqid': parts[0],
                        'source': parts[1],
                        'type': rrna_type,
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'score': parts[5],
                        'strand': parts[6],
                        'attributes': attributes
                    })
    
    except FileNotFoundError:
        print(f"警告: rRNA结果文件 {rrna_file} 未找到")
        return []
    
    return rrnas

def parse_small_rna_results(srna_file):
    """
    解析Infernal小RNA预测结果
    """
    srnas = []
    
    try:
        with open(srna_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split()
                if len(parts) >= 15:
                    srnas.append({
                        'target_name': parts[0],
                        'accession': parts[1],
                        'query_name': parts[2],
                        'query_accession': parts[3],
                        'mdl': parts[4],
                        'mdl_from': int(parts[5]),
                        'mdl_to': int(parts[6]),
                        'seq_from': int(parts[7]),
                        'seq_to': int(parts[8]),
                        'strand': parts[9],
                        'trunc': parts[10],
                        'pass': parts[11],
                        'gc': float(parts[12]),
                        'bias': float(parts[13]),
                        'score': float(parts[14]),
                        'evalue': float(parts[15]) if len(parts) > 15 else 1.0
                    })
    
    except FileNotFoundError:
        print(f"警告: 小RNA结果文件 {srna_file} 未找到")
        return []
    
    return srnas

def classify_small_rnas(srnas):
    """
    对小RNA进行功能分类
    """
    classification = defaultdict(list)
    
    for srna in srnas:
        rna_name = srna['target_name'].lower()
        
        # 基于Rfam家族名称进行分类
        if any(x in rna_name for x in ['riboswitch', 'switch']):
            classification['核糖开关'].append(srna)
        elif any(x in rna_name for x in ['srp', 'signal']):
            classification['信号识别颗粒RNA'].append(srna)
        elif any(x in rna_name for x in ['rnase', 'ribonuclease']):
            classification['RNase RNA'].append(srna)
        elif any(x in rna_name for x in ['tmrna', 'ssra']):
            classification['tmRNA'].append(srna)
        elif any(x in rna_name for x in ['antisense', 'asrna']):
            classification['反义RNA'].append(srna)
        elif any(x in rna_name for x in ['mirna', 'microrna']):
            classification['microRNA'].append(srna)
        else:
            classification['其他调控RNA'].append(srna)
    
    return classification

def analyze_ncrna():
    """
    主要分析函数
    """
    print("开始非编码RNA分析...")
    
    # 解析各类非编码RNA结果
    print("解析tRNA预测结果...")
    trnas = parse_trna_results("trna_results.txt")
    
    print("解析rRNA预测结果...")
    rrnas = parse_rrna_results("rrna_results.gff")
    
    print("解析小RNA预测结果...")
    srnas = parse_small_rna_results("small_rna_filtered.tbl")
    
    # 统计分析
    print("进行统计分析...")
    
    # tRNA统计
    trna_types = Counter([trna['aa_type'] for trna in trnas])
    trna_with_introns = len([trna for trna in trnas if trna['intron_start']])
    
    # rRNA统计
    rrna_types = Counter([rrna['type'] for rrna in rrnas])
    
    # 小RNA分类
    srna_classification = classify_small_rnas(srnas)
    
    # 生成结果目录
    os.makedirs('results', exist_ok=True)
    
    # 生成汇总报告
    with open('results/ncrna_summary.txt', 'w') as f:
        f.write("非编码RNA预测结果汇总报告\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("1. tRNA分析结果\n")
        f.write("-" * 30 + "\n")
        f.write(f"总tRNA数量: {len(trnas)}\n")
        f.write(f"含内含子tRNA: {trna_with_introns}\n")
        f.write("氨基酸类型分布:\n")
        for aa_type, count in trna_types.most_common():
            f.write(f"  {aa_type}: {count}\n")
        
        f.write(f"\n2. rRNA分析结果\n")
        f.write("-" * 30 + "\n")
        f.write(f"总rRNA数量: {len(rrnas)}\n")
        f.write("rRNA类型分布:\n")
        for rrna_type, count in rrna_types.items():
            f.write(f"  {rrna_type}: {count}\n")
        
        f.write(f"\n3. 小RNA分析结果\n")
        f.write("-" * 30 + "\n")
        f.write(f"总小RNA数量: {len(srnas)}\n")
        f.write("功能分类:\n")
        for category, rnas in srna_classification.items():
            f.write(f"  {category}: {len(rnas)}\n")
        
        f.write(f"\n4. 总体统计\n")
        f.write("-" * 20 + "\n")
        total_ncrnas = len(trnas) + len(rrnas) + len(srnas)
        f.write(f"非编码RNA总数: {total_ncrnas}\n")
        f.write(f"  tRNA: {len(trnas)} ({len(trnas)/total_ncrnas*100:.1f}%)\n")
        f.write(f"  rRNA: {len(rrnas)} ({len(rrnas)/total_ncrnas*100:.1f}%)\n")
        f.write(f"  小RNA: {len(srnas)} ({len(srnas)/total_ncrnas*100:.1f}%)\n")
    
    # 生成详细的tRNA列表
    with open('results/ncrna_trna_details.txt', 'w') as f:
        f.write("tRNA详细信息\n")
        f.write("=" * 60 + "\n")
        f.write("序列名\t位置\t氨基酸\t反密码子\t得分\t内含子\n")
        f.write("-" * 60 + "\n")
        
        for trna in sorted(trnas, key=lambda x: (x['sequence'], x['start'])):
            intron_info = "是" if trna['intron_start'] else "否"
            f.write(f"{trna['sequence']}\t{trna['start']}-{trna['end']}\t"
                   f"{trna['aa_type']}\t{trna['anticodon']}\t"
                   f"{trna['score']:.1f}\t{intron_info}\n")
    
    # 生成可视化图表
    create_ncrna_plots(trnas, rrnas, srnas, trna_types, rrna_types, srna_classification)
    
    print(f"\n非编码RNA分析完成！")
    print(f"汇总报告: results/ncrna_summary.txt")
    print(f"tRNA详情: results/ncrna_trna_details.txt")
    print(f"可视化图表: figures/目录")

def create_ncrna_plots(trnas, rrnas, srnas, trna_types, rrna_types, srna_classification):
    """
    创建非编码RNA可视化图表
    """
    os.makedirs('figures', exist_ok=True)
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 1. 非编码RNA类型分布饼图
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # 总体分布
    ncrna_counts = [len(trnas), len(rrnas), len(srnas)]
    ncrna_labels = ['tRNA', 'rRNA', '小RNA']
    colors1 = ['#3498db', '#e74c3c', '#2ecc71']
    
    wedges, texts, autotexts = ax1.pie(ncrna_counts, labels=ncrna_labels, 
                                      colors=colors1, autopct='%1.1f%%', 
                                      startangle=90)
    ax1.set_title('非编码RNA类型分布')
    
    # tRNA氨基酸类型分布（前10个）
    top_trna_types = dict(trna_types.most_common(10))
    ax2.bar(range(len(top_trna_types)), list(top_trna_types.values()), 
           color='#3498db', alpha=0.7)
    ax2.set_xticks(range(len(top_trna_types)))
    ax2.set_xticklabels(list(top_trna_types.keys()), rotation=45)
    ax2.set_ylabel('数量')
    ax2.set_title('tRNA氨基酸类型分布（前10）')
    ax2.grid(axis='y', alpha=0.3)
    
    # rRNA类型分布
    if rrna_types:
        ax3.bar(rrna_types.keys(), rrna_types.values(), 
               color='#e74c3c', alpha=0.7)
        ax3.set_ylabel('数量')
        ax3.set_title('rRNA类型分布')
        ax3.grid(axis='y', alpha=0.3)
    else:
        ax3.text(0.5, 0.5, '未检测到rRNA', ha='center', va='center', 
                transform=ax3.transAxes, fontsize=14)
        ax3.set_title('rRNA类型分布')
    
    # 小RNA功能分类
    if srna_classification:
        srna_counts = [len(rnas) for rnas in srna_classification.values()]
        srna_labels = list(srna_classification.keys())
        
        # 只显示前8个类别，其余合并为"其他"
        if len(srna_labels) > 8:
            top_8_counts = srna_counts[:8]
            other_count = sum(srna_counts[8:])
            srna_counts = top_8_counts + [other_count]
            srna_labels = srna_labels[:8] + ['其他']
        
        colors4 = plt.cm.Set3(range(len(srna_labels)))
        wedges, texts, autotexts = ax4.pie(srna_counts, labels=srna_labels, 
                                          colors=colors4, autopct='%1.1f%%', 
                                          startangle=90)
        ax4.set_title('小RNA功能分类')
    else:
        ax4.text(0.5, 0.5, '未检测到小RNA', ha='center', va='center', 
                transform=ax4.transAxes, fontsize=14)
        ax4.set_title('小RNA功能分类')
    
    plt.tight_layout()
    plt.savefig('figures/ncrna_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. tRNA基因组分布图
    if trnas:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # 按位置排序tRNA
        sorted_trnas = sorted(trnas, key=lambda x: x['start'])
        positions = [trna['start'] for trna in sorted_trnas]
        aa_types = [trna['aa_type'] for trna in sorted_trnas]
        
        # 为不同氨基酸类型分配颜色
        unique_aa = list(set(aa_types))
        colors = plt.cm.tab20(range(len(unique_aa)))
        aa_color_map = dict(zip(unique_aa, colors))
        
        # 绘制散点图
        for i, (pos, aa) in enumerate(zip(positions, aa_types)):
            ax.scatter(pos, 1, c=[aa_color_map[aa]], s=50, alpha=0.7)
        
        ax.set_xlabel('基因组位置 (bp)')
        ax.set_ylabel('')
        ax.set_title('tRNA在基因组上的分布')
        ax.set_yticks([])
        ax.grid(axis='x', alpha=0.3)
        
        # 添加图例（只显示前10个氨基酸类型）
        legend_elements = []
        for aa in unique_aa[:10]:
            legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                            markerfacecolor=aa_color_map[aa], 
                                            markersize=8, label=aa))
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        plt.savefig('figures/trna_genome_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    analyze_ncrna()