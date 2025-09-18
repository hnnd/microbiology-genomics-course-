#!/usr/bin/env python3
"""
基因组序列特征分析
"""

from Bio import SeqIO
from Bio.SeqUtils import GC, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def analyze_genome_sequence(fasta_file, genome_name):
    """分析基因组序列特征"""
    print(f"\n=== 分析 {genome_name} ===")
    
    sequences = []
    total_length = 0
    gc_contents = []
    
    # 读取序列
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_len = len(record.seq)
        gc_content = GC(record.seq)
        
        sequences.append({
            'id': record.id,
            'length': seq_len,
            'gc_content': gc_content
        })
        
        total_length += seq_len
        gc_contents.append(gc_content)
    
    # 基本统计
    print(f"序列数量: {len(sequences)}")
    print(f"总长度: {total_length:,} bp ({total_length/1e6:.2f} Mb)")
    print(f"平均GC含量: {np.mean(gc_contents):.2f}%")
    print(f"GC含量标准差: {np.std(gc_contents):.2f}%")
    
    # 序列长度分布
    lengths = [s['length'] for s in sequences]
    if len(lengths) > 1:
        print(f"最长序列: {max(lengths):,} bp")
        print(f"最短序列: {min(lengths):,} bp")
    
    return {
        'name': genome_name,
        'total_length': total_length,
        'num_sequences': len(sequences),
        'mean_gc': np.mean(gc_contents),
        'std_gc': np.std(gc_contents),
        'sequences': sequences
    }

def calculate_coding_density(genome_info):
    """估算编码密度（简化版本）"""
    # 这里使用经验公式估算
    # 实际应用中需要基因注释信息
    total_length = genome_info['total_length']
    
    # 原核生物典型编码密度约90%
    estimated_coding_bp = total_length * 0.9
    estimated_genes = estimated_coding_bp / 1000  # 假设平均基因长度1kb
    
    print(f"估算编码区长度: {estimated_coding_bp:,.0f} bp")
    print(f"估算基因数量: {estimated_genes:,.0f}")
    print(f"估算编码密度: 90%")
    
    return estimated_coding_bp, estimated_genes

def compare_genomes(genome_files, genome_names):
    """比较多个基因组"""
    results = []
    
    for fasta_file, name in zip(genome_files, genome_names):
        if os.path.exists(fasta_file):
            result = analyze_genome_sequence(fasta_file, name)
            coding_bp, genes = calculate_coding_density(result)
            result['estimated_genes'] = genes
            result['estimated_coding_density'] = 0.9
            results.append(result)
        else:
            print(f"文件不存在: {fasta_file}")
    
    return results

def create_comparison_plot(results):
    """创建比较图表"""
    if len(results) < 2:
        print("需要至少2个基因组进行比较")
        return
    
    # 准备数据
    names = [r['name'] for r in results]
    sizes = [r['total_length']/1e6 for r in results]  # 转换为Mb
    gc_contents = [r['mean_gc'] for r in results]
    gene_counts = [r['estimated_genes'] for r in results]
    
    # 创建子图
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # 基因组大小比较
    bars1 = ax1.bar(names, sizes, color=['#1f77b4', '#ff7f0e', '#2ca02c'])
    ax1.set_ylabel('基因组大小 (Mb)')
    ax1.set_title('基因组大小比较')
    ax1.tick_params(axis='x', rotation=45)
    
    # 添加数值标签
    for bar, size in zip(bars1, sizes):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                f'{size:.2f}', ha='center', va='bottom')
    
    # GC含量比较
    bars2 = ax2.bar(names, gc_contents, color=['#d62728', '#9467bd', '#8c564b'])
    ax2.set_ylabel('GC含量 (%)')
    ax2.set_title('GC含量比较')
    ax2.tick_params(axis='x', rotation=45)
    
    for bar, gc in zip(bars2, gc_contents):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{gc:.1f}%', ha='center', va='bottom')
    
    # 基因数量比较
    bars3 = ax3.bar(names, gene_counts, color=['#17becf', '#bcbd22', '#e377c2'])
    ax3.set_ylabel('估算基因数量')
    ax3.set_title('基因数量比较')
    ax3.tick_params(axis='x', rotation=45)
    
    for bar, genes in zip(bars3, gene_counts):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50,
                f'{genes:.0f}', ha='center', va='bottom')
    
    # 大小vs GC含量散点图
    ax4.scatter(sizes, gc_contents, s=100, alpha=0.7)
    for i, name in enumerate(names):
        ax4.annotate(name, (sizes[i], gc_contents[i]), 
                    xytext=(5, 5), textcoords='offset points')
    ax4.set_xlabel('基因组大小 (Mb)')
    ax4.set_ylabel('GC含量 (%)')
    ax4.set_title('基因组大小 vs GC含量')
    
    plt.tight_layout()
    plt.savefig('figures/genome_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("比较图表已保存到 figures/genome_comparison.png")

if __name__ == "__main__":
    # 定义要分析的基因组
    genome_files = [
        'data/ecoli_k12.fna',
        'data/buchnera.fna', 
        'data/mtb_h37rv.fna'
    ]
    
    genome_names = [
        'E. coli K-12 (自由生活)',
        'Buchnera (专性共生)',
        'M. tuberculosis (病原菌)'
    ]
    
    # 比较分析
    results = compare_genomes(genome_files, genome_names)
    
    # 创建比较图表
    if results:
        create_comparison_plot(results)
        
        # 保存结果到CSV
        df_results = pd.DataFrame([
            {
                'Genome': r['name'],
                'Size_Mb': r['total_length']/1e6,
                'GC_Content': r['mean_gc'],
                'Estimated_Genes': r['estimated_genes'],
                'Sequences': r['num_sequences']
            }
            for r in results
        ])
        
        df_results.to_csv('results/genome_comparison.csv', index=False)
        print("详细结果已保存到 results/genome_comparison.csv")