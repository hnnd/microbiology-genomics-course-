#!/usr/bin/env python3
"""
ANI Calculator Script for Microbial Phylogenomics Course
计算基因组间平均核苷酸一致性(ANI)的脚本

Author: Course Development Team
Date: 2025
Version: 1.0.0
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

def run_fastani(genome_list, output_dir):
    """
    使用FastANI计算基因组间ANI值
    
    Args:
        genome_list (str): 基因组文件列表
        output_dir (str): 输出目录
    
    Returns:
        str: ANI结果文件路径
    """
    print("正在运行FastANI计算...")
    
    ani_output = os.path.join(output_dir, "fastani_results.txt")
    
    # 构建FastANI命令
    cmd = [
        "fastANI",
        "--ql", genome_list,
        "--rl", genome_list,
        "-o", ani_output,
        "--matrix"
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"FastANI计算完成，结果保存至: {ani_output}")
        return ani_output
    except subprocess.CalledProcessError as e:
        print(f"FastANI运行失败: {e}")
        print(f"错误输出: {e.stderr}")
        return None

def parse_fastani_output(fastani_file, output_dir):
    """
    解析FastANI输出并生成相似性矩阵
    
    Args:
        fastani_file (str): FastANI输出文件
        output_dir (str): 输出目录
    
    Returns:
        pd.DataFrame: ANI相似性矩阵
    """
    print("正在解析FastANI结果...")
    
    # 读取FastANI结果
    try:
        data = pd.read_csv(fastani_file, sep='\t', header=None,
                          names=['query', 'reference', 'ani', 'fragments', 'total_fragments'])
    except Exception as e:
        print(f"读取FastANI结果失败: {e}")
        return None
    
    # 提取基因组名称
    genomes = sorted(list(set(data['query'].tolist() + data['reference'].tolist())))
    genome_names = [os.path.basename(g).replace('.fna', '') for g in genomes]
    
    # 创建ANI矩阵
    ani_matrix = pd.DataFrame(index=genome_names, columns=genome_names, dtype=float)
    
    # 填充矩阵
    for _, row in data.iterrows():
        query_name = os.path.basename(row['query']).replace('.fna', '')
        ref_name = os.path.basename(row['reference']).replace('.fna', '')
        ani_matrix.loc[query_name, ref_name] = row['ani']
    
    # 填充对角线（自身比较为100%）
    for genome in genome_names:
        ani_matrix.loc[genome, genome] = 100.0
    
    # 对称填充矩阵
    for i in range(len(genome_names)):
        for j in range(len(genome_names)):
            if pd.isna(ani_matrix.iloc[i, j]) and not pd.isna(ani_matrix.iloc[j, i]):
                ani_matrix.iloc[i, j] = ani_matrix.iloc[j, i]
    
    # 保存矩阵
    matrix_file = os.path.join(output_dir, "ani_matrix.txt")
    ani_matrix.to_csv(matrix_file, sep='\t', float_format='%.2f')
    print(f"ANI矩阵保存至: {matrix_file}")
    
    return ani_matrix

def create_ani_heatmap(ani_matrix, output_dir):
    """
    创建ANI热图
    
    Args:
        ani_matrix (pd.DataFrame): ANI相似性矩阵
        output_dir (str): 输出目录
    """
    print("正在生成ANI热图...")
    
    plt.figure(figsize=(10, 8))
    
    # 创建热图
    mask = ani_matrix.isna()
    sns.heatmap(ani_matrix, 
                annot=True, 
                fmt='.1f',
                cmap='RdYlBu_r',
                vmin=80, 
                vmax=100,
                mask=mask,
                square=True,
                cbar_kws={'label': 'ANI (%)'})
    
    plt.title('Average Nucleotide Identity (ANI) Matrix', fontsize=14, fontweight='bold')
    plt.xlabel('Reference Genomes', fontsize=12)
    plt.ylabel('Query Genomes', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    # 保存图片
    heatmap_file = os.path.join(output_dir, "ani_heatmap.png")
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"ANI热图保存至: {heatmap_file}")

def perform_clustering(ani_matrix, output_dir):
    """
    基于ANI值进行层次聚类分析
    
    Args:
        ani_matrix (pd.DataFrame): ANI相似性矩阵
        output_dir (str): 输出目录
    """
    print("正在进行聚类分析...")
    
    # 转换为距离矩阵（1-ANI/100）
    distance_matrix = 1 - (ani_matrix / 100)
    distance_matrix = distance_matrix.fillna(1.0)  # 缺失值设为最大距离
    
    # 确保矩阵对称
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    
    # 转换为压缩距离矩阵
    condensed_distances = squareform(distance_matrix.values)
    
    # 执行层次聚类
    linkage_matrix = linkage(condensed_distances, method='average')
    
    # 绘制聚类树
    plt.figure(figsize=(12, 8))
    dendrogram(linkage_matrix, 
               labels=ani_matrix.index.tolist(),
               orientation='top',
               leaf_rotation=45,
               leaf_font_size=10)
    
    plt.title('Hierarchical Clustering Based on ANI Values', fontsize=14, fontweight='bold')
    plt.xlabel('Genomes', fontsize=12)
    plt.ylabel('Distance (1 - ANI/100)', fontsize=12)
    plt.tight_layout()
    
    # 保存聚类图
    cluster_file = os.path.join(output_dir, "ani_clustering.png")
    plt.savefig(cluster_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # 保存聚类结果
    cluster_result_file = os.path.join(output_dir, "ani_clustering.txt")
    with open(cluster_result_file, 'w') as f:
        f.write("# ANI-based Hierarchical Clustering Results\n")
        f.write("# Linkage Matrix (ward linkage method)\n")
        f.write("# Columns: cluster1, cluster2, distance, cluster_size\n")
        for i, row in enumerate(linkage_matrix):
            f.write(f"{i+1}\t{row[0]:.0f}\t{row[1]:.0f}\t{row[2]:.4f}\t{row[3]:.0f}\n")
    
    print(f"聚类分析完成，结果保存至: {cluster_file}")
    print(f"聚类数据保存至: {cluster_result_file}")

def generate_summary_report(ani_matrix, output_dir):
    """
    生成ANI分析总结报告
    
    Args:
        ani_matrix (pd.DataFrame): ANI相似性矩阵
        output_dir (str): 输出目录
    """
    print("正在生成总结报告...")
    
    report_file = os.path.join(output_dir, "ani_summary_report.txt")
    
    with open(report_file, 'w') as f:
        f.write("# ANI Analysis Summary Report\n")
        f.write("=" * 50 + "\n\n")
        
        # 基本统计信息
        f.write("## Basic Statistics\n")
        f.write(f"Number of genomes analyzed: {len(ani_matrix)}\n")
        f.write(f"Total pairwise comparisons: {len(ani_matrix) * (len(ani_matrix) - 1) // 2}\n\n")
        
        # ANI值分布
        ani_values = []
        for i in range(len(ani_matrix)):
            for j in range(i+1, len(ani_matrix)):
                if not pd.isna(ani_matrix.iloc[i, j]):
                    ani_values.append(ani_matrix.iloc[i, j])
        
        if ani_values:
            f.write("## ANI Value Distribution\n")
            f.write(f"Mean ANI: {np.mean(ani_values):.2f}%\n")
            f.write(f"Median ANI: {np.median(ani_values):.2f}%\n")
            f.write(f"Min ANI: {np.min(ani_values):.2f}%\n")
            f.write(f"Max ANI: {np.max(ani_values):.2f}%\n")
            f.write(f"Standard deviation: {np.std(ani_values):.2f}%\n\n")
        
        # 分类学解释
        f.write("## Taxonomic Interpretation\n")
        same_species = sum(1 for v in ani_values if v >= 95)
        related_species = sum(1 for v in ani_values if 85 <= v < 95)
        distant_species = sum(1 for v in ani_values if v < 85)
        
        f.write(f"Same species pairs (ANI ≥ 95%): {same_species}\n")
        f.write(f"Related species pairs (85% ≤ ANI < 95%): {related_species}\n")
        f.write(f"Distant species pairs (ANI < 85%): {distant_species}\n\n")
        
        # 详细比较结果
        f.write("## Detailed Pairwise Comparisons\n")
        f.write("Genome1\tGenome2\tANI(%)\tInterpretation\n")
        
        genomes = ani_matrix.index.tolist()
        for i in range(len(genomes)):
            for j in range(i+1, len(genomes)):
                ani_val = ani_matrix.iloc[i, j]
                if not pd.isna(ani_val):
                    if ani_val >= 95:
                        interpretation = "Same species"
                    elif ani_val >= 85:
                        interpretation = "Related species"
                    else:
                        interpretation = "Distant species"
                    
                    f.write(f"{genomes[i]}\t{genomes[j]}\t{ani_val:.2f}\t{interpretation}\n")
    
    print(f"总结报告保存至: {report_file}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='ANI Calculator for Microbial Genomes')
    parser.add_argument('--input', '-i', default='genome_list.txt',
                       help='Input file containing list of genome files')
    parser.add_argument('--output', '-o', default='results',
                       help='Output directory')
    parser.add_argument('--cluster', action='store_true',
                       help='Perform hierarchical clustering analysis')
    
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    
    print("开始ANI分析...")
    print(f"输入文件: {args.input}")
    print(f"输出目录: {args.output}")
    
    # 检查输入文件
    if not os.path.exists(args.input):
        print(f"错误: 输入文件 {args.input} 不存在")
        sys.exit(1)
    
    # 运行FastANI
    fastani_output = run_fastani(args.input, args.output)
    if not fastani_output:
        print("FastANI计算失败，程序退出")
        sys.exit(1)
    
    # 解析结果
    ani_matrix = parse_fastani_output(fastani_output, args.output)
    if ani_matrix is None:
        print("ANI结果解析失败，程序退出")
        sys.exit(1)
    
    # 生成热图
    create_ani_heatmap(ani_matrix, 'figures')
    
    # 执行聚类分析（如果指定）
    if args.cluster:
        perform_clustering(ani_matrix, 'figures')
    
    # 生成总结报告
    generate_summary_report(ani_matrix, args.output)
    
    print("\nANI分析完成！")
    print("主要输出文件:")
    print(f"  - ANI矩阵: {args.output}/ani_matrix.txt")
    print(f"  - ANI热图: figures/ani_heatmap.png")
    print(f"  - 总结报告: {args.output}/ani_summary_report.txt")
    if args.cluster:
        print(f"  - 聚类分析: figures/ani_clustering.png")

if __name__ == "__main__":
    main()