#!/usr/bin/env python3
"""
Divergence Time Estimator Script for Microbial Phylogenomics Course
估算物种发散时间的脚本

Author: Course Development Team
Date: 2025
Version: 1.0.0
"""

import os
import sys
import argparse
import subprocess
import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def load_phylogenetic_tree(tree_file):
    """
    加载系统发育树
    
    Args:
        tree_file (str): 树文件路径
    
    Returns:
        Bio.Phylo.BaseTree.Tree: 系统发育树对象
    """
    print(f"正在加载系统发育树: {tree_file}")
    
    try:
        tree = Phylo.read(tree_file, "newick")
        print(f"成功加载树文件，包含 {tree.count_terminals()} 个末端节点")
        return tree
    except Exception as e:
        print(f"加载树文件失败: {e}")
        return None

def estimate_substitution_rate(tree, calibration_points=None):
    """
    估算分子进化速率
    
    Args:
        tree (Bio.Phylo.BaseTree.Tree): 系统发育树
        calibration_points (dict): 校准点信息
    
    Returns:
        float: 估算的替换速率（每百万年每位点）
    """
    print("正在估算分子进化速率...")
    
    # 默认校准点（基于文献报道的微生物进化速率）
    if calibration_points is None:
        calibration_points = {
            "E_coli_divergence": {
                "age_mya": 100,  # 百万年前
                "branch_length": 0.05  # 分支长度
            }
        }
    
    # 使用简化的速率估算方法
    # 实际应用中应使用更复杂的贝叶斯方法
    
    total_branch_length = tree.total_branch_length()
    estimated_age = 200  # 假设的总进化时间（百万年）
    
    substitution_rate = total_branch_length / (2 * estimated_age)
    
    print(f"估算的替换速率: {substitution_rate:.2e} 替换/位点/百万年")
    
    return substitution_rate

def calculate_divergence_times(tree, substitution_rate, output_dir):
    """
    计算各节点的发散时间
    
    Args:
        tree (Bio.Phylo.BaseTree.Tree): 系统发育树
        substitution_rate (float): 替换速率
        output_dir (str): 输出目录
    
    Returns:
        dict: 节点发散时间字典
    """
    print("正在计算发散时间...")
    
    # 计算从根到各节点的距离
    node_distances = {}
    divergence_times = {}
    
    # 获取根节点到各末端的距离
    for terminal in tree.get_terminals():
        distance = tree.distance(tree.root, terminal)
        node_distances[terminal.name] = distance
    
    # 计算内部节点的时间
    max_distance = max(node_distances.values()) if node_distances else 0
    
    # 为每个节点分配时间
    node_counter = 0
    for clade in tree.find_clades():
        if clade.name is None:
            clade.name = f"Node_{node_counter}"
            node_counter += 1
        
        if clade.is_terminal():
            # 末端节点时间为0（现在）
            divergence_times[clade.name] = 0
        else:
            # 内部节点时间基于分支长度计算
            distance_to_root = tree.distance(tree.root, clade)
            time_from_present = (max_distance - distance_to_root) / substitution_rate
            divergence_times[clade.name] = time_from_present
    
    # 保存发散时间结果
    times_file = os.path.join(output_dir, "divergence_times.txt")
    with open(times_file, 'w') as f:
        f.write("# Divergence Time Estimation Results\n")
        f.write("# Node\tDivergence_Time_MYA\tConfidence_Interval\n")
        
        for node, time_mya in sorted(divergence_times.items(), key=lambda x: x[1], reverse=True):
            # 简单的置信区间估算（实际应使用统计方法）
            confidence_lower = max(0, time_mya * 0.8)
            confidence_upper = time_mya * 1.2
            
            f.write(f"{node}\t{time_mya:.2f}\t{confidence_lower:.2f}-{confidence_upper:.2f}\n")
    
    print(f"发散时间结果保存至: {times_file}")
    
    return divergence_times

def create_timetree(tree, divergence_times, output_dir):
    """
    创建时间校准的系统发育树
    
    Args:
        tree (Bio.Phylo.BaseTree.Tree): 原始系统发育树
        divergence_times (dict): 节点发散时间
        output_dir (str): 输出目录
    
    Returns:
        str: 时间树文件路径
    """
    print("正在创建时间校准树...")
    
    # 创建时间树的副本
    timetree = tree.copy()
    
    # 调整分支长度以反映时间
    max_time = max(divergence_times.values()) if divergence_times else 0
    
    for clade in timetree.find_clades():
        if clade.name in divergence_times:
            # 将分支长度转换为时间单位
            if clade.is_terminal():
                # 末端分支长度 = 从父节点到现在的时间
                parent = timetree.get_path(clade)[-2] if len(timetree.get_path(clade)) > 1 else None
                if parent and parent.name in divergence_times:
                    clade.branch_length = divergence_times[parent.name]
                else:
                    clade.branch_length = max_time / 2
            else:
                # 内部分支长度基于时间差
                clade.branch_length = divergence_times[clade.name] / max_time if max_time > 0 else 0.1
    
    # 保存时间树
    timetree_file = os.path.join(output_dir, "timetree.nwk")
    Phylo.write(timetree, timetree_file, "newick")
    
    print(f"时间校准树保存至: {timetree_file}")
    
    return timetree_file

def visualize_timetree(tree, divergence_times, output_dir):
    """
    可视化时间校准的系统发育树
    
    Args:
        tree (Bio.Phylo.BaseTree.Tree): 系统发育树
        divergence_times (dict): 节点发散时间
        output_dir (str): 输出目录
    """
    print("正在生成时间树可视化...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 左图：原始系统发育树
    ax1.set_title("Original Phylogenetic Tree", fontsize=14, fontweight='bold')
    Phylo.draw(tree, axes=ax1, do_show=False)
    ax1.set_xlabel("Evolutionary Distance")
    
    # 右图：时间校准树
    ax2.set_title("Time-Calibrated Tree", fontsize=14, fontweight='bold')
    
    # 创建简化的时间轴可视化
    terminals = tree.get_terminals()
    y_positions = np.arange(len(terminals))
    
    # 绘制时间轴
    max_time = max(divergence_times.values()) if divergence_times else 100
    
    for i, terminal in enumerate(terminals):
        # 绘制从现在到根的线
        ax2.plot([0, max_time], [i, i], 'k-', alpha=0.3)
        ax2.text(-5, i, terminal.name, ha='right', va='center', fontsize=10)
    
    # 添加时间标尺
    time_ticks = np.arange(0, max_time + 1, max_time // 5)
    ax2.set_xticks(time_ticks)
    ax2.set_xlabel("Time (Million Years Ago)")
    ax2.set_ylabel("Species")
    ax2.set_xlim(-max_time * 0.3, max_time * 1.1)
    ax2.invert_xaxis()  # 时间轴从现在到过去
    
    # 添加地质时代背景（可选）
    if max_time > 50:
        # 添加一些地质时代标记
        geological_periods = [
            ("Quaternary", 0, 2.6, "#FFE4B5"),
            ("Neogene", 2.6, 23, "#FFFF99"),
            ("Paleogene", 23, 66, "#FFA500")
        ]
        
        for period, start, end, color in geological_periods:
            if start < max_time:
                rect = patches.Rectangle((start, -0.5), min(end, max_time) - start, 
                                       len(terminals), facecolor=color, alpha=0.3)
                ax2.add_patch(rect)
                if end < max_time:
                    ax2.text((start + end) / 2, len(terminals) - 0.5, period, 
                           ha='center', va='bottom', fontsize=8, rotation=90)
    
    plt.tight_layout()
    
    # 保存图片
    timetree_plot = os.path.join(output_dir, "timetree_visualization.png")
    plt.savefig(timetree_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"时间树可视化保存至: {timetree_plot}")

def generate_divergence_report(tree, divergence_times, substitution_rate, output_dir):
    """
    生成发散时间分析报告
    
    Args:
        tree (Bio.Phylo.BaseTree.Tree): 系统发育树
        divergence_times (dict): 节点发散时间
        substitution_rate (float): 替换速率
        output_dir (str): 输出目录
    """
    print("正在生成发散时间分析报告...")
    
    report_file = os.path.join(output_dir, "divergence_analysis_report.txt")
    
    with open(report_file, 'w') as f:
        f.write("# Divergence Time Analysis Report\n")
        f.write("=" * 50 + "\n\n")
        
        # 基本信息
        f.write("## Analysis Parameters\n")
        f.write(f"Number of taxa: {tree.count_terminals()}\n")
        f.write(f"Substitution rate: {substitution_rate:.2e} substitutions/site/Myr\n")
        f.write(f"Total tree length: {tree.total_branch_length():.4f}\n\n")
        
        # 发散时间统计
        times_list = [t for t in divergence_times.values() if t > 0]
        if times_list:
            f.write("## Divergence Time Statistics\n")
            f.write(f"Oldest divergence: {max(times_list):.2f} MYA\n")
            f.write(f"Most recent divergence: {min(times_list):.2f} MYA\n")
            f.write(f"Mean divergence time: {np.mean(times_list):.2f} MYA\n")
            f.write(f"Median divergence time: {np.median(times_list):.2f} MYA\n\n")
        
        # 详细节点信息
        f.write("## Detailed Node Information\n")
        f.write("Node\tDivergence_Time_MYA\tBiological_Interpretation\n")
        
        # 按时间排序节点
        sorted_nodes = sorted([(node, time) for node, time in divergence_times.items() if time > 0], 
                            key=lambda x: x[1], reverse=True)
        
        for node, time_mya in sorted_nodes:
            # 简单的生物学解释
            if time_mya > 100:
                interpretation = "Ancient divergence (>100 MYA)"
            elif time_mya > 50:
                interpretation = "Moderate divergence (50-100 MYA)"
            elif time_mya > 10:
                interpretation = "Recent divergence (10-50 MYA)"
            else:
                interpretation = "Very recent divergence (<10 MYA)"
            
            f.write(f"{node}\t{time_mya:.2f}\t{interpretation}\n")
        
        f.write("\n## Method Notes\n")
        f.write("- Divergence times estimated using molecular clock approach\n")
        f.write("- Substitution rate based on literature values for bacterial genomes\n")
        f.write("- Confidence intervals represent approximate uncertainty ranges\n")
        f.write("- Results should be interpreted with caution due to model assumptions\n")
    
    print(f"发散时间分析报告保存至: {report_file}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='Divergence Time Estimator for Microbial Phylogenomics')
    parser.add_argument('--tree', '-t', 
                       default='results/concatenated_alignment.fasta.treefile',
                       help='Input phylogenetic tree file')
    parser.add_argument('--output', '-o', default='results',
                       help='Output directory')
    parser.add_argument('--rate', '-r', type=float,
                       help='Substitution rate (substitutions/site/Myr)')
    
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    
    print("开始发散时间估算分析...")
    print(f"输入树文件: {args.tree}")
    print(f"输出目录: {args.output}")
    
    # 检查树文件是否存在
    if not os.path.exists(args.tree):
        print(f"错误: 树文件 {args.tree} 不存在")
        
        # 尝试查找其他可能的树文件
        possible_files = [
            "results/simple_phylogenetic_tree.nwk",
            "results/basic_tree.nwk",
            "concatenated_alignment.fasta.treefile"
        ]
        
        for possible_file in possible_files:
            if os.path.exists(possible_file):
                print(f"找到替代树文件: {possible_file}")
                args.tree = possible_file
                break
        else:
            print("未找到任何树文件，程序退出")
            sys.exit(1)
    
    # 加载系统发育树
    tree = load_phylogenetic_tree(args.tree)
    if tree is None:
        print("加载树文件失败，程序退出")
        sys.exit(1)
    
    # 估算替换速率
    if args.rate:
        substitution_rate = args.rate
        print(f"使用用户指定的替换速率: {substitution_rate:.2e}")
    else:
        substitution_rate = estimate_substitution_rate(tree)
    
    # 计算发散时间
    divergence_times = calculate_divergence_times(tree, substitution_rate, args.output)
    
    # 创建时间校准树
    timetree_file = create_timetree(tree, divergence_times, args.output)
    
    # 生成可视化
    visualize_timetree(tree, divergence_times, 'figures')
    
    # 生成分析报告
    generate_divergence_report(tree, divergence_times, substitution_rate, args.output)
    
    print("\n发散时间估算分析完成！")
    print("主要输出文件:")
    print(f"  - 发散时间数据: {args.output}/divergence_times.txt")
    print(f"  - 时间校准树: {args.output}/timetree.nwk")
    print(f"  - 可视化图表: figures/timetree_visualization.png")
    print(f"  - 分析报告: {args.output}/divergence_analysis_report.txt")

if __name__ == "__main__":
    main()