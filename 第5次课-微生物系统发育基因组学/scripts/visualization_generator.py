#!/usr/bin/env python3
"""
Visualization Generator Script for Microbial Phylogenomics Course
生成系统发育基因组学分析结果可视化图表的脚本

Author: Course Development Team
Date: 2025
Version: 1.0.0
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

def load_ani_matrix(ani_file):
    """
    加载ANI相似性矩阵
    
    Args:
        ani_file (str): ANI矩阵文件路径
    
    Returns:
        pd.DataFrame: ANI矩阵
    """
    try:
        ani_matrix = pd.read_csv(ani_file, sep='\t', index_col=0)
        print(f"成功加载ANI矩阵: {ani_matrix.shape}")
        return ani_matrix
    except Exception as e:
        print(f"加载ANI矩阵失败: {e}")
        return None

def create_comprehensive_ani_heatmap(ani_matrix, output_dir):
    """
    创建综合ANI热图
    
    Args:
        ani_matrix (pd.DataFrame): ANI相似性矩阵
        output_dir (str): 输出目录
    """
    print("正在生成综合ANI热图...")
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # 子图1: 基本ANI热图
    mask = ani_matrix.isna()
    im1 = sns.heatmap(ani_matrix, 
                      annot=True, 
                      fmt='.1f',
                      cmap='RdYlBu_r',
                      vmin=80, 
                      vmax=100,
                      mask=mask,
                      square=True,
                      ax=ax1,
                      cbar_kws={'label': 'ANI (%)'})
    ax1.set_title('ANI Similarity Matrix', fontweight='bold')
    ax1.set_xlabel('Reference Genomes')
    ax1.set_ylabel('Query Genomes')
    
    # 子图2: ANI分布直方图
    ani_values = []
    for i in range(len(ani_matrix)):
        for j in range(i+1, len(ani_matrix)):
            if not pd.isna(ani_matrix.iloc[i, j]):
                ani_values.append(ani_matrix.iloc[i, j])
    
    if ani_values:
        ax2.hist(ani_values, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        ax2.axvline(95, color='red', linestyle='--', label='Species boundary (95%)')
        ax2.axvline(85, color='orange', linestyle='--', label='Genus boundary (85%)')
        ax2.set_xlabel('ANI Value (%)')
        ax2.set_ylabel('Frequency')
        ax2.set_title('ANI Value Distribution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    
    # 子图3: 层次聚类树状图
    if len(ani_matrix) > 2:
        distance_matrix = 1 - (ani_matrix / 100)
        distance_matrix = distance_matrix.fillna(1.0)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
        
        try:
            condensed_distances = squareform(distance_matrix.values)
            linkage_matrix = linkage(condensed_distances, method='average')
            
            dendrogram(linkage_matrix, 
                      labels=ani_matrix.index.tolist(),
                      orientation='top',
                      leaf_rotation=45,
                      leaf_font_size=10,
                      ax=ax3)
            ax3.set_title('Hierarchical Clustering (ANI-based)')
            ax3.set_ylabel('Distance (1 - ANI/100)')
        except Exception as e:
            ax3.text(0.5, 0.5, f'Clustering failed:\n{str(e)}', 
                    ha='center', va='center', transform=ax3.transAxes)
            ax3.set_title('Hierarchical Clustering (Failed)')
    
    # 子图4: ANI vs 分类学关系
    if ani_values:
        same_species = [v for v in ani_values if v >= 95]
        related_species = [v for v in ani_values if 85 <= v < 95]
        distant_species = [v for v in ani_values if v < 85]
        
        categories = ['Same Species\n(≥95%)', 'Related Species\n(85-95%)', 'Distant Species\n(<85%)']
        counts = [len(same_species), len(related_species), len(distant_species)]
        colors = ['#2E8B57', '#FF8C00', '#DC143C']
        
        bars = ax4.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
        ax4.set_ylabel('Number of Pairs')
        ax4.set_title('Taxonomic Relationship Distribution')
        ax4.grid(True, alpha=0.3, axis='y')
        
        # 添加数值标签
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{count}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    
    # 保存图片
    heatmap_file = os.path.join(output_dir, "comprehensive_ani_analysis.png")
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"综合ANI分析图保存至: {heatmap_file}")

def visualize_phylogenetic_tree(tree_file, output_dir):
    """
    可视化系统发育树
    
    Args:
        tree_file (str): 树文件路径
        output_dir (str): 输出目录
    """
    print("正在生成系统发育树可视化...")
    
    try:
        tree = Phylo.read(tree_file, "newick")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # 左图: 标准系统发育树
        Phylo.draw(tree, axes=ax1, do_show=False)
        ax1.set_title('Phylogenetic Tree (Standard View)', fontweight='bold')
        ax1.set_xlabel('Evolutionary Distance')
        
        # 右图: 圆形系统发育树
        try:
            Phylo.draw(tree, axes=ax2, do_show=False, branch_labels=lambda c: c.branch_length)
            ax2.set_title('Phylogenetic Tree (Circular View)', fontweight='bold')
        except:
            # 如果圆形图失败，绘制简化版本
            ax2.text(0.5, 0.5, 'Circular tree view\nnot available', 
                    ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('Phylogenetic Tree (Alternative View)')
        
        plt.tight_layout()
        
        # 保存图片
        tree_plot = os.path.join(output_dir, "phylogenetic_tree_visualization.png")
        plt.savefig(tree_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"系统发育树可视化保存至: {tree_plot}")
        
    except Exception as e:
        print(f"系统发育树可视化失败: {e}")

def create_divergence_timeline(divergence_file, output_dir):
    """
    创建发散时间时间轴
    
    Args:
        divergence_file (str): 发散时间文件路径
        output_dir (str): 输出目录
    """
    print("正在生成发散时间时间轴...")
    
    try:
        # 读取发散时间数据
        divergence_data = pd.read_csv(divergence_file, sep='\t', comment='#')
        
        # 过滤掉末端节点（时间为0的节点）
        internal_nodes = divergence_data[divergence_data['Divergence_Time_MYA'] > 0]
        
        if len(internal_nodes) == 0:
            print("未找到内部节点数据，跳过时间轴生成")
            return
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
        
        # 上图: 时间轴
        times = internal_nodes['Divergence_Time_MYA'].values
        nodes = internal_nodes['Node'].values
        
        # 创建时间轴
        y_pos = np.arange(len(nodes))
        
        # 绘制时间点
        scatter = ax1.scatter(times, y_pos, s=100, c=times, cmap='viridis', 
                            alpha=0.7, edgecolors='black')
        
        # 添加节点标签
        for i, (time, node) in enumerate(zip(times, nodes)):
            ax1.text(time + max(times) * 0.02, i, node, 
                    va='center', fontsize=9)
        
        ax1.set_xlabel('Time (Million Years Ago)')
        ax1.set_ylabel('Divergence Events')
        ax1.set_title('Divergence Time Timeline', fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.invert_xaxis()  # 时间从现在到过去
        
        # 添加颜色条
        cbar = plt.colorbar(scatter, ax=ax1)
        cbar.set_label('Divergence Time (MYA)')
        
        # 下图: 发散时间分布
        ax2.hist(times, bins=min(10, len(times)), alpha=0.7, 
                color='lightblue', edgecolor='black')
        ax2.set_xlabel('Divergence Time (Million Years Ago)')
        ax2.set_ylabel('Number of Events')
        ax2.set_title('Distribution of Divergence Times')
        ax2.grid(True, alpha=0.3)
        
        # 添加统计信息
        stats_text = f'Mean: {np.mean(times):.1f} MYA\n'
        stats_text += f'Median: {np.median(times):.1f} MYA\n'
        stats_text += f'Range: {np.min(times):.1f} - {np.max(times):.1f} MYA'
        
        ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # 保存图片
        timeline_plot = os.path.join(output_dir, "divergence_timeline.png")
        plt.savefig(timeline_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"发散时间时间轴保存至: {timeline_plot}")
        
    except Exception as e:
        print(f"发散时间时间轴生成失败: {e}")

def create_summary_dashboard(ani_file, tree_file, divergence_file, output_dir):
    """
    创建综合分析仪表板
    
    Args:
        ani_file (str): ANI矩阵文件
        tree_file (str): 系统发育树文件
        divergence_file (str): 发散时间文件
        output_dir (str): 输出目录
    """
    print("正在生成综合分析仪表板...")
    
    fig = plt.figure(figsize=(20, 12))
    
    # 创建网格布局
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
    
    # 1. ANI热图 (左上)
    ax1 = fig.add_subplot(gs[0, :2])
    if os.path.exists(ani_file):
        ani_matrix = load_ani_matrix(ani_file)
        if ani_matrix is not None:
            mask = ani_matrix.isna()
            sns.heatmap(ani_matrix, annot=True, fmt='.1f', cmap='RdYlBu_r',
                       vmin=80, vmax=100, mask=mask, square=True, ax=ax1,
                       cbar_kws={'label': 'ANI (%)'})
            ax1.set_title('ANI Similarity Matrix', fontweight='bold')
    
    # 2. 系统发育树 (右上)
    ax2 = fig.add_subplot(gs[0, 2:])
    if os.path.exists(tree_file):
        try:
            tree = Phylo.read(tree_file, "newick")
            Phylo.draw(tree, axes=ax2, do_show=False)
            ax2.set_title('Phylogenetic Tree', fontweight='bold')
        except:
            ax2.text(0.5, 0.5, 'Tree visualization\nnot available', 
                    ha='center', va='center', transform=ax2.transAxes)
    
    # 3. ANI分布 (左中)
    ax3 = fig.add_subplot(gs[1, :2])
    if os.path.exists(ani_file) and ani_matrix is not None:
        ani_values = []
        for i in range(len(ani_matrix)):
            for j in range(i+1, len(ani_matrix)):
                if not pd.isna(ani_matrix.iloc[i, j]):
                    ani_values.append(ani_matrix.iloc[i, j])
        
        if ani_values:
            ax3.hist(ani_values, bins=15, alpha=0.7, color='skyblue', edgecolor='black')
            ax3.axvline(95, color='red', linestyle='--', label='Species (95%)')
            ax3.axvline(85, color='orange', linestyle='--', label='Genus (85%)')
            ax3.set_xlabel('ANI Value (%)')
            ax3.set_ylabel('Frequency')
            ax3.set_title('ANI Distribution')
            ax3.legend()
            ax3.grid(True, alpha=0.3)
    
    # 4. 发散时间 (右中)
    ax4 = fig.add_subplot(gs[1, 2:])
    if os.path.exists(divergence_file):
        try:
            divergence_data = pd.read_csv(divergence_file, sep='\t', comment='#')
            internal_nodes = divergence_data[divergence_data['Divergence_Time_MYA'] > 0]
            
            if len(internal_nodes) > 0:
                times = internal_nodes['Divergence_Time_MYA'].values
                ax4.hist(times, bins=min(8, len(times)), alpha=0.7, 
                        color='lightgreen', edgecolor='black')
                ax4.set_xlabel('Divergence Time (MYA)')
                ax4.set_ylabel('Frequency')
                ax4.set_title('Divergence Time Distribution')
                ax4.grid(True, alpha=0.3)
        except:
            ax4.text(0.5, 0.5, 'Divergence time data\nnot available', 
                    ha='center', va='center', transform=ax4.transAxes)
    
    # 5. 分类学关系统计 (左下)
    ax5 = fig.add_subplot(gs[2, :2])
    if os.path.exists(ani_file) and ani_matrix is not None and ani_values:
        same_species = len([v for v in ani_values if v >= 95])
        related_species = len([v for v in ani_values if 85 <= v < 95])
        distant_species = len([v for v in ani_values if v < 85])
        
        categories = ['Same Species\n(≥95%)', 'Related Species\n(85-95%)', 'Distant Species\n(<85%)']
        counts = [same_species, related_species, distant_species]
        colors = ['#2E8B57', '#FF8C00', '#DC143C']
        
        bars = ax5.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
        ax5.set_ylabel('Number of Pairs')
        ax5.set_title('Taxonomic Relationships')
        ax5.grid(True, alpha=0.3, axis='y')
        
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            ax5.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{count}', ha='center', va='bottom', fontweight='bold')
    
    # 6. 分析总结 (右下)
    ax6 = fig.add_subplot(gs[2, 2:])
    ax6.axis('off')
    
    # 生成总结文本
    summary_text = "Analysis Summary\n" + "="*20 + "\n\n"
    
    if os.path.exists(ani_file) and ani_matrix is not None:
        summary_text += f"• Genomes analyzed: {len(ani_matrix)}\n"
        if ani_values:
            summary_text += f"• Mean ANI: {np.mean(ani_values):.1f}%\n"
            summary_text += f"• ANI range: {np.min(ani_values):.1f}-{np.max(ani_values):.1f}%\n"
    
    if os.path.exists(tree_file):
        summary_text += f"• Phylogenetic tree: Available\n"
    
    if os.path.exists(divergence_file):
        try:
            divergence_data = pd.read_csv(divergence_file, sep='\t', comment='#')
            internal_nodes = divergence_data[divergence_data['Divergence_Time_MYA'] > 0]
            if len(internal_nodes) > 0:
                times = internal_nodes['Divergence_Time_MYA'].values
                summary_text += f"• Oldest divergence: {np.max(times):.1f} MYA\n"
                summary_text += f"• Most recent: {np.min(times):.1f} MYA\n"
        except:
            pass
    
    summary_text += "\nKey Findings:\n"
    summary_text += "• Phylogenomic analysis completed\n"
    summary_text += "• ANI-based relationships determined\n"
    summary_text += "• Divergence times estimated\n"
    
    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, 
            verticalalignment='top', fontsize=11, fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8))
    
    # 添加主标题
    fig.suptitle('Microbial Phylogenomics Analysis Dashboard', 
                fontsize=16, fontweight='bold', y=0.98)
    
    # 保存仪表板
    dashboard_file = os.path.join(output_dir, "phylogenomics_dashboard.png")
    plt.savefig(dashboard_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"综合分析仪表板保存至: {dashboard_file}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='Visualization Generator for Phylogenomics Analysis')
    parser.add_argument('--ani', default='results/ani_matrix.txt',
                       help='ANI matrix file')
    parser.add_argument('--tree', default='results/concatenated_alignment.fasta.treefile',
                       help='Phylogenetic tree file')
    parser.add_argument('--divergence', default='results/divergence_times.txt',
                       help='Divergence times file')
    parser.add_argument('--output', '-o', default='figures',
                       help='Output directory')
    
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    
    print("开始生成可视化图表...")
    print(f"输出目录: {args.output}")
    
    # 检查输入文件
    files_to_check = [
        ('ANI matrix', args.ani),
        ('Phylogenetic tree', args.tree),
        ('Divergence times', args.divergence)
    ]
    
    available_files = []
    for file_type, file_path in files_to_check:
        if os.path.exists(file_path):
            print(f"✓ {file_type}: {file_path}")
            available_files.append(file_path)
        else:
            print(f"✗ {file_type}: {file_path} (not found)")
    
    if not available_files:
        print("未找到任何输入文件，程序退出")
        sys.exit(1)
    
    # 生成各种可视化
    if os.path.exists(args.ani):
        ani_matrix = load_ani_matrix(args.ani)
        if ani_matrix is not None:
            create_comprehensive_ani_heatmap(ani_matrix, args.output)
    
    if os.path.exists(args.tree):
        visualize_phylogenetic_tree(args.tree, args.output)
    
    if os.path.exists(args.divergence):
        create_divergence_timeline(args.divergence, args.output)
    
    # 生成综合仪表板
    create_summary_dashboard(args.ani, args.tree, args.divergence, args.output)
    
    print("\n可视化生成完成！")
    print("生成的图表文件:")
    
    output_files = [
        "comprehensive_ani_analysis.png",
        "phylogenetic_tree_visualization.png", 
        "divergence_timeline.png",
        "phylogenomics_dashboard.png"
    ]
    
    for output_file in output_files:
        file_path = os.path.join(args.output, output_file)
        if os.path.exists(file_path):
            print(f"  ✓ {file_path}")
        else:
            print(f"  ✗ {file_path} (not generated)")

if __name__ == "__main__":
    main()