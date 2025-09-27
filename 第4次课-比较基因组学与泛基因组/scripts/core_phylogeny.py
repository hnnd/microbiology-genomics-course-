#!/usr/bin/env python3
"""
核心基因组系统发育分析脚本
用于处理核心基因组比对文件，构建系统发育树

作者: 微生物基因组学课程组
版本: 1.0.0
日期: 2025
"""

import os
import sys
import argparse
import subprocess
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='核心基因组系统发育分析')
    parser.add_argument('--input', required=True,
                       help='核心基因组比对文件 (FASTA格式)')
    parser.add_argument('--output', default='results',
                       help='输出目录 (默认: results)')
    parser.add_argument('--method', choices=['nj', 'upgma'], default='nj',
                       help='建树方法: nj=邻接法, upgma=UPGMA (默认: nj)')
    parser.add_argument('--bootstrap', type=int, default=100,
                       help='Bootstrap重复次数 (默认: 100)')
    return parser.parse_args()

def check_alignment_file(alignment_file):
    """检查比对文件质量"""
    print(f"检查比对文件: {alignment_file}")
    
    if not os.path.exists(alignment_file):
        raise FileNotFoundError(f"比对文件不存在: {alignment_file}")
    
    # 读取序列
    sequences = list(SeqIO.parse(alignment_file, "fasta"))
    
    if len(sequences) < 2:
        raise ValueError("序列数量少于2条，无法构建系统发育树")
    
    # 检查序列长度一致性
    seq_lengths = [len(seq) for seq in sequences]
    if len(set(seq_lengths)) > 1:
        raise ValueError("序列长度不一致，请检查比对质量")
    
    print(f"✅ 序列数量: {len(sequences)}")
    print(f"✅ 序列长度: {seq_lengths[0]} bp")
    print(f"✅ 序列ID: {[seq.id for seq in sequences]}")
    
    return sequences

def calculate_distances(sequences, output_dir):
    """计算序列间距离"""
    print("计算序列间进化距离...")
    
    # 创建多序列比对对象
    alignment = MultipleSeqAlignment(sequences)
    
    # 计算距离矩阵
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    
    # 保存距离矩阵
    distance_file = os.path.join(output_dir, 'distance_matrix.txt')
    with open(distance_file, 'w') as f:
        f.write("序列间进化距离矩阵\n")
        f.write("=" * 40 + "\n\n")
        
        # 写入矩阵头
        seq_names = [seq.id for seq in sequences]
        f.write("\t" + "\t".join(seq_names) + "\n")
        
        # 写入距离数据
        for i, name1 in enumerate(seq_names):
            row = [name1]
            for j, name2 in enumerate(seq_names):
                if i <= j:
                    distance = distance_matrix[i, j]
                    row.append(f"{distance:.6f}")
                else:
                    row.append("-")
            f.write("\t".join(row) + "\n")
    
    print(f"距离矩阵已保存到: {distance_file}")
    return distance_matrix

def build_phylogenetic_tree(alignment, distance_matrix, method, output_dir):
    """构建系统发育树"""
    print(f"使用{method.upper()}方法构建系统发育树...")
    
    # 构建系统发育树
    constructor = DistanceTreeConstructor()
    
    if method == 'nj':
        tree = constructor.nj(distance_matrix)
        tree_method = "Neighbor-Joining"
    elif method == 'upgma':
        tree = constructor.upgma(distance_matrix)
        tree_method = "UPGMA"
    
    # 保存树文件
    tree_file = os.path.join(output_dir, f'phylogenetic_tree_{method}.newick')
    Phylo.write(tree, tree_file, 'newick')
    print(f"系统发育树已保存到: {tree_file}")
    
    # 生成树的可视化
    plt.figure(figsize=(12, 8))
    Phylo.draw(tree, do_show=False)
    plt.title(f'核心基因组系统发育树 ({tree_method})', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    tree_image = os.path.join(output_dir, f'phylogenetic_tree_{method}.png')
    plt.savefig(tree_image, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"系统发育树图片已保存到: {tree_image}")
    
    return tree

def run_fasttree_analysis(alignment_file, output_dir):
    """使用FastTree进行系统发育分析"""
    print("使用FastTree构建最大似然树...")
    
    try:
        # 检查FastTree是否可用
        subprocess.run(['FastTree', '-help'], capture_output=True, check=True)
        
        # 运行FastTree
        fasttree_output = os.path.join(output_dir, 'fasttree_phylogeny.newick')
        
        with open(fasttree_output, 'w') as f:
            result = subprocess.run([
                'FastTree', '-nt', '-gtr', alignment_file
            ], stdout=f, stderr=subprocess.PIPE, text=True)
        
        if result.returncode == 0:
            print(f"FastTree系统发育树已保存到: {fasttree_output}")
            
            # 可视化FastTree结果
            tree = Phylo.read(fasttree_output, 'newick')
            plt.figure(figsize=(12, 8))
            Phylo.draw(tree, do_show=False)
            plt.title('核心基因组系统发育树 (FastTree ML)', fontsize=14, fontweight='bold')
            plt.tight_layout()
            
            fasttree_image = os.path.join(output_dir, 'fasttree_phylogeny.png')
            plt.savefig(fasttree_image, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"FastTree系统发育树图片已保存到: {fasttree_image}")
            
        else:
            print(f"FastTree运行失败: {result.stderr}")
            
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("⚠️  FastTree未安装或不可用，跳过最大似然分析")

def generate_phylogeny_report(sequences, distance_matrix, tree, method, output_dir):
    """生成系统发育分析报告"""
    print("生成系统发育分析报告...")
    
    report = []
    report.append("=== 核心基因组系统发育分析报告 ===\n")
    report.append(f"分析时间: {pd.Timestamp.now()}")
    report.append(f"建树方法: {method.upper()}")
    report.append(f"序列数量: {len(sequences)}")
    report.append(f"比对长度: {len(sequences[0])} bp\n")
    
    # 序列信息
    report.append("=== 序列信息 ===")
    for i, seq in enumerate(sequences):
        gc_content = (seq.seq.count('G') + seq.seq.count('C')) / len(seq.seq) * 100
        report.append(f"{i+1}. {seq.id}")
        report.append(f"   长度: {len(seq.seq)} bp")
        report.append(f"   GC含量: {gc_content:.1f}%")
    
    # 距离矩阵摘要
    report.append(f"\n=== 进化距离摘要 ===")
    distances = []
    seq_names = [seq.id for seq in sequences]
    
    for i in range(len(seq_names)):
        for j in range(i+1, len(seq_names)):
            distance = distance_matrix[i, j]
            distances.append(distance)
            report.append(f"{seq_names[i]} <-> {seq_names[j]}: {distance:.6f}")
    
    if distances:
        report.append(f"\n平均距离: {sum(distances)/len(distances):.6f}")
        report.append(f"最小距离: {min(distances):.6f}")
        report.append(f"最大距离: {max(distances):.6f}")
    
    # 系统发育树信息
    report.append(f"\n=== 系统发育树信息 ===")
    report.append(f"建树方法: {method.upper()}")
    report.append(f"分支数: {tree.count_terminals()}")
    report.append(f"内部节点数: {len(tree.get_nonterminals())}")
    
    # 保存报告
    report_file = os.path.join(output_dir, 'phylogeny_report.txt')
    with open(report_file, 'w') as f:
        f.write('\n'.join(report))
    
    print(f"系统发育分析报告已保存到: {report_file}")

def main():
    """主函数"""
    args = parse_arguments()
    
    print(f"核心基因组系统发育分析脚本 v1.0.0")
    print(f"输入文件: {args.input}")
    print(f"输出目录: {args.output}")
    print(f"建树方法: {args.method}")
    print("-" * 50)
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    
    try:
        # 检查比对文件
        sequences = check_alignment_file(args.input)
        
        # 计算距离矩阵
        distance_matrix = calculate_distances(sequences, args.output)
        
        # 构建系统发育树
        alignment = MultipleSeqAlignment(sequences)
        tree = build_phylogenetic_tree(alignment, distance_matrix, args.method, args.output)
        
        # 尝试使用FastTree
        run_fasttree_analysis(args.input, args.output)
        
        # 生成分析报告
        import pandas as pd
        generate_phylogeny_report(sequences, distance_matrix, tree, args.method, args.output)
        
        print("\n✅ 系统发育分析完成！")
        print(f"结果文件保存在: {args.output}")
        
    except Exception as e:
        print(f"❌ 错误: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()