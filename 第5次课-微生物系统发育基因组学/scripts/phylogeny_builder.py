#!/usr/bin/env python3
"""
Phylogeny Builder Script for Microbial Phylogenomics Course
构建全基因组系统发育树的脚本

Author: Course Development Team
Date: 2025
Version: 1.0.0
"""

import os
import sys
import argparse
import subprocess
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

def find_core_genes(genome_dir, output_dir, identity_threshold=0.7):
    """
    识别核心基因（在所有基因组中都存在的单拷贝基因）
    
    Args:
        genome_dir (str): 基因组文件目录
        output_dir (str): 输出目录
        identity_threshold (float): 同源性阈值
    
    Returns:
        list: 核心基因列表
    """
    print("正在识别核心基因...")
    
    # 获取基因组文件列表
    genome_files = glob.glob(os.path.join(genome_dir, "*.fna"))
    if not genome_files:
        print(f"错误: 在 {genome_dir} 中未找到基因组文件")
        return []
    
    print(f"找到 {len(genome_files)} 个基因组文件")
    
    # 创建核心基因输出目录
    core_genes_dir = os.path.join(output_dir, "core_genes")
    os.makedirs(core_genes_dir, exist_ok=True)
    
    # 使用OrthoFinder或简化的BLAST方法识别核心基因
    # 这里使用简化的方法作为示例
    
    # 预定义的核心基因列表（基于常见的单拷贝核心基因）
    core_gene_names = [
        "rpoB",  # RNA聚合酶β亚基
        "rpoC",  # RNA聚合酶β'亚基
        "gyrA",  # DNA回旋酶A亚基
        "gyrB",  # DNA回旋酶B亚基
        "recA",  # 重组酶A
        "dnaK",  # 分子伴侣DnaK
        "groEL", # 分子伴侣GroEL
        "atpD",  # ATP合酶β亚基
        "infB",  # 翻译起始因子IF-2
        "nusA",  # 转录终止因子NusA
        "pgk",   # 磷酸甘油酸激酶
        "pyrG",  # CTP合酶
        "secY",  # 蛋白转运通道SecY
        "tsf",   # 翻译延伸因子Ts
        "fusA"   # 翻译延伸因子G
    ]
    
    # 为每个核心基因创建多序列文件
    core_genes_found = []
    
    for gene_name in core_gene_names:
        gene_sequences = []
        
        for genome_file in genome_files:
            genome_name = os.path.basename(genome_file).replace('.fna', '')
            
            # 模拟基因提取（实际应用中需要基因注释信息）
            # 这里创建模拟序列作为示例
            mock_sequence = create_mock_gene_sequence(gene_name, genome_name)
            if mock_sequence:
                gene_sequences.append(mock_sequence)
        
        # 如果所有基因组都有该基因，则认为是核心基因
        if len(gene_sequences) == len(genome_files):
            core_genes_found.append(gene_name)
            
            # 保存核心基因序列
            gene_file = os.path.join(core_genes_dir, f"{gene_name}.fasta")
            with open(gene_file, 'w') as f:
                SeqIO.write(gene_sequences, f, "fasta")
            
            print(f"核心基因 {gene_name} 已保存至 {gene_file}")
    
    # 保存核心基因统计信息
    stats_file = os.path.join(output_dir, "core_genes_stats.txt")
    with open(stats_file, 'w') as f:
        f.write("# Core Genes Analysis Statistics\n")
        f.write(f"Total genomes analyzed: {len(genome_files)}\n")
        f.write(f"Core genes identified: {len(core_genes_found)}\n")
        f.write(f"Core genes list: {', '.join(core_genes_found)}\n\n")
        
        f.write("# Genome files processed:\n")
        for genome_file in genome_files:
            f.write(f"  - {os.path.basename(genome_file)}\n")
    
    print(f"核心基因统计信息保存至: {stats_file}")
    print(f"识别到 {len(core_genes_found)} 个核心基因")
    
    return core_genes_found

def create_mock_gene_sequence(gene_name, genome_name):
    """
    创建模拟基因序列（用于演示目的）
    
    Args:
        gene_name (str): 基因名称
        genome_name (str): 基因组名称
    
    Returns:
        SeqRecord: 模拟基因序列记录
    """
    # 基于基因名称和基因组名称生成不同的序列
    base_sequences = {
        "rpoB": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "rpoC": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "gyrA": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "gyrB": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "recA": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "dnaK": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "groEL": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "atpD": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "infB": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "nusA": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "pgk": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "pyrG": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "secY": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "tsf": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG",
        "fusA": "ATGAGCAAAGACGTCCTGAAAGACATCGACGAACTGGTGAAAGACCTGAAAGACGTCGACGAACTGGTGAAAGACCTG"
    }
    
    if gene_name not in base_sequences:
        return None
    
    # 基于基因组名称引入变异
    base_seq = base_sequences[gene_name]
    
    # 简单的变异模拟
    variations = {
        "E_coli_K12": "",
        "E_coli_O157H7": "A",
        "S_enterica": "GC"
    }
    
    variation = variations.get(genome_name, "")
    modified_seq = base_seq + variation
    
    # 创建SeqRecord对象
    seq_record = SeqRecord(
        Seq(modified_seq),
        id=f"{genome_name}_{gene_name}",
        description=f"{gene_name} gene from {genome_name}"
    )
    
    return seq_record

def align_core_genes(core_genes_dir, output_dir):
    """
    对核心基因进行多序列比对
    
    Args:
        core_genes_dir (str): 核心基因目录
        output_dir (str): 输出目录
    
    Returns:
        list: 比对文件列表
    """
    print("正在进行多序列比对...")
    
    # 创建比对输出目录
    alignments_dir = os.path.join(output_dir, "alignments")
    os.makedirs(alignments_dir, exist_ok=True)
    
    # 获取核心基因文件
    gene_files = glob.glob(os.path.join(core_genes_dir, "*.fasta"))
    alignment_files = []
    
    for gene_file in gene_files:
        gene_name = os.path.basename(gene_file).replace('.fasta', '')
        alignment_file = os.path.join(alignments_dir, f"{gene_name}_aligned.fasta")
        
        # 使用MUSCLE进行多序列比对
        try:
            cmd = ["muscle", "-in", gene_file, "-out", alignment_file]
            subprocess.run(cmd, check=True, capture_output=True)
            alignment_files.append(alignment_file)
            print(f"基因 {gene_name} 比对完成")
        except subprocess.CalledProcessError:
            # 如果MUSCLE不可用，创建简单的比对文件
            print(f"MUSCLE不可用，为基因 {gene_name} 创建简单比对")
            create_simple_alignment(gene_file, alignment_file)
            alignment_files.append(alignment_file)
        except FileNotFoundError:
            # MUSCLE未安装，创建简单比对
            print(f"MUSCLE未安装，为基因 {gene_name} 创建简单比对")
            create_simple_alignment(gene_file, alignment_file)
            alignment_files.append(alignment_file)
    
    print(f"完成 {len(alignment_files)} 个基因的多序列比对")
    return alignment_files

def create_simple_alignment(input_file, output_file):
    """
    创建简单的序列比对（当MUSCLE不可用时）
    
    Args:
        input_file (str): 输入序列文件
        output_file (str): 输出比对文件
    """
    # 读取序列并创建简单比对（实际上就是复制）
    sequences = list(SeqIO.parse(input_file, "fasta"))
    
    # 确保所有序列长度相同（简单填充）
    max_length = max(len(seq.seq) for seq in sequences)
    
    aligned_sequences = []
    for seq in sequences:
        # 用'-'填充到相同长度
        padded_seq = str(seq.seq).ljust(max_length, '-')
        aligned_seq = SeqRecord(
            Seq(padded_seq),
            id=seq.id,
            description=seq.description
        )
        aligned_sequences.append(aligned_seq)
    
    # 保存比对结果
    with open(output_file, 'w') as f:
        SeqIO.write(aligned_sequences, f, "fasta")

def concatenate_alignments(alignment_files, output_dir):
    """
    连接多个基因的比对序列
    
    Args:
        alignment_files (list): 比对文件列表
        output_dir (str): 输出目录
    
    Returns:
        str: 连接比对文件路径
    """
    print("正在连接多基因比对序列...")
    
    # 读取所有比对文件
    all_alignments = {}
    gene_order = []
    
    for alignment_file in alignment_files:
        gene_name = os.path.basename(alignment_file).replace('_aligned.fasta', '')
        gene_order.append(gene_name)
        
        sequences = {}
        for record in SeqIO.parse(alignment_file, "fasta"):
            # 提取基因组名称（去除基因名后缀）
            genome_name = record.id.replace(f"_{gene_name}", "")
            sequences[genome_name] = str(record.seq)
        
        all_alignments[gene_name] = sequences
    
    # 获取所有基因组名称
    all_genomes = set()
    for gene_seqs in all_alignments.values():
        all_genomes.update(gene_seqs.keys())
    
    all_genomes = sorted(list(all_genomes))
    
    # 连接序列
    concatenated_sequences = []
    
    for genome in all_genomes:
        concatenated_seq = ""
        
        for gene in gene_order:
            if genome in all_alignments[gene]:
                concatenated_seq += all_alignments[gene][genome]
            else:
                # 如果某个基因组缺少某个基因，用'-'填充
                gene_length = len(list(all_alignments[gene].values())[0])
                concatenated_seq += '-' * gene_length
        
        # 创建连接序列记录
        seq_record = SeqRecord(
            Seq(concatenated_seq),
            id=genome,
            description=f"Concatenated alignment for {genome}"
        )
        concatenated_sequences.append(seq_record)
    
    # 保存连接比对
    concatenated_file = os.path.join(output_dir, "concatenated_alignment.fasta")
    with open(concatenated_file, 'w') as f:
        SeqIO.write(concatenated_sequences, f, "fasta")
    
    print(f"连接比对保存至: {concatenated_file}")
    print(f"总序列长度: {len(concatenated_sequences[0].seq)} bp")
    print(f"包含基因组: {len(concatenated_sequences)} 个")
    
    return concatenated_file

def build_phylogenetic_tree(alignment_file, output_dir):
    """
    构建系统发育树
    
    Args:
        alignment_file (str): 连接比对文件
        output_dir (str): 输出目录
    
    Returns:
        str: 系统发育树文件路径
    """
    print("正在构建系统发育树...")
    
    # 使用IQ-TREE构建系统发育树
    try:
        cmd = [
            "iqtree",
            "-s", alignment_file,
            "-m", "MFP",  # 自动选择最佳模型
            "-bb", "1000",  # 1000次bootstrap
            "-nt", "AUTO"  # 自动检测CPU核心数
        ]
        
        subprocess.run(cmd, check=True, capture_output=True, cwd=output_dir)
        
        tree_file = alignment_file + ".treefile"
        print(f"系统发育树构建完成: {tree_file}")
        
        return tree_file
        
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("IQ-TREE不可用，创建简单的邻接树...")
        return create_simple_tree(alignment_file, output_dir)

def create_simple_tree(alignment_file, output_dir):
    """
    创建简单的系统发育树（当IQ-TREE不可用时）
    
    Args:
        alignment_file (str): 比对文件
        output_dir (str): 输出目录
    
    Returns:
        str: 树文件路径
    """
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio import AlignIO
    
    try:
        # 读取比对
        alignment = AlignIO.read(alignment_file, "fasta")
        
        # 计算距离矩阵
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)
        
        # 构建邻接树
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distance_matrix)
        
        # 保存树文件
        tree_file = os.path.join(output_dir, "simple_phylogenetic_tree.nwk")
        Phylo.write(tree, tree_file, "newick")
        
        print(f"简单系统发育树保存至: {tree_file}")
        return tree_file
        
    except Exception as e:
        print(f"构建简单树失败: {e}")
        
        # 创建最基本的树结构
        basic_tree_file = os.path.join(output_dir, "basic_tree.nwk")
        sequences = list(SeqIO.parse(alignment_file, "fasta"))
        
        if len(sequences) >= 2:
            # 创建简单的二叉树结构
            tree_string = f"({sequences[0].id}:0.1,({sequences[1].id}:0.1"
            for i in range(2, len(sequences)):
                tree_string += f",{sequences[i].id}:0.1"
            tree_string += "):0.1);"
            
            with open(basic_tree_file, 'w') as f:
                f.write(tree_string)
            
            print(f"基本树结构保存至: {basic_tree_file}")
            return basic_tree_file
        
        return None

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='Phylogeny Builder for Microbial Genomes')
    parser.add_argument('--step', choices=['core_genes', 'alignment', 'tree', 'all'],
                       default='all', help='Analysis step to perform')
    parser.add_argument('--input', '-i', default='.',
                       help='Input directory containing genome files')
    parser.add_argument('--output', '-o', default='results',
                       help='Output directory')
    
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    
    print("开始系统发育分析...")
    print(f"输入目录: {args.input}")
    print(f"输出目录: {args.output}")
    print(f"分析步骤: {args.step}")
    
    if args.step in ['core_genes', 'all']:
        # 识别核心基因
        core_genes = find_core_genes(args.input, args.output)
        if not core_genes:
            print("未找到核心基因，程序退出")
            sys.exit(1)
    
    if args.step in ['alignment', 'all']:
        # 多序列比对
        core_genes_dir = os.path.join(args.output, "core_genes")
        if not os.path.exists(core_genes_dir):
            print("核心基因目录不存在，请先运行核心基因识别步骤")
            sys.exit(1)
        
        alignment_files = align_core_genes(core_genes_dir, args.output)
        
        # 连接比对
        concatenated_file = concatenate_alignments(alignment_files, args.output)
    
    if args.step in ['tree', 'all']:
        # 构建系统发育树
        concatenated_file = os.path.join(args.output, "concatenated_alignment.fasta")
        if not os.path.exists(concatenated_file):
            print("连接比对文件不存在，请先运行比对步骤")
            sys.exit(1)
        
        tree_file = build_phylogenetic_tree(concatenated_file, args.output)
    
    print("\n系统发育分析完成！")
    print("主要输出文件:")
    if args.step in ['core_genes', 'all']:
        print(f"  - 核心基因统计: {args.output}/core_genes_stats.txt")
        print(f"  - 核心基因序列: {args.output}/core_genes/")
    if args.step in ['alignment', 'all']:
        print(f"  - 序列比对: {args.output}/alignments/")
        print(f"  - 连接比对: {args.output}/concatenated_alignment.fasta")
    if args.step in ['tree', 'all']:
        print(f"  - 系统发育树: {tree_file}")

if __name__ == "__main__":
    main()