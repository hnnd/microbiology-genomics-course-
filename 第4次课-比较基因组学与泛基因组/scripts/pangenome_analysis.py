#!/usr/bin/env python3
"""
泛基因组统计分析脚本
用于分析Roary输出结果，生成泛基因组统计报告

作者: 微生物基因组学课程组
版本: 1.0.0
日期: 2025
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='泛基因组统计分析')
    parser.add_argument('--mode', choices=['qc', 'stats'], required=True,
                       help='分析模式: qc=质量控制, stats=统计分析')
    parser.add_argument('--input', required=True,
                       help='输入目录或文件路径')
    parser.add_argument('--output', default='results',
                       help='输出目录 (默认: results)')
    return parser.parse_args()

def quality_control(input_dir, output_dir):
    """质量控制检查"""
    print("执行质量控制检查...")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 检查GFF文件
    gff_files = [f for f in os.listdir(input_dir) if f.endswith('.gff')]
    
    qc_report = []
    qc_report.append("=== 泛基因组分析质量控制报告 ===\n")
    qc_report.append(f"检查时间: {pd.Timestamp.now()}\n")
    qc_report.append(f"输入目录: {input_dir}\n")
    qc_report.append(f"发现GFF文件数量: {len(gff_files)}\n\n")
    
    gene_counts = {}
    
    for gff_file in gff_files:
        filepath = os.path.join(input_dir, gff_file)
        
        # 统计基因数量
        cds_count = 0
        total_lines = 0
        
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    total_lines += 1
                    if not line.startswith('#') and '\tCDS\t' in line:
                        cds_count += 1
            
            gene_counts[gff_file] = cds_count
            
            qc_report.append(f"文件: {gff_file}")
            qc_report.append(f"  总行数: {total_lines}")
            qc_report.append(f"  CDS数量: {cds_count}")
            
            # 质量评估
            if cds_count < 1000:
                qc_report.append(f"  ⚠️  警告: CDS数量过少 (<1000)")
            elif cds_count > 10000:
                qc_report.append(f"  ⚠️  警告: CDS数量过多 (>10000)")
            else:
                qc_report.append(f"  ✅ CDS数量正常")
            
            qc_report.append("")
            
        except Exception as e:
            qc_report.append(f"  ❌ 错误: 无法读取文件 - {str(e)}\n")
    
    # 统计摘要
    if gene_counts:
        qc_report.append("=== 统计摘要 ===")
        qc_report.append(f"平均基因数: {np.mean(list(gene_counts.values())):.0f}")
        qc_report.append(f"基因数范围: {min(gene_counts.values())} - {max(gene_counts.values())}")
        qc_report.append(f"标准差: {np.std(list(gene_counts.values())):.0f}")
    
    # 保存报告
    report_file = os.path.join(output_dir, 'quality_report.txt')
    with open(report_file, 'w') as f:
        f.write('\n'.join(qc_report))
    
    print(f"质量控制报告已保存到: {report_file}")
    return gene_counts

def analyze_pangenome_stats(input_dir, output_dir):
    """分析泛基因组统计数据"""
    print("执行泛基因组统计分析...")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 查找Roary输出文件
    gene_presence_file = os.path.join(input_dir, 'gene_presence_absence.csv')
    summary_stats_file = os.path.join(input_dir, 'summary_statistics.txt')
    
    stats_report = []
    stats_report.append("=== 泛基因组统计分析报告 ===\n")
    stats_report.append(f"分析时间: {pd.Timestamp.now()}\n")
    
    if os.path.exists(gene_presence_file):
        # 读取基因存在/缺失矩阵
        print("分析基因存在/缺失矩阵...")
        df = pd.read_csv(gene_presence_file, low_memory=False)
        
        # 获取菌株列（跳过前14列的注释信息）
        strain_columns = df.columns[14:]
        num_strains = len(strain_columns)
        
        stats_report.append(f"菌株数量: {num_strains}")
        stats_report.append(f"基因家族总数: {len(df)}")
        
        # 计算核心、辅助、独特基因
        gene_presence_counts = []
        
        for idx, row in df.iterrows():
            # 计算每个基因家族在多少个菌株中存在
            presence_count = sum(1 for col in strain_columns if pd.notna(row[col]) and row[col] != '')
            gene_presence_counts.append(presence_count)
        
        # 分类基因
        core_genes = sum(1 for count in gene_presence_counts if count == num_strains)
        accessory_genes = sum(1 for count in gene_presence_counts if 1 < count < num_strains)
        unique_genes = sum(1 for count in gene_presence_counts if count == 1)
        
        stats_report.append(f"\n=== 基因组成分析 ===")
        stats_report.append(f"核心基因数量: {core_genes} ({core_genes/len(df)*100:.1f}%)")
        stats_report.append(f"辅助基因数量: {accessory_genes} ({accessory_genes/len(df)*100:.1f}%)")
        stats_report.append(f"独特基因数量: {unique_genes} ({unique_genes/len(df)*100:.1f}%)")
        
        # 每个菌株的基因数量
        strain_gene_counts = {}
        for strain in strain_columns:
            gene_count = sum(1 for val in df[strain] if pd.notna(val) and val != '')
            strain_gene_counts[strain] = gene_count
        
        stats_report.append(f"\n=== 各菌株基因数量 ===")
        for strain, count in strain_gene_counts.items():
            stats_report.append(f"{strain}: {count} 个基因")
        
        # 保存详细统计数据
        detailed_stats = {
            'total_gene_families': len(df),
            'num_strains': num_strains,
            'core_genes': core_genes,
            'accessory_genes': accessory_genes,
            'unique_genes': unique_genes,
            'strain_gene_counts': strain_gene_counts
        }
        
        # 保存为JSON格式
        import json
        with open(os.path.join(output_dir, 'pangenome_stats.json'), 'w') as f:
            json.dump(detailed_stats, f, indent=2)
        
    else:
        stats_report.append("❌ 未找到 gene_presence_absence.csv 文件")
    
    # 读取Roary摘要统计
    if os.path.exists(summary_stats_file):
        stats_report.append(f"\n=== Roary摘要统计 ===")
        with open(summary_stats_file, 'r') as f:
            roary_stats = f.read()
            stats_report.append(roary_stats)
    else:
        stats_report.append("❌ 未找到 summary_statistics.txt 文件")
    
    # 保存统计报告
    report_file = os.path.join(output_dir, 'pangenome_statistics.txt')
    with open(report_file, 'w') as f:
        f.write('\n'.join(stats_report))
    
    print(f"统计分析报告已保存到: {report_file}")

def main():
    """主函数"""
    args = parse_arguments()
    
    print(f"泛基因组分析脚本 v1.0.0")
    print(f"模式: {args.mode}")
    print(f"输入: {args.input}")
    print(f"输出: {args.output}")
    print("-" * 50)
    
    try:
        if args.mode == 'qc':
            quality_control(args.input, args.output)
        elif args.mode == 'stats':
            analyze_pangenome_stats(args.input, args.output)
        
        print("分析完成！")
        
    except Exception as e:
        print(f"❌ 错误: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()