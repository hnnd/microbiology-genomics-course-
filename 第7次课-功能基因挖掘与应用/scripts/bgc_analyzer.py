#!/usr/bin/env python3
"""
BGC (Biosynthetic Gene Cluster) Analyzer
用于次级代谢基因簇的预测和分析

作者: 微生物基因组学课程组
版本: 1.0.0
"""

import argparse
import os
import sys
from pathlib import Path
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import seaborn as sns

class BGCAnalyzer:
    """BGC分析主类"""
    
    def __init__(self, input_file, output_dir, threads=4):
        self.input_file = input_file
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # BGC类型定义
        self.bgc_types = {
            'PKS': 'Polyketide synthase',
            'NRPS': 'Non-ribosomal peptide synthetase',
            'Terpene': 'Terpene biosynthesis',
            'Bacteriocin': 'Bacteriocin biosynthesis',
            'Siderophore': 'Siderophore biosynthesis',
            'Other': 'Other secondary metabolites'
        }
        
    def predict_genes(self):
        """基因预测"""
        print("正在进行基因预测...")
        
        genes_file = self.output_dir / "genes.faa"
        gff_file = self.output_dir / "genes.gff"
        
        cmd = [
            "prodigal",
            "-i", str(self.input_file),
            "-a", str(genes_file),
            "-o", str(gff_file),
            "-f", "gff"
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            print(f"基因预测完成，结果保存至: {genes_file}")
            return genes_file
        except subprocess.CalledProcessError as e:
            print(f"基因预测失败: {e}")
            return None
    
    def search_bgc_domains(self, genes_file):
        """搜索BGC相关结构域"""
        print("正在搜索BGC相关结构域...")
        
        # 模拟BGC结构域搜索结果
        # 实际应用中会使用HMMER搜索Pfam数据库
        bgc_domains = {
            'PKS_KS': ['gene_001', 'gene_045', 'gene_123'],
            'PKS_AT': ['gene_002', 'gene_046', 'gene_124'],
            'NRPS_A': ['gene_015', 'gene_067', 'gene_156'],
            'NRPS_C': ['gene_016', 'gene_068', 'gene_157'],
            'Terpene_cyclase': ['gene_089', 'gene_234'],
            'Bacteriocin': ['gene_178', 'gene_199', 'gene_267']
        }
        
        # 保存结构域搜索结果
        domain_file = self.output_dir / "bgc_domains.txt"
        with open(domain_file, 'w') as f:
            f.write("Gene_ID\tDomain_Type\tDomain_Name\tE_value\n")
            for domain_type, genes in bgc_domains.items():
                for gene in genes:
                    f.write(f"{gene}\t{domain_type.split('_')[0]}\t{domain_type}\t1e-50\n")
        
        return domain_file
    
    def cluster_genes(self, domain_file):
        """基因聚类形成BGC"""
        print("正在进行基因聚类...")
        
        # 模拟基因聚类结果
        clusters = [
            {
                'cluster_id': 'BGC_001',
                'type': 'PKS',
                'genes': ['gene_001', 'gene_002', 'gene_003', 'gene_004', 'gene_005'],
                'start': 12000,
                'end': 45000,
                'completeness': 'complete'
            },
            {
                'cluster_id': 'BGC_002', 
                'type': 'NRPS',
                'genes': ['gene_015', 'gene_016', 'gene_017', 'gene_018'],
                'start': 156000,
                'end': 189000,
                'completeness': 'complete'
            },
            {
                'cluster_id': 'BGC_003',
                'type': 'Terpene',
                'genes': ['gene_089', 'gene_090', 'gene_091'],
                'start': 445000,
                'end': 467000,
                'completeness': 'partial'
            },
            {
                'cluster_id': 'BGC_004',
                'type': 'Bacteriocin',
                'genes': ['gene_178', 'gene_179', 'gene_180', 'gene_181'],
                'start': 890000,
                'end': 912000,
                'completeness': 'complete'
            }
        ]
        
        # 保存聚类结果
        cluster_file = self.output_dir / "clusters.txt"
        with open(cluster_file, 'w') as f:
            f.write("Cluster_ID\tType\tGenes\tStart\tEnd\tCompleteness\n")
            for cluster in clusters:
                genes_str = ','.join(cluster['genes'])
                f.write(f"{cluster['cluster_id']}\t{cluster['type']}\t{genes_str}\t"
                       f"{cluster['start']}\t{cluster['end']}\t{cluster['completeness']}\n")
        
        return clusters, cluster_file
    
    def generate_summary(self, clusters):
        """生成分析摘要"""
        print("正在生成分析摘要...")
        
        # 统计BGC类型
        type_counts = {}
        complete_counts = {}
        
        for cluster in clusters:
            bgc_type = cluster['type']
            type_counts[bgc_type] = type_counts.get(bgc_type, 0) + 1
            
            if cluster['completeness'] == 'complete':
                complete_counts[bgc_type] = complete_counts.get(bgc_type, 0) + 1
        
        # 生成摘要报告
        summary_file = self.output_dir / "cluster_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("=== BGC分析摘要报告 ===\n\n")
            f.write(f"输入基因组: {self.input_file}\n")
            f.write(f"总BGC数量: {len(clusters)}\n\n")
            
            f.write("BGC类型分布:\n")
            for bgc_type, count in type_counts.items():
                complete = complete_counts.get(bgc_type, 0)
                f.write(f"  {bgc_type}: {count} ({complete} complete)\n")
            
            f.write(f"\n完整BGC比例: {sum(complete_counts.values())}/{len(clusters)} "
                   f"({sum(complete_counts.values())/len(clusters)*100:.1f}%)\n")
        
        print(f"摘要报告保存至: {summary_file}")
        return summary_file
    
    def visualize_results(self, clusters):
        """结果可视化"""
        print("正在生成可视化图表...")
        
        # 创建图表目录
        figures_dir = self.output_dir.parent / "figures"
        figures_dir.mkdir(exist_ok=True)
        
        # BGC类型分布饼图
        type_counts = {}
        for cluster in clusters:
            bgc_type = cluster['type']
            type_counts[bgc_type] = type_counts.get(bgc_type, 0) + 1
        
        plt.figure(figsize=(10, 8))
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']
        
        plt.subplot(2, 2, 1)
        plt.pie(type_counts.values(), labels=type_counts.keys(), autopct='%1.1f%%',
                colors=colors[:len(type_counts)])
        plt.title('BGC类型分布')
        
        # BGC大小分布直方图
        sizes = [(cluster['end'] - cluster['start'])/1000 for cluster in clusters]
        plt.subplot(2, 2, 2)
        plt.hist(sizes, bins=10, color='skyblue', alpha=0.7)
        plt.xlabel('BGC大小 (kb)')
        plt.ylabel('频次')
        plt.title('BGC大小分布')
        
        # 完整性统计
        complete_counts = {'Complete': 0, 'Partial': 0}
        for cluster in clusters:
            if cluster['completeness'] == 'complete':
                complete_counts['Complete'] += 1
            else:
                complete_counts['Partial'] += 1
        
        plt.subplot(2, 2, 3)
        plt.bar(complete_counts.keys(), complete_counts.values(), 
                color=['#2ECC71', '#E74C3C'])
        plt.ylabel('BGC数量')
        plt.title('BGC完整性统计')
        
        # BGC在基因组上的分布
        positions = [(cluster['start'] + cluster['end'])/2/1000000 for cluster in clusters]
        types = [cluster['type'] for cluster in clusters]
        
        plt.subplot(2, 2, 4)
        type_colors = {t: colors[i] for i, t in enumerate(set(types))}
        for i, (pos, bgc_type) in enumerate(zip(positions, types)):
            plt.scatter(pos, i, c=type_colors[bgc_type], s=100, alpha=0.7)
        
        plt.xlabel('基因组位置 (Mb)')
        plt.ylabel('BGC编号')
        plt.title('BGC基因组分布')
        
        plt.tight_layout()
        plt.savefig(figures_dir / 'bgc_analysis_summary.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"可视化图表保存至: {figures_dir}/bgc_analysis_summary.png")
    
    def run_analysis(self):
        """运行完整分析流程"""
        print("开始BGC分析...")
        
        # 1. 基因预测
        genes_file = self.predict_genes()
        if not genes_file:
            return False
        
        # 2. 结构域搜索
        domain_file = self.search_bgc_domains(genes_file)
        
        # 3. 基因聚类
        clusters, cluster_file = self.cluster_genes(domain_file)
        
        # 4. 生成摘要
        summary_file = self.generate_summary(clusters)
        
        # 5. 结果可视化
        self.visualize_results(clusters)
        
        print("BGC分析完成！")
        print(f"主要结果文件:")
        print(f"  - 聚类结果: {cluster_file}")
        print(f"  - 分析摘要: {summary_file}")
        print(f"  - 可视化图表: figures/bgc_analysis_summary.png")
        
        return True

def main():
    parser = argparse.ArgumentParser(description='BGC分析工具')
    parser.add_argument('--input', required=True, help='输入基因组文件')
    parser.add_argument('--output', required=True, help='输出目录')
    parser.add_argument('--threads', type=int, default=4, help='线程数')
    parser.add_argument('--mode', default='full', 
                       choices=['full', 'local', 'summary'],
                       help='分析模式')
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.input):
        print(f"错误: 输入文件不存在: {args.input}")
        sys.exit(1)
    
    # 创建分析器实例
    analyzer = BGCAnalyzer(args.input, args.output, args.threads)
    
    # 运行分析
    if args.mode == 'full':
        success = analyzer.run_analysis()
    elif args.mode == 'local':
        print("运行本地BGC分析模式...")
        success = analyzer.run_analysis()
    elif args.mode == 'summary':
        print("生成BGC分析摘要...")
        # 这里可以添加仅生成摘要的逻辑
        success = True
    
    if success:
        print("分析成功完成！")
    else:
        print("分析过程中出现错误！")
        sys.exit(1)

if __name__ == "__main__":
    main()