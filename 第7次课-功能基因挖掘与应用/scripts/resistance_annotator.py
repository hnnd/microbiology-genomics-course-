#!/usr/bin/env python3
"""
Resistance Gene Annotator
用于抗生素抗性基因的注释和分析

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
import matplotlib.pyplot as plt
import seaborn as sns

class ResistanceAnnotator:
    """抗性基因注释主类"""
    
    def __init__(self, proteins_file, genome_file, output_dir, database='CARD'):
        self.proteins_file = proteins_file
        self.genome_file = genome_file
        self.output_dir = Path(output_dir)
        self.database = database
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 抗性机制分类
        self.resistance_mechanisms = {
            'enzymatic_inactivation': '酶修饰失活',
            'target_alteration': '靶点改变',
            'efflux_pump': '主动外排',
            'reduced_permeability': '通透性降低',
            'target_protection': '靶点保护',
            'metabolic_bypass': '代谢旁路'
        }
        
        # 抗生素类别
        self.antibiotic_classes = {
            'beta_lactam': 'β-内酰胺类',
            'aminoglycoside': '氨基糖苷类',
            'macrolide': '大环内酯类',
            'tetracycline': '四环素类',
            'quinolone': '喹诺酮类',
            'glycopeptide': '糖肽类',
            'chloramphenicol': '氯霉素类',
            'sulfonamide': '磺胺类'
        }
    
    def search_resistance_genes(self):
        """搜索抗性基因"""
        print("正在搜索抗性基因...")
        
        # 模拟抗性基因搜索结果
        # 实际应用中会使用BLAST搜索CARD数据库
        resistance_genes = [
            {
                'gene_id': 'gene_0045',
                'resistance_gene': 'mecA',
                'mechanism': 'target_alteration',
                'antibiotic_class': 'beta_lactam',
                'identity': 98.5,
                'coverage': 95.2,
                'evalue': 1e-120,
                'description': 'methicillin resistance protein'
            },
            {
                'gene_id': 'gene_0123',
                'resistance_gene': 'blaZ',
                'mechanism': 'enzymatic_inactivation',
                'antibiotic_class': 'beta_lactam',
                'identity': 96.8,
                'coverage': 98.1,
                'evalue': 1e-95,
                'description': 'beta-lactamase'
            },
            {
                'gene_id': 'gene_0234',
                'resistance_gene': 'aac(6\')-Ie-aph(2\'\')-Ia',
                'mechanism': 'enzymatic_inactivation',
                'antibiotic_class': 'aminoglycoside',
                'identity': 99.2,
                'coverage': 100.0,
                'evalue': 1e-150,
                'description': 'aminoglycoside acetyltransferase/phosphotransferase'
            },
            {
                'gene_id': 'gene_0345',
                'resistance_gene': 'ermA',
                'mechanism': 'target_alteration',
                'antibiotic_class': 'macrolide',
                'identity': 97.3,
                'coverage': 94.5,
                'evalue': 1e-88,
                'description': 'rRNA methyltransferase'
            },
            {
                'gene_id': 'gene_0456',
                'resistance_gene': 'tetK',
                'mechanism': 'efflux_pump',
                'antibiotic_class': 'tetracycline',
                'identity': 95.7,
                'coverage': 92.3,
                'evalue': 1e-76,
                'description': 'tetracycline efflux protein'
            },
            {
                'gene_id': 'gene_0567',
                'resistance_gene': 'norA',
                'mechanism': 'efflux_pump',
                'antibiotic_class': 'quinolone',
                'identity': 94.2,
                'coverage': 89.7,
                'evalue': 1e-65,
                'description': 'quinolone efflux pump'
            }
        ]
        
        # 保存搜索结果
        results_file = self.output_dir / "resistance_genes.txt"
        with open(results_file, 'w') as f:
            f.write("Gene_ID\tResistance_Gene\tMechanism\tAntibiotic_Class\t"
                   "Identity\tCoverage\tE_value\tDescription\n")
            for gene in resistance_genes:
                f.write(f"{gene['gene_id']}\t{gene['resistance_gene']}\t"
                       f"{gene['mechanism']}\t{gene['antibiotic_class']}\t"
                       f"{gene['identity']}\t{gene['coverage']}\t"
                       f"{gene['evalue']}\t{gene['description']}\n")
        
        print(f"发现 {len(resistance_genes)} 个抗性基因")
        print(f"结果保存至: {results_file}")
        return resistance_genes, results_file
    
    def classify_mechanisms(self, resistance_genes):
        """按抗性机制分类"""
        print("正在按抗性机制分类...")
        
        mechanism_stats = {}
        for gene in resistance_genes:
            mechanism = gene['mechanism']
            if mechanism not in mechanism_stats:
                mechanism_stats[mechanism] = []
            mechanism_stats[mechanism].append(gene)
        
        # 保存分类结果
        classify_file = self.output_dir / "mechanism_summary.txt"
        with open(classify_file, 'w') as f:
            f.write("=== 抗性机制分类统计 ===\n\n")
            for mechanism, genes in mechanism_stats.items():
                mechanism_name = self.resistance_mechanisms.get(mechanism, mechanism)
                f.write(f"{mechanism_name} ({mechanism}): {len(genes)} 个基因\n")
                for gene in genes:
                    f.write(f"  - {gene['resistance_gene']} ({gene['gene_id']})\n")
                f.write("\n")
        
        print(f"机制分类结果保存至: {classify_file}")
        return mechanism_stats
    
    def predict_resistance_profile(self, resistance_genes):
        """预测抗生素抗性谱"""
        print("正在预测抗生素抗性谱...")
        
        # 统计各类抗生素的抗性基因
        antibiotic_resistance = {}
        for gene in resistance_genes:
            antibiotic_class = gene['antibiotic_class']
            if antibiotic_class not in antibiotic_resistance:
                antibiotic_resistance[antibiotic_class] = []
            antibiotic_resistance[antibiotic_class].append(gene)
        
        # 生成抗性谱预测
        profile_file = self.output_dir / "resistance_profile.txt"
        with open(profile_file, 'w') as f:
            f.write("=== 抗生素抗性谱预测 ===\n\n")
            f.write("抗生素类别\t抗性状态\t相关基因\t置信度\n")
            
            for antibiotic_class in self.antibiotic_classes:
                class_name = self.antibiotic_classes[antibiotic_class]
                if antibiotic_class in antibiotic_resistance:
                    genes = antibiotic_resistance[antibiotic_class]
                    gene_names = [g['resistance_gene'] for g in genes]
                    avg_identity = sum(g['identity'] for g in genes) / len(genes)
                    
                    if avg_identity >= 95:
                        confidence = "高"
                        status = "抗性"
                    elif avg_identity >= 80:
                        confidence = "中"
                        status = "可能抗性"
                    else:
                        confidence = "低"
                        status = "不确定"
                    
                    f.write(f"{class_name}\t{status}\t{','.join(gene_names)}\t{confidence}\n")
                else:
                    f.write(f"{class_name}\t敏感\t-\t-\n")
        
        print(f"抗性谱预测保存至: {profile_file}")
        return antibiotic_resistance
    
    def analyze_mobile_elements(self, resistance_genes):
        """分析可移动遗传元件"""
        print("正在分析可移动遗传元件...")
        
        # 模拟可移动元件分析
        mobile_elements = [
            {
                'element_type': 'plasmid',
                'element_name': 'pSA01',
                'resistance_genes': ['mecA', 'blaZ'],
                'size': 45000,
                'mobility': 'high'
            },
            {
                'element_type': 'transposon',
                'element_name': 'Tn554',
                'resistance_genes': ['ermA'],
                'size': 6800,
                'mobility': 'medium'
            }
        ]
        
        mobile_file = self.output_dir / "mobile_elements.txt"
        with open(mobile_file, 'w') as f:
            f.write("=== 可移动遗传元件分析 ===\n\n")
            f.write("元件类型\t元件名称\t携带抗性基因\t大小(bp)\t可移动性\n")
            for element in mobile_elements:
                genes_str = ','.join(element['resistance_genes'])
                f.write(f"{element['element_type']}\t{element['element_name']}\t"
                       f"{genes_str}\t{element['size']}\t{element['mobility']}\n")
        
        print(f"可移动元件分析保存至: {mobile_file}")
        return mobile_elements
    
    def visualize_results(self, resistance_genes, mechanism_stats):
        """结果可视化"""
        print("正在生成可视化图表...")
        
        # 创建图表目录
        figures_dir = self.output_dir.parent / "figures"
        figures_dir.mkdir(exist_ok=True)
        
        plt.figure(figsize=(15, 10))
        
        # 1. 抗性机制分布饼图
        plt.subplot(2, 3, 1)
        mechanism_counts = {self.resistance_mechanisms.get(k, k): len(v) 
                           for k, v in mechanism_stats.items()}
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']
        plt.pie(mechanism_counts.values(), labels=mechanism_counts.keys(), 
                autopct='%1.1f%%', colors=colors)
        plt.title('抗性机制分布')
        
        # 2. 抗生素类别分布
        plt.subplot(2, 3, 2)
        antibiotic_counts = {}
        for gene in resistance_genes:
            antibiotic_class = self.antibiotic_classes.get(gene['antibiotic_class'], 
                                                          gene['antibiotic_class'])
            antibiotic_counts[antibiotic_class] = antibiotic_counts.get(antibiotic_class, 0) + 1
        
        plt.bar(range(len(antibiotic_counts)), list(antibiotic_counts.values()),
                color='skyblue')
        plt.xticks(range(len(antibiotic_counts)), list(antibiotic_counts.keys()),
                  rotation=45, ha='right')
        plt.ylabel('基因数量')
        plt.title('抗生素类别分布')
        
        # 3. 序列相似性分布
        plt.subplot(2, 3, 3)
        identities = [gene['identity'] for gene in resistance_genes]
        plt.hist(identities, bins=10, color='lightgreen', alpha=0.7)
        plt.xlabel('序列一致性 (%)')
        plt.ylabel('基因数量')
        plt.title('序列相似性分布')
        
        # 4. 覆盖度分布
        plt.subplot(2, 3, 4)
        coverages = [gene['coverage'] for gene in resistance_genes]
        plt.hist(coverages, bins=10, color='orange', alpha=0.7)
        plt.xlabel('覆盖度 (%)')
        plt.ylabel('基因数量')
        plt.title('序列覆盖度分布')
        
        # 5. E值分布
        plt.subplot(2, 3, 5)
        evalues = [-np.log10(gene['evalue']) for gene in resistance_genes]
        plt.hist(evalues, bins=10, color='purple', alpha=0.7)
        plt.xlabel('-log10(E-value)')
        plt.ylabel('基因数量')
        plt.title('E值分布')
        
        # 6. 抗性基因热图
        plt.subplot(2, 3, 6)
        # 创建抗性基因-机制矩阵
        genes = [gene['resistance_gene'] for gene in resistance_genes]
        mechanisms = [self.resistance_mechanisms.get(gene['mechanism'], gene['mechanism']) 
                     for gene in resistance_genes]
        
        # 简化的热图数据
        unique_genes = list(set(genes))[:10]  # 取前10个基因
        unique_mechanisms = list(set(mechanisms))
        
        heatmap_data = np.random.rand(len(unique_genes), len(unique_mechanisms))
        
        sns.heatmap(heatmap_data, xticklabels=unique_mechanisms, 
                   yticklabels=unique_genes, cmap='YlOrRd', cbar=True)
        plt.title('抗性基因-机制关联')
        plt.xticks(rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(figures_dir / 'resistance_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"可视化图表保存至: {figures_dir}/resistance_analysis.png")

def main():
    parser = argparse.ArgumentParser(description='抗性基因注释工具')
    parser.add_argument('--proteins', required=True, help='蛋白质序列文件')
    parser.add_argument('--genome', required=True, help='基因组序列文件')
    parser.add_argument('--output', required=True, help='输出目录')
    parser.add_argument('--database', default='CARD', 
                       choices=['CARD', 'ResFinder', 'ARG-ANNOT'],
                       help='抗性基因数据库')
    parser.add_argument('--identity', type=float, default=80.0,
                       help='序列一致性阈值')
    parser.add_argument('--coverage', type=float, default=70.0,
                       help='覆盖度阈值')
    parser.add_argument('--evalue', type=float, default=1e-10,
                       help='E值阈值')
    parser.add_argument('--mode', default='annotate',
                       choices=['annotate', 'classify', 'predict'],
                       help='分析模式')
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.proteins):
        print(f"错误: 蛋白质文件不存在: {args.proteins}")
        sys.exit(1)
    
    if not os.path.exists(args.genome):
        print(f"错误: 基因组文件不存在: {args.genome}")
        sys.exit(1)
    
    # 创建注释器实例
    annotator = ResistanceAnnotator(args.proteins, args.genome, 
                                   args.output, args.database)
    
    # 运行分析
    if args.mode == 'annotate':
        # 完整注释流程
        resistance_genes, results_file = annotator.search_resistance_genes()
        mechanism_stats = annotator.classify_mechanisms(resistance_genes)
        antibiotic_resistance = annotator.predict_resistance_profile(resistance_genes)
        mobile_elements = annotator.analyze_mobile_elements(resistance_genes)
        annotator.visualize_results(resistance_genes, mechanism_stats)
        
        print("抗性基因注释完成！")
        print(f"主要结果文件:")
        print(f"  - 抗性基因列表: {results_file}")
        print(f"  - 机制分类: {annotator.output_dir}/mechanism_summary.txt")
        print(f"  - 抗性谱预测: {annotator.output_dir}/resistance_profile.txt")
        
    elif args.mode == 'classify':
        # 仅分类模式
        resistance_genes, _ = annotator.search_resistance_genes()
        mechanism_stats = annotator.classify_mechanisms(resistance_genes)
        print("抗性机制分类完成！")
        
    elif args.mode == 'predict':
        # 仅预测模式
        resistance_genes, _ = annotator.search_resistance_genes()
        antibiotic_resistance = annotator.predict_resistance_profile(resistance_genes)
        print("抗性谱预测完成！")

if __name__ == "__main__":
    # 添加必要的导入
    import numpy as np
    main()