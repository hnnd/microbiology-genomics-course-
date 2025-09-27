#!/usr/bin/env python3
"""
Virulence Factor Scanner
用于毒力因子的扫描和分析

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
import numpy as np

class VirulenceScanner:
    """毒力因子扫描主类"""
    
    def __init__(self, genome_file, proteins_file, output_dir, database='VFDB'):
        self.genome_file = genome_file
        self.proteins_file = proteins_file
        self.output_dir = Path(output_dir)
        self.database = database
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 毒力因子分类
        self.virulence_categories = {
            'adhesion': '粘附因子',
            'invasion': '侵袭因子',
            'toxin': '毒素',
            'immune_evasion': '免疫逃逸',
            'secretion_system': '分泌系统',
            'regulation': '毒力调控',
            'iron_uptake': '铁摄取',
            'stress_survival': '应激生存'
        }
        
        # 毒力因子功能描述
        self.virulence_functions = {
            'hla': 'α-溶血素，细胞膜孔道形成',
            'hlb': 'β-溶血素，磷脂酶活性',
            'tst': '中毒性休克综合征毒素',
            'sea': '肠毒素A，食物中毒',
            'spa': '蛋白A，免疫球蛋白结合',
            'clfA': '凝集因子A，纤维蛋白原结合',
            'fnbA': '纤连蛋白结合蛋白A',
            'icaA': '胞间粘附蛋白A，生物膜形成',
            'sarA': '毒力调控蛋白',
            'agr': '群体感应调控系统'
        }
    
    def scan_virulence_factors(self):
        """扫描毒力因子"""
        print("正在扫描毒力因子...")
        
        # 模拟毒力因子搜索结果
        # 实际应用中会使用BLAST搜索VFDB数据库
        virulence_factors = [
            {
                'gene_id': 'gene_0078',
                'vf_gene': 'hla',
                'category': 'toxin',
                'function': 'alpha-hemolysin',
                'identity': 98.7,
                'coverage': 99.2,
                'evalue': 1e-145,
                'pathogen': 'Staphylococcus aureus',
                'description': 'pore-forming cytotoxin'
            },
            {
                'gene_id': 'gene_0156',
                'vf_gene': 'hlb',
                'category': 'toxin',
                'function': 'beta-hemolysin',
                'identity': 96.3,
                'coverage': 94.8,
                'evalue': 1e-98,
                'pathogen': 'Staphylococcus aureus',
                'description': 'sphingomyelinase C'
            },
            {
                'gene_id': 'gene_0234',
                'vf_gene': 'tst',
                'category': 'toxin',
                'function': 'toxic shock syndrome toxin',
                'identity': 99.1,
                'coverage': 98.5,
                'evalue': 1e-156,
                'pathogen': 'Staphylococcus aureus',
                'description': 'superantigen toxin'
            },
            {
                'gene_id': 'gene_0345',
                'vf_gene': 'spa',
                'category': 'immune_evasion',
                'function': 'protein A',
                'identity': 97.8,
                'coverage': 96.7,
                'evalue': 1e-112,
                'pathogen': 'Staphylococcus aureus',
                'description': 'immunoglobulin-binding protein'
            },
            {
                'gene_id': 'gene_0456',
                'vf_gene': 'clfA',
                'category': 'adhesion',
                'function': 'clumping factor A',
                'identity': 95.4,
                'coverage': 92.1,
                'evalue': 1e-87,
                'pathogen': 'Staphylococcus aureus',
                'description': 'fibrinogen-binding protein'
            },
            {
                'gene_id': 'gene_0567',
                'vf_gene': 'fnbA',
                'category': 'adhesion',
                'function': 'fibronectin-binding protein A',
                'identity': 94.2,
                'coverage': 89.3,
                'evalue': 1e-76,
                'pathogen': 'Staphylococcus aureus',
                'description': 'extracellular matrix binding'
            },
            {
                'gene_id': 'gene_0678',
                'vf_gene': 'icaA',
                'category': 'adhesion',
                'function': 'intercellular adhesin A',
                'identity': 96.8,
                'coverage': 95.4,
                'evalue': 1e-102,
                'pathogen': 'Staphylococcus aureus',
                'description': 'biofilm formation'
            },
            {
                'gene_id': 'gene_0789',
                'vf_gene': 'sarA',
                'category': 'regulation',
                'function': 'staphylococcal accessory regulator A',
                'identity': 98.2,
                'coverage': 97.6,
                'evalue': 1e-89,
                'pathogen': 'Staphylococcus aureus',
                'description': 'global virulence regulator'
            }
        ]
        
        # 保存扫描结果
        results_file = self.output_dir / "virulence_factors.txt"
        with open(results_file, 'w') as f:
            f.write("Gene_ID\tVF_Gene\tCategory\tFunction\tIdentity\tCoverage\t"
                   "E_value\tPathogen\tDescription\n")
            for vf in virulence_factors:
                f.write(f"{vf['gene_id']}\t{vf['vf_gene']}\t{vf['category']}\t"
                       f"{vf['function']}\t{vf['identity']}\t{vf['coverage']}\t"
                       f"{vf['evalue']}\t{vf['pathogen']}\t{vf['description']}\n")
        
        print(f"发现 {len(virulence_factors)} 个毒力因子")
        print(f"结果保存至: {results_file}")
        return virulence_factors, results_file
    
    def predict_pathogenicity_islands(self):
        """预测病原性岛"""
        print("正在预测病原性岛...")
        
        # 模拟病原性岛预测结果
        pathogenicity_islands = [
            {
                'island_id': 'PAI_1',
                'start': 1234567,
                'end': 1267890,
                'size': 33323,
                'gc_content': 31.2,
                'virulence_genes': ['tst', 'sea', 'seb'],
                'mobility_genes': ['int', 'tra'],
                'confidence': 'high'
            },
            {
                'island_id': 'PAI_2',
                'start': 2145678,
                'end': 2167890,
                'size': 22212,
                'gc_content': 32.8,
                'virulence_genes': ['hla', 'hlb'],
                'mobility_genes': ['int'],
                'confidence': 'medium'
            }
        ]
        
        # 保存预测结果
        islands_file = self.output_dir / "pathogenicity_islands.txt"
        with open(islands_file, 'w') as f:
            f.write("Island_ID\tStart\tEnd\tSize\tGC_Content\tVirulence_Genes\t"
                   "Mobility_Genes\tConfidence\n")
            for island in pathogenicity_islands:
                vf_genes = ','.join(island['virulence_genes'])
                mob_genes = ','.join(island['mobility_genes'])
                f.write(f"{island['island_id']}\t{island['start']}\t{island['end']}\t"
                       f"{island['size']}\t{island['gc_content']}\t{vf_genes}\t"
                       f"{mob_genes}\t{island['confidence']}\n")
        
        print(f"预测到 {len(pathogenicity_islands)} 个病原性岛")
        print(f"结果保存至: {islands_file}")
        return pathogenicity_islands
    
    def analyze_gene_distribution(self, virulence_factors):
        """分析毒力基因分布"""
        print("正在分析毒力基因分布...")
        
        # 模拟基因组位置信息
        gene_positions = {}
        genome_length = 2800000  # 假设基因组长度
        
        for i, vf in enumerate(virulence_factors):
            # 模拟基因位置
            position = (i + 1) * (genome_length // len(virulence_factors)) + \
                      np.random.randint(-50000, 50000)
            gene_positions[vf['gene_id']] = {
                'position': max(1, position),
                'vf_gene': vf['vf_gene'],
                'category': vf['category']
            }
        
        # 保存分布分析
        distribution_file = self.output_dir / "gene_distribution.txt"
        with open(distribution_file, 'w') as f:
            f.write("Gene_ID\tVF_Gene\tCategory\tPosition\tRegion\n")
            for gene_id, info in gene_positions.items():
                region = f"Region_{info['position'] // 500000 + 1}"
                f.write(f"{gene_id}\t{info['vf_gene']}\t{info['category']}\t"
                       f"{info['position']}\t{region}\n")
        
        print(f"基因分布分析保存至: {distribution_file}")
        return gene_positions
    
    def assess_pathogenicity(self, virulence_factors, resistance_genes_file=None):
        """综合致病性评估"""
        print("正在进行综合致病性评估...")
        
        # 毒力因子评分
        category_scores = {
            'toxin': 3.0,
            'adhesion': 2.0,
            'invasion': 2.5,
            'immune_evasion': 2.0,
            'secretion_system': 1.5,
            'regulation': 1.0,
            'iron_uptake': 1.0,
            'stress_survival': 0.5
        }
        
        # 计算毒力评分
        virulence_score = 0
        category_counts = {}
        
        for vf in virulence_factors:
            category = vf['category']
            score = category_scores.get(category, 1.0)
            identity_factor = vf['identity'] / 100.0
            virulence_score += score * identity_factor
            
            category_counts[category] = category_counts.get(category, 0) + 1
        
        # 读取抗性基因信息（如果提供）
        resistance_score = 0
        if resistance_genes_file and os.path.exists(resistance_genes_file):
            with open(resistance_genes_file, 'r') as f:
                resistance_lines = f.readlines()[1:]  # 跳过标题行
                resistance_score = len(resistance_lines) * 0.5  # 每个抗性基因0.5分
        
        # 综合评分
        total_score = virulence_score + resistance_score
        
        # 致病性等级
        if total_score >= 15:
            pathogenicity_level = "高致病性"
        elif total_score >= 10:
            pathogenicity_level = "中等致病性"
        elif total_score >= 5:
            pathogenicity_level = "低致病性"
        else:
            pathogenicity_level = "非致病性"
        
        # 保存评估结果
        assessment_file = self.output_dir / "pathogenicity_assessment.txt"
        with open(assessment_file, 'w') as f:
            f.write("=== 综合致病性评估报告 ===\n\n")
            f.write(f"毒力因子总数: {len(virulence_factors)}\n")
            f.write(f"毒力评分: {virulence_score:.2f}\n")
            f.write(f"抗性评分: {resistance_score:.2f}\n")
            f.write(f"综合评分: {total_score:.2f}\n")
            f.write(f"致病性等级: {pathogenicity_level}\n\n")
            
            f.write("毒力因子分类统计:\n")
            for category, count in category_counts.items():
                category_name = self.virulence_categories.get(category, category)
                f.write(f"  {category_name}: {count} 个\n")
            
            f.write("\n风险评估:\n")
            if total_score >= 15:
                f.write("- 具有高致病性潜力，需要特别关注\n")
                f.write("- 建议采取严格的生物安全措施\n")
            elif total_score >= 10:
                f.write("- 具有中等致病性，需要适当防护\n")
                f.write("- 建议监测其传播和进化\n")
            else:
                f.write("- 致病性相对较低\n")
                f.write("- 仍需关注其潜在风险\n")
        
        print(f"致病性评估完成，等级: {pathogenicity_level}")
        print(f"评估报告保存至: {assessment_file}")
        return total_score, pathogenicity_level
    
    def visualize_results(self, virulence_factors, gene_positions):
        """结果可视化"""
        print("正在生成可视化图表...")
        
        # 创建图表目录
        figures_dir = self.output_dir.parent / "figures"
        figures_dir.mkdir(exist_ok=True)
        
        plt.figure(figsize=(16, 12))
        
        # 1. 毒力因子分类分布
        plt.subplot(2, 3, 1)
        category_counts = {}
        for vf in virulence_factors:
            category = self.virulence_categories.get(vf['category'], vf['category'])
            category_counts[category] = category_counts.get(category, 0) + 1
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD', '#FFB6C1', '#98FB98']
        plt.pie(category_counts.values(), labels=category_counts.keys(), 
                autopct='%1.1f%%', colors=colors[:len(category_counts)])
        plt.title('毒力因子分类分布')
        
        # 2. 序列相似性分布
        plt.subplot(2, 3, 2)
        identities = [vf['identity'] for vf in virulence_factors]
        plt.hist(identities, bins=8, color='lightcoral', alpha=0.7)
        plt.xlabel('序列一致性 (%)')
        plt.ylabel('基因数量')
        plt.title('序列相似性分布')
        
        # 3. 基因组分布图
        plt.subplot(2, 3, 3)
        positions = [info['position']/1000000 for info in gene_positions.values()]
        categories = [info['category'] for info in gene_positions.values()]
        
        category_colors = {cat: colors[i] for i, cat in enumerate(set(categories))}
        for pos, cat in zip(positions, categories):
            plt.scatter(pos, 0, c=category_colors[cat], s=100, alpha=0.7)
        
        plt.xlabel('基因组位置 (Mb)')
        plt.ylabel('')
        plt.title('毒力基因基因组分布')
        plt.yticks([])
        
        # 4. 毒力因子功能热图
        plt.subplot(2, 3, 4)
        # 创建功能-类别矩阵
        functions = [vf['vf_gene'] for vf in virulence_factors]
        categories = [vf['category'] for vf in virulence_factors]
        
        unique_functions = list(set(functions))[:8]  # 取前8个功能
        unique_categories = list(set(categories))
        
        # 创建热图数据
        heatmap_data = np.zeros((len(unique_functions), len(unique_categories)))
        for i, func in enumerate(unique_functions):
            for j, cat in enumerate(unique_categories):
                # 计算该功能在该类别中的出现次数
                count = sum(1 for vf in virulence_factors 
                           if vf['vf_gene'] == func and vf['category'] == cat)
                heatmap_data[i, j] = count
        
        sns.heatmap(heatmap_data, xticklabels=unique_categories, 
                   yticklabels=unique_functions, cmap='Reds', annot=True, fmt='g')
        plt.title('毒力基因-功能关联')
        plt.xticks(rotation=45, ha='right')
        
        # 5. 致病性评分雷达图
        plt.subplot(2, 3, 5)
        category_scores = {}
        for vf in virulence_factors:
            category = vf['category']
            if category not in category_scores:
                category_scores[category] = []
            category_scores[category].append(vf['identity'])
        
        # 计算各类别平均得分
        avg_scores = {cat: np.mean(scores) for cat, scores in category_scores.items()}
        
        categories_radar = list(avg_scores.keys())
        scores_radar = list(avg_scores.values())
        
        # 简化的条形图代替雷达图
        plt.bar(range(len(categories_radar)), scores_radar, color='lightblue')
        plt.xticks(range(len(categories_radar)), categories_radar, rotation=45, ha='right')
        plt.ylabel('平均相似性 (%)')
        plt.title('各类别毒力因子质量')
        
        # 6. E值分布
        plt.subplot(2, 3, 6)
        evalues = [-np.log10(vf['evalue']) for vf in virulence_factors]
        plt.hist(evalues, bins=8, color='gold', alpha=0.7)
        plt.xlabel('-log10(E-value)')
        plt.ylabel('基因数量')
        plt.title('E值分布')
        
        plt.tight_layout()
        plt.savefig(figures_dir / 'virulence_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"可视化图表保存至: {figures_dir}/virulence_analysis.png")

def main():
    parser = argparse.ArgumentParser(description='毒力因子扫描工具')
    parser.add_argument('--genome', required=True, help='基因组序列文件')
    parser.add_argument('--proteins', required=True, help='蛋白质序列文件')
    parser.add_argument('--output', required=True, help='输出目录')
    parser.add_argument('--database', default='VFDB',
                       choices=['VFDB', 'MVirDB', 'VirulenceFinder'],
                       help='毒力因子数据库')
    parser.add_argument('--identity', type=float, default=70.0,
                       help='序列一致性阈值')
    parser.add_argument('--coverage', type=float, default=60.0,
                       help='覆盖度阈值')
    parser.add_argument('--evalue', type=float, default=1e-5,
                       help='E值阈值')
    parser.add_argument('--mode', default='scan',
                       choices=['scan', 'island', 'distribution', 'assessment'],
                       help='分析模式')
    parser.add_argument('--virulence', help='毒力因子结果文件（用于评估模式）')
    parser.add_argument('--resistance', help='抗性基因结果文件（用于评估模式）')
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.genome):
        print(f"错误: 基因组文件不存在: {args.genome}")
        sys.exit(1)
    
    if not os.path.exists(args.proteins):
        print(f"错误: 蛋白质文件不存在: {args.proteins}")
        sys.exit(1)
    
    # 创建扫描器实例
    scanner = VirulenceScanner(args.genome, args.proteins, args.output, args.database)
    
    # 运行分析
    if args.mode == 'scan':
        # 完整扫描流程
        virulence_factors, results_file = scanner.scan_virulence_factors()
        pathogenicity_islands = scanner.predict_pathogenicity_islands()
        gene_positions = scanner.analyze_gene_distribution(virulence_factors)
        total_score, pathogenicity_level = scanner.assess_pathogenicity(
            virulence_factors, args.resistance)
        scanner.visualize_results(virulence_factors, gene_positions)
        
        print("毒力因子扫描完成！")
        print(f"主要结果文件:")
        print(f"  - 毒力因子列表: {results_file}")
        print(f"  - 病原性岛: {scanner.output_dir}/pathogenicity_islands.txt")
        print(f"  - 致病性评估: {scanner.output_dir}/pathogenicity_assessment.txt")
        
    elif args.mode == 'island':
        # 仅预测病原性岛
        pathogenicity_islands = scanner.predict_pathogenicity_islands()
        print("病原性岛预测完成！")
        
    elif args.mode == 'distribution':
        # 仅分析基因分布
        if args.virulence and os.path.exists(args.virulence):
            # 从文件读取毒力因子信息
            virulence_factors = []
            with open(args.virulence, 'r') as f:
                lines = f.readlines()[1:]  # 跳过标题行
                for line in lines:
                    parts = line.strip().split('\t')
                    if len(parts) >= 9:
                        vf = {
                            'gene_id': parts[0],
                            'vf_gene': parts[1],
                            'category': parts[2],
                            'identity': float(parts[4]),
                            'evalue': float(parts[6])
                        }
                        virulence_factors.append(vf)
            
            gene_positions = scanner.analyze_gene_distribution(virulence_factors)
            print("基因分布分析完成！")
        else:
            print("错误: 需要提供毒力因子结果文件")
            sys.exit(1)
            
    elif args.mode == 'assessment':
        # 仅进行致病性评估
        if args.virulence and os.path.exists(args.virulence):
            # 从文件读取毒力因子信息
            virulence_factors = []
            with open(args.virulence, 'r') as f:
                lines = f.readlines()[1:]  # 跳过标题行
                for line in lines:
                    parts = line.strip().split('\t')
                    if len(parts) >= 9:
                        vf = {
                            'gene_id': parts[0],
                            'vf_gene': parts[1],
                            'category': parts[2],
                            'identity': float(parts[4])
                        }
                        virulence_factors.append(vf)
            
            total_score, pathogenicity_level = scanner.assess_pathogenicity(
                virulence_factors, args.resistance)
            print("致病性评估完成！")
        else:
            print("错误: 需要提供毒力因子结果文件")
            sys.exit(1)

if __name__ == "__main__":
    main()