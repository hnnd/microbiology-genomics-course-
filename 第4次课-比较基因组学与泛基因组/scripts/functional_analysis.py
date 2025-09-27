#!/usr/bin/env python3
"""
功能分类分析脚本
用于分析泛基因组中基因的功能分类和注释

作者: 微生物基因组学课程组
版本: 1.0.0
日期: 2025
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='泛基因组功能分类分析')
    parser.add_argument('--input', required=True,
                       help='基因存在/缺失矩阵文件 (gene_presence_absence.csv)')
    parser.add_argument('--mode', choices=['basic', 'accessory'], default='basic',
                       help='分析模式: basic=基础功能分析, accessory=辅助基因组分析')
    parser.add_argument('--output', default='results',
                       help='输出目录 (默认: results)')
    return parser.parse_args()

def load_gene_presence_data(input_file):
    """加载基因存在/缺失数据"""
    print(f"加载基因存在/缺失数据: {input_file}")
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"输入文件不存在: {input_file}")
    
    # 读取数据
    df = pd.read_csv(input_file, low_memory=False)
    
    print(f"✅ 数据加载完成")
    print(f"   基因家族数量: {len(df)}")
    print(f"   数据列数: {len(df.columns)}")
    
    # 获取菌株列（跳过前14列的注释信息）
    strain_columns = df.columns[14:]
    print(f"   菌株数量: {len(strain_columns)}")
    
    return df, strain_columns

def classify_genes_by_presence(df, strain_columns):
    """根据存在模式分类基因"""
    print("根据存在模式分类基因...")
    
    num_strains = len(strain_columns)
    gene_classification = {
        'core': [],
        'accessory': [],
        'unique': []
    }
    
    for idx, row in df.iterrows():
        # 计算基因在多少个菌株中存在
        presence_count = sum(1 for col in strain_columns if pd.notna(row[col]) and row[col] != '')
        
        gene_info = {
            'gene_id': row['Gene'],
            'annotation': row['Annotation'] if 'Annotation' in row else 'Unknown',
            'presence_count': presence_count,
            'frequency': presence_count / num_strains
        }
        
        if presence_count == num_strains:
            gene_classification['core'].append(gene_info)
        elif presence_count == 1:
            gene_classification['unique'].append(gene_info)
        else:
            gene_classification['accessory'].append(gene_info)
    
    print(f"✅ 基因分类完成:")
    print(f"   核心基因: {len(gene_classification['core'])}")
    print(f"   辅助基因: {len(gene_classification['accessory'])}")
    print(f"   独特基因: {len(gene_classification['unique'])}")
    
    return gene_classification

def analyze_functional_categories(gene_classification, output_dir):
    """分析功能分类"""
    print("分析基因功能分类...")
    
    # 简化的功能分类（基于注释关键词）
    functional_categories = {
        'Metabolism': ['metabol', 'synthase', 'dehydrogenase', 'oxidase', 'reductase', 'kinase'],
        'Transport': ['transport', 'permease', 'channel', 'pump', 'porter'],
        'Transcription': ['transcription', 'RNA polymerase', 'sigma', 'regulator'],
        'Translation': ['ribosom', 'tRNA', 'translation', 'elongation', 'initiation'],
        'DNA Replication': ['DNA', 'replication', 'helicase', 'polymerase', 'primase'],
        'Cell Wall': ['cell wall', 'peptidoglycan', 'murein', 'penicillin'],
        'Virulence': ['virulence', 'toxin', 'pathogen', 'invasion'],
        'Resistance': ['resistance', 'antibiotic', 'drug', 'efflux'],
        'Stress Response': ['stress', 'heat shock', 'cold shock', 'SOS'],
        'Unknown': ['hypothetical', 'unknown', 'uncharacterized']
    }
    
    # 分析每类基因的功能分布
    results = {}
    
    for gene_type in ['core', 'accessory', 'unique']:
        genes = gene_classification[gene_type]
        category_counts = Counter()
        
        for gene in genes:
            annotation = gene['annotation'].lower()
            categorized = False
            
            for category, keywords in functional_categories.items():
                if any(keyword in annotation for keyword in keywords):
                    category_counts[category] += 1
                    categorized = True
                    break
            
            if not categorized:
                category_counts['Other'] += 1
        
        results[gene_type] = category_counts
    
    # 保存功能分类结果
    func_file = os.path.join(output_dir, 'functional_classification.txt')
    with open(func_file, 'w') as f:
        f.write("=== 泛基因组功能分类分析 ===\n\n")
        
        for gene_type in ['core', 'accessory', 'unique']:
            f.write(f"=== {gene_type.upper()}基因功能分类 ===\n")
            total = sum(results[gene_type].values())
            
            for category, count in results[gene_type].most_common():
                percentage = count / total * 100 if total > 0 else 0
                f.write(f"{category}: {count} ({percentage:.1f}%)\n")
            f.write(f"总计: {total}\n\n")
    
    print(f"功能分类结果已保存到: {func_file}")
    
    # 生成功能分类图表
    create_functional_plots(results, output_dir)
    
    return results

def create_functional_plots(functional_results, output_dir):
    """创建功能分类图表"""
    print("生成功能分类图表...")
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 1. 饼图显示各类基因的功能分布
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    gene_types = ['core', 'accessory', 'unique']
    titles = ['核心基因功能分类', '辅助基因功能分类', '独特基因功能分类']
    
    for i, (gene_type, title) in enumerate(zip(gene_types, titles)):
        if functional_results[gene_type]:
            categories = list(functional_results[gene_type].keys())
            counts = list(functional_results[gene_type].values())
            
            # 只显示前8个最多的类别
            if len(categories) > 8:
                top_categories = categories[:7]
                top_counts = counts[:7]
                other_count = sum(counts[7:])
                top_categories.append('其他')
                top_counts.append(other_count)
            else:
                top_categories = categories
                top_counts = counts
            
            axes[i].pie(top_counts, labels=top_categories, autopct='%1.1f%%', startangle=90)
            axes[i].set_title(title, fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    pie_file = os.path.join(output_dir, 'functional_pie_charts.png')
    plt.savefig(pie_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"功能分类饼图已保存到: {pie_file}")
    
    # 2. 柱状图比较不同基因类型的功能分布
    all_categories = set()
    for results in functional_results.values():
        all_categories.update(results.keys())
    
    all_categories = sorted(list(all_categories))
    
    # 准备数据
    data_for_plot = []
    for category in all_categories:
        for gene_type in gene_types:
            count = functional_results[gene_type].get(category, 0)
            data_for_plot.append({
                'Category': category,
                'Gene_Type': gene_type.capitalize(),
                'Count': count
            })
    
    plot_df = pd.DataFrame(data_for_plot)
    
    # 创建柱状图
    plt.figure(figsize=(14, 8))
    sns.barplot(data=plot_df, x='Category', y='Count', hue='Gene_Type')
    plt.title('不同基因类型的功能分类比较', fontsize=14, fontweight='bold')
    plt.xlabel('功能类别', fontsize=12)
    plt.ylabel('基因数量', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='基因类型')
    plt.tight_layout()
    
    bar_file = os.path.join(output_dir, 'functional_comparison.png')
    plt.savefig(bar_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"功能分类比较图已保存到: {bar_file}")

def analyze_strain_specific_genes(df, strain_columns, output_dir):
    """分析菌株特异性基因"""
    print("分析菌株特异性基因...")
    
    strain_specific = {}
    
    for strain in strain_columns:
        specific_genes = []
        
        for idx, row in df.iterrows():
            # 检查是否只在当前菌株中存在
            presence_in_strains = []
            for col in strain_columns:
                if pd.notna(row[col]) and row[col] != '':
                    presence_in_strains.append(col)
            
            if len(presence_in_strains) == 1 and presence_in_strains[0] == strain:
                specific_genes.append({
                    'gene_id': row['Gene'],
                    'annotation': row['Annotation'] if 'Annotation' in row else 'Unknown'
                })
        
        strain_specific[strain] = specific_genes
    
    # 保存菌株特异性基因
    specific_file = os.path.join(output_dir, 'strain_specific_genes.txt')
    with open(specific_file, 'w') as f:
        f.write("=== 菌株特异性基因分析 ===\n\n")
        
        for strain, genes in strain_specific.items():
            f.write(f"=== {strain} 特异性基因 ({len(genes)}个) ===\n")
            
            for gene in genes[:20]:  # 只显示前20个
                f.write(f"- {gene['gene_id']}: {gene['annotation']}\n")
            
            if len(genes) > 20:
                f.write(f"... 还有 {len(genes) - 20} 个基因\n")
            f.write("\n")
    
    print(f"菌株特异性基因分析已保存到: {specific_file}")
    return strain_specific

def analyze_accessory_genome(df, strain_columns, output_dir):
    """深入分析辅助基因组"""
    print("深入分析辅助基因组...")
    
    num_strains = len(strain_columns)
    accessory_genes = []
    
    # 收集辅助基因信息
    for idx, row in df.iterrows():
        presence_count = sum(1 for col in strain_columns if pd.notna(row[col]) and row[col] != '')
        
        if 1 < presence_count < num_strains:
            # 找出存在该基因的菌株
            present_strains = []
            for col in strain_columns:
                if pd.notna(row[col]) and row[col] != '':
                    present_strains.append(col)
            
            accessory_genes.append({
                'gene_id': row['Gene'],
                'annotation': row['Annotation'] if 'Annotation' in row else 'Unknown',
                'presence_count': presence_count,
                'frequency': presence_count / num_strains,
                'present_strains': present_strains
            })
    
    # 按频率分组
    frequency_groups = {
        'rare': [],      # 2-25%
        'intermediate': [], # 25-75%
        'common': []     # 75-99%
    }
    
    for gene in accessory_genes:
        freq = gene['frequency']
        if freq <= 0.25:
            frequency_groups['rare'].append(gene)
        elif freq <= 0.75:
            frequency_groups['intermediate'].append(gene)
        else:
            frequency_groups['common'].append(gene)
    
    # 保存辅助基因组分析
    accessory_file = os.path.join(output_dir, 'accessory_functions.txt')
    with open(accessory_file, 'w') as f:
        f.write("=== 辅助基因组深入分析 ===\n\n")
        
        f.write(f"辅助基因总数: {len(accessory_genes)}\n\n")
        
        for group_name, genes in frequency_groups.items():
            f.write(f"=== {group_name.upper()}辅助基因 ({len(genes)}个) ===\n")
            
            # 显示前10个基因的详细信息
            for gene in genes[:10]:
                f.write(f"基因ID: {gene['gene_id']}\n")
                f.write(f"注释: {gene['annotation']}\n")
                f.write(f"频率: {gene['frequency']:.2f} ({gene['presence_count']}/{num_strains})\n")
                f.write(f"存在菌株: {', '.join(gene['present_strains'])}\n")
                f.write("-" * 50 + "\n")
            
            if len(genes) > 10:
                f.write(f"... 还有 {len(genes) - 10} 个基因\n")
            f.write("\n")
    
    print(f"辅助基因组分析已保存到: {accessory_file}")
    return frequency_groups

def main():
    """主函数"""
    args = parse_arguments()
    
    print(f"功能分类分析脚本 v1.0.0")
    print(f"输入文件: {args.input}")
    print(f"分析模式: {args.mode}")
    print(f"输出目录: {args.output}")
    print("-" * 50)
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    
    try:
        # 加载数据
        df, strain_columns = load_gene_presence_data(args.input)
        
        # 基因分类
        gene_classification = classify_genes_by_presence(df, strain_columns)
        
        if args.mode == 'basic':
            # 基础功能分析
            functional_results = analyze_functional_categories(gene_classification, args.output)
            
            # 菌株特异性基因分析
            strain_specific = analyze_strain_specific_genes(df, strain_columns, args.output)
            
        elif args.mode == 'accessory':
            # 辅助基因组深入分析
            frequency_groups = analyze_accessory_genome(df, strain_columns, args.output)
        
        print("\n✅ 功能分析完成！")
        print(f"结果文件保存在: {args.output}")
        
    except Exception as e:
        print(f"❌ 错误: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()