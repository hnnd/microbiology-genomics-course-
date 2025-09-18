#!/usr/bin/env python3
"""
微生物生活方式与基因组特征关联分析
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def categorize_lifestyle(df):
    """根据基因组特征推断生活方式"""
    df = df.copy()
    
    # 基于基因组大小和已知信息分类
    conditions = [
        (df['Size (Mb)'] < 1.0),  # 极小基因组，可能是共生菌
        (df['Size (Mb)'] >= 1.0) & (df['Size (Mb)'] < 3.0),  # 小基因组
        (df['Size (Mb)'] >= 3.0) & (df['Size (Mb)'] < 6.0),  # 中等基因组
        (df['Size (Mb)'] >= 6.0)  # 大基因组
    ]
    
    choices = ['可能共生菌', '小基因组菌', '中等基因组菌', '大基因组菌']
    
    df['Lifestyle_Category'] = np.select(conditions, choices, default='未分类')
    
    return df

def analyze_lifestyle_patterns(df):
    """分析生活方式模式"""
    print("\n=== 生活方式与基因组特征关联分析 ===")
    
    # 按生活方式分类统计
    lifestyle_stats = df.groupby('Lifestyle_Category').agg({
        'Size (Mb)': ['count', 'mean', 'std', 'min', 'max'],
        'GC%': ['mean', 'std', 'min', 'max'],
        'Genes': ['mean', 'std']
    }).round(2)
    
    print("\n各生活方式类别统计:")
    print(lifestyle_stats)
    
    return lifestyle_stats

def create_lifestyle_plots(df):
    """创建生活方式相关图表"""
    # 设置图表样式
    plt.style.use('default')
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 1. 基因组大小分布
    ax1 = axes[0, 0]
    df['Size (Mb)'].hist(bins=50, alpha=0.7, ax=ax1)
    ax1.axvline(df['Size (Mb)'].mean(), color='red', linestyle='--', 
                label=f'平均值: {df["Size (Mb)"].mean():.2f} Mb')
    ax1.set_xlabel('基因组大小 (Mb)')
    ax1.set_ylabel('频数')
    ax1.set_title('基因组大小分布')
    ax1.legend()
    
    # 2. GC含量分布
    ax2 = axes[0, 1]
    df['GC%'].hist(bins=50, alpha=0.7, ax=ax2, color='orange')
    ax2.axvline(df['GC%'].mean(), color='red', linestyle='--',
                label=f'平均值: {df["GC%"].mean():.1f}%')
    ax2.set_xlabel('GC含量 (%)')
    ax2.set_ylabel('频数')
    ax2.set_title('GC含量分布')
    ax2.legend()
    
    # 3. 基因组大小 vs GC含量散点图
    ax3 = axes[0, 2]
    scatter = ax3.scatter(df['Size (Mb)'], df['GC%'], alpha=0.6, s=20)
    ax3.set_xlabel('基因组大小 (Mb)')
    ax3.set_ylabel('GC含量 (%)')
    ax3.set_title('基因组大小 vs GC含量')
    
    # 添加相关系数
    correlation = df['Size (Mb)'].corr(df['GC%'])
    ax3.text(0.05, 0.95, f'相关系数: {correlation:.3f}', 
             transform=ax3.transAxes, bbox=dict(boxstyle="round", facecolor='wheat'))
    
    # 4. 按生活方式分类的基因组大小箱线图
    ax4 = axes[1, 0]
    df_categorized = categorize_lifestyle(df)
    df_categorized.boxplot(column='Size (Mb)', by='Lifestyle_Category', ax=ax4)
    ax4.set_xlabel('生活方式类别')
    ax4.set_ylabel('基因组大小 (Mb)')
    ax4.set_title('不同生活方式的基因组大小分布')
    plt.setp(ax4.xaxis.get_majorticklabels(), rotation=45)
    
    # 5. 按生活方式分类的GC含量箱线图
    ax5 = axes[1, 1]
    df_categorized.boxplot(column='GC%', by='Lifestyle_Category', ax=ax5)
    ax5.set_xlabel('生活方式类别')
    ax5.set_ylabel('GC含量 (%)')
    ax5.set_title('不同生活方式的GC含量分布')
    plt.setp(ax5.xaxis.get_majorticklabels(), rotation=45)
    
    # 6. 基因数量与基因组大小关系
    ax6 = axes[1, 2]
    if 'Genes' in df.columns:
        ax6.scatter(df['Size (Mb)'], df['Genes'], alpha=0.6, s=20, color='green')
        ax6.set_xlabel('基因组大小 (Mb)')
        ax6.set_ylabel('基因数量')
        ax6.set_title('基因组大小 vs 基因数量')
        
        # 添加趋势线
        z = np.polyfit(df['Size (Mb)'].dropna(), df['Genes'].dropna(), 1)
        p = np.poly1d(z)
        ax6.plot(df['Size (Mb)'], p(df['Size (Mb)']), "r--", alpha=0.8)
        
        # 计算基因密度
        gene_density = df['Genes'] / df['Size (Mb)']
        ax6.text(0.05, 0.95, f'平均基因密度: {gene_density.mean():.0f} 基因/Mb', 
                 transform=ax6.transAxes, bbox=dict(boxstyle="round", facecolor='lightgreen'))
    
    plt.tight_layout()
    plt.savefig('figures/lifestyle_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("生活方式分析图表已保存到 figures/lifestyle_analysis.png")
    
    return df_categorized

def generate_summary_report(df, df_categorized):
    """生成分析总结报告"""
    report = []
    report.append("# 微生物基因组特征与生活方式关联分析报告\n")
    
    report.append("## 数据概况")
    report.append(f"- 分析基因组数量: {len(df):,}")
    report.append(f"- 基因组大小范围: {df['Size (Mb)'].min():.3f} - {df['Size (Mb)'].max():.1f} Mb")
    report.append(f"- GC含量范围: {df['GC%'].min():.1f}% - {df['GC%'].max():.1f}%")
    
    if 'Genes' in df.columns:
        report.append(f"- 基因数量范围: {df['Genes'].min():,.0f} - {df['Genes'].max():,.0f}")
    
    report.append("\n## 主要发现")
    
    # 基因组大小分析
    small_genomes = df[df['Size (Mb)'] < 1.0]
    large_genomes = df[df['Size (Mb)'] > 8.0]
    
    report.append(f"- 极小基因组(<1Mb): {len(small_genomes)} 个 ({len(small_genomes)/len(df)*100:.1f}%)")
    report.append(f"- 大基因组(>8Mb): {len(large_genomes)} 个 ({len(large_genomes)/len(df)*100:.1f}%)")
    
    # GC含量分析
    low_gc = df[df['GC%'] < 40]
    high_gc = df[df['GC%'] > 60]
    
    report.append(f"- 低GC含量(<40%): {len(low_gc)} 个 ({len(low_gc)/len(df)*100:.1f}%)")
    report.append(f"- 高GC含量(>60%): {len(high_gc)} 个 ({len(high_gc)/len(df)*100:.1f}%)")
    
    # 相关性分析
    size_gc_corr = df['Size (Mb)'].corr(df['GC%'])
    report.append(f"- 基因组大小与GC含量相关系数: {size_gc_corr:.3f}")
    
    if 'Genes' in df.columns:
        size_genes_corr = df['Size (Mb)'].corr(df['Genes'])
        report.append(f"- 基因组大小与基因数量相关系数: {size_genes_corr:.3f}")
    
    report.append("\n## 生活方式分类统计")
    lifestyle_counts = df_categorized['Lifestyle_Category'].value_counts()
    for category, count in lifestyle_counts.items():
        percentage = count / len(df_categorized) * 100
        report.append(f"- {category}: {count} 个 ({percentage:.1f}%)")
    
    report.append("\n## 结论")
    report.append("1. 微生物基因组展现出极大的多样性，基因组大小变异范围超过100倍")
    report.append("2. 极小基因组通常与专性共生或寄生生活方式相关")
    report.append("3. 基因组大小与基因数量呈强正相关，反映编码密度相对稳定")
    report.append("4. GC含量变异显著，可能反映不同的环境适应策略")
    
    # 保存报告
    with open('results/analysis_report.md', 'w', encoding='utf-8') as f:
        f.write('\n'.join(report))
    
    print("分析报告已保存到 results/analysis_report.md")

if __name__ == "__main__":
    # 加载清理后的数据
    try:
        df = pd.read_csv('results/cleaned_genome_data.csv')
        print(f"加载数据: {len(df)} 个基因组记录")
        
        # 生活方式分析
        df_categorized = create_lifestyle_plots(df)
        
        # 统计分析
        analyze_lifestyle_patterns(df_categorized)
        
        # 生成报告
        generate_summary_report(df, df_categorized)
        
        # 保存分类结果
        df_categorized.to_csv('results/lifestyle_categorized_data.csv', index=False)
        
    except FileNotFoundError:
        print("请先运行 genome_analysis.py 生成清理后的数据")