#!/usr/bin/env python3
"""
微生物基因组特征分析脚本
用于分析NCBI原核生物基因组数据库信息
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import re
import os

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

def load_prokaryote_data(file_path):
    """加载NCBI原核生物数据"""
    try:
        df = pd.read_csv(file_path, sep='\t', low_memory=False)
        print(f"成功加载数据：{len(df)} 个基因组记录")
        print(f"数据列数：{len(df.columns)}")
        return df
    except Exception as e:
        print(f"数据加载失败：{e}")
        return None

def clean_genome_data(df):
    """清理和预处理基因组数据"""
    # 选择关键列
    key_columns = [
        '#Organism Name', 'Size (Mb)', 'GC%', 'Genes', 'Proteins', 
        'Status', 'Group', 'SubGroup', 'Assembly Accession'
    ]
    
    # 检查列是否存在
    available_columns = [col for col in key_columns if col in df.columns]
    df_clean = df[available_columns].copy()
    
    # 数据类型转换
    if 'Size (Mb)' in df_clean.columns:
        df_clean['Size (Mb)'] = pd.to_numeric(df_clean['Size (Mb)'], errors='coerce')
    if 'GC%' in df_clean.columns:
        df_clean['GC%'] = pd.to_numeric(df_clean['GC%'], errors='coerce')
    if 'Genes' in df_clean.columns:
        df_clean['Genes'] = pd.to_numeric(df_clean['Genes'], errors='coerce')
    
    # 移除缺失值
    df_clean = df_clean.dropna(subset=['Size (Mb)', 'GC%'])
    
    print(f"清理后数据：{len(df_clean)} 个有效记录")
    return df_clean

def basic_statistics(df):
    """计算基本统计信息"""
    print("\n=== 基因组基本统计信息 ===")
    
    if 'Size (Mb)' in df.columns:
        size_stats = df['Size (Mb)'].describe()
        print(f"\n基因组大小统计 (Mb):")
        print(size_stats)
    
    if 'GC%' in df.columns:
        gc_stats = df['GC%'].describe()
        print(f"\nGC含量统计 (%):")
        print(gc_stats)
    
    if 'Genes' in df.columns:
        gene_stats = df['Genes'].describe()
        print(f"\n基因数量统计:")
        print(gene_stats)
    
    return df

def analyze_by_group(df):
    """按分类群分析基因组特征"""
    if 'Group' not in df.columns:
        print("缺少分类群信息")
        return
    
    print("\n=== 按分类群统计 ===")
    group_stats = df.groupby('Group').agg({
        'Size (Mb)': ['count', 'mean', 'std'],
        'GC%': ['mean', 'std'],
        'Genes': ['mean', 'std']
    }).round(2)
    
    print(group_stats)
    return group_stats

if __name__ == "__main__":
    # 加载数据
    df = load_prokaryote_data('data/prokaryotes_summary.txt')
    if df is not None:
        # 清理数据
        df_clean = clean_genome_data(df)
        
        # 基本统计
        basic_statistics(df_clean)
        
        # 分组分析
        analyze_by_group(df_clean)
        
        # 保存清理后的数据
        df_clean.to_csv('results/cleaned_genome_data.csv', index=False)
        print("\n清理后的数据已保存到 results/cleaned_genome_data.csv")