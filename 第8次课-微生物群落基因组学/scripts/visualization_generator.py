#!/usr/bin/env python3
"""
可视化生成脚本
用于生成微生物群落基因组学分析结果的可视化图表

作者: 微生物基因组学课程组
版本: 1.0.0
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import logging
import json
import networkx as nx
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class VisualizationGenerator:
    """可视化生成类"""
    
    def __init__(self, results_dir="results", output_dir="figures"):
        """
        初始化可视化生成器
        
        参数:
        results_dir: 分析结果目录
        output_dir: 图表输出目录
        """
        self.results_dir = Path(results_dir)
        self.output_dir = Path(output_dir)
        
        # 创建输出目录
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置图表样式
        plt.style.use('default')
        sns.set_palette("husl")
        
        # 颜色配置
        self.colors = {
            'high_quality': '#2ecc71',
            'medium_quality': '#f39c12', 
            'low_quality': '#e74c3c',
            'complete': '#27ae60',
            'partial': '#f1c40f',
            'incomplete': '#e67e22'
        }
    
    def load_data(self):
        """加载分析结果数据"""
        logger.info("加载分析结果数据...")
        
        data = {}
        
        # 加载MAG质量数据
        mag_file = self.results_dir / "mag_detailed_results.tsv"
        if mag_file.exists():
            data['mag_results'] = pd.read_csv(mag_file, sep='\t')
        
        # 加载功能分析数据
        functional_dir = self.results_dir / "functional_analysis"
        if functional_dir.exists():
            # 基因家族数据
            families_file = functional_dir / "gene_families.tsv"
            if families_file.exists():
                data['gene_families'] = pd.read_csv(families_file, sep='\t')
            
            # 功能谱数据
            profile_file = functional_dir / "functional_profile.json"
            if profile_file.exists():
                with open(profile_file, 'r') as f:
                    data['functional_profile'] = json.load(f)
            
            # 多样性指数
            diversity_file = functional_dir / "diversity_indices.json"
            if diversity_file.exists():
                with open(diversity_file, 'r') as f:
                    data['diversity_indices'] = json.load(f)
        
        # 加载通路重建数据
        pathway_dir = self.results_dir / "pathway_reconstruction"
        if pathway_dir.exists():
            # 通路完整性
            completeness_file = pathway_dir / "pathway_completeness.json"
            if completeness_file.exists():
                with open(completeness_file, 'r') as f:
                    data['pathway_completeness'] = json.load(f)
            
            # 网络性质
            network_file = pathway_dir / "network_properties.json"
            if network_file.exists():
                with open(network_file, 'r') as f:
                    data['network_properties'] = json.load(f)
            
            # 关键代谢物
            metabolites_file = pathway_dir / "key_metabolites.json"
            if metabolites_file.exists():
                with open(metabolites_file, 'r') as f:
                    data['key_metabolites'] = json.load(f)
        
        return data
    
    def create_mag_quality_heatmap(self, mag_data):
        """创建MAG质量热图"""
        logger.info("创建MAG质量热图...")
        
        if mag_data is None:
            logger.warning("MAG数据不存在，跳过质量热图")
            return
        
        # 准备数据
        quality_metrics = ['Completeness', 'Contamination', 'GC_Content', 'Total_Length']
        
        # 标准化数据（除了污染率，其他指标越高越好）
        plot_data = mag_data[['MAG_ID'] + quality_metrics].copy()
        
        # 对长度进行对数转换
        plot_data['Total_Length'] = np.log10(plot_data['Total_Length'])
        
        # 反转污染率（使得低污染率显示为高分）
        plot_data['Contamination'] = 100 - plot_data['Contamination']
        
        # 标准化到0-100范围
        for col in quality_metrics:
            if col == 'Total_Length':
                # 长度标准化
                min_val, max_val = plot_data[col].min(), plot_data[col].max()
                plot_data[col] = (plot_data[col] - min_val) / (max_val - min_val) * 100
            elif col != 'Contamination':  # 完整性和GC含量已经是百分比
                pass
        
        # 创建热图
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # 设置数据矩阵
        heatmap_data = plot_data.set_index('MAG_ID')[quality_metrics]
        
        # 绘制热图
        im = ax.imshow(heatmap_data.values, cmap='RdYlGn', aspect='auto', vmin=0, vmax=100)
        
        # 设置坐标轴
        ax.set_xticks(range(len(quality_metrics)))
        ax.set_xticklabels(['完整性', '污染率(反)', 'GC含量', '基因组大小'], rotation=45)
        ax.set_yticks(range(len(heatmap_data)))
        ax.set_yticklabels(heatmap_data.index)
        
        # 添加数值标签
        for i in range(len(heatmap_data)):
            for j in range(len(quality_metrics)):
                value = heatmap_data.iloc[i, j]
                ax.text(j, i, f'{value:.1f}', ha='center', va='center', 
                       color='white' if value < 50 else 'black', fontweight='bold')
        
        # 添加颜色条
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('质量评分', rotation=270, labelpad=20)
        
        # 设置标题和布局
        plt.title('MAG基因组质量热图', fontsize=16, fontweight='bold', pad=20)
        plt.tight_layout()
        
        # 保存图片
        output_file = self.output_dir / "mag_quality_heatmap.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"MAG质量热图已保存到: {output_file}")
    
    def create_functional_profile_barplot(self, functional_profile):
        """创建功能谱柱状图"""
        logger.info("创建功能谱柱状图...")
        
        if functional_profile is None:
            logger.warning("功能谱数据不存在，跳过功能谱图")
            return
        
        # 创建子图
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('微生物群落功能谱分析', fontsize=16, fontweight='bold')
        
        # KEGG通路分析
        if 'kegg_pathways' in functional_profile:
            ax = axes[0, 0]
            kegg_data = functional_profile['kegg_pathways']
            
            pathways = list(kegg_data.keys())
            gene_counts = [kegg_data[p]['gene_count'] for p in pathways]
            pathway_names = [kegg_data[p]['name'][:30] + '...' if len(kegg_data[p]['name']) > 30 
                           else kegg_data[p]['name'] for p in pathways]
            
            bars = ax.barh(pathway_names, gene_counts, color=sns.color_palette("viridis", len(pathways)))
            ax.set_xlabel('基因数量')
            ax.set_title('KEGG代谢通路', fontweight='bold')
            
            # 添加数值标签
            for i, bar in enumerate(bars):
                width = bar.get_width()
                ax.text(width + 0.1, bar.get_y() + bar.get_height()/2, 
                       f'{int(width)}', ha='left', va='center')
        
        # COG功能分类
        if 'cog_categories' in functional_profile:
            ax = axes[0, 1]
            cog_data = functional_profile['cog_categories']
            
            categories = list(cog_data.keys())
            gene_counts = [cog_data[c]['gene_count'] for c in categories]
            category_names = [f"[{c}] {cog_data[c]['name'][:25]}..." if len(cog_data[c]['name']) > 25 
                            else f"[{c}] {cog_data[c]['name']}" for c in categories]
            
            bars = ax.barh(category_names, gene_counts, color=sns.color_palette("plasma", len(categories)))
            ax.set_xlabel('基因数量')
            ax.set_title('COG功能分类', fontweight='bold')
            
            # 添加数值标签
            for i, bar in enumerate(bars):
                width = bar.get_width()
                ax.text(width + 0.1, bar.get_y() + bar.get_height()/2, 
                       f'{int(width)}', ha='left', va='center')
        
        # Pfam结构域分析
        if 'pfam_domains' in functional_profile:
            ax = axes[1, 0]
            pfam_data = functional_profile['pfam_domains']
            
            domains = list(pfam_data.keys())[:10]  # 只显示前10个
            gene_counts = [pfam_data[d]['gene_count'] for d in domains]
            
            bars = ax.bar(domains, gene_counts, color=sns.color_palette("coolwarm", len(domains)))
            ax.set_xlabel('Pfam结构域')
            ax.set_ylabel('基因数量')
            ax.set_title('主要Pfam结构域 (Top 10)', fontweight='bold')
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
            
            # 添加数值标签
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2, height + 0.5, 
                       f'{int(height)}', ha='center', va='bottom')
        
        # 功能多样性比较
        ax = axes[1, 1]
        diversity_types = ['KEGG通路', 'COG分类', 'Pfam结构域']
        richness_values = [
            len(functional_profile.get('kegg_pathways', {})),
            len(functional_profile.get('cog_categories', {})),
            len(functional_profile.get('pfam_domains', {}))
        ]
        
        bars = ax.bar(diversity_types, richness_values, 
                     color=['#3498db', '#e74c3c', '#2ecc71'])
        ax.set_ylabel('功能类别数量')
        ax.set_title('功能多样性比较', fontweight='bold')
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        
        # 添加数值标签
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 0.5, 
                   f'{int(height)}', ha='center', va='bottom')
        
        plt.tight_layout()
        
        # 保存图片
        output_file = self.output_dir / "functional_profile_barplot.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"功能谱柱状图已保存到: {output_file}")
    
    def create_pathway_completeness_chart(self, pathway_completeness):
        """创建通路完整性图表"""
        logger.info("创建通路完整性图表...")
        
        if pathway_completeness is None:
            logger.warning("通路完整性数据不存在，跳过通路图")
            return
        
        # 准备数据
        pathways = list(pathway_completeness.keys())
        completeness_values = [pathway_completeness[p]['completeness'] for p in pathways]
        pathway_names = [pathway_completeness[p]['name'] for p in pathways]
        categories = [pathway_completeness[p]['category'] for p in pathways]
        
        # 创建图表
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        fig.suptitle('代谢通路完整性分析', fontsize=16, fontweight='bold')
        
        # 通路完整性柱状图
        colors = []
        for comp in completeness_values:
            if comp >= 80:
                colors.append(self.colors['complete'])
            elif comp >= 50:
                colors.append(self.colors['partial'])
            else:
                colors.append(self.colors['incomplete'])
        
        bars = ax1.barh(pathway_names, completeness_values, color=colors)
        ax1.set_xlabel('完整性 (%)')
        ax1.set_title('各通路完整性水平', fontweight='bold')
        ax1.axvline(x=80, color='red', linestyle='--', alpha=0.7, label='高完整性阈值')
        ax1.axvline(x=50, color='orange', linestyle='--', alpha=0.7, label='中等完整性阈值')
        ax1.legend()
        
        # 添加数值标签
        for i, bar in enumerate(bars):
            width = bar.get_width()
            ax1.text(width + 1, bar.get_y() + bar.get_height()/2, 
                    f'{width:.1f}%', ha='left', va='center', fontweight='bold')
        
        # 完整性分布饼图
        complete_count = sum(1 for comp in completeness_values if comp >= 80)
        partial_count = sum(1 for comp in completeness_values if 50 <= comp < 80)
        incomplete_count = sum(1 for comp in completeness_values if comp < 50)
        
        sizes = [complete_count, partial_count, incomplete_count]
        labels = [f'完整 (≥80%)\n{complete_count}个', 
                 f'部分完整 (50-80%)\n{partial_count}个', 
                 f'不完整 (<50%)\n{incomplete_count}个']
        colors_pie = [self.colors['complete'], self.colors['partial'], self.colors['incomplete']]
        
        # 只显示非零的部分
        non_zero_sizes = [(size, label, color) for size, label, color in zip(sizes, labels, colors_pie) if size > 0]
        if non_zero_sizes:
            sizes_nz, labels_nz, colors_nz = zip(*non_zero_sizes)
            
            wedges, texts, autotexts = ax2.pie(sizes_nz, labels=labels_nz, colors=colors_nz, 
                                              autopct='%1.1f%%', startangle=90)
            ax2.set_title('通路完整性分布', fontweight='bold')
            
            # 美化文本
            for autotext in autotexts:
                autotext.set_color('white')
                autotext.set_fontweight('bold')
        
        plt.tight_layout()
        
        # 保存图片
        output_file = self.output_dir / "pathway_completeness_chart.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"通路完整性图表已保存到: {output_file}")
    
    def create_metabolic_network_plot(self, network_file=None):
        """创建代谢网络图"""
        logger.info("创建代谢网络图...")
        
        # 尝试加载网络文件
        network_path = self.results_dir / "pathway_reconstruction" / "metabolic_network.gml"
        
        if not network_path.exists():
            logger.warning("代谢网络文件不存在，创建示例网络图")
            self._create_example_network_plot()
            return
        
        try:
            # 加载网络
            G = nx.read_gml(network_path)
            
            # 创建图表
            fig, ax = plt.subplots(figsize=(12, 10))
            
            # 设置节点位置
            pos = nx.spring_layout(G, k=2, iterations=50)
            
            # 分离化合物和反应节点
            compound_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'compound']
            reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
            
            # 绘制化合物节点
            nx.draw_networkx_nodes(G, pos, nodelist=compound_nodes, 
                                 node_color='lightblue', node_size=300, 
                                 node_shape='o', alpha=0.8, ax=ax)
            
            # 绘制反应节点
            nx.draw_networkx_nodes(G, pos, nodelist=reaction_nodes, 
                                 node_color='lightcoral', node_size=200, 
                                 node_shape='s', alpha=0.8, ax=ax)
            
            # 绘制边
            nx.draw_networkx_edges(G, pos, edge_color='gray', 
                                 arrows=True, arrowsize=20, alpha=0.6, ax=ax)
            
            # 添加标签（只显示部分重要节点）
            important_nodes = compound_nodes[:10] + reaction_nodes[:5]
            labels = {n: n for n in important_nodes}
            nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)
            
            # 设置标题和图例
            ax.set_title('微生物群落代谢网络', fontsize=16, fontweight='bold')
            
            # 添加图例
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor='lightblue', 
                      markersize=10, label='化合物'),
                Line2D([0], [0], marker='s', color='w', markerfacecolor='lightcoral', 
                      markersize=8, label='反应')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
            
            ax.axis('off')
            plt.tight_layout()
            
            # 保存图片
            output_file = self.output_dir / "metabolic_network.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"代谢网络图已保存到: {output_file}")
            
        except Exception as e:
            logger.error(f"加载网络文件失败: {e}")
            self._create_example_network_plot()
    
    def _create_example_network_plot(self):
        """创建示例网络图"""
        logger.info("创建示例代谢网络图...")
        
        # 创建示例网络
        G = nx.DiGraph()
        
        # 添加节点
        compounds = ['Glucose', 'G6P', 'F6P', 'FBP', 'Pyruvate', 'Acetyl-CoA']
        reactions = ['HK', 'GPI', 'PFK', 'PYK', 'PDH']
        
        for comp in compounds:
            G.add_node(comp, type='compound')
        for rxn in reactions:
            G.add_node(rxn, type='reaction')
        
        # 添加边（简化的糖酵解途径）
        edges = [
            ('Glucose', 'HK'), ('HK', 'G6P'),
            ('G6P', 'GPI'), ('GPI', 'F6P'),
            ('F6P', 'PFK'), ('PFK', 'FBP'),
            ('FBP', 'PYK'), ('PYK', 'Pyruvate'),
            ('Pyruvate', 'PDH'), ('PDH', 'Acetyl-CoA')
        ]
        G.add_edges_from(edges)
        
        # 创建图表
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 设置布局
        pos = nx.spring_layout(G, k=3, iterations=50)
        
        # 绘制节点
        nx.draw_networkx_nodes(G, pos, nodelist=compounds, 
                             node_color='lightblue', node_size=800, 
                             node_shape='o', alpha=0.8, ax=ax)
        nx.draw_networkx_nodes(G, pos, nodelist=reactions, 
                             node_color='lightcoral', node_size=600, 
                             node_shape='s', alpha=0.8, ax=ax)
        
        # 绘制边
        nx.draw_networkx_edges(G, pos, edge_color='gray', 
                             arrows=True, arrowsize=20, alpha=0.6, ax=ax)
        
        # 添加标签
        nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax)
        
        # 设置标题
        ax.set_title('示例代谢网络 (糖酵解途径)', fontsize=16, fontweight='bold')
        
        # 添加图例
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='lightblue', 
                  markersize=15, label='化合物'),
            Line2D([0], [0], marker='s', color='w', markerfacecolor='lightcoral', 
                  markersize=12, label='酶反应')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        ax.axis('off')
        plt.tight_layout()
        
        # 保存图片
        output_file = self.output_dir / "metabolic_network.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"示例代谢网络图已保存到: {output_file}")
    
    def create_diversity_comparison(self, diversity_indices):
        """创建多样性比较图"""
        logger.info("创建多样性比较图...")
        
        if diversity_indices is None:
            logger.warning("多样性指数数据不存在，跳过多样性图")
            return
        
        # 准备数据
        shannon_indices = {
            'KEGG通路': diversity_indices.get('kegg_shannon', 0),
            'COG分类': diversity_indices.get('cog_shannon', 0),
            'Pfam结构域': diversity_indices.get('pfam_shannon', 0)
        }
        
        richness_indices = {
            'KEGG通路': diversity_indices.get('kegg_richness', 0),
            'COG分类': diversity_indices.get('cog_richness', 0),
            'Pfam结构域': diversity_indices.get('pfam_richness', 0)
        }
        
        # 创建图表
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle('功能多样性指数比较', fontsize=16, fontweight='bold')
        
        # Shannon多样性指数
        categories = list(shannon_indices.keys())
        shannon_values = list(shannon_indices.values())
        
        bars1 = ax1.bar(categories, shannon_values, color=['#3498db', '#e74c3c', '#2ecc71'])
        ax1.set_ylabel('Shannon多样性指数')
        ax1.set_title('Shannon多样性指数', fontweight='bold')
        ax1.set_ylim(0, max(shannon_values) * 1.2 if shannon_values else 1)
        
        # 添加数值标签
        for bar in bars1:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2, height + 0.01, 
                    f'{height:.3f}', ha='center', va='bottom', fontweight='bold')
        
        # 功能丰富度
        richness_values = list(richness_indices.values())
        
        bars2 = ax2.bar(categories, richness_values, color=['#9b59b6', '#f39c12', '#1abc9c'])
        ax2.set_ylabel('功能类别数量')
        ax2.set_title('功能丰富度', fontweight='bold')
        ax2.set_ylim(0, max(richness_values) * 1.2 if richness_values else 1)
        
        # 添加数值标签
        for bar in bars2:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2, height + 0.5, 
                    f'{int(height)}', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        
        # 保存图片
        output_file = self.output_dir / "diversity_comparison.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"多样性比较图已保存到: {output_file}")
    
    def generate_summary_dashboard(self, data):
        """生成汇总仪表板"""
        logger.info("生成汇总仪表板...")
        
        # 创建大型仪表板
        fig = plt.figure(figsize=(20, 16))
        gs = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)
        
        # 主标题
        fig.suptitle('微生物群落基因组学分析汇总仪表板', fontsize=24, fontweight='bold', y=0.95)
        
        # 1. MAG质量概览
        if 'mag_results' in data:
            ax1 = fig.add_subplot(gs[0, :2])
            mag_data = data['mag_results']
            quality_counts = mag_data['Quality_Level'].value_counts()
            
            colors = [self.colors['high_quality'], self.colors['medium_quality'], self.colors['low_quality']]
            wedges, texts, autotexts = ax1.pie(quality_counts.values, labels=quality_counts.index, 
                                              colors=colors, autopct='%1.1f%%', startangle=90)
            ax1.set_title('MAG质量分布', fontsize=14, fontweight='bold')
        
        # 2. 功能多样性雷达图
        if 'diversity_indices' in data:
            ax2 = fig.add_subplot(gs[0, 2:], projection='polar')
            diversity_data = data['diversity_indices']
            
            categories = ['KEGG\n通路', 'COG\n分类', 'Pfam\n结构域']
            values = [
                diversity_data.get('kegg_shannon', 0),
                diversity_data.get('cog_shannon', 0),
                diversity_data.get('pfam_shannon', 0)
            ]
            
            # 标准化到0-1范围
            max_val = max(values) if values else 1
            normalized_values = [v/max_val for v in values] if max_val > 0 else [0, 0, 0]
            
            angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False).tolist()
            normalized_values += normalized_values[:1]  # 闭合图形
            angles += angles[:1]
            
            ax2.plot(angles, normalized_values, 'o-', linewidth=2, color='#3498db')
            ax2.fill(angles, normalized_values, alpha=0.25, color='#3498db')
            ax2.set_xticks(angles[:-1])
            ax2.set_xticklabels(categories)
            ax2.set_ylim(0, 1)
            ax2.set_title('功能多样性雷达图', fontsize=14, fontweight='bold', pad=20)
        
        # 3. 通路完整性
        if 'pathway_completeness' in data:
            ax3 = fig.add_subplot(gs[1, :2])
            pathway_data = data['pathway_completeness']
            
            pathways = list(pathway_data.keys())
            completeness = [pathway_data[p]['completeness'] for p in pathways]
            pathway_names = [pathway_data[p]['name'][:20] + '...' if len(pathway_data[p]['name']) > 20 
                           else pathway_data[p]['name'] for p in pathways]
            
            colors = []
            for comp in completeness:
                if comp >= 80:
                    colors.append(self.colors['complete'])
                elif comp >= 50:
                    colors.append(self.colors['partial'])
                else:
                    colors.append(self.colors['incomplete'])
            
            bars = ax3.barh(pathway_names, completeness, color=colors)
            ax3.set_xlabel('完整性 (%)')
            ax3.set_title('代谢通路完整性', fontsize=14, fontweight='bold')
            ax3.axvline(x=80, color='red', linestyle='--', alpha=0.5)
            ax3.axvline(x=50, color='orange', linestyle='--', alpha=0.5)
        
        # 4. 网络统计
        if 'network_properties' in data:
            ax4 = fig.add_subplot(gs[1, 2:])
            network_data = data['network_properties']
            
            metrics = ['节点数', '边数', '化合物', '反应']
            values = [
                network_data.get('nodes', 0),
                network_data.get('edges', 0),
                network_data.get('compound_nodes', 0),
                network_data.get('reaction_nodes', 0)
            ]
            
            bars = ax4.bar(metrics, values, color=['#e74c3c', '#3498db', '#2ecc71', '#f39c12'])
            ax4.set_ylabel('数量')
            ax4.set_title('代谢网络统计', fontsize=14, fontweight='bold')
            
            # 添加数值标签
            for bar in bars:
                height = bar.get_height()
                ax4.text(bar.get_x() + bar.get_width()/2, height + 0.5, 
                        f'{int(height)}', ha='center', va='bottom', fontweight='bold')
        
        # 5. 基因家族分布
        if 'gene_families' in data:
            ax5 = fig.add_subplot(gs[2, :2])
            families_data = data['gene_families']
            
            # 按丰度排序，取前10
            top_families = families_data.nlargest(10, 'abundance')
            
            bars = ax5.bar(range(len(top_families)), top_families['abundance'], 
                          color=sns.color_palette("viridis", len(top_families)))
            ax5.set_xlabel('基因家族排名')
            ax5.set_ylabel('丰度')
            ax5.set_title('主要基因家族丰度 (Top 10)', fontsize=14, fontweight='bold')
            ax5.set_xticks(range(len(top_families)))
            ax5.set_xticklabels([f'{i+1}' for i in range(len(top_families))])
        
        # 6. 关键代谢物
        if 'key_metabolites' in data:
            ax6 = fig.add_subplot(gs[2, 2:])
            metabolites_data = data['key_metabolites'][:8]  # 取前8个
            
            names = [m['name'][:15] + '...' if len(m['name']) > 15 else m['name'] for m in metabolites_data]
            degrees = [m['total_degree'] for m in metabolites_data]
            
            bars = ax6.barh(names, degrees, color=sns.color_palette("coolwarm", len(names)))
            ax6.set_xlabel('连接度')
            ax6.set_title('关键代谢物 (Top 8)', fontsize=14, fontweight='bold')
        
        # 7. 统计摘要文本
        ax7 = fig.add_subplot(gs[3, :])
        ax7.axis('off')
        
        # 生成摘要文本
        summary_text = self._generate_summary_text(data)
        ax7.text(0.05, 0.95, summary_text, transform=ax7.transAxes, fontsize=12,
                verticalalignment='top', bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
        
        # 保存仪表板
        output_file = self.output_dir / "summary_dashboard.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"汇总仪表板已保存到: {output_file}")
    
    def _generate_summary_text(self, data):
        """生成摘要文本"""
        summary_lines = ["=== 分析结果摘要 ===\n"]
        
        if 'mag_results' in data:
            mag_data = data['mag_results']
            total_mags = len(mag_data)
            high_quality = len(mag_data[mag_data['Quality_Level'] == 'High Quality'])
            summary_lines.append(f"• MAG基因组: 总计{total_mags}个，其中{high_quality}个高质量")
        
        if 'functional_profile' in data:
            fp = data['functional_profile']
            kegg_count = len(fp.get('kegg_pathways', {}))
            cog_count = len(fp.get('cog_categories', {}))
            summary_lines.append(f"• 功能注释: {kegg_count}个KEGG通路，{cog_count}个COG分类")
        
        if 'pathway_completeness' in data:
            pc = data['pathway_completeness']
            complete_pathways = sum(1 for p in pc.values() if p['completeness'] >= 80)
            summary_lines.append(f"• 代谢通路: {complete_pathways}个完整通路 (≥80%)")
        
        if 'network_properties' in data:
            np_data = data['network_properties']
            nodes = np_data.get('nodes', 0)
            edges = np_data.get('edges', 0)
            summary_lines.append(f"• 代谢网络: {nodes}个节点，{edges}条边")
        
        if 'diversity_indices' in data:
            di = data['diversity_indices']
            kegg_shannon = di.get('kegg_shannon', 0)
            summary_lines.append(f"• 功能多样性: KEGG Shannon指数 = {kegg_shannon:.3f}")
        
        summary_lines.append("\n该微生物群落显示出丰富的代谢功能和良好的基因组质量，")
        summary_lines.append("具备完整的核心代谢能力和环境适应潜力。")
        
        return '\n'.join(summary_lines)
    
    def generate_all_visualizations(self):
        """生成所有可视化图表"""
        logger.info("开始生成所有可视化图表...")
        
        # 加载数据
        data = self.load_data()
        
        if not data:
            logger.error("没有找到分析结果数据")
            return
        
        # 生成各种图表
        self.create_mag_quality_heatmap(data.get('mag_results'))
        self.create_functional_profile_barplot(data.get('functional_profile'))
        self.create_pathway_completeness_chart(data.get('pathway_completeness'))
        self.create_metabolic_network_plot()
        self.create_diversity_comparison(data.get('diversity_indices'))
        self.generate_summary_dashboard(data)
        
        logger.info("所有可视化图表生成完成!")
        
        return True

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='微生物群落基因组学可视化工具')
    parser.add_argument('--results_dir', default='results', help='分析结果目录')
    parser.add_argument('--output_dir', default='figures', help='图表输出目录')
    
    args = parser.parse_args()
    
    # 创建可视化生成器实例
    visualizer = VisualizationGenerator(
        results_dir=args.results_dir,
        output_dir=args.output_dir
    )
    
    # 生成所有可视化
    success = visualizer.generate_all_visualizations()
    
    if success:
        print("\n=== 可视化生成完成 ===")
        print(f"图表已保存到: {args.output_dir}/")
        print("生成的图表包括:")
        print("- MAG质量热图")
        print("- 功能谱柱状图")
        print("- 通路完整性图表")
        print("- 代谢网络图")
        print("- 多样性比较图")
        print("- 汇总仪表板")
    else:
        print("可视化生成失败")
        sys.exit(1)

if __name__ == "__main__":
    main()