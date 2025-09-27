#!/usr/bin/env python3
"""
代谢网络可视化器
用于生成代谢网络的各种可视化图表

作者: 微生物基因组学课程组
版本: 1.0.0
"""

import argparse
import pandas as pd
import json
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

class NetworkVisualizer:
    def __init__(self, output_dir="results", figure_dir="figures"):
        self.output_dir = Path(output_dir)
        self.figure_dir = Path(figure_dir)
        self.figure_dir.mkdir(exist_ok=True)
        
        # 设置matplotlib中文字体
        plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False
        
        # 加载网络数据
        self.metabolic_network = None
        self.nutrition_data = None
        self._load_data()
        
        # 颜色配置
        self.colors = {
            'reaction': '#FF6B6B',
            'metabolite': '#4ECDC4',
            'glycolysis': '#FFE66D',
            'TCA_cycle': '#95E1D3',
            'amino_acid': '#A8E6CF',
            'nucleotide': '#C7CEEA'
        }
    
    def _load_data(self):
        """加载网络和营养数据"""
        # 加载代谢网络
        network_file = self.output_dir / "metabolic_network.json"
        if network_file.exists():
            with open(network_file, 'r') as f:
                network_data = json.load(f)
            self.metabolic_network = nx.node_link_graph(network_data)
            print(f"已加载代谢网络: {self.metabolic_network.number_of_nodes()} 个节点")
        
        # 加载营养预测数据
        nutrition_file = self.output_dir / "nutrition_prediction_detailed.json"
        if nutrition_file.exists():
            with open(nutrition_file, 'r') as f:
                self.nutrition_data = json.load(f)
            print("已加载营养预测数据")
    
    def visualize_network(self):
        """生成完整代谢网络可视化"""
        if self.metabolic_network is None:
            print("代谢网络数据不存在，请先进行网络重建")
            return
        
        print("正在生成代谢网络可视化...")
        
        # 创建图形
        plt.figure(figsize=(16, 12))
        
        # 设置布局
        pos = nx.spring_layout(self.metabolic_network, k=3, iterations=50)
        
        # 分离反应和代谢物节点
        reaction_nodes = [n for n, d in self.metabolic_network.nodes(data=True) 
                         if d.get('type') == 'reaction']
        metabolite_nodes = [n for n, d in self.metabolic_network.nodes(data=True) 
                           if d.get('type') == 'metabolite']
        
        # 绘制代谢物节点
        nx.draw_networkx_nodes(self.metabolic_network, pos,
                              nodelist=metabolite_nodes,
                              node_color=self.colors['metabolite'],
                              node_size=300,
                              alpha=0.8)
        
        # 绘制反应节点
        nx.draw_networkx_nodes(self.metabolic_network, pos,
                              nodelist=reaction_nodes,
                              node_color=self.colors['reaction'],
                              node_size=200,
                              node_shape='s',
                              alpha=0.8)
        
        # 绘制边
        nx.draw_networkx_edges(self.metabolic_network, pos,
                              edge_color='gray',
                              alpha=0.6,
                              arrows=True,
                              arrowsize=10)
        
        # 添加标签（仅对重要节点）
        important_nodes = [n for n in self.metabolic_network.nodes() 
                          if self.metabolic_network.degree(n) > 3]
        labels = {n: n.replace('-', '\\n') for n in important_nodes}
        nx.draw_networkx_labels(self.metabolic_network, pos,
                               labels=labels,
                               font_size=8)
        
        plt.title('微生物代谢网络', fontsize=16, fontweight='bold')
        plt.axis('off')
        
        # 添加图例
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor=self.colors['metabolite'], 
                      markersize=10, label='代谢物'),
            plt.Line2D([0], [0], marker='s', color='w', 
                      markerfacecolor=self.colors['reaction'], 
                      markersize=10, label='反应')
        ]
        plt.legend(handles=legend_elements, loc='upper right')
        
        # 保存图片
        output_file = self.figure_dir / "metabolic_network.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"代谢网络图已保存到: {output_file}")
    
    def pathway_heatmap(self):
        """生成通路完整性热图"""
        if self.nutrition_data is None:
            print("营养预测数据不存在，请先运行营养需求预测")
            return
        
        print("正在生成通路完整性热图...")
        
        # 准备数据
        pathway_data = []
        
        # 氨基酸合成数据
        if 'amino_acid_synthesis' in self.nutrition_data:
            for aa, info in self.nutrition_data['amino_acid_synthesis'].items():
                pathway_data.append({
                    'pathway': f'AA_{aa}',
                    'category': '氨基酸合成',
                    'completeness': info['completeness']
                })
        
        # 碳源利用数据
        if 'carbon_sources' in self.nutrition_data:
            for carbon, info in self.nutrition_data['carbon_sources'].items():
                pathway_data.append({
                    'pathway': f'Carbon_{carbon}',
                    'category': '碳源利用',
                    'completeness': info['completeness']
                })
        
        # 维生素合成数据
        if 'vitamin_requirements' in self.nutrition_data:
            for vitamin, info in self.nutrition_data['vitamin_requirements'].items():
                pathway_data.append({
                    'pathway': f'Vitamin_{vitamin}',
                    'category': '维生素合成',
                    'completeness': info['completeness']
                })
        
        if not pathway_data:
            print("没有可用的通路数据")
            return
        
        # 转换为DataFrame
        df = pd.DataFrame(pathway_data)
        
        # 创建透视表
        pivot_df = df.pivot(index='pathway', columns='category', values='completeness')
        pivot_df = pivot_df.fillna(0)
        
        # 创建热图
        plt.figure(figsize=(12, 10))
        sns.heatmap(pivot_df, 
                   annot=True, 
                   fmt='.1f',
                   cmap='RdYlGn',
                   vmin=0, 
                   vmax=100,
                   cbar_kws={'label': '完整性 (%)'})
        
        plt.title('代谢通路完整性热图', fontsize=16, fontweight='bold')
        plt.xlabel('通路类别', fontsize=12)
        plt.ylabel('具体通路', fontsize=12)
        plt.xticks(rotation=45)
        plt.yticks(rotation=0)
        
        # 保存图片
        output_file = self.figure_dir / "pathway_heatmap.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"通路热图已保存到: {output_file}")
    
    def nutrition_summary(self):
        """生成营养需求总结图"""
        if self.nutrition_data is None:
            print("营养预测数据不存在")
            return
        
        print("正在生成营养需求总结图...")
        
        # 创建子图
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('微生物营养需求分析总结', fontsize=16, fontweight='bold')
        
        # 1. 氨基酸合成能力饼图
        if 'amino_acid_synthesis' in self.nutrition_data:
            aa_data = self.nutrition_data['amino_acid_synthesis']
            synthesis_status = defaultdict(int)
            
            for aa, info in aa_data.items():
                if info['completeness'] >= 90:
                    synthesis_status['完全合成'] += 1
                elif info['completeness'] >= 70:
                    synthesis_status['可能合成'] += 1
                elif info['completeness'] >= 30:
                    synthesis_status['部分合成'] += 1
                else:
                    synthesis_status['需要外源'] += 1
            
            axes[0, 0].pie(synthesis_status.values(), 
                          labels=synthesis_status.keys(),
                          autopct='%1.1f%%',
                          colors=['#2ECC71', '#F39C12', '#E74C3C', '#95A5A6'])
            axes[0, 0].set_title('氨基酸合成能力分布')
        
        # 2. 碳源利用能力条形图
        if 'carbon_sources' in self.nutrition_data:
            carbon_data = self.nutrition_data['carbon_sources']
            carbon_names = list(carbon_data.keys())
            carbon_completeness = [info['completeness'] for info in carbon_data.values()]
            
            bars = axes[0, 1].bar(range(len(carbon_names)), carbon_completeness,
                                 color=['#3498DB' if x >= 70 else '#E74C3C' for x in carbon_completeness])
            axes[0, 1].set_xlabel('碳源类型')
            axes[0, 1].set_ylabel('利用能力 (%)')
            axes[0, 1].set_title('碳源利用能力')
            axes[0, 1].set_xticks(range(len(carbon_names)))
            axes[0, 1].set_xticklabels(carbon_names, rotation=45)
            axes[0, 1].axhline(y=70, color='red', linestyle='--', alpha=0.7)
        
        # 3. 维生素需求状态
        if 'vitamin_requirements' in self.nutrition_data:
            vitamin_data = self.nutrition_data['vitamin_requirements']
            vitamin_status = defaultdict(int)
            
            for vitamin, info in vitamin_data.items():
                if info['completeness'] >= 90:
                    vitamin_status['可自主合成'] += 1
                elif info['completeness'] >= 30:
                    vitamin_status['部分合成'] += 1
                else:
                    vitamin_status['需要外源'] += 1
            
            axes[1, 0].bar(vitamin_status.keys(), vitamin_status.values(),
                          color=['#27AE60', '#F39C12', '#E74C3C'])
            axes[1, 0].set_xlabel('合成状态')
            axes[1, 0].set_ylabel('维生素数量')
            axes[1, 0].set_title('维生素合成能力')
            axes[1, 0].tick_params(axis='x', rotation=45)
        
        # 4. 氨基酸家族分析
        if 'amino_acid_synthesis' in self.nutrition_data:
            aa_data = self.nutrition_data['amino_acid_synthesis']
            family_completeness = defaultdict(list)
            
            for aa, info in aa_data.items():
                family_completeness[info['family']].append(info['completeness'])
            
            family_avg = {family: np.mean(completeness) 
                         for family, completeness in family_completeness.items()}
            
            families = list(family_avg.keys())
            avg_completeness = list(family_avg.values())
            
            bars = axes[1, 1].bar(range(len(families)), avg_completeness,
                                 color=['#9B59B6', '#E67E22', '#1ABC9C', '#34495E', '#F1C40F'])
            axes[1, 1].set_xlabel('氨基酸家族')
            axes[1, 1].set_ylabel('平均完整性 (%)')
            axes[1, 1].set_title('氨基酸家族合成能力')
            axes[1, 1].set_xticks(range(len(families)))
            axes[1, 1].set_xticklabels([f.replace('_', '\\n') for f in families], rotation=45)
        
        plt.tight_layout()
        
        # 保存图片
        output_file = self.figure_dir / "nutrition_summary.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"营养需求总结图已保存到: {output_file}")
    
    def pathway_specific(self, pathway='glycolysis'):
        """生成特定通路的详细图"""
        if self.metabolic_network is None:
            print("代谢网络数据不存在")
            return
        
        print(f"正在生成 {pathway} 通路详细图...")
        
        # 筛选特定通路的节点
        pathway_nodes = []
        for node, data in self.metabolic_network.nodes(data=True):
            if data.get('pathway') == pathway or pathway.lower() in str(data.get('name', '')).lower():
                pathway_nodes.append(node)
        
        if not pathway_nodes:
            print(f"未找到 {pathway} 通路的节点")
            return
        
        # 创建子图
        subgraph = self.metabolic_network.subgraph(pathway_nodes)
        
        plt.figure(figsize=(14, 10))
        
        # 设置布局
        pos = nx.spring_layout(subgraph, k=2, iterations=50)
        
        # 绘制节点
        reaction_nodes = [n for n, d in subgraph.nodes(data=True) 
                         if d.get('type') == 'reaction']
        metabolite_nodes = [n for n, d in subgraph.nodes(data=True) 
                           if d.get('type') == 'metabolite']
        
        # 绘制代谢物
        nx.draw_networkx_nodes(subgraph, pos,
                              nodelist=metabolite_nodes,
                              node_color=self.colors['metabolite'],
                              node_size=500,
                              alpha=0.8)
        
        # 绘制反应
        nx.draw_networkx_nodes(subgraph, pos,
                              nodelist=reaction_nodes,
                              node_color=self.colors['reaction'],
                              node_size=300,
                              node_shape='s',
                              alpha=0.8)
        
        # 绘制边
        nx.draw_networkx_edges(subgraph, pos,
                              edge_color='gray',
                              alpha=0.7,
                              arrows=True,
                              arrowsize=15)
        
        # 添加标签
        labels = {n: n.replace('-', '\\n') for n in subgraph.nodes()}
        nx.draw_networkx_labels(subgraph, pos,
                               labels=labels,
                               font_size=10)
        
        plt.title(f'{pathway.upper()} 代谢通路', fontsize=16, fontweight='bold')
        plt.axis('off')
        
        # 保存图片
        output_file = self.figure_dir / f"pathway_{pathway}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"{pathway} 通路图已保存到: {output_file}")
    
    def detailed_pathway(self, pathway='tca_cycle'):
        """生成通路的详细分析图"""
        print(f"正在生成 {pathway} 详细分析图...")
        
        # 这里可以添加更详细的通路分析
        # 包括酶活性、通量分布等信息
        
        plt.figure(figsize=(12, 8))
        
        # 示例：创建一个简单的通路流程图
        pathway_steps = {
            'tca_cycle': [
                'Acetyl-CoA', 'Citrate', 'Isocitrate', 
                'α-Ketoglutarate', 'Succinate', 'Fumarate', 
                'Malate', 'Oxaloacetate'
            ],
            'glycolysis': [
                'Glucose', 'G6P', 'F6P', 'FBP', 
                'G3P', 'PEP', 'Pyruvate'
            ]
        }
        
        if pathway in pathway_steps:
            steps = pathway_steps[pathway]
            
            # 创建线性布局
            x_positions = np.linspace(0, 10, len(steps))
            y_positions = [0] * len(steps)
            
            # 绘制步骤
            for i, (x, y, step) in enumerate(zip(x_positions, y_positions, steps)):
                plt.scatter(x, y, s=500, c=self.colors['metabolite'], alpha=0.8)
                plt.text(x, y-0.3, step, ha='center', va='top', fontsize=10)
                
                # 绘制箭头
                if i < len(steps) - 1:
                    plt.arrow(x+0.3, y, 0.4, 0, head_width=0.1, 
                             head_length=0.1, fc='gray', ec='gray')
            
            plt.title(f'{pathway.upper()} 详细通路', fontsize=16, fontweight='bold')
            plt.xlim(-1, 11)
            plt.ylim(-1, 1)
            plt.axis('off')
        else:
            plt.text(0.5, 0.5, f'通路 {pathway} 的详细信息暂未实现', 
                    ha='center', va='center', transform=plt.gca().transAxes,
                    fontsize=14)
            plt.axis('off')
        
        # 保存图片
        output_file = self.figure_dir / f"detailed_{pathway}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"{pathway} 详细图已保存到: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='代谢网络可视化器')
    parser.add_argument('--visualize-network', action='store_true',
                       help='生成完整代谢网络图')
    parser.add_argument('--pathway-heatmap', action='store_true',
                       help='生成通路完整性热图')
    parser.add_argument('--nutrition-summary', action='store_true',
                       help='生成营养需求总结图')
    parser.add_argument('--pathway-specific', action='store_true',
                       help='生成特定通路图')
    parser.add_argument('--detailed-pathway', action='store_true',
                       help='生成详细通路图')
    parser.add_argument('--pathway', type=str, default='glycolysis',
                       help='指定通路名称')
    parser.add_argument('--output-dir', type=str, default='results',
                       help='数据输入目录')
    parser.add_argument('--figure-dir', type=str, default='figures',
                       help='图片输出目录')
    
    args = parser.parse_args()
    
    # 创建可视化器实例
    visualizer = NetworkVisualizer(output_dir=args.output_dir, 
                                  figure_dir=args.figure_dir)
    
    if args.visualize_network:
        visualizer.visualize_network()
    
    if args.pathway_heatmap:
        visualizer.pathway_heatmap()
    
    if args.nutrition_summary:
        visualizer.nutrition_summary()
    
    if args.pathway_specific:
        visualizer.pathway_specific(args.pathway)
    
    if args.detailed_pathway:
        visualizer.detailed_pathway(args.pathway)
    
    if not any([args.visualize_network, args.pathway_heatmap, 
                args.nutrition_summary, args.pathway_specific,
                args.detailed_pathway]):
        # 如果没有指定参数，生成所有图表
        print("生成所有可视化图表...")
        visualizer.visualize_network()
        visualizer.pathway_heatmap()
        visualizer.nutrition_summary()
        visualizer.pathway_specific('glycolysis')
        visualizer.detailed_pathway('tca_cycle')

if __name__ == "__main__":
    main()