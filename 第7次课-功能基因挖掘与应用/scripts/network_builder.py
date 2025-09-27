#!/usr/bin/env python3
"""
Functional Network Builder
用于构建功能基因关联网络

作者: 微生物基因组学课程组
版本: 1.0.0
"""

import argparse
import os
import sys
from pathlib import Path
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict

class NetworkBuilder:
    """功能网络构建主类"""
    
    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 功能类别颜色映射
        self.category_colors = {
            'BGC': '#FF6B6B',
            'Resistance': '#4ECDC4', 
            'Virulence': '#45B7D1',
            'Metabolism': '#96CEB4',
            'Regulation': '#FFEAA7',
            'Transport': '#DDA0DD'
        }
    
    def integrate_data(self, bgc_file, resistance_file, virulence_file):
        """整合多种分析结果"""
        print("正在整合分析数据...")
        
        integrated_data = []
        
        # 读取BGC数据
        if bgc_file and os.path.exists(bgc_file):
            print(f"读取BGC数据: {bgc_file}")
            try:
                with open(bgc_file, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if 'BGC_' in line and 'Type:' in line:
                            # 解析BGC信息
                            parts = line.strip().split()
                            bgc_id = parts[0] if parts else f"BGC_{len(integrated_data)+1}"
                            bgc_type = 'PKS' if 'PKS' in line else 'NRPS' if 'NRPS' in line else 'Other'
                            
                            integrated_data.append({
                                'gene_id': bgc_id,
                                'gene_name': bgc_id,
                                'category': 'BGC',
                                'subcategory': bgc_type,
                                'function': f'{bgc_type} biosynthetic gene cluster',
                                'score': 90.0
                            })
            except Exception as e:
                print(f"读取BGC文件出错: {e}")
        
        # 读取抗性基因数据
        if resistance_file and os.path.exists(resistance_file):
            print(f"读取抗性基因数据: {resistance_file}")
            try:
                with open(resistance_file, 'r') as f:
                    lines = f.readlines()[1:]  # 跳过标题行
                    for line in lines:
                        parts = line.strip().split('\t')
                        if len(parts) >= 8:
                            integrated_data.append({
                                'gene_id': parts[0],
                                'gene_name': parts[1],
                                'category': 'Resistance',
                                'subcategory': parts[2],
                                'function': parts[7],
                                'score': float(parts[4]) if parts[4].replace('.','').isdigit() else 85.0
                            })
            except Exception as e:
                print(f"读取抗性基因文件出错: {e}")
        
        # 读取毒力因子数据
        if virulence_file and os.path.exists(virulence_file):
            print(f"读取毒力因子数据: {virulence_file}")
            try:
                with open(virulence_file, 'r') as f:
                    lines = f.readlines()[1:]  # 跳过标题行
                    for line in lines:
                        parts = line.strip().split('\t')
                        if len(parts) >= 9:
                            integrated_data.append({
                                'gene_id': parts[0],
                                'gene_name': parts[1],
                                'category': 'Virulence',
                                'subcategory': parts[2],
                                'function': parts[8],
                                'score': float(parts[4]) if parts[4].replace('.','').isdigit() else 80.0
                            })
            except Exception as e:
                print(f"读取毒力因子文件出错: {e}")
        
        # 如果没有真实数据，生成模拟数据
        if not integrated_data:
            print("未找到输入数据，生成模拟数据...")
            integrated_data = self._generate_mock_data()
        
        # 保存整合数据
        integrated_file = self.output_dir / "integrated_analysis.txt"
        with open(integrated_file, 'w') as f:
            f.write("Gene_ID\tGene_Name\tCategory\tSubcategory\tFunction\tScore\n")
            for data in integrated_data:
                f.write(f"{data['gene_id']}\t{data['gene_name']}\t{data['category']}\t"
                       f"{data['subcategory']}\t{data['function']}\t{data['score']}\n")
        
        print(f"数据整合完成，共 {len(integrated_data)} 个功能基因")
        print(f"整合数据保存至: {integrated_file}")
        return integrated_data
    
    def _generate_mock_data(self):
        """生成模拟数据"""
        mock_data = [
            # BGC数据
            {'gene_id': 'BGC_001', 'gene_name': 'pksA', 'category': 'BGC', 
             'subcategory': 'PKS', 'function': 'polyketide synthase', 'score': 95.0},
            {'gene_id': 'BGC_002', 'gene_name': 'nrpsB', 'category': 'BGC',
             'subcategory': 'NRPS', 'function': 'nonribosomal peptide synthetase', 'score': 92.0},
            {'gene_id': 'BGC_003', 'gene_name': 'terC', 'category': 'BGC',
             'subcategory': 'Terpene', 'function': 'terpene cyclase', 'score': 88.0},
            
            # 抗性基因数据
            {'gene_id': 'RES_001', 'gene_name': 'mecA', 'category': 'Resistance',
             'subcategory': 'target_alteration', 'function': 'methicillin resistance', 'score': 98.0},
            {'gene_id': 'RES_002', 'gene_name': 'blaZ', 'category': 'Resistance',
             'subcategory': 'enzymatic_inactivation', 'function': 'beta-lactamase', 'score': 96.0},
            {'gene_id': 'RES_003', 'gene_name': 'tetK', 'category': 'Resistance',
             'subcategory': 'efflux_pump', 'function': 'tetracycline efflux', 'score': 94.0},
            
            # 毒力因子数据
            {'gene_id': 'VIR_001', 'gene_name': 'hla', 'category': 'Virulence',
             'subcategory': 'toxin', 'function': 'alpha-hemolysin', 'score': 99.0},
            {'gene_id': 'VIR_002', 'gene_name': 'spa', 'category': 'Virulence',
             'subcategory': 'immune_evasion', 'function': 'protein A', 'score': 97.0},
            {'gene_id': 'VIR_003', 'gene_name': 'clfA', 'category': 'Virulence',
             'subcategory': 'adhesion', 'function': 'clumping factor A', 'score': 95.0}
        ]
        return mock_data
    
    def build_network(self, integrated_data, output_format='graphml'):
        """构建功能关联网络"""
        print("正在构建功能关联网络...")
        
        # 创建网络图
        G = nx.Graph()
        
        # 添加节点
        for data in integrated_data:
            G.add_node(data['gene_id'], 
                      name=data['gene_name'],
                      category=data['category'],
                      subcategory=data['subcategory'],
                      function=data['function'],
                      score=data['score'])
        
        # 添加边（基于功能相关性）
        self._add_functional_edges(G, integrated_data)
        self._add_colocalization_edges(G, integrated_data)
        self._add_regulatory_edges(G, integrated_data)
        
        # 保存网络文件
        if output_format == 'graphml':
            network_file = self.output_dir / "functional_network.graphml"
            nx.write_graphml(G, network_file)
        elif output_format == 'gexf':
            network_file = self.output_dir / "functional_network.gexf"
            nx.write_gexf(G, network_file)
        else:
            network_file = self.output_dir / "functional_network.txt"
            nx.write_edgelist(G, network_file, data=True)
        
        print(f"网络构建完成: {len(G.nodes())} 个节点, {len(G.edges())} 条边")
        print(f"网络文件保存至: {network_file}")
        return G, network_file
    
    def _add_functional_edges(self, G, integrated_data):
        """添加功能相关边"""
        # 同类功能基因之间的连接
        category_groups = defaultdict(list)
        for data in integrated_data:
            category_groups[data['category']].append(data['gene_id'])
        
        for category, genes in category_groups.items():
            for i in range(len(genes)):
                for j in range(i+1, len(genes)):
                    G.add_edge(genes[i], genes[j], 
                              edge_type='functional_similarity',
                              weight=0.7)
    
    def _add_colocalization_edges(self, G, integrated_data):
        """添加共定位边"""
        # 模拟基因组共定位关系
        colocalized_pairs = [
            ('BGC_001', 'BGC_002'),  # PKS和NRPS常共定位
            ('RES_001', 'RES_002'),  # 抗性基因常聚集
            ('VIR_001', 'VIR_002'),  # 毒力因子常聚集
        ]
        
        for gene1, gene2 in colocalized_pairs:
            if G.has_node(gene1) and G.has_node(gene2):
                G.add_edge(gene1, gene2, 
                          edge_type='colocalization',
                          weight=0.8)
    
    def _add_regulatory_edges(self, G, integrated_data):
        """添加调控关系边"""
        # 模拟调控关系
        regulatory_pairs = [
            ('VIR_002', 'VIR_001'),  # spa调控hla
            ('BGC_001', 'RES_001'),  # 次级代谢与抗性的关联
        ]
        
        for regulator, target in regulatory_pairs:
            if G.has_node(regulator) and G.has_node(target):
                G.add_edge(regulator, target,
                          edge_type='regulation',
                          weight=0.6)
    
    def calculate_network_statistics(self, G):
        """计算网络统计指标"""
        print("正在计算网络统计指标...")
        
        stats = {
            'nodes': len(G.nodes()),
            'edges': len(G.edges()),
            'density': nx.density(G),
            'clustering_coefficient': nx.average_clustering(G),
            'average_path_length': 0,
            'diameter': 0,
            'centrality_measures': {}
        }
        
        # 计算路径长度（仅对连通图）
        if nx.is_connected(G):
            stats['average_path_length'] = nx.average_shortest_path_length(G)
            stats['diameter'] = nx.diameter(G)
        else:
            # 对于非连通图，计算最大连通分量的指标
            largest_cc = max(nx.connected_components(G), key=len)
            subgraph = G.subgraph(largest_cc)
            if len(subgraph) > 1:
                stats['average_path_length'] = nx.average_shortest_path_length(subgraph)
                stats['diameter'] = nx.diameter(subgraph)
        
        # 计算中心性指标
        stats['centrality_measures'] = {
            'degree_centrality': nx.degree_centrality(G),
            'betweenness_centrality': nx.betweenness_centrality(G),
            'closeness_centrality': nx.closeness_centrality(G),
            'eigenvector_centrality': nx.eigenvector_centrality(G, max_iter=1000)
        }
        
        # 保存统计结果
        stats_file = self.output_dir / "network_statistics.txt"
        with open(stats_file, 'w') as f:
            f.write("=== 功能网络统计分析 ===\n\n")
            f.write(f"节点数量: {stats['nodes']}\n")
            f.write(f"边数量: {stats['edges']}\n")
            f.write(f"网络密度: {stats['density']:.4f}\n")
            f.write(f"聚类系数: {stats['clustering_coefficient']:.4f}\n")
            f.write(f"平均路径长度: {stats['average_path_length']:.4f}\n")
            f.write(f"网络直径: {stats['diameter']}\n\n")
            
            f.write("中心性分析 (Top 5):\n")
            for centrality_type, centrality_values in stats['centrality_measures'].items():
                f.write(f"\n{centrality_type}:\n")
                sorted_nodes = sorted(centrality_values.items(), 
                                    key=lambda x: x[1], reverse=True)[:5]
                for node, value in sorted_nodes:
                    node_name = G.nodes[node].get('name', node)
                    f.write(f"  {node_name} ({node}): {value:.4f}\n")
        
        print(f"网络统计完成，结果保存至: {stats_file}")
        return stats
    
    def visualize_network(self, G, output_dir):
        """网络可视化"""
        print("正在生成网络可视化...")
        
        # 创建图表目录
        figures_dir = Path(output_dir) / "figures"
        figures_dir.mkdir(exist_ok=True)
        
        plt.figure(figsize=(16, 12))
        
        # 1. 主网络图
        plt.subplot(2, 2, 1)
        
        # 设置节点位置
        pos = nx.spring_layout(G, k=1, iterations=50)
        
        # 按类别着色节点
        node_colors = []
        for node in G.nodes():
            category = G.nodes[node].get('category', 'Unknown')
            node_colors.append(self.category_colors.get(category, '#CCCCCC'))
        
        # 绘制网络
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, 
                              node_size=500, alpha=0.8)
        nx.draw_networkx_edges(G, pos, alpha=0.5, width=1)
        nx.draw_networkx_labels(G, pos, 
                               {node: G.nodes[node].get('name', node)[:6] 
                                for node in G.nodes()}, 
                               font_size=8)
        
        plt.title('功能基因关联网络')
        plt.axis('off')
        
        # 2. 度分布
        plt.subplot(2, 2, 2)
        degrees = [G.degree(n) for n in G.nodes()]
        plt.hist(degrees, bins=max(1, len(set(degrees))), color='skyblue', alpha=0.7)
        plt.xlabel('节点度数')
        plt.ylabel('频次')
        plt.title('度分布')
        
        # 3. 类别分布饼图
        plt.subplot(2, 2, 3)
        category_counts = defaultdict(int)
        for node in G.nodes():
            category = G.nodes[node].get('category', 'Unknown')
            category_counts[category] += 1
        
        plt.pie(category_counts.values(), labels=category_counts.keys(),
                autopct='%1.1f%%', colors=[self.category_colors.get(cat, '#CCCCCC') 
                                          for cat in category_counts.keys()])
        plt.title('功能类别分布')
        
        # 4. 中心性分析
        plt.subplot(2, 2, 4)
        centrality = nx.degree_centrality(G)
        nodes = list(centrality.keys())
        values = list(centrality.values())
        
        # 取前10个节点
        sorted_indices = sorted(range(len(values)), key=lambda i: values[i], reverse=True)[:10]
        top_nodes = [G.nodes[nodes[i]].get('name', nodes[i])[:8] for i in sorted_indices]
        top_values = [values[i] for i in sorted_indices]
        
        plt.barh(range(len(top_nodes)), top_values, color='lightcoral')
        plt.yticks(range(len(top_nodes)), top_nodes)
        plt.xlabel('度中心性')
        plt.title('关键节点分析')
        
        plt.tight_layout()
        plt.savefig(figures_dir / 'functional_network.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 生成单独的大网络图
        plt.figure(figsize=(12, 10))
        pos = nx.spring_layout(G, k=2, iterations=100)
        
        # 绘制边
        edge_colors = []
        edge_widths = []
        for edge in G.edges(data=True):
            edge_type = edge[2].get('edge_type', 'unknown')
            weight = edge[2].get('weight', 0.5)
            
            if edge_type == 'functional_similarity':
                edge_colors.append('#FF6B6B')
            elif edge_type == 'colocalization':
                edge_colors.append('#4ECDC4')
            elif edge_type == 'regulation':
                edge_colors.append('#45B7D1')
            else:
                edge_colors.append('#CCCCCC')
            
            edge_widths.append(weight * 3)
        
        nx.draw_networkx_edges(G, pos, edge_color=edge_colors, 
                              width=edge_widths, alpha=0.6)
        
        # 绘制节点
        for category in self.category_colors:
            category_nodes = [node for node in G.nodes() 
                            if G.nodes[node].get('category') == category]
            if category_nodes:
                nx.draw_networkx_nodes(G, pos, nodelist=category_nodes,
                                     node_color=self.category_colors[category],
                                     node_size=800, alpha=0.8, label=category)
        
        # 添加标签
        labels = {node: G.nodes[node].get('name', node) for node in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=10, font_weight='bold')
        
        plt.title('功能基因关联网络详图', fontsize=16, fontweight='bold')
        plt.legend(loc='upper right')
        plt.axis('off')
        
        plt.savefig(figures_dir / 'network_detailed.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"网络可视化完成，图表保存至: {figures_dir}/")

def main():
    parser = argparse.ArgumentParser(description='功能网络构建工具')
    parser.add_argument('--mode', required=True,
                       choices=['integrate', 'network', 'stats', 'visualize'],
                       help='运行模式')
    parser.add_argument('--bgc', help='BGC分析结果文件')
    parser.add_argument('--resistance', help='抗性基因分析结果文件')
    parser.add_argument('--virulence', help='毒力因子分析结果文件')
    parser.add_argument('--input', help='整合数据文件或网络文件')
    parser.add_argument('--network', help='网络文件')
    parser.add_argument('--output', required=True, help='输出目录')
    parser.add_argument('--format', default='graphml',
                       choices=['graphml', 'gexf', 'txt'],
                       help='网络文件格式')
    
    args = parser.parse_args()
    
    # 创建网络构建器实例
    builder = NetworkBuilder(args.output)
    
    if args.mode == 'integrate':
        # 数据整合模式
        integrated_data = builder.integrate_data(args.bgc, args.resistance, args.virulence)
        print("数据整合完成！")
        
    elif args.mode == 'network':
        # 网络构建模式
        if args.input and os.path.exists(args.input):
            # 从整合数据文件读取
            integrated_data = []
            with open(args.input, 'r') as f:
                lines = f.readlines()[1:]  # 跳过标题行
                for line in lines:
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        integrated_data.append({
                            'gene_id': parts[0],
                            'gene_name': parts[1],
                            'category': parts[2],
                            'subcategory': parts[3],
                            'function': parts[4],
                            'score': float(parts[5])
                        })
        else:
            # 直接整合数据
            integrated_data = builder.integrate_data(args.bgc, args.resistance, args.virulence)
        
        G, network_file = builder.build_network(integrated_data, args.format)
        print("网络构建完成！")
        
    elif args.mode == 'stats':
        # 网络统计模式
        if args.network and os.path.exists(args.network):
            if args.network.endswith('.graphml'):
                G = nx.read_graphml(args.network)
            elif args.network.endswith('.gexf'):
                G = nx.read_gexf(args.network)
            else:
                G = nx.read_edgelist(args.network)
            
            stats = builder.calculate_network_statistics(G)
            print("网络统计完成！")
        else:
            print("错误: 需要提供网络文件")
            sys.exit(1)
            
    elif args.mode == 'visualize':
        # 可视化模式
        if args.network and os.path.exists(args.network):
            if args.network.endswith('.graphml'):
                G = nx.read_graphml(args.network)
            elif args.network.endswith('.gexf'):
                G = nx.read_gexf(args.network)
            else:
                G = nx.read_edgelist(args.network)
            
            builder.visualize_network(G, args.output)
            print("网络可视化完成！")
        else:
            print("错误: 需要提供网络文件")
            sys.exit(1)

if __name__ == "__main__":
    main()