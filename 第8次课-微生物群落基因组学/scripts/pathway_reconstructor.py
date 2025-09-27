#!/usr/bin/env python3
"""
代谢通路重建脚本
用于重建微生物群落的代谢网络和通路

作者: 微生物基因组学课程组
版本: 1.0.0
"""

import os
import sys
import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
import argparse
import logging
import json
from collections import defaultdict, Counter

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PathwayReconstructor:
    """代谢通路重建类"""
    
    def __init__(self, input_dir="results/functional_analysis", output_dir="results", threads=8):
        """
        初始化代谢通路重建器
        
        参数:
        input_dir: 功能分析结果目录
        output_dir: 输出结果目录
        threads: 使用的线程数
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.threads = threads
        
        # 创建输出目录
        self.pathway_dir = self.output_dir / "pathway_reconstruction"
        self.pathway_dir.mkdir(parents=True, exist_ok=True)
        
        # KEGG通路数据库（简化版）
        self.kegg_pathways = self._load_kegg_pathways()
        self.kegg_reactions = self._load_kegg_reactions()
        self.kegg_compounds = self._load_kegg_compounds()
    
    def _load_kegg_pathways(self):
        """加载KEGG通路信息"""
        pathways = {
            'ko00010': {
                'name': 'Glycolysis / Gluconeogenesis',
                'category': 'Carbohydrate metabolism',
                'genes': ['K00001', 'K00002', 'K00003', 'K00844', 'K01810'],
                'reactions': ['R00200', 'R00299', 'R00658', 'R01015', 'R01512']
            },
            'ko00020': {
                'name': 'Citrate cycle (TCA cycle)',
                'category': 'Carbohydrate metabolism',
                'genes': ['K00024', 'K00025', 'K00026', 'K00030', 'K00031'],
                'reactions': ['R00351', 'R00405', 'R00709', 'R01082', 'R01900']
            },
            'ko00260': {
                'name': 'Glycine, serine and threonine metabolism',
                'category': 'Amino acid metabolism',
                'genes': ['K00003', 'K00058', 'K00281', 'K00600', 'K01695'],
                'reactions': ['R00945', 'R01001', 'R01775', 'R02291', 'R03260']
            },
            'ko00561': {
                'name': 'Glycerolipid metabolism',
                'category': 'Lipid metabolism',
                'genes': ['K00010', 'K00111', 'K00631', 'K01054', 'K02446'],
                'reactions': ['R00847', 'R01011', 'R01369', 'R02250', 'R02687']
            },
            'ko00620': {
                'name': 'Pyruvate metabolism',
                'category': 'Carbohydrate metabolism',
                'genes': ['K00016', 'K00161', 'K00169', 'K00382', 'K01006'],
                'reactions': ['R00014', 'R00200', 'R00344', 'R00709', 'R01196']
            }
        }
        return pathways
    
    def _load_kegg_reactions(self):
        """加载KEGG反应信息"""
        reactions = {
            'R00200': {
                'name': 'Aldehyde:NAD+ oxidoreductase',
                'substrates': ['C00118'],  # D-Glucose
                'products': ['C00221'],    # D-Glucose 6-phosphate
                'enzymes': ['K00001']
            },
            'R00299': {
                'name': 'ATP:D-glucose 6-phosphotransferase',
                'substrates': ['C00221'],  # D-Glucose 6-phosphate
                'products': ['C00085'],    # D-Fructose 6-phosphate
                'enzymes': ['K00844']
            },
            'R00658': {
                'name': 'ATP:D-fructose-6-phosphate 1-phosphotransferase',
                'substrates': ['C00085'],  # D-Fructose 6-phosphate
                'products': ['C00354'],    # D-Fructose 1,6-bisphosphate
                'enzymes': ['K00850']
            }
        }
        return reactions
    
    def _load_kegg_compounds(self):
        """加载KEGG化合物信息"""
        compounds = {
            'C00118': {'name': 'D-Glucose', 'formula': 'C6H12O6'},
            'C00221': {'name': 'D-Glucose 6-phosphate', 'formula': 'C6H13O9P'},
            'C00085': {'name': 'D-Fructose 6-phosphate', 'formula': 'C6H13O9P'},
            'C00354': {'name': 'D-Fructose 1,6-bisphosphate', 'formula': 'C6H14O12P2'},
            'C00022': {'name': 'Pyruvate', 'formula': 'C3H4O3'},
            'C00024': {'name': 'Acetyl-CoA', 'formula': 'C23H38N7O17P3S'},
            'C00036': {'name': 'Oxaloacetate', 'formula': 'C4H4O5'},
            'C00158': {'name': 'Citrate', 'formula': 'C6H8O7'}
        }
        return compounds
    
    def load_functional_data(self):
        """加载功能分析数据"""
        logger.info("加载功能分析数据...")
        
        # 加载KEGG注释结果
        kegg_file = self.input_dir / "kegg_annotation.tsv"
        if kegg_file.exists():
            kegg_df = pd.read_csv(kegg_file, sep='\t')
        else:
            logger.warning("KEGG注释文件不存在，创建模拟数据")
            kegg_df = self._create_mock_kegg_data()
        
        # 加载基因家族丰度数据
        families_file = self.input_dir / "gene_families.tsv"
        if families_file.exists():
            families_df = pd.read_csv(families_file, sep='\t')
        else:
            logger.warning("基因家族文件不存在，创建模拟数据")
            families_df = self._create_mock_families_data()
        
        return kegg_df, families_df
    
    def _create_mock_kegg_data(self):
        """创建模拟KEGG数据"""
        mock_data = [
            {'gene_id': 'gene_001', 'kegg_id': 'K00001', 'description': 'alcohol dehydrogenase', 'pathway': 'ko00010'},
            {'gene_id': 'gene_002', 'kegg_id': 'K00002', 'description': 'aldehyde dehydrogenase', 'pathway': 'ko00010'},
            {'gene_id': 'gene_003', 'kegg_id': 'K00003', 'description': 'homoserine dehydrogenase', 'pathway': 'ko00260'},
            {'gene_id': 'gene_004', 'kegg_id': 'K00010', 'description': 'glycerol kinase', 'pathway': 'ko00561'},
            {'gene_id': 'gene_005', 'kegg_id': 'K00016', 'description': 'L-lactate dehydrogenase', 'pathway': 'ko00620'},
            {'gene_id': 'gene_006', 'kegg_id': 'K00024', 'description': 'malate dehydrogenase', 'pathway': 'ko00020'},
            {'gene_id': 'gene_007', 'kegg_id': 'K00844', 'description': 'hexokinase', 'pathway': 'ko00010'},
            {'gene_id': 'gene_008', 'kegg_id': 'K01810', 'description': 'glucose-6-phosphate isomerase', 'pathway': 'ko00010'}
        ]
        return pd.DataFrame(mock_data)
    
    def _create_mock_families_data(self):
        """创建模拟基因家族数据"""
        mock_data = [
            {'family_id': 'K00001', 'family_type': 'KEGG', 'description': 'alcohol dehydrogenase', 'abundance': 150},
            {'family_id': 'K00002', 'family_type': 'KEGG', 'description': 'aldehyde dehydrogenase', 'abundance': 120},
            {'family_id': 'K00003', 'family_type': 'KEGG', 'description': 'homoserine dehydrogenase', 'abundance': 80},
            {'family_id': 'K00010', 'family_type': 'KEGG', 'description': 'glycerol kinase', 'abundance': 200},
            {'family_id': 'K00016', 'family_type': 'KEGG', 'description': 'L-lactate dehydrogenase', 'abundance': 180},
            {'family_id': 'K00024', 'family_type': 'KEGG', 'description': 'malate dehydrogenase', 'abundance': 160},
            {'family_id': 'K00844', 'family_type': 'KEGG', 'description': 'hexokinase', 'abundance': 140},
            {'family_id': 'K01810', 'family_type': 'KEGG', 'description': 'glucose-6-phosphate isomerase', 'abundance': 110}
        ]
        return pd.DataFrame(mock_data)
    
    def analyze_pathway_completeness(self, kegg_df):
        """分析代谢通路完整性"""
        logger.info("分析代谢通路完整性...")
        
        pathway_completeness = {}
        
        # 获取检测到的KEGG基因
        detected_genes = set(kegg_df['kegg_id'].unique())
        
        for pathway_id, pathway_info in self.kegg_pathways.items():
            required_genes = set(pathway_info['genes'])
            present_genes = detected_genes.intersection(required_genes)
            
            completeness = len(present_genes) / len(required_genes) * 100
            
            pathway_completeness[pathway_id] = {
                'name': pathway_info['name'],
                'category': pathway_info['category'],
                'total_genes': len(required_genes),
                'present_genes': len(present_genes),
                'completeness': round(completeness, 2),
                'missing_genes': list(required_genes - present_genes),
                'present_gene_list': list(present_genes)
            }
        
        # 保存通路完整性分析结果
        completeness_file = self.pathway_dir / "pathway_completeness.json"
        with open(completeness_file, 'w') as f:
            json.dump(pathway_completeness, f, indent=2)
        
        logger.info(f"通路完整性分析结果保存到: {completeness_file}")
        return pathway_completeness
    
    def build_metabolic_network(self, kegg_df, families_df):
        """构建代谢网络"""
        logger.info("构建代谢网络...")
        
        # 创建网络图
        G = nx.DiGraph()
        
        # 添加化合物节点
        for compound_id, compound_info in self.kegg_compounds.items():
            G.add_node(compound_id, 
                      type='compound',
                      name=compound_info['name'],
                      formula=compound_info['formula'])
        
        # 添加反应节点和边
        detected_genes = set(kegg_df['kegg_id'].unique())
        gene_abundance = dict(zip(families_df['family_id'], families_df['abundance']))
        
        for reaction_id, reaction_info in self.kegg_reactions.items():
            # 检查反应所需的酶是否存在
            reaction_enzymes = set(reaction_info['enzymes'])
            present_enzymes = detected_genes.intersection(reaction_enzymes)
            
            if present_enzymes:
                # 添加反应节点
                enzyme_abundance = sum(gene_abundance.get(enzyme, 0) for enzyme in present_enzymes)
                
                G.add_node(reaction_id,
                          type='reaction',
                          name=reaction_info['name'],
                          enzymes=list(present_enzymes),
                          abundance=enzyme_abundance)
                
                # 添加底物到反应的边
                for substrate in reaction_info['substrates']:
                    if substrate in G:
                        G.add_edge(substrate, reaction_id, type='substrate')
                
                # 添加反应到产物的边
                for product in reaction_info['products']:
                    if product in G:
                        G.add_edge(reaction_id, product, type='product')
        
        # 保存网络
        network_file = self.pathway_dir / "metabolic_network.gml"
        nx.write_gml(G, network_file)
        
        logger.info(f"代谢网络已保存到: {network_file}")
        return G
    
    def analyze_network_properties(self, G):
        """分析网络拓扑性质"""
        logger.info("分析网络拓扑性质...")
        
        # 基本网络统计
        network_stats = {
            'nodes': G.number_of_nodes(),
            'edges': G.number_of_edges(),
            'compound_nodes': len([n for n, d in G.nodes(data=True) if d.get('type') == 'compound']),
            'reaction_nodes': len([n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']),
            'density': nx.density(G),
            'is_connected': nx.is_weakly_connected(G)
        }
        
        # 节点度分布
        degrees = dict(G.degree())
        network_stats['average_degree'] = np.mean(list(degrees.values()))
        network_stats['max_degree'] = max(degrees.values())
        
        # 中心性分析
        try:
            betweenness = nx.betweenness_centrality(G)
            closeness = nx.closeness_centrality(G)
            
            # 找出中心性最高的节点
            top_betweenness = sorted(betweenness.items(), key=lambda x: x[1], reverse=True)[:5]
            top_closeness = sorted(closeness.items(), key=lambda x: x[1], reverse=True)[:5]
            
            network_stats['top_betweenness_nodes'] = [(node, round(score, 4)) for node, score in top_betweenness]
            network_stats['top_closeness_nodes'] = [(node, round(score, 4)) for node, score in top_closeness]
            
        except:
            logger.warning("无法计算中心性指标（可能由于网络不连通）")
        
        # 保存网络性质分析
        stats_file = self.pathway_dir / "network_properties.json"
        with open(stats_file, 'w') as f:
            json.dump(network_stats, f, indent=2)
        
        logger.info(f"网络性质分析结果保存到: {stats_file}")
        return network_stats
    
    def identify_key_metabolites(self, G):
        """识别关键代谢物"""
        logger.info("识别关键代谢物...")
        
        key_metabolites = []
        
        # 基于度中心性识别关键代谢物
        compound_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'compound']
        
        for node in compound_nodes:
            degree = G.degree(node)
            in_degree = G.in_degree(node)
            out_degree = G.out_degree(node)
            
            node_info = G.nodes[node]
            
            key_metabolites.append({
                'compound_id': node,
                'name': node_info.get('name', 'Unknown'),
                'formula': node_info.get('formula', 'Unknown'),
                'total_degree': degree,
                'in_degree': in_degree,
                'out_degree': out_degree,
                'hub_score': degree  # 简单的hub评分
            })
        
        # 按度排序
        key_metabolites.sort(key=lambda x: x['hub_score'], reverse=True)
        
        # 保存关键代谢物
        metabolites_file = self.pathway_dir / "key_metabolites.json"
        with open(metabolites_file, 'w') as f:
            json.dump(key_metabolites, f, indent=2)
        
        logger.info(f"关键代谢物分析结果保存到: {metabolites_file}")
        return key_metabolites
    
    def generate_pathway_report(self, pathway_completeness, network_stats, key_metabolites):
        """生成代谢通路重建报告"""
        logger.info("生成代谢通路重建报告...")
        
        report_file = self.pathway_dir / "pathway_reconstruction_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("代谢通路重建分析报告\n")
            f.write("=" * 50 + "\n\n")
            
            # 通路完整性总结
            f.write("1. 代谢通路完整性分析\n")
            f.write("-" * 30 + "\n")
            
            complete_pathways = [p for p in pathway_completeness.values() if p['completeness'] >= 80]
            partial_pathways = [p for p in pathway_completeness.values() if 50 <= p['completeness'] < 80]
            incomplete_pathways = [p for p in pathway_completeness.values() if p['completeness'] < 50]
            
            f.write(f"完整通路 (≥80%): {len(complete_pathways)}\n")
            f.write(f"部分完整通路 (50-80%): {len(partial_pathways)}\n")
            f.write(f"不完整通路 (<50%): {len(incomplete_pathways)}\n\n")
            
            # 详细通路信息
            f.write("主要代谢通路状态:\n")
            for pathway_id, info in pathway_completeness.items():
                status = "完整" if info['completeness'] >= 80 else "部分" if info['completeness'] >= 50 else "缺失"
                f.write(f"- {info['name']}: {info['completeness']}% ({status})\n")
            f.write("\n")
            
            # 网络拓扑性质
            f.write("2. 代谢网络拓扑性质\n")
            f.write("-" * 30 + "\n")
            f.write(f"网络节点数: {network_stats['nodes']}\n")
            f.write(f"网络边数: {network_stats['edges']}\n")
            f.write(f"化合物节点: {network_stats['compound_nodes']}\n")
            f.write(f"反应节点: {network_stats['reaction_nodes']}\n")
            f.write(f"网络密度: {network_stats['density']:.4f}\n")
            f.write(f"网络连通性: {'是' if network_stats['is_connected'] else '否'}\n")
            f.write(f"平均度: {network_stats['average_degree']:.2f}\n\n")
            
            # 关键代谢物
            f.write("3. 关键代谢物 (Top 5)\n")
            f.write("-" * 30 + "\n")
            for i, metabolite in enumerate(key_metabolites[:5], 1):
                f.write(f"{i}. {metabolite['name']} ({metabolite['compound_id']})\n")
                f.write(f"   连接度: {metabolite['total_degree']}, 分子式: {metabolite['formula']}\n")
            f.write("\n")
            
            # 功能预测
            f.write("4. 群落代谢功能预测\n")
            f.write("-" * 30 + "\n")
            
            if len(complete_pathways) >= 3:
                f.write("该微生物群落具有较完整的核心代谢能力:\n")
                f.write("- 能够进行基本的碳水化合物代谢\n")
                f.write("- 具备氨基酸合成和分解能力\n")
                f.write("- 可能具有完整的能量代谢途径\n")
            elif len(complete_pathways) >= 1:
                f.write("该微生物群落具有部分代谢功能:\n")
                f.write("- 某些核心代谢通路完整\n")
                f.write("- 可能依赖环境提供部分营养物质\n")
                f.write("- 群落内可能存在代谢互补\n")
            else:
                f.write("该微生物群落代谢功能有限:\n")
                f.write("- 大多数代谢通路不完整\n")
                f.write("- 可能高度依赖环境营养供应\n")
                f.write("- 需要进一步分析特殊代谢途径\n")
        
        logger.info(f"代谢通路重建报告保存到: {report_file}")
    
    def run_reconstruction(self):
        """运行完整的代谢通路重建流程"""
        logger.info("开始代谢通路重建...")
        
        # 加载功能数据
        kegg_df, families_df = self.load_functional_data()
        
        # 分析通路完整性
        pathway_completeness = self.analyze_pathway_completeness(kegg_df)
        
        # 构建代谢网络
        metabolic_network = self.build_metabolic_network(kegg_df, families_df)
        
        # 分析网络性质
        network_stats = self.analyze_network_properties(metabolic_network)
        
        # 识别关键代谢物
        key_metabolites = self.identify_key_metabolites(metabolic_network)
        
        # 生成报告
        self.generate_pathway_report(pathway_completeness, network_stats, key_metabolites)
        
        logger.info("代谢通路重建完成!")
        
        return {
            'pathway_completeness': pathway_completeness,
            'network_stats': network_stats,
            'key_metabolites': key_metabolites,
            'metabolic_network': metabolic_network
        }

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='代谢通路重建工具')
    parser.add_argument('--input_dir', default='results/functional_analysis', help='功能分析结果目录')
    parser.add_argument('--output_dir', default='results', help='输出结果目录')
    parser.add_argument('--threads', type=int, default=8, help='使用的线程数')
    
    args = parser.parse_args()
    
    # 创建重建器实例
    reconstructor = PathwayReconstructor(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        threads=args.threads
    )
    
    # 运行重建
    results = reconstructor.run_reconstruction()
    
    if results:
        print("\n=== 代谢通路重建完成 ===")
        complete_pathways = [p for p in results['pathway_completeness'].values() if p['completeness'] >= 80]
        print(f"完整代谢通路数: {len(complete_pathways)}")
        print(f"网络节点数: {results['network_stats']['nodes']}")
        print(f"网络边数: {results['network_stats']['edges']}")
        print(f"关键代谢物数: {len(results['key_metabolites'])}")
        print("\n详细结果请查看 results/pathway_reconstruction/ 目录")
    else:
        print("代谢通路重建失败")
        sys.exit(1)

if __name__ == "__main__":
    main()