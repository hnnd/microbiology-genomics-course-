#!/usr/bin/env python3
"""
代谢网络重建器
基于酶基因信息重建微生物代谢网络

作者: 微生物基因组学课程组
版本: 1.0.0
"""

import argparse
import pandas as pd
import json
import networkx as nx
from pathlib import Path
import xml.etree.ElementTree as ET
from collections import defaultdict
import logging

class MetabolicReconstructor:
    def __init__(self, output_dir="results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # 设置日志
        logging.basicConfig(
            filename=self.output_dir / "reconstruction.log",
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        # 代谢网络
        self.metabolic_network = nx.DiGraph()
        
        # 反应数据库（简化版本）
        self.reaction_database = self._load_reaction_database()
        
        # 代谢物数据库
        self.metabolite_database = self._load_metabolite_database()
        
    def _load_reaction_database(self):
        """加载反应数据库（简化版本）"""
        # 这里是一些核心代谢反应的简化表示
        reactions = {
            '2.7.1.1': {  # hexokinase
                'name': 'hexokinase',
                'equation': 'glucose + ATP -> glucose-6-phosphate + ADP',
                'substrates': ['glucose', 'ATP'],
                'products': ['glucose-6-phosphate', 'ADP'],
                'pathway': 'glycolysis'
            },
            '5.3.1.9': {  # glucose-6-phosphate isomerase
                'name': 'glucose-6-phosphate isomerase',
                'equation': 'glucose-6-phosphate -> fructose-6-phosphate',
                'substrates': ['glucose-6-phosphate'],
                'products': ['fructose-6-phosphate'],
                'pathway': 'glycolysis'
            },
            '2.7.1.11': {  # 6-phosphofructokinase
                'name': '6-phosphofructokinase',
                'equation': 'fructose-6-phosphate + ATP -> fructose-1,6-bisphosphate + ADP',
                'substrates': ['fructose-6-phosphate', 'ATP'],
                'products': ['fructose-1,6-bisphosphate', 'ADP'],
                'pathway': 'glycolysis'
            },
            '4.1.2.13': {  # fructose-bisphosphate aldolase
                'name': 'fructose-bisphosphate aldolase',
                'equation': 'fructose-1,6-bisphosphate -> glyceraldehyde-3-phosphate + dihydroxyacetone-phosphate',
                'substrates': ['fructose-1,6-bisphosphate'],
                'products': ['glyceraldehyde-3-phosphate', 'dihydroxyacetone-phosphate'],
                'pathway': 'glycolysis'
            },
            '1.2.1.12': {  # glyceraldehyde-3-phosphate dehydrogenase
                'name': 'glyceraldehyde-3-phosphate dehydrogenase',
                'equation': 'glyceraldehyde-3-phosphate + NAD+ + Pi -> 1,3-bisphosphoglycerate + NADH',
                'substrates': ['glyceraldehyde-3-phosphate', 'NAD+', 'phosphate'],
                'products': ['1,3-bisphosphoglycerate', 'NADH'],
                'pathway': 'glycolysis'
            },
            '2.7.2.3': {  # phosphoglycerate kinase
                'name': 'phosphoglycerate kinase',
                'equation': '1,3-bisphosphoglycerate + ADP -> 3-phosphoglycerate + ATP',
                'substrates': ['1,3-bisphosphoglycerate', 'ADP'],
                'products': ['3-phosphoglycerate', 'ATP'],
                'pathway': 'glycolysis'
            },
            '2.7.1.40': {  # pyruvate kinase
                'name': 'pyruvate kinase',
                'equation': 'phosphoenolpyruvate + ADP -> pyruvate + ATP',
                'substrates': ['phosphoenolpyruvate', 'ADP'],
                'products': ['pyruvate', 'ATP'],
                'pathway': 'glycolysis'
            },
            # TCA循环反应
            '2.3.3.1': {  # citrate synthase
                'name': 'citrate synthase',
                'equation': 'acetyl-CoA + oxaloacetate -> citrate + CoA',
                'substrates': ['acetyl-CoA', 'oxaloacetate'],
                'products': ['citrate', 'CoA'],
                'pathway': 'TCA_cycle'
            },
            '4.2.1.3': {  # aconitase
                'name': 'aconitase',
                'equation': 'citrate -> isocitrate',
                'substrates': ['citrate'],
                'products': ['isocitrate'],
                'pathway': 'TCA_cycle'
            },
            '1.1.1.41': {  # isocitrate dehydrogenase
                'name': 'isocitrate dehydrogenase',
                'equation': 'isocitrate + NAD+ -> alpha-ketoglutarate + NADH + CO2',
                'substrates': ['isocitrate', 'NAD+'],
                'products': ['alpha-ketoglutarate', 'NADH', 'CO2'],
                'pathway': 'TCA_cycle'
            }
        }
        return reactions
    
    def _load_metabolite_database(self):
        """加载代谢物数据库"""
        metabolites = {
            'glucose': {'name': 'Glucose', 'formula': 'C6H12O6', 'charge': 0},
            'ATP': {'name': 'Adenosine triphosphate', 'formula': 'C10H16N5O13P3', 'charge': -4},
            'ADP': {'name': 'Adenosine diphosphate', 'formula': 'C10H15N5O10P2', 'charge': -3},
            'glucose-6-phosphate': {'name': 'Glucose 6-phosphate', 'formula': 'C6H13O9P', 'charge': -2},
            'fructose-6-phosphate': {'name': 'Fructose 6-phosphate', 'formula': 'C6H13O9P', 'charge': -2},
            'fructose-1,6-bisphosphate': {'name': 'Fructose 1,6-bisphosphate', 'formula': 'C6H14O12P2', 'charge': -4},
            'pyruvate': {'name': 'Pyruvate', 'formula': 'C3H4O3', 'charge': -1},
            'acetyl-CoA': {'name': 'Acetyl-CoA', 'formula': 'C23H38N7O17P3S', 'charge': -4},
            'citrate': {'name': 'Citrate', 'formula': 'C6H8O7', 'charge': -3},
            'oxaloacetate': {'name': 'Oxaloacetate', 'formula': 'C4H4O5', 'charge': -2},
            'NAD+': {'name': 'Nicotinamide adenine dinucleotide', 'formula': 'C21H27N7O14P2', 'charge': -1},
            'NADH': {'name': 'Nicotinamide adenine dinucleotide (reduced)', 'formula': 'C21H29N7O14P2', 'charge': -2}
        }
        return metabolites
    
    def auto_reconstruct(self, enzyme_file="results/enzyme_genes.txt"):
        """自动重建代谢网络"""
        self.logger.info("开始自动代谢网络重建")
        print("正在进行自动代谢网络重建...")
        
        try:
            # 读取酶基因文件
            if not Path(enzyme_file).exists():
                print("酶基因文件不存在，请先运行酶基因提取")
                return
            
            enzyme_df = pd.read_csv(enzyme_file, sep='\\t')
            ec_numbers = set(enzyme_df['EC_Number'].unique())
            
            self.logger.info(f"发现 {len(ec_numbers)} 个不同的EC编号")
            
            # 构建网络
            reactions_added = 0
            metabolites_added = set()
            
            for ec_number in ec_numbers:
                if ec_number in self.reaction_database:
                    reaction = self.reaction_database[ec_number]
                    
                    # 添加反应节点
                    reaction_id = f"R_{ec_number}"
                    self.metabolic_network.add_node(
                        reaction_id,
                        type='reaction',
                        ec_number=ec_number,
                        name=reaction['name'],
                        pathway=reaction['pathway']
                    )
                    
                    # 添加代谢物和连接
                    for substrate in reaction['substrates']:
                        if substrate not in metabolites_added:
                            self.metabolic_network.add_node(
                                substrate,
                                type='metabolite',
                                **self.metabolite_database.get(substrate, {'name': substrate})
                            )
                            metabolites_added.add(substrate)
                        
                        # 底物到反应的边
                        self.metabolic_network.add_edge(substrate, reaction_id, type='substrate')
                    
                    for product in reaction['products']:
                        if product not in metabolites_added:
                            self.metabolic_network.add_node(
                                product,
                                type='metabolite',
                                **self.metabolite_database.get(product, {'name': product})
                            )
                            metabolites_added.add(product)
                        
                        # 反应到产物的边
                        self.metabolic_network.add_edge(reaction_id, product, type='product')
                    
                    reactions_added += 1
            
            self.logger.info(f"重建完成: {reactions_added} 个反应, {len(metabolites_added)} 个代谢物")
            print(f"网络重建完成:")
            print(f"  反应数: {reactions_added}")
            print(f"  代谢物数: {len(metabolites_added)}")
            print(f"  总节点数: {self.metabolic_network.number_of_nodes()}")
            print(f"  总边数: {self.metabolic_network.number_of_edges()}")
            
            # 保存网络
            self._save_network()
            
        except Exception as e:
            self.logger.error(f"自动重建失败: {e}")
            print(f"重建过程中出错: {e}")
    
    def validate_network(self):
        """验证重建网络的质量"""
        print("正在验证网络质量...")
        
        if self.metabolic_network.number_of_nodes() == 0:
            print("网络为空，请先进行重建")
            return
        
        validation_results = {}
        
        # 基本统计
        validation_results['total_nodes'] = self.metabolic_network.number_of_nodes()
        validation_results['total_edges'] = self.metabolic_network.number_of_edges()
        
        # 节点类型统计
        reaction_nodes = [n for n, d in self.metabolic_network.nodes(data=True) if d.get('type') == 'reaction']
        metabolite_nodes = [n for n, d in self.metabolic_network.nodes(data=True) if d.get('type') == 'metabolite']
        
        validation_results['reaction_count'] = len(reaction_nodes)
        validation_results['metabolite_count'] = len(metabolite_nodes)
        
        # 连通性检查
        validation_results['is_connected'] = nx.is_weakly_connected(self.metabolic_network)
        validation_results['connected_components'] = nx.number_weakly_connected_components(self.metabolic_network)
        
        # 度分布
        degrees = dict(self.metabolic_network.degree())
        validation_results['avg_degree'] = sum(degrees.values()) / len(degrees) if degrees else 0
        validation_results['max_degree'] = max(degrees.values()) if degrees else 0
        
        # 孤立节点
        isolated_nodes = list(nx.isolates(self.metabolic_network))
        validation_results['isolated_nodes'] = len(isolated_nodes)
        
        # 生成验证报告
        report = self._generate_validation_report(validation_results)
        
        # 保存验证结果
        validation_file = self.output_dir / "network_validation.txt"
        with open(validation_file, 'w') as f:
            f.write(report)
        
        print(f"网络验证完成，报告保存到: {validation_file}")
        
        return validation_results
    
    def _generate_validation_report(self, results):
        """生成验证报告"""
        report = "代谢网络质量验证报告\\n"
        report += "=" * 40 + "\\n\\n"
        
        report += f"网络规模:\\n"
        report += f"  总节点数: {results['total_nodes']}\\n"
        report += f"  总边数: {results['total_edges']}\\n"
        report += f"  反应数: {results['reaction_count']}\\n"
        report += f"  代谢物数: {results['metabolite_count']}\\n\\n"
        
        report += f"网络拓扑:\\n"
        report += f"  连通性: {'是' if results['is_connected'] else '否'}\\n"
        report += f"  连通组件数: {results['connected_components']}\\n"
        report += f"  平均度: {results['avg_degree']:.2f}\\n"
        report += f"  最大度: {results['max_degree']}\\n"
        report += f"  孤立节点数: {results['isolated_nodes']}\\n\\n"
        
        # 质量评估
        quality_score = self._calculate_quality_score(results)
        report += f"质量评分: {quality_score:.1f}/100\\n"
        
        return report
    
    def _calculate_quality_score(self, results):
        """计算网络质量评分"""
        score = 0
        
        # 规模评分 (30分)
        if results['reaction_count'] >= 100:
            score += 30
        elif results['reaction_count'] >= 50:
            score += 20
        elif results['reaction_count'] >= 20:
            score += 10
        
        # 连通性评分 (40分)
        if results['is_connected']:
            score += 40
        elif results['connected_components'] <= 3:
            score += 20
        
        # 完整性评分 (30分)
        if results['isolated_nodes'] == 0:
            score += 30
        elif results['isolated_nodes'] <= 5:
            score += 20
        elif results['isolated_nodes'] <= 10:
            score += 10
        
        return score
    
    def identify_gaps(self):
        """识别网络中的代谢间隙"""
        print("正在识别代谢间隙...")
        
        gaps = []
        
        # 查找死端代谢物（只有输入或只有输出的代谢物）
        for node in self.metabolic_network.nodes():
            if self.metabolic_network.nodes[node].get('type') == 'metabolite':
                in_degree = self.metabolic_network.in_degree(node)
                out_degree = self.metabolic_network.out_degree(node)
                
                if in_degree == 0 and out_degree > 0:
                    gaps.append({
                        'type': 'source_metabolite',
                        'metabolite': node,
                        'description': f'{node} 只有输出，可能需要转运反应'
                    })
                elif out_degree == 0 and in_degree > 0:
                    gaps.append({
                        'type': 'sink_metabolite',
                        'metabolite': node,
                        'description': f'{node} 只有输入，可能是终产物或需要消耗反应'
                    })
        
        # 保存间隙分析结果
        gaps_file = self.output_dir / "metabolic_gaps.json"
        with open(gaps_file, 'w') as f:
            json.dump(gaps, f, indent=2)
        
        print(f"发现 {len(gaps)} 个潜在的代谢间隙")
        print(f"间隙分析结果保存到: {gaps_file}")
        
        return gaps
    
    def add_reactions(self, reaction_file="data/additional_reactions.txt"):
        """添加额外的反应"""
        print("正在添加额外反应...")
        
        if not Path(reaction_file).exists():
            print(f"反应文件 {reaction_file} 不存在")
            return
        
        # 这里应该实现从文件读取反应并添加到网络的逻辑
        print("额外反应添加功能需要进一步实现")
    
    def finalize_model(self):
        """生成最终的代谢模型"""
        print("正在生成最终代谢模型...")
        
        # 验证网络
        validation_results = self.validate_network()
        
        # 生成SBML格式的模型文件
        self._export_to_sbml()
        
        # 生成网络统计报告
        self._generate_network_statistics()
        
        print("最终模型生成完成")
    
    def _save_network(self):
        """保存网络到文件"""
        # 保存为GraphML格式
        graphml_file = self.output_dir / "metabolic_network.graphml"
        nx.write_graphml(self.metabolic_network, graphml_file)
        
        # 保存为JSON格式
        json_file = self.output_dir / "metabolic_network.json"
        network_data = nx.node_link_data(self.metabolic_network)
        with open(json_file, 'w') as f:
            json.dump(network_data, f, indent=2)
        
        print(f"网络已保存到: {graphml_file} 和 {json_file}")
    
    def _export_to_sbml(self):
        """导出为SBML格式"""
        # 创建简化的SBML文件
        root = ET.Element("sbml")
        model = ET.SubElement(root, "model", id="metabolic_model")
        
        # 添加代谢物
        species_list = ET.SubElement(model, "listOfSpecies")
        for node in self.metabolic_network.nodes():
            if self.metabolic_network.nodes[node].get('type') == 'metabolite':
                species = ET.SubElement(species_list, "species")
                species.set("id", node.replace('-', '_').replace(' ', '_'))
                species.set("name", node)
        
        # 添加反应
        reactions_list = ET.SubElement(model, "listOfReactions")
        for node in self.metabolic_network.nodes():
            if self.metabolic_network.nodes[node].get('type') == 'reaction':
                reaction = ET.SubElement(reactions_list, "reaction")
                reaction.set("id", node)
                reaction.set("name", self.metabolic_network.nodes[node].get('name', ''))
        
        # 保存SBML文件
        sbml_file = self.output_dir / "metabolic_model.xml"
        tree = ET.ElementTree(root)
        tree.write(sbml_file, encoding='utf-8', xml_declaration=True)
        
        print(f"SBML模型已保存到: {sbml_file}")
    
    def _generate_network_statistics(self):
        """生成网络统计信息"""
        stats = {}
        
        # 基本统计
        stats['nodes'] = self.metabolic_network.number_of_nodes()
        stats['edges'] = self.metabolic_network.number_of_edges()
        
        # 节点类型统计
        node_types = defaultdict(int)
        for node, data in self.metabolic_network.nodes(data=True):
            node_types[data.get('type', 'unknown')] += 1
        stats['node_types'] = dict(node_types)
        
        # 通路统计
        pathway_counts = defaultdict(int)
        for node, data in self.metabolic_network.nodes(data=True):
            if data.get('type') == 'reaction':
                pathway = data.get('pathway', 'unknown')
                pathway_counts[pathway] += 1
        stats['pathway_distribution'] = dict(pathway_counts)
        
        # 保存统计信息
        stats_file = self.output_dir / "network_statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        print(f"网络统计信息已保存到: {stats_file}")

def main():
    parser = argparse.ArgumentParser(description='代谢网络重建器')
    parser.add_argument('--auto-reconstruct', action='store_true',
                       help='自动重建代谢网络')
    parser.add_argument('--validate-network', action='store_true',
                       help='验证网络质量')
    parser.add_argument('--identify-gaps', action='store_true',
                       help='识别代谢间隙')
    parser.add_argument('--add-reactions', action='store_true',
                       help='添加额外反应')
    parser.add_argument('--finalize-model', action='store_true',
                       help='生成最终模型')
    parser.add_argument('--reaction-file', type=str,
                       help='额外反应文件路径')
    parser.add_argument('--output-dir', type=str, default='results',
                       help='输出目录')
    
    args = parser.parse_args()
    
    # 创建重建器实例
    reconstructor = MetabolicReconstructor(output_dir=args.output_dir)
    
    if args.auto_reconstruct:
        reconstructor.auto_reconstruct()
    
    if args.validate_network:
        reconstructor.validate_network()
    
    if args.identify_gaps:
        reconstructor.identify_gaps()
    
    if args.add_reactions:
        reaction_file = args.reaction_file or "data/additional_reactions.txt"
        reconstructor.add_reactions(reaction_file)
    
    if args.finalize_model:
        reconstructor.finalize_model()
    
    if not any([args.auto_reconstruct, args.validate_network, args.identify_gaps,
                args.add_reactions, args.finalize_model]):
        parser.print_help()

if __name__ == "__main__":
    main()