#!/usr/bin/env python3
"""
营养需求预测器
基于代谢网络预测微生物的营养需求和代谢能力

作者: 微生物基因组学课程组
版本: 1.0.0
"""

import argparse
import pandas as pd
import json
import networkx as nx
from pathlib import Path
from collections import defaultdict
import logging

class NutritionPredictor:
    def __init__(self, output_dir="results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # 设置日志
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # 加载代谢网络
        self.metabolic_network = None
        self._load_network()
        
        # 氨基酸合成通路数据
        self.amino_acid_pathways = self._load_amino_acid_pathways()
        
        # 维生素和辅因子合成数据
        self.vitamin_pathways = self._load_vitamin_pathways()
        
        # 碳源利用相关酶
        self.carbon_utilization_enzymes = self._load_carbon_enzymes()
    
    def _load_network(self):
        """加载代谢网络"""
        network_file = self.output_dir / "metabolic_network.json"
        if network_file.exists():
            with open(network_file, 'r') as f:
                network_data = json.load(f)
            self.metabolic_network = nx.node_link_graph(network_data)
            self.logger.info(f"已加载代谢网络: {self.metabolic_network.number_of_nodes()} 个节点")
        else:
            self.logger.warning("代谢网络文件不存在，请先进行网络重建")
    
    def _load_amino_acid_pathways(self):
        """加载氨基酸合成通路信息"""
        pathways = {
            # 丝氨酸家族
            'serine': {
                'enzymes': ['2.1.2.1', '4.3.1.17', '1.1.1.95'],
                'precursor': '3-phosphoglycerate',
                'family': 'serine_family'
            },
            'glycine': {
                'enzymes': ['2.1.2.1', '4.3.1.17'],
                'precursor': 'serine',
                'family': 'serine_family'
            },
            'cysteine': {
                'enzymes': ['2.5.1.47', '4.2.99.8'],
                'precursor': 'serine',
                'family': 'serine_family'
            },
            
            # 丙酮酸家族
            'alanine': {
                'enzymes': ['2.6.1.2'],
                'precursor': 'pyruvate',
                'family': 'pyruvate_family'
            },
            'valine': {
                'enzymes': ['2.2.1.6', '4.2.1.9', '2.6.1.42'],
                'precursor': 'pyruvate',
                'family': 'pyruvate_family'
            },
            'leucine': {
                'enzymes': ['2.3.3.13', '4.2.1.33', '2.6.1.6'],
                'precursor': 'pyruvate',
                'family': 'pyruvate_family'
            },
            'isoleucine': {
                'enzymes': ['2.3.3.13', '4.2.1.9', '2.6.1.42'],
                'precursor': 'threonine',
                'family': 'pyruvate_family'
            },
            
            # 天冬氨酸家族
            'aspartate': {
                'enzymes': ['2.6.1.1'],
                'precursor': 'oxaloacetate',
                'family': 'aspartate_family'
            },
            'asparagine': {
                'enzymes': ['6.3.5.4'],
                'precursor': 'aspartate',
                'family': 'aspartate_family'
            },
            'threonine': {
                'enzymes': ['2.7.2.4', '1.1.1.3', '4.3.1.19'],
                'precursor': 'aspartate',
                'family': 'aspartate_family'
            },
            'lysine': {
                'enzymes': ['2.3.1.117', '1.5.1.7', '1.4.1.16'],
                'precursor': 'aspartate',
                'family': 'aspartate_family'
            },
            'methionine': {
                'enzymes': ['2.5.1.6', '4.4.1.5'],
                'precursor': 'aspartate',
                'family': 'aspartate_family'
            },
            
            # 谷氨酸家族
            'glutamate': {
                'enzymes': ['1.4.1.3', '2.6.1.2'],
                'precursor': 'alpha-ketoglutarate',
                'family': 'glutamate_family'
            },
            'glutamine': {
                'enzymes': ['6.3.1.2'],
                'precursor': 'glutamate',
                'family': 'glutamate_family'
            },
            'proline': {
                'enzymes': ['2.7.2.11', '1.2.1.41', '1.5.1.2'],
                'precursor': 'glutamate',
                'family': 'glutamate_family'
            },
            'arginine': {
                'enzymes': ['2.1.3.3', '3.5.3.6', '4.3.2.1'],
                'precursor': 'glutamate',
                'family': 'glutamate_family'
            },
            
            # 芳香族氨基酸
            'phenylalanine': {
                'enzymes': ['4.2.1.51', '2.6.1.57'],
                'precursor': 'chorismate',
                'family': 'aromatic_family'
            },
            'tyrosine': {
                'enzymes': ['4.2.1.51', '2.6.1.57'],
                'precursor': 'chorismate',
                'family': 'aromatic_family'
            },
            'tryptophan': {
                'enzymes': ['4.1.3.27', '2.4.2.18', '4.2.1.20'],
                'precursor': 'chorismate',
                'family': 'aromatic_family'
            },
            
            # 组氨酸
            'histidine': {
                'enzymes': ['2.4.2.17', '3.5.4.12', '1.1.1.23'],
                'precursor': 'ribose-5-phosphate',
                'family': 'histidine_family'
            }
        }
        return pathways
    
    def _load_vitamin_pathways(self):
        """加载维生素合成通路信息"""
        pathways = {
            'thiamine': {  # 维生素B1
                'enzymes': ['2.5.1.3', '2.7.6.2'],
                'name': '硫胺素',
                'essential': True
            },
            'riboflavin': {  # 维生素B2
                'enzymes': ['2.5.1.9', '2.7.1.26'],
                'name': '核黄素',
                'essential': True
            },
            'pyridoxine': {  # 维生素B6
                'enzymes': ['4.3.3.1', '1.4.3.5'],
                'name': '吡哆醇',
                'essential': True
            },
            'cobalamin': {  # 维生素B12
                'enzymes': ['2.1.1.107', '6.3.5.10'],
                'name': '钴胺素',
                'essential': False
            },
            'biotin': {  # 生物素
                'enzymes': ['2.3.1.47', '6.3.3.3'],
                'name': '生物素',
                'essential': True
            },
            'folate': {  # 叶酸
                'enzymes': ['2.5.1.15', '6.3.2.12'],
                'name': '叶酸',
                'essential': True
            }
        }
        return pathways
    
    def _load_carbon_enzymes(self):
        """加载碳源利用相关酶"""
        enzymes = {
            'glucose': ['2.7.1.1', '2.7.1.2'],  # hexokinase, glucokinase
            'fructose': ['2.7.1.4'],  # fructokinase
            'galactose': ['2.7.1.6'],  # galactokinase
            'lactose': ['3.2.1.108'],  # beta-galactosidase
            'sucrose': ['3.2.1.48'],  # sucrase
            'maltose': ['3.2.1.20'],  # alpha-glucosidase
            'starch': ['3.2.1.1'],  # alpha-amylase
            'cellulose': ['3.2.1.4'],  # cellulase
            'xylose': ['5.3.1.5'],  # xylose isomerase
            'arabinose': ['5.3.1.4'],  # L-arabinose isomerase
            'glycerol': ['2.7.1.30'],  # glycerol kinase
            'acetate': ['6.2.1.1'],  # acetyl-CoA synthetase
        }
        return enzymes
    
    def predict_requirements(self, enzyme_file="results/enzyme_genes.txt"):
        """预测营养需求"""
        print("正在预测营养需求...")
        
        if not Path(enzyme_file).exists():
            print("酶基因文件不存在，请先运行酶基因提取")
            return
        
        # 读取酶基因
        enzyme_df = pd.read_csv(enzyme_file, sep='\\t')
        available_enzymes = set(enzyme_df['EC_Number'].unique())
        
        # 预测结果
        predictions = {
            'carbon_sources': self._predict_carbon_utilization(available_enzymes),
            'amino_acid_synthesis': self._predict_amino_acid_synthesis(available_enzymes),
            'vitamin_requirements': self._predict_vitamin_requirements(available_enzymes),
            'growth_factors': self._predict_growth_factors(available_enzymes)
        }
        
        # 生成预测报告
        report = self._generate_nutrition_report(predictions)
        
        # 保存结果
        prediction_file = self.output_dir / "nutrition_prediction.txt"
        with open(prediction_file, 'w') as f:
            f.write(report)
        
        # 保存详细数据
        detailed_file = self.output_dir / "nutrition_prediction_detailed.json"
        with open(detailed_file, 'w') as f:
            json.dump(predictions, f, indent=2)
        
        print(f"营养需求预测完成，结果保存到: {prediction_file}")
        
        return predictions
    
    def _predict_carbon_utilization(self, available_enzymes):
        """预测碳源利用能力"""
        carbon_utilization = {}
        
        for carbon_source, required_enzymes in self.carbon_utilization_enzymes.items():
            available_count = sum(1 for enzyme in required_enzymes if enzyme in available_enzymes)
            total_count = len(required_enzymes)
            
            if available_count == total_count:
                status = "可利用"
            elif available_count > 0:
                status = "部分利用"
            else:
                status = "不可利用"
            
            carbon_utilization[carbon_source] = {
                'status': status,
                'available_enzymes': available_count,
                'total_enzymes': total_count,
                'completeness': (available_count / total_count) * 100
            }
        
        return carbon_utilization
    
    def _predict_amino_acid_synthesis(self, available_enzymes):
        """预测氨基酸合成能力"""
        aa_synthesis = {}
        
        for aa, pathway_info in self.amino_acid_pathways.items():
            required_enzymes = pathway_info['enzymes']
            available_count = sum(1 for enzyme in required_enzymes if enzyme in available_enzymes)
            total_count = len(required_enzymes)
            
            if available_count == total_count:
                status = "可合成"
            elif available_count >= total_count * 0.7:
                status = "可能合成"
            else:
                status = "需要外源"
            
            aa_synthesis[aa] = {
                'status': status,
                'family': pathway_info['family'],
                'precursor': pathway_info['precursor'],
                'available_enzymes': available_count,
                'total_enzymes': total_count,
                'completeness': (available_count / total_count) * 100
            }
        
        return aa_synthesis
    
    def _predict_vitamin_requirements(self, available_enzymes):
        """预测维生素需求"""
        vitamin_requirements = {}
        
        for vitamin, pathway_info in self.vitamin_pathways.items():
            required_enzymes = pathway_info['enzymes']
            available_count = sum(1 for enzyme in required_enzymes if enzyme in available_enzymes)
            total_count = len(required_enzymes)
            
            if available_count == total_count:
                requirement = "可自主合成"
            elif available_count > 0:
                requirement = "部分合成能力"
            else:
                requirement = "需要外源补充"
            
            vitamin_requirements[vitamin] = {
                'name': pathway_info['name'],
                'requirement': requirement,
                'essential': pathway_info['essential'],
                'available_enzymes': available_count,
                'total_enzymes': total_count,
                'completeness': (available_count / total_count) * 100
            }
        
        return vitamin_requirements
    
    def _predict_growth_factors(self, available_enzymes):
        """预测其他生长因子需求"""
        # 简化的生长因子预测
        growth_factors = {
            'heme': {
                'enzymes': ['2.3.1.37', '4.99.1.1'],  # 血红素合成酶
                'requirement': 'unknown'
            },
            'coenzyme_a': {
                'enzymes': ['2.7.1.33', '6.3.2.5'],  # 辅酶A合成酶
                'requirement': 'unknown'
            }
        }
        
        for factor, info in growth_factors.items():
            required_enzymes = info['enzymes']
            available_count = sum(1 for enzyme in required_enzymes if enzyme in available_enzymes)
            total_count = len(required_enzymes)
            
            if available_count == total_count:
                requirement = "可自主合成"
            elif available_count > 0:
                requirement = "部分合成能力"
            else:
                requirement = "可能需要外源"
            
            growth_factors[factor]['requirement'] = requirement
            growth_factors[factor]['completeness'] = (available_count / total_count) * 100
        
        return growth_factors
    
    def _generate_nutrition_report(self, predictions):
        """生成营养需求报告"""
        report = "微生物营养需求预测报告\\n"
        report += "=" * 50 + "\\n\\n"
        
        # 碳源利用
        report += "1. 碳源利用能力\\n"
        report += "-" * 20 + "\\n"
        for carbon, info in predictions['carbon_sources'].items():
            report += f"{carbon}: {info['status']} ({info['completeness']:.1f}%)\\n"
        report += "\\n"
        
        # 氨基酸合成
        report += "2. 氨基酸合成能力\\n"
        report += "-" * 20 + "\\n"
        
        # 按家族分组
        families = defaultdict(list)
        for aa, info in predictions['amino_acid_synthesis'].items():
            families[info['family']].append((aa, info))
        
        for family, aa_list in families.items():
            report += f"\\n{family}:\\n"
            for aa, info in aa_list:
                report += f"  {aa}: {info['status']} ({info['completeness']:.1f}%)\\n"
        
        report += "\\n"
        
        # 维生素需求
        report += "3. 维生素需求\\n"
        report += "-" * 20 + "\\n"
        for vitamin, info in predictions['vitamin_requirements'].items():
            essential_mark = " *" if info['essential'] else ""
            report += f"{info['name']} ({vitamin}): {info['requirement']}{essential_mark}\\n"
        report += "\\n* 表示必需维生素\\n\\n"
        
        # 其他生长因子
        report += "4. 其他生长因子\\n"
        report += "-" * 20 + "\\n"
        for factor, info in predictions['growth_factors'].items():
            report += f"{factor}: {info['requirement']} ({info['completeness']:.1f}%)\\n"
        
        return report
    
    def amino_acid_synthesis(self, enzyme_file="results/enzyme_genes.txt"):
        """详细分析氨基酸合成通路"""
        print("正在分析氨基酸合成通路...")
        
        if not Path(enzyme_file).exists():
            print("酶基因文件不存在")
            return
        
        enzyme_df = pd.read_csv(enzyme_file, sep='\\t')
        available_enzymes = set(enzyme_df['EC_Number'].unique())
        
        # 详细分析每个氨基酸
        detailed_analysis = {}
        
        for aa, pathway_info in self.amino_acid_pathways.items():
            analysis = {
                'family': pathway_info['family'],
                'precursor': pathway_info['precursor'],
                'required_enzymes': pathway_info['enzymes'],
                'available_enzymes': [],
                'missing_enzymes': [],
                'synthesis_capability': 'unknown'
            }
            
            for enzyme in pathway_info['enzymes']:
                if enzyme in available_enzymes:
                    analysis['available_enzymes'].append(enzyme)
                else:
                    analysis['missing_enzymes'].append(enzyme)
            
            # 判断合成能力
            completeness = len(analysis['available_enzymes']) / len(pathway_info['enzymes'])
            if completeness == 1.0:
                analysis['synthesis_capability'] = 'complete'
            elif completeness >= 0.7:
                analysis['synthesis_capability'] = 'likely'
            elif completeness >= 0.3:
                analysis['synthesis_capability'] = 'partial'
            else:
                analysis['synthesis_capability'] = 'unlikely'
            
            detailed_analysis[aa] = analysis
        
        # 保存详细分析结果
        analysis_file = self.output_dir / "amino_acid_synthesis_analysis.json"
        with open(analysis_file, 'w') as f:
            json.dump(detailed_analysis, f, indent=2)
        
        print(f"氨基酸合成分析完成，结果保存到: {analysis_file}")
        
        return detailed_analysis
    
    def generate_aa_report(self):
        """生成氨基酸需求报告"""
        analysis_file = self.output_dir / "amino_acid_synthesis_analysis.json"
        
        if not analysis_file.exists():
            print("请先运行氨基酸合成分析")
            return
        
        with open(analysis_file, 'r') as f:
            analysis = json.load(f)
        
        # 生成报告
        report = "氨基酸合成能力详细报告\\n"
        report += "=" * 40 + "\\n\\n"
        
        # 按合成能力分类
        categories = {
            'complete': '完全合成',
            'likely': '可能合成',
            'partial': '部分合成',
            'unlikely': '需要外源'
        }
        
        for category, description in categories.items():
            aa_list = [aa for aa, info in analysis.items() 
                      if info['synthesis_capability'] == category]
            
            if aa_list:
                report += f"{description} ({len(aa_list)}个):\\n"
                for aa in aa_list:
                    info = analysis[aa]
                    completeness = len(info['available_enzymes']) / len(info['required_enzymes']) * 100
                    report += f"  {aa}: {completeness:.1f}% ({info['family']})\\n"
                report += "\\n"
        
        # 按家族统计
        report += "按氨基酸家族统计:\\n"
        report += "-" * 20 + "\\n"
        
        families = defaultdict(list)
        for aa, info in analysis.items():
            families[info['family']].append(aa)
        
        for family, aa_list in families.items():
            complete_count = sum(1 for aa in aa_list 
                               if analysis[aa]['synthesis_capability'] == 'complete')
            report += f"{family}: {complete_count}/{len(aa_list)} 完全合成\\n"
        
        # 保存报告
        report_file = self.output_dir / "amino_acid_synthesis_report.txt"
        with open(report_file, 'w') as f:
            f.write(report)
        
        print(f"氨基酸报告生成完成: {report_file}")
    
    def compare_metabolism(self, reference_file="data/reference_metabolism.json"):
        """与参考菌株比较代谢能力"""
        print("正在进行代谢能力比较...")
        
        if not Path(reference_file).exists():
            print(f"参考文件 {reference_file} 不存在")
            return
        
        # 这里应该实现与参考数据的比较逻辑
        print("代谢比较功能需要进一步实现")
    
    def comparison_report(self):
        """生成比较报告"""
        print("正在生成比较报告...")
        print("比较报告功能需要进一步实现")

def main():
    parser = argparse.ArgumentParser(description='营养需求预测器')
    parser.add_argument('--predict-requirements', action='store_true',
                       help='预测营养需求')
    parser.add_argument('--amino-acid-synthesis', action='store_true',
                       help='分析氨基酸合成能力')
    parser.add_argument('--generate-aa-report', action='store_true',
                       help='生成氨基酸报告')
    parser.add_argument('--compare-metabolism', action='store_true',
                       help='比较代谢能力')
    parser.add_argument('--comparison-report', action='store_true',
                       help='生成比较报告')
    parser.add_argument('--reference', type=str,
                       help='参考数据文件')
    parser.add_argument('--output-dir', type=str, default='results',
                       help='输出目录')
    
    args = parser.parse_args()
    
    # 创建预测器实例
    predictor = NutritionPredictor(output_dir=args.output_dir)
    
    if args.predict_requirements:
        predictor.predict_requirements()
    
    if args.amino_acid_synthesis:
        predictor.amino_acid_synthesis()
    
    if args.generate_aa_report:
        predictor.generate_aa_report()
    
    if args.compare_metabolism:
        reference_file = args.reference or "data/reference_metabolism.json"
        predictor.compare_metabolism(reference_file)
    
    if args.comparison_report:
        predictor.comparison_report()
    
    if not any([args.predict_requirements, args.amino_acid_synthesis,
                args.generate_aa_report, args.compare_metabolism,
                args.comparison_report]):
        parser.print_help()

if __name__ == "__main__":
    main()