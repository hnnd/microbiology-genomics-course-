#!/usr/bin/env python3
"""
KEGG通路分析器
用于分析微生物基因组的KEGG通路完整性和酶基因分布

作者: 微生物基因组学课程组
版本: 1.0.0
"""

import argparse
import requests
import time
import pandas as pd
import json
from pathlib import Path
from Bio import SeqIO
import re
import sys

class KEGGPathwayAnalyzer:
    def __init__(self, output_dir="results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # KEGG API基础URL
        self.kegg_base_url = "https://rest.kegg.jp"
        
        # 核心代谢通路列表
        self.core_pathways = {
            'glycolysis': 'map00010',
            'tca_cycle': 'map00020', 
            'pentose_phosphate': 'map00030',
            'amino_acid_synthesis': ['map00260', 'map00270', 'map00280', 'map00290'],
            'nucleotide_synthesis': ['map00230', 'map00240'],
            'fatty_acid_synthesis': 'map00061'
        }
        
        # 延时设置（避免API限制）
        self.api_delay = 1.0
        
    def extract_enzymes(self, annotation_file="ecoli_features.txt"):
        """从基因组注释文件中提取酶基因信息"""
        print("正在提取酶基因信息...")
        
        enzyme_genes = []
        
        try:
            # 读取注释文件
            df = pd.read_csv(annotation_file, sep='\t', low_memory=False)
            
            # 查找包含EC编号的基因
            for idx, row in df.iterrows():
                if pd.notna(row.get('product', '')):
                    product = str(row['product'])
                    
                    # 使用正则表达式查找EC编号
                    ec_pattern = r'EC[:\s]?(\d+\.\d+\.\d+\.\d+)'
                    ec_matches = re.findall(ec_pattern, product, re.IGNORECASE)
                    
                    if ec_matches:
                        for ec_number in ec_matches:
                            enzyme_genes.append({
                                'Gene_ID': row.get('locus_tag', ''),
                                'EC_Number': ec_number,
                                'Gene_Name': row.get('gene', ''),
                                'Function': product
                            })
            
            # 保存结果
            enzyme_df = pd.DataFrame(enzyme_genes)
            output_file = self.output_dir / "enzyme_genes.txt"
            enzyme_df.to_csv(output_file, sep='\t', index=False)
            
            print(f"提取到 {len(enzyme_genes)} 个酶基因")
            print(f"结果保存到: {output_file}")
            
            return enzyme_df
            
        except Exception as e:
            print(f"提取酶基因时出错: {e}")
            return pd.DataFrame()
    
    def get_kegg_pathway_info(self, pathway_id):
        """获取KEGG通路信息"""
        try:
            url = f"{self.kegg_base_url}/get/{pathway_id}"
            response = requests.get(url)
            time.sleep(self.api_delay)
            
            if response.status_code == 200:
                return response.text
            else:
                print(f"无法获取通路 {pathway_id} 的信息")
                return None
                
        except Exception as e:
            print(f"获取KEGG通路信息时出错: {e}")
            return None
    
    def pathway_mapping(self, enzyme_file="results/enzyme_genes.txt"):
        """将酶基因映射到KEGG通路"""
        print("正在进行KEGG通路映射...")
        
        # 创建输出目录
        pathway_dir = self.output_dir / "pathway_analysis"
        pathway_dir.mkdir(exist_ok=True)
        
        try:
            # 读取酶基因文件
            if not Path(enzyme_file).exists():
                print("酶基因文件不存在，请先运行 --extract-enzymes")
                return
                
            enzyme_df = pd.read_csv(enzyme_file, sep='\t')
            ec_numbers = enzyme_df['EC_Number'].unique()
            
            # 通路完整性分析
            pathway_completeness = {}
            
            for pathway_name, pathway_id in self.core_pathways.items():
                if isinstance(pathway_id, list):
                    # 处理多个通路ID的情况
                    total_enzymes = 0
                    found_enzymes = 0
                    
                    for pid in pathway_id:
                        pathway_info = self.get_kegg_pathway_info(pid)
                        if pathway_info:
                            # 简化的通路分析（实际应用中需要更复杂的解析）
                            pathway_enzymes = self._extract_pathway_enzymes(pathway_info)
                            total_enzymes += len(pathway_enzymes)
                            found_enzymes += len([ec for ec in pathway_enzymes if ec in ec_numbers])
                else:
                    pathway_info = self.get_kegg_pathway_info(pathway_id)
                    if pathway_info:
                        pathway_enzymes = self._extract_pathway_enzymes(pathway_info)
                        total_enzymes = len(pathway_enzymes)
                        found_enzymes = len([ec for ec in pathway_enzymes if ec in ec_numbers])
                    else:
                        total_enzymes = 0
                        found_enzymes = 0
                
                if total_enzymes > 0:
                    completeness = (found_enzymes / total_enzymes) * 100
                else:
                    completeness = 0
                    
                pathway_completeness[pathway_name] = {
                    'total_enzymes': total_enzymes,
                    'found_enzymes': found_enzymes,
                    'completeness_percent': completeness
                }
            
            # 保存通路完整性结果
            completeness_df = pd.DataFrame.from_dict(pathway_completeness, orient='index')
            completeness_file = pathway_dir / "pathway_completeness.txt"
            completeness_df.to_csv(completeness_file, sep='\t')
            
            # 生成缺失酶基因列表（示例）
            missing_enzymes = self._generate_missing_enzymes_list(ec_numbers)
            missing_file = pathway_dir / "missing_enzymes.txt"
            with open(missing_file, 'w') as f:
                f.write("\\n".join(missing_enzymes))
            
            # 生成通路覆盖度报告
            coverage_report = self._generate_coverage_report(pathway_completeness)
            coverage_file = pathway_dir / "pathway_coverage.txt"
            with open(coverage_file, 'w') as f:
                f.write(coverage_report)
            
            print(f"通路映射完成，结果保存到: {pathway_dir}")
            
        except Exception as e:
            print(f"通路映射时出错: {e}")
    
    def _extract_pathway_enzymes(self, pathway_info):
        """从通路信息中提取EC编号（简化版本）"""
        # 这是一个简化的实现，实际应用中需要更复杂的解析
        ec_pattern = r'(\d+\.\d+\.\d+\.\d+)'
        ec_numbers = re.findall(ec_pattern, pathway_info)
        return list(set(ec_numbers))
    
    def _generate_missing_enzymes_list(self, found_ec_numbers):
        """生成缺失酶基因列表（示例）"""
        # 这里是一些常见的核心代谢酶
        essential_enzymes = [
            '2.7.1.1',   # hexokinase
            '5.3.1.9',   # glucose-6-phosphate isomerase
            '2.7.1.11',  # 6-phosphofructokinase
            '4.1.2.13',  # fructose-bisphosphate aldolase
            '1.2.1.12',  # glyceraldehyde-3-phosphate dehydrogenase
            '2.7.2.3',   # phosphoglycerate kinase
            '5.4.2.11',  # phosphoglycerate mutase
            '4.2.1.11',  # enolase
            '2.7.1.40',  # pyruvate kinase
        ]
        
        missing = [ec for ec in essential_enzymes if ec not in found_ec_numbers]
        return missing
    
    def _generate_coverage_report(self, pathway_completeness):
        """生成通路覆盖度报告"""
        report = "KEGG通路覆盖度分析报告\\n"
        report += "=" * 50 + "\\n\\n"
        
        for pathway, data in pathway_completeness.items():
            report += f"通路: {pathway}\\n"
            report += f"  总酶数: {data['total_enzymes']}\\n"
            report += f"  发现酶数: {data['found_enzymes']}\\n"
            report += f"  完整性: {data['completeness_percent']:.1f}%\\n\\n"
        
        return report
    
    def assess_completeness(self, pathways=None):
        """评估指定通路的完整性"""
        if pathways is None:
            pathways = ['glycolysis', 'tca_cycle', 'amino_acid_synthesis']
        
        print(f"正在评估通路完整性: {', '.join(pathways)}")
        
        # 读取已有的通路分析结果
        completeness_file = self.output_dir / "pathway_analysis" / "pathway_completeness.txt"
        
        if completeness_file.exists():
            df = pd.read_csv(completeness_file, sep='\\t', index_col=0)
            
            report = "通路完整性评估报告\\n"
            report += "=" * 40 + "\\n\\n"
            
            for pathway in pathways:
                if pathway in df.index:
                    completeness = df.loc[pathway, 'completeness_percent']
                    status = self._get_completeness_status(completeness)
                    
                    report += f"{pathway}:\\n"
                    report += f"  完整性: {completeness:.1f}%\\n"
                    report += f"  状态: {status}\\n\\n"
                else:
                    report += f"{pathway}: 数据不可用\\n\\n"
            
            # 保存报告
            report_file = self.output_dir / "pathway_completeness_report.txt"
            with open(report_file, 'w') as f:
                f.write(report)
            
            print(f"完整性评估完成，报告保存到: {report_file}")
        else:
            print("请先运行通路映射分析 (--pathway-mapping)")
    
    def _get_completeness_status(self, completeness):
        """根据完整性百分比返回状态"""
        if completeness >= 90:
            return "完整"
        elif completeness >= 70:
            return "基本完整"
        elif completeness >= 50:
            return "部分完整"
        else:
            return "不完整"

def main():
    parser = argparse.ArgumentParser(description='KEGG通路分析器')
    parser.add_argument('--extract-enzymes', action='store_true',
                       help='从注释文件中提取酶基因')
    parser.add_argument('--pathway-mapping', action='store_true',
                       help='进行KEGG通路映射分析')
    parser.add_argument('--assess-completeness', action='store_true',
                       help='评估通路完整性')
    parser.add_argument('--pathways', type=str,
                       help='指定要分析的通路（逗号分隔）')
    parser.add_argument('--slow-mode', action='store_true',
                       help='使用慢速模式（增加API延时）')
    parser.add_argument('--output-dir', type=str, default='results',
                       help='输出目录')
    
    args = parser.parse_args()
    
    # 创建分析器实例
    analyzer = KEGGPathwayAnalyzer(output_dir=args.output_dir)
    
    if args.slow_mode:
        analyzer.api_delay = 3.0
        print("使用慢速模式，API延时增加到3秒")
    
    if args.extract_enzymes:
        analyzer.extract_enzymes()
    
    if args.pathway_mapping:
        analyzer.pathway_mapping()
    
    if args.assess_completeness:
        pathways = None
        if args.pathways:
            pathways = [p.strip() for p in args.pathways.split(',')]
        analyzer.assess_completeness(pathways)
    
    if not any([args.extract_enzymes, args.pathway_mapping, args.assess_completeness]):
        parser.print_help()

if __name__ == "__main__":
    main()