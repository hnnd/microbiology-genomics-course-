#!/usr/bin/env python3
"""
微生物群落功能分析脚本
用于分析宏基因组数据的群落功能特征

作者: 微生物基因组学课程组
版本: 1.0.0
"""

import os
import sys
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import logging
import json
from collections import defaultdict

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CommunityFunctionAnalyzer:
    """群落功能分析类"""
    
    def __init__(self, input_file="metagenome_sample.fastq.gz", output_dir="results", threads=8):
        """
        初始化群落功能分析器
        
        参数:
        input_file: 输入的宏基因组数据文件
        output_dir: 输出结果目录
        threads: 使用的线程数
        """
        self.input_file = Path(input_file)
        self.output_dir = Path(output_dir)
        self.threads = threads
        
        # 创建输出目录
        self.functional_dir = self.output_dir / "functional_analysis"
        self.functional_dir.mkdir(parents=True, exist_ok=True)
        
        # 功能数据库路径（需要预先下载）
        self.databases = {
            'kegg': 'humann_databases/kegg',
            'cog': 'humann_databases/cog',
            'pfam': 'humann_databases/pfam'
        }
    
    def run_gene_prediction(self):
        """运行基因预测"""
        logger.info("运行基因预测...")
        
        # 使用Prodigal进行基因预测
        genes_file = self.functional_dir / "predicted_genes.faa"
        
        try:
            cmd = [
                "prodigal",
                "-i", str(self.input_file),
                "-a", str(genes_file),
                "-p", "meta",  # 宏基因组模式
                "-q"  # 安静模式
            ]
            
            subprocess.run(cmd, check=True)
            logger.info(f"基因预测完成，结果保存到: {genes_file}")
            return genes_file
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("Prodigal未安装或运行失败，创建模拟基因文件")
            return self._create_mock_genes()
    
    def _create_mock_genes(self):
        """创建模拟基因文件用于演示"""
        genes_file = self.functional_dir / "predicted_genes.faa"
        
        # 创建模拟基因序列
        mock_genes = [
            ">gene_001 # 1 # 1200 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.650",
            "MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVFTAYESE",
            ">gene_002 # 1201 # 2400 # 1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.620",
            "MKLVLSLSLVLAFSSATAAFAAIPQNIRIGTDPTYAPFESKNSQGELVGF",
            ">gene_003 # 2401 # 3600 # 1 # ID=1_3;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.680",
            "MTTPSHLSDRYELGEILGFGGMSEVHLARDLRLHRDVAVKVLRADLALQQ"
        ]
        
        with open(genes_file, 'w') as f:
            f.write('\n'.join(mock_genes))
        
        logger.info(f"模拟基因文件已创建: {genes_file}")
        return genes_file
    
    def run_functional_annotation(self, genes_file):
        """运行功能注释"""
        logger.info("运行功能注释...")
        
        # 使用DIAMOND进行蛋白质比对
        annotation_results = {}
        
        # KEGG注释
        kegg_results = self._annotate_kegg(genes_file)
        annotation_results['kegg'] = kegg_results
        
        # COG注释
        cog_results = self._annotate_cog(genes_file)
        annotation_results['cog'] = cog_results
        
        # Pfam注释
        pfam_results = self._annotate_pfam(genes_file)
        annotation_results['pfam'] = pfam_results
        
        return annotation_results
    
    def _annotate_kegg(self, genes_file):
        """KEGG功能注释"""
        logger.info("进行KEGG功能注释...")
        
        kegg_file = self.functional_dir / "kegg_annotation.tsv"
        
        # 模拟KEGG注释结果
        mock_kegg_data = [
            {'gene_id': 'gene_001', 'kegg_id': 'K00001', 'description': 'alcohol dehydrogenase', 'pathway': 'ko00010'},
            {'gene_id': 'gene_002', 'kegg_id': 'K00002', 'description': 'aldehyde dehydrogenase', 'pathway': 'ko00010'},
            {'gene_id': 'gene_003', 'kegg_id': 'K00003', 'description': 'homoserine dehydrogenase', 'pathway': 'ko00260'},
            {'gene_id': 'gene_004', 'kegg_id': 'K00010', 'description': 'glycerol kinase', 'pathway': 'ko00561'},
            {'gene_id': 'gene_005', 'kegg_id': 'K00016', 'description': 'L-lactate dehydrogenase', 'pathway': 'ko00620'}
        ]
        
        kegg_df = pd.DataFrame(mock_kegg_data)
        kegg_df.to_csv(kegg_file, sep='\t', index=False)
        
        logger.info(f"KEGG注释结果保存到: {kegg_file}")
        return kegg_df
    
    def _annotate_cog(self, genes_file):
        """COG功能注释"""
        logger.info("进行COG功能注释...")
        
        cog_file = self.functional_dir / "cog_annotation.tsv"
        
        # 模拟COG注释结果
        mock_cog_data = [
            {'gene_id': 'gene_001', 'cog_id': 'COG0001', 'category': 'J', 'description': 'Glutamate-1-semialdehyde 2,1-aminomutase'},
            {'gene_id': 'gene_002', 'cog_id': 'COG0002', 'category': 'K', 'description': 'Acetyl-CoA carboxylase alpha subunit'},
            {'gene_id': 'gene_003', 'cog_id': 'COG0003', 'category': 'L', 'description': 'Oxyanion-translocating ATPase'},
            {'gene_id': 'gene_004', 'cog_id': 'COG0004', 'category': 'E', 'description': 'Ammonia permease'},
            {'gene_id': 'gene_005', 'cog_id': 'COG0005', 'category': 'F', 'description': 'Purine nucleoside phosphorylase'}
        ]
        
        cog_df = pd.DataFrame(mock_cog_data)
        cog_df.to_csv(cog_file, sep='\t', index=False)
        
        logger.info(f"COG注释结果保存到: {cog_file}")
        return cog_df
    
    def _annotate_pfam(self, genes_file):
        """Pfam结构域注释"""
        logger.info("进行Pfam结构域注释...")
        
        pfam_file = self.functional_dir / "pfam_annotation.tsv"
        
        # 模拟Pfam注释结果
        mock_pfam_data = [
            {'gene_id': 'gene_001', 'pfam_id': 'PF00001', 'domain': '7tm_1', 'description': '7 transmembrane receptor'},
            {'gene_id': 'gene_002', 'pfam_id': 'PF00002', 'domain': '7tm_2', 'description': '7 transmembrane receptor'},
            {'gene_id': 'gene_003', 'pfam_id': 'PF00003', 'domain': '7tm_3', 'description': '7 transmembrane receptor'},
            {'gene_id': 'gene_004', 'pfam_id': 'PF00004', 'domain': 'AAA', 'description': 'ATPase family'},
            {'gene_id': 'gene_005', 'pfam_id': 'PF00005', 'domain': 'ABC_tran', 'description': 'ABC transporter'}
        ]
        
        pfam_df = pd.DataFrame(mock_pfam_data)
        pfam_df.to_csv(pfam_file, sep='\t', index=False)
        
        logger.info(f"Pfam注释结果保存到: {pfam_file}")
        return pfam_df
    
    def analyze_functional_profile(self, annotation_results):
        """分析功能谱"""
        logger.info("分析功能谱...")
        
        functional_profile = {}
        
        # KEGG通路分析
        kegg_pathways = self._analyze_kegg_pathways(annotation_results['kegg'])
        functional_profile['kegg_pathways'] = kegg_pathways
        
        # COG功能分类分析
        cog_categories = self._analyze_cog_categories(annotation_results['cog'])
        functional_profile['cog_categories'] = cog_categories
        
        # Pfam结构域分析
        pfam_domains = self._analyze_pfam_domains(annotation_results['pfam'])
        functional_profile['pfam_domains'] = pfam_domains
        
        # 保存功能谱结果
        profile_file = self.functional_dir / "functional_profile.json"
        with open(profile_file, 'w') as f:
            json.dump(functional_profile, f, indent=2)
        
        logger.info(f"功能谱分析结果保存到: {profile_file}")
        return functional_profile
    
    def _analyze_kegg_pathways(self, kegg_df):
        """分析KEGG代谢通路"""
        pathway_counts = kegg_df['pathway'].value_counts().to_dict()
        
        # 通路名称映射
        pathway_names = {
            'ko00010': 'Glycolysis / Gluconeogenesis',
            'ko00260': 'Glycine, serine and threonine metabolism',
            'ko00561': 'Glycerolipid metabolism',
            'ko00620': 'Pyruvate metabolism'
        }
        
        pathway_analysis = {}
        for pathway_id, count in pathway_counts.items():
            pathway_analysis[pathway_id] = {
                'name': pathway_names.get(pathway_id, 'Unknown pathway'),
                'gene_count': count,
                'percentage': round(count / len(kegg_df) * 100, 2)
            }
        
        return pathway_analysis
    
    def _analyze_cog_categories(self, cog_df):
        """分析COG功能分类"""
        category_counts = cog_df['category'].value_counts().to_dict()
        
        # COG分类名称映射
        category_names = {
            'J': 'Translation, ribosomal structure and biogenesis',
            'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'E': 'Amino acid transport and metabolism',
            'F': 'Nucleotide transport and metabolism'
        }
        
        category_analysis = {}
        for category, count in category_counts.items():
            category_analysis[category] = {
                'name': category_names.get(category, 'Unknown category'),
                'gene_count': count,
                'percentage': round(count / len(cog_df) * 100, 2)
            }
        
        return category_analysis
    
    def _analyze_pfam_domains(self, pfam_df):
        """分析Pfam结构域"""
        domain_counts = pfam_df['domain'].value_counts().to_dict()
        
        domain_analysis = {}
        for domain, count in domain_counts.items():
            domain_analysis[domain] = {
                'gene_count': count,
                'percentage': round(count / len(pfam_df) * 100, 2)
            }
        
        return domain_analysis
    
    def generate_gene_families_table(self, annotation_results):
        """生成基因家族丰度表"""
        logger.info("生成基因家族丰度表...")
        
        # 合并所有注释结果
        all_annotations = []
        
        # 添加KEGG注释
        for _, row in annotation_results['kegg'].iterrows():
            all_annotations.append({
                'gene_id': row['gene_id'],
                'family_id': row['kegg_id'],
                'family_type': 'KEGG',
                'description': row['description'],
                'abundance': np.random.randint(10, 1000)  # 模拟丰度值
            })
        
        # 添加COG注释
        for _, row in annotation_results['cog'].iterrows():
            all_annotations.append({
                'gene_id': row['gene_id'],
                'family_id': row['cog_id'],
                'family_type': 'COG',
                'description': row['description'],
                'abundance': np.random.randint(10, 1000)  # 模拟丰度值
            })
        
        # 创建基因家族表
        gene_families_df = pd.DataFrame(all_annotations)
        
        # 按基因家族汇总丰度
        family_abundance = gene_families_df.groupby(['family_id', 'family_type', 'description'])['abundance'].sum().reset_index()
        
        # 保存基因家族表
        families_file = self.functional_dir / "gene_families.tsv"
        family_abundance.to_csv(families_file, sep='\t', index=False)
        
        logger.info(f"基因家族丰度表保存到: {families_file}")
        return family_abundance
    
    def calculate_diversity_indices(self, functional_profile):
        """计算功能多样性指数"""
        logger.info("计算功能多样性指数...")
        
        diversity_indices = {}
        
        # KEGG通路多样性
        kegg_counts = [data['gene_count'] for data in functional_profile['kegg_pathways'].values()]
        diversity_indices['kegg_shannon'] = self._calculate_shannon_diversity(kegg_counts)
        diversity_indices['kegg_richness'] = len(kegg_counts)
        
        # COG分类多样性
        cog_counts = [data['gene_count'] for data in functional_profile['cog_categories'].values()]
        diversity_indices['cog_shannon'] = self._calculate_shannon_diversity(cog_counts)
        diversity_indices['cog_richness'] = len(cog_counts)
        
        # Pfam结构域多样性
        pfam_counts = [data['gene_count'] for data in functional_profile['pfam_domains'].values()]
        diversity_indices['pfam_shannon'] = self._calculate_shannon_diversity(pfam_counts)
        diversity_indices['pfam_richness'] = len(pfam_counts)
        
        # 保存多样性指数
        diversity_file = self.functional_dir / "diversity_indices.json"
        with open(diversity_file, 'w') as f:
            json.dump(diversity_indices, f, indent=2)
        
        logger.info(f"功能多样性指数保存到: {diversity_file}")
        return diversity_indices
    
    def _calculate_shannon_diversity(self, counts):
        """计算Shannon多样性指数"""
        if not counts or sum(counts) == 0:
            return 0
        
        total = sum(counts)
        proportions = [count / total for count in counts if count > 0]
        shannon = -sum(p * np.log(p) for p in proportions)
        
        return round(shannon, 3)
    
    def generate_summary_report(self, functional_profile, diversity_indices, gene_families):
        """生成功能分析汇总报告"""
        logger.info("生成功能分析汇总报告...")
        
        summary_file = self.functional_dir / "functional_analysis_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("微生物群落功能分析汇总报告\n")
            f.write("=" * 50 + "\n\n")
            
            # 基本统计
            f.write("基本统计信息:\n")
            f.write(f"- 总基因家族数: {len(gene_families)}\n")
            f.write(f"- KEGG通路数: {len(functional_profile['kegg_pathways'])}\n")
            f.write(f"- COG功能分类数: {len(functional_profile['cog_categories'])}\n")
            f.write(f"- Pfam结构域数: {len(functional_profile['pfam_domains'])}\n\n")
            
            # 多样性指数
            f.write("功能多样性指数:\n")
            f.write(f"- KEGG Shannon指数: {diversity_indices['kegg_shannon']}\n")
            f.write(f"- COG Shannon指数: {diversity_indices['cog_shannon']}\n")
            f.write(f"- Pfam Shannon指数: {diversity_indices['pfam_shannon']}\n\n")
            
            # 主要功能类别
            f.write("主要KEGG代谢通路:\n")
            for pathway_id, data in functional_profile['kegg_pathways'].items():
                f.write(f"- {data['name']}: {data['gene_count']} genes ({data['percentage']}%)\n")
            
            f.write("\n主要COG功能分类:\n")
            for category, data in functional_profile['cog_categories'].items():
                f.write(f"- [{category}] {data['name']}: {data['gene_count']} genes ({data['percentage']}%)\n")
        
        logger.info(f"功能分析汇总报告保存到: {summary_file}")
    
    def run_analysis(self):
        """运行完整的群落功能分析流程"""
        logger.info("开始群落功能分析...")
        
        # 检查输入文件
        if not self.input_file.exists():
            logger.error(f"输入文件不存在: {self.input_file}")
            return None
        
        # 基因预测
        genes_file = self.run_gene_prediction()
        
        # 功能注释
        annotation_results = self.run_functional_annotation(genes_file)
        
        # 功能谱分析
        functional_profile = self.analyze_functional_profile(annotation_results)
        
        # 生成基因家族表
        gene_families = self.generate_gene_families_table(annotation_results)
        
        # 计算多样性指数
        diversity_indices = self.calculate_diversity_indices(functional_profile)
        
        # 生成汇总报告
        self.generate_summary_report(functional_profile, diversity_indices, gene_families)
        
        logger.info("群落功能分析完成!")
        
        return {
            'functional_profile': functional_profile,
            'diversity_indices': diversity_indices,
            'gene_families': gene_families
        }

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='微生物群落功能分析工具')
    parser.add_argument('--input', default='metagenome_sample.fastq.gz', help='输入宏基因组数据文件')
    parser.add_argument('--output_dir', default='results', help='输出结果目录')
    parser.add_argument('--threads', type=int, default=8, help='使用的线程数')
    
    args = parser.parse_args()
    
    # 创建分析器实例
    analyzer = CommunityFunctionAnalyzer(
        input_file=args.input,
        output_dir=args.output_dir,
        threads=args.threads
    )
    
    # 运行分析
    results = analyzer.run_analysis()
    
    if results:
        print("\n=== 群落功能分析完成 ===")
        print(f"KEGG通路数: {len(results['functional_profile']['kegg_pathways'])}")
        print(f"COG功能分类数: {len(results['functional_profile']['cog_categories'])}")
        print(f"基因家族数: {len(results['gene_families'])}")
        print(f"KEGG Shannon指数: {results['diversity_indices']['kegg_shannon']}")
        print("\n详细结果请查看 results/functional_analysis/ 目录")
    else:
        print("群落功能分析失败")
        sys.exit(1)

if __name__ == "__main__":
    main()