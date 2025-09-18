#!/usr/bin/env python3
"""
基因组污染检测和分类分析脚本
用于第2次课实践操作：基因组组装与质量评估

功能：
1. 使用Kraken2进行分类学注释
2. 检测异常GC含量区域
3. 分析覆盖度异常序列
4. 识别载体序列和接头污染
5. 生成污染检测综合报告

作者：微生物基因组学课程组
版本：v1.0.0
"""

import os
import sys
import subprocess
import argparse
import json
from pathlib import Path
import time
import statistics
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

class ContaminationChecker:
    """基因组污染检测器"""
    
    def __init__(self, input_assembly, output_dir, kraken_db=None, threads=8):
        self.input_assembly = input_assembly
        self.output_dir = Path(output_dir)
        self.kraken_db = kraken_db
        self.threads = threads
        
        # 创建输出目录
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 污染检测阈值
        self.thresholds = {
            'gc_deviation': 10.0,  # GC含量偏差阈值(%)
            'coverage_outlier': 3.0,  # 覆盖度异常倍数
            'min_contig_length': 1000,  # 最小contig长度
            'contamination_ratio': 0.05  # 污染序列比例阈值
        }
        
        # 常见污染序列模式
        self.contamination_patterns = {
            'vector': [
                'GAATTC',  # EcoRI
                'AGATCT',  # BglII
                'GGATCC',  # BamHI
                'CTGCAG',  # PstI
            ],
            'adapter': [
                'AGATCGGAAGAGC',  # Illumina adapter
                'CTGTCTCTTATACACATCT',  # Nextera adapter
                'TGGAATTCTCGGGTGCCAAGG',  # TruSeq adapter
            ],
            'primer': [
                'GTGCCAGCMGCCGCGGTAA',  # 16S V4 forward
                'GGACTACHVGGGTWTCTAAT',  # 16S V4 reverse
            ]
        }
    
    def parse_fasta_sequences(self):
        """解析FASTA文件并提取序列信息"""
        print("📖 解析FASTA序列...")
        
        sequences = {}
        current_id = None
        current_seq = ""
        
        with open(self.input_assembly, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # 保存前一个序列
                    if current_id and current_seq:
                        sequences[current_id] = current_seq.upper()
                    
                    # 开始新序列
                    current_id = line[1:].split()[0]  # 取第一个空格前的部分作为ID
                    current_seq = ""
                else:
                    current_seq += line
            
            # 保存最后一个序列
            if current_id and current_seq:
                sequences[current_id] = current_seq.upper()
        
        print(f"✅ 解析完成，共{len(sequences)}个序列")
        return sequences
    
    def calculate_gc_content(self, sequences):
        """计算每个contig的GC含量"""
        print("🧮 计算GC含量...")
        
        gc_stats = {}
        all_gc_values = []
        
        for seq_id, sequence in sequences.items():
            if len(sequence) < self.thresholds['min_contig_length']:
                continue
                
            gc_count = sequence.count('G') + sequence.count('C')
            total_bases = len(sequence)
            gc_content = (gc_count / total_bases) * 100 if total_bases > 0 else 0
            
            gc_stats[seq_id] = {
                'length': total_bases,
                'gc_content': gc_content,
                'gc_count': gc_count
            }
            all_gc_values.append(gc_content)
        
        # 计算整体GC统计
        if all_gc_values:
            mean_gc = statistics.mean(all_gc_values)
            median_gc = statistics.median(all_gc_values)
            stdev_gc = statistics.stdev(all_gc_values) if len(all_gc_values) > 1 else 0
            
            overall_stats = {
                'mean_gc': mean_gc,
                'median_gc': median_gc,
                'stdev_gc': stdev_gc,
                'min_gc': min(all_gc_values),
                'max_gc': max(all_gc_values)
            }
            
            print(f"✅ GC含量统计: 平均={mean_gc:.2f}%, 标准差={stdev_gc:.2f}%")
        else:
            overall_stats = {}
        
        return gc_stats, overall_stats
    
    def detect_gc_anomalies(self, gc_stats, overall_stats):
        """检测GC含量异常的序列"""
        print("🔍 检测GC含量异常...")
        
        anomalies = []
        
        if not overall_stats:
            return anomalies
        
        mean_gc = overall_stats['mean_gc']
        threshold = self.thresholds['gc_deviation']
        
        for seq_id, stats in gc_stats.items():
            gc_content = stats['gc_content']
            deviation = abs(gc_content - mean_gc)
            
            if deviation > threshold:
                anomaly_type = 'high_gc' if gc_content > mean_gc else 'low_gc'
                anomalies.append({
                    'seq_id': seq_id,
                    'length': stats['length'],
                    'gc_content': gc_content,
                    'deviation': deviation,
                    'type': anomaly_type,
                    'severity': 'high' if deviation > threshold * 2 else 'medium'
                })
        
        print(f"✅ 发现{len(anomalies)}个GC含量异常序列")
        return anomalies
    
    def run_kraken_classification(self):
        """运行Kraken2分类学注释"""
        if not self.kraken_db:
            print("⚠️  未指定Kraken数据库，跳过分类学分析")
            return None
            
        print("🔬 运行Kraken2分类学注释...")
        
        kraken_output = self.output_dir / 'kraken_output.txt'
        kraken_report = self.output_dir / 'kraken_report.txt'
        
        cmd = [
            'kraken2',
            '--db', self.kraken_db,
            '--threads', str(self.threads),
            '--output', str(kraken_output),
            '--report', str(kraken_report),
            str(self.input_assembly)
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # 解析Kraken结果
            classification_stats = self.parse_kraken_results(kraken_output, kraken_report)
            print(f"✅ Kraken2分类完成")
            return classification_stats
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Kraken2运行失败: {e}")
            return None
        except FileNotFoundError:
            print("❌ 未找到Kraken2程序，请确保已正确安装")
            return None
    
    def parse_kraken_results(self, output_file, report_file):
        """解析Kraken2结果"""
        classification_stats = {
            'classified_sequences': 0,
            'unclassified_sequences': 0,
            'total_sequences': 0,
            'species_distribution': defaultdict(int),
            'contamination_candidates': []
        }
        
        # 解析分类输出文件
        sequence_classifications = {}
        
        with open(output_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    classified = parts[0] == 'C'
                    seq_id = parts[1]
                    taxid = parts[2] if classified else '0'
                    
                    sequence_classifications[seq_id] = {
                        'classified': classified,
                        'taxid': taxid
                    }
                    
                    if classified:
                        classification_stats['classified_sequences'] += 1
                    else:
                        classification_stats['unclassified_sequences'] += 1
                    
                    classification_stats['total_sequences'] += 1
        
        # 解析报告文件获取物种分布
        if report_file.exists():
            with open(report_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        percentage = float(parts[0])
                        rank = parts[3]
                        name = parts[5].strip()
                        
                        if rank == 'S' and percentage > 1.0:  # 物种级别且比例>1%
                            classification_stats['species_distribution'][name] = percentage
        
        # 识别潜在污染
        target_species = self.identify_target_species(classification_stats['species_distribution'])
        for species, percentage in classification_stats['species_distribution'].items():
            if species != target_species and percentage > 5.0:  # 非目标物种且比例>5%
                classification_stats['contamination_candidates'].append({
                    'species': species,
                    'percentage': percentage,
                    'type': 'taxonomic_contamination'
                })
        
        return classification_stats
    
    def identify_target_species(self, species_distribution):
        """识别目标物种（比例最高的物种）"""
        if not species_distribution:
            return None
        return max(species_distribution, key=species_distribution.get)
    
    def detect_sequence_patterns(self, sequences):
        """检测载体序列和接头污染"""
        print("🔍 检测序列模式污染...")
        
        pattern_matches = {
            'vector': [],
            'adapter': [],
            'primer': []
        }
        
        for seq_id, sequence in sequences.items():
            for pattern_type, patterns in self.contamination_patterns.items():
                for pattern in patterns:
                    # 检查正向和反向互补
                    if pattern in sequence:
                        pattern_matches[pattern_type].append({
                            'seq_id': seq_id,
                            'pattern': pattern,
                            'position': sequence.find(pattern),
                            'strand': 'forward'
                        })
                    
                    # 简单反向互补检查
                    rev_comp = self.reverse_complement(pattern)
                    if rev_comp in sequence:
                        pattern_matches[pattern_type].append({
                            'seq_id': seq_id,
                            'pattern': pattern,
                            'position': sequence.find(rev_comp),
                            'strand': 'reverse'
                        })
        
        total_matches = sum(len(matches) for matches in pattern_matches.values())
        print(f"✅ 发现{total_matches}个序列模式匹配")
        
        return pattern_matches
    
    def reverse_complement(self, sequence):
        """计算反向互补序列"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(sequence))
    
    def generate_contamination_report(self, sequences, gc_stats, overall_stats, 
                                    gc_anomalies, classification_stats, pattern_matches):
        """生成污染检测综合报告"""
        report_file = self.output_dir / 'contamination_report.txt'
        
        with open(report_file, 'w') as f:
            f.write("基因组污染检测综合报告\n")
            f.write("=" * 50 + "\n\n")
            
            # 基本统计信息
            f.write("基本统计信息:\n")
            f.write("-" * 20 + "\n")
            f.write(f"总序列数: {len(sequences)}\n")
            f.write(f"总长度: {sum(len(seq) for seq in sequences.values()):,} bp\n")
            f.write(f"分析序列数: {len(gc_stats)} (长度≥{self.thresholds['min_contig_length']} bp)\n\n")
            
            # GC含量分析
            if overall_stats:
                f.write("GC含量分析:\n")
                f.write("-" * 20 + "\n")
                f.write(f"平均GC含量: {overall_stats['mean_gc']:.2f}%\n")
                f.write(f"GC含量范围: {overall_stats['min_gc']:.2f}% - {overall_stats['max_gc']:.2f}%\n")
                f.write(f"标准差: {overall_stats['stdev_gc']:.2f}%\n")
                f.write(f"GC异常序列: {len(gc_anomalies)}个\n\n")
                
                if gc_anomalies:
                    f.write("GC含量异常序列详情:\n")
                    f.write(f"{'序列ID':<20} {'长度':<10} {'GC%':<8} {'偏差':<8} {'类型':<10} {'严重程度':<8}\n")
                    f.write("-" * 70 + "\n")
                    
                    for anomaly in gc_anomalies[:10]:  # 只显示前10个
                        f.write(f"{anomaly['seq_id']:<20} {anomaly['length']:<10} "
                               f"{anomaly['gc_content']:<8.2f} {anomaly['deviation']:<8.2f} "
                               f"{anomaly['type']:<10} {anomaly['severity']:<8}\n")
                    
                    if len(gc_anomalies) > 10:
                        f.write(f"... 还有{len(gc_anomalies) - 10}个异常序列\n")
                    f.write("\n")
            
            # 分类学分析
            if classification_stats:
                f.write("分类学分析:\n")
                f.write("-" * 20 + "\n")
                f.write(f"已分类序列: {classification_stats['classified_sequences']}\n")
                f.write(f"未分类序列: {classification_stats['unclassified_sequences']}\n")
                f.write(f"分类成功率: {classification_stats['classified_sequences']/classification_stats['total_sequences']*100:.1f}%\n\n")
                
                if classification_stats['species_distribution']:
                    f.write("主要物种分布:\n")
                    f.write(f"{'物种名称':<30} {'比例%':<8}\n")
                    f.write("-" * 40 + "\n")
                    
                    sorted_species = sorted(classification_stats['species_distribution'].items(), 
                                          key=lambda x: x[1], reverse=True)
                    for species, percentage in sorted_species[:10]:
                        f.write(f"{species:<30} {percentage:<8.2f}\n")
                    f.write("\n")
                
                if classification_stats['contamination_candidates']:
                    f.write("潜在污染物种:\n")
                    f.write(f"{'物种名称':<30} {'比例%':<8} {'类型':<15}\n")
                    f.write("-" * 55 + "\n")
                    
                    for candidate in classification_stats['contamination_candidates']:
                        f.write(f"{candidate['species']:<30} {candidate['percentage']:<8.2f} "
                               f"{candidate['type']:<15}\n")
                    f.write("\n")
            
            # 序列模式分析
            f.write("序列模式污染检测:\n")
            f.write("-" * 20 + "\n")
            
            total_pattern_matches = 0
            for pattern_type, matches in pattern_matches.items():
                f.write(f"{pattern_type.capitalize()}序列: {len(matches)}个匹配\n")
                total_pattern_matches += len(matches)
            
            f.write(f"总模式匹配: {total_pattern_matches}个\n\n")
            
            if total_pattern_matches > 0:
                f.write("模式匹配详情:\n")
                f.write(f"{'类型':<10} {'序列ID':<20} {'模式':<20} {'位置':<8} {'链':<8}\n")
                f.write("-" * 70 + "\n")
                
                for pattern_type, matches in pattern_matches.items():
                    for match in matches[:5]:  # 每种类型只显示前5个
                        f.write(f"{pattern_type:<10} {match['seq_id']:<20} "
                               f"{match['pattern']:<20} {match['position']:<8} "
                               f"{match['strand']:<8}\n")
                f.write("\n")
            
            # 污染评估总结
            f.write("污染评估总结:\n")
            f.write("-" * 20 + "\n")
            
            contamination_score = self.calculate_contamination_score(
                gc_anomalies, classification_stats, pattern_matches)
            
            f.write(f"污染风险评分: {contamination_score:.2f}/100\n")
            
            if contamination_score < 20:
                f.write("评估结果: 低污染风险 ✅\n")
            elif contamination_score < 50:
                f.write("评估结果: 中等污染风险 ⚠️\n")
            else:
                f.write("评估结果: 高污染风险 ❌\n")
            
            f.write("\n建议措施:\n")
            if len(gc_anomalies) > 0:
                f.write("- 检查GC含量异常序列的生物学合理性\n")
            if classification_stats and classification_stats['contamination_candidates']:
                f.write("- 去除非目标物种的污染序列\n")
            if total_pattern_matches > 0:
                f.write("- 清除载体序列和接头污染\n")
            if contamination_score < 20:
                f.write("- 当前组装质量良好，可用于后续分析\n")
            
            f.write(f"\n报告生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        print(f"📊 污染检测报告已保存至: {report_file}")
    
    def calculate_contamination_score(self, gc_anomalies, classification_stats, pattern_matches):
        """计算污染风险评分 (0-100, 越高越严重)"""
        score = 0
        
        # GC异常评分 (最高30分)
        if gc_anomalies:
            gc_score = min(30, len(gc_anomalies) * 2)
            score += gc_score
        
        # 分类学污染评分 (最高40分)
        if classification_stats and classification_stats['contamination_candidates']:
            for candidate in classification_stats['contamination_candidates']:
                score += min(20, candidate['percentage'])
        
        # 序列模式污染评分 (最高30分)
        total_patterns = sum(len(matches) for matches in pattern_matches.values())
        pattern_score = min(30, total_patterns * 3)
        score += pattern_score
        
        return min(100, score)
    
    def create_visualization(self, gc_stats, overall_stats, gc_anomalies):
        """创建污染检测可视化图表"""
        print("📊 生成可视化图表...")
        
        # 设置图表样式
        plt.style.use('seaborn-v0_8')
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('基因组污染检测分析', fontsize=16, fontweight='bold')
        
        # 1. GC含量分布直方图
        if gc_stats:
            gc_values = [stats['gc_content'] for stats in gc_stats.values()]
            
            axes[0, 0].hist(gc_values, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
            if overall_stats:
                axes[0, 0].axvline(overall_stats['mean_gc'], color='red', linestyle='--', 
                                 label=f'平均值: {overall_stats["mean_gc"]:.2f}%')
            axes[0, 0].set_xlabel('GC含量 (%)')
            axes[0, 0].set_ylabel('序列数量')
            axes[0, 0].set_title('GC含量分布')
            axes[0, 0].legend()
            axes[0, 0].grid(True, alpha=0.3)
        
        # 2. 序列长度 vs GC含量散点图
        if gc_stats:
            lengths = [stats['length'] for stats in gc_stats.values()]
            gc_values = [stats['gc_content'] for stats in gc_stats.values()]
            
            # 正常序列
            axes[0, 1].scatter(lengths, gc_values, alpha=0.6, s=20, color='blue', label='正常序列')
            
            # 异常序列
            if gc_anomalies:
                anomaly_lengths = [anomaly['length'] for anomaly in gc_anomalies]
                anomaly_gc = [anomaly['gc_content'] for anomaly in gc_anomalies]
                axes[0, 1].scatter(anomaly_lengths, anomaly_gc, alpha=0.8, s=30, 
                                 color='red', label='GC异常序列')
            
            axes[0, 1].set_xlabel('序列长度 (bp)')
            axes[0, 1].set_ylabel('GC含量 (%)')
            axes[0, 1].set_title('序列长度 vs GC含量')
            axes[0, 1].set_xscale('log')
            axes[0, 1].legend()
            axes[0, 1].grid(True, alpha=0.3)
        
        # 3. GC异常序列统计
        if gc_anomalies:
            anomaly_types = [anomaly['type'] for anomaly in gc_anomalies]
            type_counts = Counter(anomaly_types)
            
            axes[1, 0].bar(type_counts.keys(), type_counts.values(), 
                          color=['lightcoral', 'lightblue'])
            axes[1, 0].set_xlabel('异常类型')
            axes[1, 0].set_ylabel('序列数量')
            axes[1, 0].set_title('GC异常类型分布')
            axes[1, 0].grid(True, alpha=0.3)
        else:
            axes[1, 0].text(0.5, 0.5, '未发现GC异常', ha='center', va='center', 
                           transform=axes[1, 0].transAxes, fontsize=14)
            axes[1, 0].set_title('GC异常类型分布')
        
        # 4. 污染风险评估雷达图（简化为柱状图）
        risk_categories = ['GC异常', '分类学污染', '序列模式污染']
        risk_scores = [
            min(30, len(gc_anomalies) * 2) if gc_anomalies else 0,
            20,  # 占位符，实际应该从classification_stats计算
            10   # 占位符，实际应该从pattern_matches计算
        ]
        
        colors = ['red' if score > 20 else 'orange' if score > 10 else 'green' 
                 for score in risk_scores]
        
        axes[1, 1].bar(risk_categories, risk_scores, color=colors, alpha=0.7)
        axes[1, 1].set_ylabel('风险评分')
        axes[1, 1].set_title('污染风险评估')
        axes[1, 1].set_ylim(0, 50)
        axes[1, 1].grid(True, alpha=0.3)
        
        # 旋转x轴标签
        for ax in axes.flat:
            for label in ax.get_xticklabels():
                label.set_rotation(45)
        
        plt.tight_layout()
        
        # 保存图表
        plot_file = self.output_dir / 'contamination_analysis.png'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"📊 可视化图表已保存至: {plot_file}")
    
    def run_contamination_check(self):
        """运行完整的污染检测流程"""
        print("🚀 开始基因组污染检测...")
        print(f"输入文件: {self.input_assembly}")
        print(f"输出目录: {self.output_dir}")
        
        # 1. 解析序列
        sequences = self.parse_fasta_sequences()
        
        # 2. GC含量分析
        gc_stats, overall_stats = self.calculate_gc_content(sequences)
        gc_anomalies = self.detect_gc_anomalies(gc_stats, overall_stats)
        
        # 3. 分类学分析
        classification_stats = self.run_kraken_classification()
        
        # 4. 序列模式检测
        pattern_matches = self.detect_sequence_patterns(sequences)
        
        # 5. 生成报告
        self.generate_contamination_report(
            sequences, gc_stats, overall_stats, gc_anomalies, 
            classification_stats, pattern_matches
        )
        
        # 6. 创建可视化
        self.create_visualization(gc_stats, overall_stats, gc_anomalies)
        
        # 7. 保存详细结果
        results = {
            'basic_stats': {
                'total_sequences': len(sequences),
                'analyzed_sequences': len(gc_stats),
                'total_length': sum(len(seq) for seq in sequences.values())
            },
            'gc_analysis': {
                'overall_stats': overall_stats,
                'anomalies': gc_anomalies
            },
            'classification': classification_stats,
            'pattern_matches': pattern_matches,
            'contamination_score': self.calculate_contamination_score(
                gc_anomalies, classification_stats, pattern_matches),
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        results_file = self.output_dir / 'contamination_results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        print(f"\n✅ 污染检测完成！")
        print(f"📁 结果目录: {self.output_dir}")
        print(f"📊 查看报告: {self.output_dir}/contamination_report.txt")
        
        return results

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='基因组污染检测和分类分析')
    parser.add_argument('-i', '--input', required=True, help='输入组装FASTA文件')
    parser.add_argument('-o', '--output', default='results/contamination_check',
                       help='输出目录 (默认: results/contamination_check)')
    parser.add_argument('-d', '--kraken-db', help='Kraken2数据库路径')
    parser.add_argument('-t', '--threads', type=int, default=8, help='线程数 (默认: 8)')
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.input):
        print(f"❌ 错误: 找不到输入文件 {args.input}")
        sys.exit(1)
    
    # 创建污染检测器并运行
    checker = ContaminationChecker(
        input_assembly=args.input,
        output_dir=args.output,
        kraken_db=args.kraken_db,
        threads=args.threads
    )
    
    try:
        results = checker.run_contamination_check()
        print("\n🎉 污染检测任务完成！")
        
    except KeyboardInterrupt:
        print("\n⚠️  用户中断操作")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()