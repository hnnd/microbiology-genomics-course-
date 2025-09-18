#!/usr/bin/env python3
"""
基因组组装质量评估脚本
用于第2次课实践操作：基因组组装与质量评估

功能：
1. 运行BUSCO基因组完整性评估
2. 运行CheckM质量分级评估
3. 运行QUAST组装统计分析
4. 整合多种评估结果生成综合报告

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
import re

class QualityAssessor:
    """基因组组装质量评估器"""
    
    def __init__(self, input_assemblies, output_dir, threads=8):
        self.input_assemblies = input_assemblies  # 字典格式: {name: path}
        self.output_dir = Path(output_dir)
        self.threads = threads
        
        # 创建输出目录
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 评估工具配置
        self.tools_config = {
            'busco': {
                'lineage': 'bacteria_odb10',
                'mode': 'genome'
            },
            'checkm': {
                'rank': 'domain',
                'taxon': 'Bacteria'
            },
            'quast': {
                'min_contig': 500,
                'threads': threads
            }
        }
    
    def run_busco_assessment(self, assembly_name, assembly_path):
        """运行BUSCO基因组完整性评估"""
        print(f"\n🔍 运行BUSCO评估: {assembly_name}")
        
        busco_dir = self.output_dir / 'busco' / assembly_name
        busco_dir.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            'busco',
            '-i', str(assembly_path),
            '-o', assembly_name,
            '--out_path', str(busco_dir.parent),
            '-l', self.tools_config['busco']['lineage'],
            '-m', self.tools_config['busco']['mode'],
            '--cpu', str(self.threads),
            '--quiet'
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # 解析BUSCO结果
            summary_file = busco_dir / 'short_summary.specific.bacteria_odb10.txt'
            if summary_file.exists():
                busco_stats = self.parse_busco_results(summary_file)
                print(f"✅ BUSCO完成: {busco_stats['summary']}")
                return busco_stats
            else:
                print(f"⚠️  BUSCO结果文件未找到")
                return None
                
        except subprocess.CalledProcessError as e:
            print(f"❌ BUSCO运行失败: {e}")
            return None
    
    def parse_busco_results(self, summary_file):
        """解析BUSCO结果文件"""
        stats = {
            'complete': 0,
            'single_copy': 0,
            'duplicated': 0,
            'fragmented': 0,
            'missing': 0,
            'total': 0,
            'summary': ''
        }
        
        with open(summary_file, 'r') as f:
            content = f.read()
            
            # 提取统计信息
            patterns = {
                'complete': r'C:(\d+\.\d+)%',
                'single_copy': r'S:(\d+\.\d+)%',
                'duplicated': r'D:(\d+\.\d+)%',
                'fragmented': r'F:(\d+\.\d+)%',
                'missing': r'M:(\d+\.\d+)%',
                'total': r'n:(\d+)'
            }
            
            for key, pattern in patterns.items():
                match = re.search(pattern, content)
                if match:
                    if key == 'total':
                        stats[key] = int(match.group(1))
                    else:
                        stats[key] = float(match.group(1))
            
            # 提取汇总行
            summary_match = re.search(r'C:[\d.]+%.*?n:\d+', content)
            if summary_match:
                stats['summary'] = summary_match.group(0)
        
        return stats
    
    def run_checkm_assessment(self, assembly_name, assembly_path):
        """运行CheckM质量评估"""
        print(f"\n🔍 运行CheckM评估: {assembly_name}")
        
        checkm_dir = self.output_dir / 'checkm' / assembly_name
        checkm_dir.mkdir(parents=True, exist_ok=True)
        
        # 创建临时目录存放单个基因组
        temp_dir = checkm_dir / 'genomes'
        temp_dir.mkdir(exist_ok=True)
        
        # 复制基因组文件到临时目录
        temp_genome = temp_dir / f"{assembly_name}.fasta"
        subprocess.run(['cp', str(assembly_path), str(temp_genome)])
        
        # 运行CheckM lineage_wf
        cmd = [
            'checkm', 'lineage_wf',
            str(temp_dir),
            str(checkm_dir / 'output'),
            '-t', str(self.threads),
            '--pplacer_threads', str(min(self.threads, 4)),
            '-x', 'fasta',
            '--tab_table',
            '-f', str(checkm_dir / 'checkm_results.txt')
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # 解析CheckM结果
            results_file = checkm_dir / 'checkm_results.txt'
            if results_file.exists():
                checkm_stats = self.parse_checkm_results(results_file)
                print(f"✅ CheckM完成: 完整性={checkm_stats['completeness']:.1f}%, "
                      f"污染={checkm_stats['contamination']:.1f}%")
                return checkm_stats
            else:
                print(f"⚠️  CheckM结果文件未找到")
                return None
                
        except subprocess.CalledProcessError as e:
            print(f"❌ CheckM运行失败: {e}")
            return None
    
    def parse_checkm_results(self, results_file):
        """解析CheckM结果文件"""
        stats = {
            'completeness': 0,
            'contamination': 0,
            'strain_heterogeneity': 0,
            'quality_score': 0,
            'grade': 'Low'
        }
        
        with open(results_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:  # 跳过标题行
                data = lines[1].strip().split('\t')
                if len(data) >= 13:
                    stats['completeness'] = float(data[11])
                    stats['contamination'] = float(data[12])
                    if len(data) > 13:
                        stats['strain_heterogeneity'] = float(data[13])
                    
                    # 计算质量分数 (CheckM质量公式)
                    stats['quality_score'] = stats['completeness'] - 5 * stats['contamination']
                    
                    # 质量分级
                    if stats['completeness'] > 90 and stats['contamination'] < 5:
                        stats['grade'] = 'High'
                    elif stats['completeness'] > 50 and stats['contamination'] < 10:
                        stats['grade'] = 'Medium'
                    else:
                        stats['grade'] = 'Low'
        
        return stats
    
    def run_quast_assessment(self):
        """运行QUAST组装统计分析"""
        print(f"\n🔍 运行QUAST统计分析")
        
        quast_dir = self.output_dir / 'quast'
        quast_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建QUAST命令
        cmd = ['quast.py']
        
        # 添加所有组装文件
        for name, path in self.input_assemblies.items():
            cmd.extend(['-l', name, str(path)])
        
        cmd.extend([
            '-o', str(quast_dir),
            '--threads', str(self.threads),
            '--min-contig', str(self.tools_config['quast']['min_contig']),
            '--no-plots'  # 跳过图表生成以节省时间
        ])
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # 解析QUAST结果
            report_file = quast_dir / 'report.txt'
            if report_file.exists():
                quast_stats = self.parse_quast_results(report_file)
                print(f"✅ QUAST完成，分析了{len(self.input_assemblies)}个组装")
                return quast_stats
            else:
                print(f"⚠️  QUAST结果文件未找到")
                return None
                
        except subprocess.CalledProcessError as e:
            print(f"❌ QUAST运行失败: {e}")
            return None
    
    def parse_quast_results(self, report_file):
        """解析QUAST结果文件"""
        stats = {}
        
        with open(report_file, 'r') as f:
            lines = f.readlines()
            
            # 找到数据开始行
            data_start = -1
            for i, line in enumerate(lines):
                if line.startswith('Assembly'):
                    data_start = i
                    break
            
            if data_start == -1:
                return stats
            
            # 解析表头
            headers = lines[data_start].strip().split('\t')[1:]  # 跳过第一列"Assembly"
            
            # 解析每个指标行
            for line in lines[data_start + 1:]:
                if line.strip() == '':
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                    
                metric = parts[0]
                values = parts[1:]
                
                # 为每个组装创建统计字典
                for i, header in enumerate(headers):
                    if header not in stats:
                        stats[header] = {}
                    
                    if i < len(values):
                        # 尝试转换为数字
                        try:
                            if '.' in values[i]:
                                stats[header][metric] = float(values[i])
                            else:
                                stats[header][metric] = int(values[i])
                        except ValueError:
                            stats[header][metric] = values[i]
        
        return stats
    
    def generate_comprehensive_report(self, busco_results, checkm_results, quast_results):
        """生成综合质量评估报告"""
        report_file = self.output_dir / 'quality_assessment_report.txt'
        
        with open(report_file, 'w') as f:
            f.write("基因组组装质量综合评估报告\n")
            f.write("=" * 60 + "\n\n")
            
            # 写入评估概览
            f.write("评估概览:\n")
            f.write("-" * 30 + "\n")
            f.write(f"评估组装数量: {len(self.input_assemblies)}\n")
            f.write(f"评估工具: BUSCO, CheckM, QUAST\n")
            f.write(f"评估时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # BUSCO结果汇总
            if busco_results:
                f.write("BUSCO基因组完整性评估:\n")
                f.write("-" * 30 + "\n")
                f.write(f"{'组装名称':<20} {'完整性%':<10} {'单拷贝%':<10} {'重复%':<8} "
                       f"{'片段化%':<10} {'缺失%':<8}\n")
                f.write("-" * 70 + "\n")
                
                for name, stats in busco_results.items():
                    if stats:
                        f.write(f"{name:<20} {stats['complete']:<10.1f} "
                               f"{stats['single_copy']:<10.1f} {stats['duplicated']:<8.1f} "
                               f"{stats['fragmented']:<10.1f} {stats['missing']:<8.1f}\n")
                f.write("\n")
            
            # CheckM结果汇总
            if checkm_results:
                f.write("CheckM质量分级评估:\n")
                f.write("-" * 30 + "\n")
                f.write(f"{'组装名称':<20} {'完整性%':<10} {'污染%':<8} {'质量分数':<10} {'等级':<8}\n")
                f.write("-" * 60 + "\n")
                
                for name, stats in checkm_results.items():
                    if stats:
                        f.write(f"{name:<20} {stats['completeness']:<10.1f} "
                               f"{stats['contamination']:<8.1f} {stats['quality_score']:<10.1f} "
                               f"{stats['grade']:<8}\n")
                f.write("\n")
            
            # QUAST结果汇总
            if quast_results:
                f.write("QUAST组装统计分析:\n")
                f.write("-" * 30 + "\n")
                
                # 选择关键指标显示
                key_metrics = ['# contigs', 'Total length', 'N50', 'L50', 'GC (%)']
                
                f.write(f"{'组装名称':<20}")
                for metric in key_metrics:
                    f.write(f"{metric:<15}")
                f.write("\n")
                f.write("-" * (20 + 15 * len(key_metrics)) + "\n")
                
                for name in self.input_assemblies.keys():
                    if name in quast_results:
                        f.write(f"{name:<20}")
                        for metric in key_metrics:
                            value = quast_results[name].get(metric, 'N/A')
                            f.write(f"{str(value):<15}")
                        f.write("\n")
                f.write("\n")
            
            # 质量评估总结
            f.write("质量评估总结:\n")
            f.write("-" * 30 + "\n")
            
            # 找出最佳组装
            best_assembly = self.find_best_assembly(busco_results, checkm_results, quast_results)
            if best_assembly:
                f.write(f"推荐的最佳组装: {best_assembly}\n")
            
            # 质量建议
            f.write("\n质量改进建议:\n")
            f.write("1. BUSCO完整性应 >90% (高质量)\n")
            f.write("2. CheckM污染水平应 <5% (高质量)\n")
            f.write("3. N50值越大表示连续性越好\n")
            f.write("4. Contig数量越少表示组装越完整\n")
        
        print(f"📊 综合评估报告已保存至: {report_file}")
    
    def find_best_assembly(self, busco_results, checkm_results, quast_results):
        """根据多个指标找出最佳组装"""
        scores = {}
        
        for name in self.input_assemblies.keys():
            score = 0
            
            # BUSCO分数 (权重: 40%)
            if busco_results and name in busco_results and busco_results[name]:
                score += busco_results[name]['complete'] * 0.4
            
            # CheckM分数 (权重: 30%)
            if checkm_results and name in checkm_results and checkm_results[name]:
                completeness = checkm_results[name]['completeness']
                contamination = checkm_results[name]['contamination']
                checkm_score = max(0, completeness - 5 * contamination)
                score += checkm_score * 0.3
            
            # QUAST N50分数 (权重: 30%)
            if quast_results and name in quast_results:
                n50 = quast_results[name].get('N50', 0)
                if isinstance(n50, (int, float)) and n50 > 0:
                    # 标准化N50分数 (假设100kb为满分)
                    n50_score = min(100, n50 / 1000)  # 转换为相对分数
                    score += n50_score * 0.3
            
            scores[name] = score
        
        # 返回得分最高的组装
        if scores:
            return max(scores, key=scores.get)
        return None
    
    def run_comprehensive_assessment(self):
        """运行完整的质量评估流程"""
        print("🚀 开始基因组组装质量评估...")
        print(f"评估组装数量: {len(self.input_assemblies)}")
        print(f"输出目录: {self.output_dir}")
        
        # 运行各种评估
        busco_results = {}
        checkm_results = {}
        quast_results = None
        
        # BUSCO评估
        print("\n" + "="*50)
        print("第1步: BUSCO基因组完整性评估")
        print("="*50)
        
        for name, path in self.input_assemblies.items():
            result = self.run_busco_assessment(name, path)
            busco_results[name] = result
        
        # CheckM评估
        print("\n" + "="*50)
        print("第2步: CheckM质量分级评估")
        print("="*50)
        
        for name, path in self.input_assemblies.items():
            result = self.run_checkm_assessment(name, path)
            checkm_results[name] = result
        
        # QUAST评估
        print("\n" + "="*50)
        print("第3步: QUAST组装统计分析")
        print("="*50)
        
        quast_results = self.run_quast_assessment()
        
        # 生成综合报告
        print("\n" + "="*50)
        print("第4步: 生成综合评估报告")
        print("="*50)
        
        self.generate_comprehensive_report(busco_results, checkm_results, quast_results)
        
        # 保存详细结果
        results = {
            'busco': busco_results,
            'checkm': checkm_results,
            'quast': quast_results,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        results_file = self.output_dir / 'detailed_results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        print(f"\n✅ 质量评估完成！")
        print(f"📁 结果目录: {self.output_dir}")
        print(f"📊 查看报告: {self.output_dir}/quality_assessment_report.txt")
        
        return results

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='基因组组装质量综合评估')
    parser.add_argument('-i', '--input', required=True, nargs='+', 
                       help='输入组装文件 (格式: name:path)')
    parser.add_argument('-o', '--output', default='results/quality_assessment',
                       help='输出目录 (默认: results/quality_assessment)')
    parser.add_argument('-t', '--threads', type=int, default=8, help='线程数 (默认: 8)')
    
    args = parser.parse_args()
    
    # 解析输入文件
    assemblies = {}
    for item in args.input:
        if ':' in item:
            name, path = item.split(':', 1)
        else:
            name = Path(item).stem
            path = item
        
        if not os.path.exists(path):
            print(f"❌ 错误: 找不到输入文件 {path}")
            sys.exit(1)
        
        assemblies[name] = path
    
    print(f"将评估以下组装:")
    for name, path in assemblies.items():
        print(f"  {name}: {path}")
    
    # 创建评估器并运行
    assessor = QualityAssessor(
        input_assemblies=assemblies,
        output_dir=args.output,
        threads=args.threads
    )
    
    try:
        results = assessor.run_comprehensive_assessment()
        print("\n🎉 所有评估任务完成！")
        
    except KeyboardInterrupt:
        print("\n⚠️  用户中断操作")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()