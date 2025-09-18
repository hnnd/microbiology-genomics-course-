#!/usr/bin/env python3
"""
基因组组装参数优化脚本
用于第2次课实践操作：基因组组装与质量评估

功能：
1. 自动化SPAdes参数优化
2. 批量测试不同k-mer组合
3. 比较不同参数设置的组装效果
4. 生成优化建议报告

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

class AssemblyOptimizer:
    """基因组组装参数优化器"""
    
    def __init__(self, input_r1, input_r2, output_dir, threads=8, memory=16):
        self.input_r1 = input_r1
        self.input_r2 = input_r2
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.memory = memory
        
        # 创建输出目录
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 参数组合配置
        self.parameter_sets = {
            'standard': {
                'k_values': '21,33,55,77',
                'careful': False,
                'cov_cutoff': 'off',
                'description': '标准参数组装'
            },
            'optimized': {
                'k_values': '21,33,55,77,99',
                'careful': True,
                'cov_cutoff': 'auto',
                'description': '优化参数组装'
            },
            'high_quality': {
                'k_values': '21,33,55,77,99,127',
                'careful': True,
                'cov_cutoff': 'auto',
                'description': '高质量参数组装'
            },
            'fast': {
                'k_values': '21,55,99',
                'careful': False,
                'cov_cutoff': 'auto',
                'description': '快速组装参数'
            }
        }
        
    def run_spades_assembly(self, param_name, params):
        """运行SPAdes组装"""
        print(f"\n开始运行 {param_name} 参数组装...")
        print(f"描述: {params['description']}")
        
        # 构建输出目录
        assembly_dir = self.output_dir / f"spades_{param_name}"
        
        # 构建SPAdes命令
        cmd = [
            'spades.py',
            '-1', str(self.input_r1),
            '-2', str(self.input_r2),
            '-o', str(assembly_dir),
            '-k', params['k_values'],
            '--threads', str(self.threads),
            '--memory', str(self.memory)
        ]
        
        # 添加可选参数
        if params['careful']:
            cmd.append('--careful')
            
        if params['cov_cutoff'] != 'off':
            cmd.extend(['--cov-cutoff', params['cov_cutoff']])
        
        print(f"执行命令: {' '.join(cmd)}")
        
        # 记录开始时间
        start_time = time.time()
        
        try:
            # 运行SPAdes
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # 记录结束时间
            end_time = time.time()
            runtime = end_time - start_time
            
            print(f"✅ {param_name} 组装完成，用时: {runtime:.2f}秒")
            
            return {
                'success': True,
                'runtime': runtime,
                'output_dir': str(assembly_dir),
                'stdout': result.stdout,
                'stderr': result.stderr
            }
            
        except subprocess.CalledProcessError as e:
            print(f"❌ {param_name} 组装失败: {e}")
            return {
                'success': False,
                'error': str(e),
                'stdout': e.stdout,
                'stderr': e.stderr
            }
    
    def analyze_assembly_stats(self, assembly_dir):
        """分析组装统计信息"""
        contigs_file = Path(assembly_dir) / 'contigs.fasta'
        
        if not contigs_file.exists():
            return None
            
        stats = {
            'total_contigs': 0,
            'total_length': 0,
            'longest_contig': 0,
            'n50': 0,
            'l50': 0,
            'gc_content': 0
        }
        
        # 读取contigs文件并计算统计信息
        contig_lengths = []
        total_bases = 0
        gc_count = 0
        
        with open(contigs_file, 'r') as f:
            current_seq = ""
            for line in f:
                if line.startswith('>'):
                    if current_seq:
                        length = len(current_seq)
                        contig_lengths.append(length)
                        total_bases += length
                        gc_count += current_seq.count('G') + current_seq.count('C')
                    current_seq = ""
                else:
                    current_seq += line.strip().upper()
            
            # 处理最后一个序列
            if current_seq:
                length = len(current_seq)
                contig_lengths.append(length)
                total_bases += length
                gc_count += current_seq.count('G') + current_seq.count('C')
        
        # 计算统计信息
        stats['total_contigs'] = len(contig_lengths)
        stats['total_length'] = total_bases
        stats['longest_contig'] = max(contig_lengths) if contig_lengths else 0
        stats['gc_content'] = (gc_count / total_bases * 100) if total_bases > 0 else 0
        
        # 计算N50和L50
        contig_lengths.sort(reverse=True)
        cumulative_length = 0
        half_genome = total_bases / 2
        
        for i, length in enumerate(contig_lengths):
            cumulative_length += length
            if cumulative_length >= half_genome:
                stats['n50'] = length
                stats['l50'] = i + 1
                break
        
        return stats
    
    def generate_comparison_report(self, results):
        """生成参数比较报告"""
        report_file = self.output_dir / 'optimization_report.txt'
        
        with open(report_file, 'w') as f:
            f.write("SPAdes参数优化比较报告\n")
            f.write("=" * 50 + "\n\n")
            
            # 写入参数设置
            f.write("参数设置:\n")
            f.write("-" * 20 + "\n")
            for param_name, params in self.parameter_sets.items():
                f.write(f"{param_name}:\n")
                f.write(f"  描述: {params['description']}\n")
                f.write(f"  k-mer值: {params['k_values']}\n")
                f.write(f"  careful模式: {params['careful']}\n")
                f.write(f"  覆盖度阈值: {params['cov_cutoff']}\n\n")
            
            # 写入组装结果比较
            f.write("组装结果比较:\n")
            f.write("-" * 20 + "\n")
            f.write(f"{'参数组':<15} {'成功':<8} {'运行时间(s)':<12} {'Contigs数':<10} "
                   f"{'总长度(bp)':<12} {'N50':<10} {'L50':<8} {'GC%':<8}\n")
            f.write("-" * 85 + "\n")
            
            for param_name, result in results.items():
                if result['success'] and result['stats']:
                    stats = result['stats']
                    f.write(f"{param_name:<15} {'是':<8} {result['runtime']:<12.1f} "
                           f"{stats['total_contigs']:<10} {stats['total_length']:<12} "
                           f"{stats['n50']:<10} {stats['l50']:<8} {stats['gc_content']:<8.2f}\n")
                else:
                    f.write(f"{param_name:<15} {'否':<8} {'N/A':<12} {'N/A':<10} "
                           f"{'N/A':<12} {'N/A':<10} {'N/A':<8} {'N/A':<8}\n")
            
            # 写入优化建议
            f.write("\n优化建议:\n")
            f.write("-" * 20 + "\n")
            
            # 找出最佳参数组合
            best_n50 = 0
            best_param = None
            
            for param_name, result in results.items():
                if result['success'] and result['stats']:
                    if result['stats']['n50'] > best_n50:
                        best_n50 = result['stats']['n50']
                        best_param = param_name
            
            if best_param:
                f.write(f"1. 最佳N50结果: {best_param} (N50={best_n50})\n")
                f.write(f"2. 推荐用于后续分析的组装结果: {best_param}\n")
                f.write(f"3. 如需平衡质量和速度，建议使用: optimized参数\n")
            
            f.write(f"\n报告生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        print(f"📊 优化报告已保存至: {report_file}")
    
    def run_optimization(self):
        """运行完整的参数优化流程"""
        print("🚀 开始SPAdes参数优化...")
        print(f"输入文件: {self.input_r1}, {self.input_r2}")
        print(f"输出目录: {self.output_dir}")
        print(f"线程数: {self.threads}, 内存: {self.memory}GB")
        
        results = {}
        
        # 运行所有参数组合
        for param_name, params in self.parameter_sets.items():
            result = self.run_spades_assembly(param_name, params)
            
            # 分析组装统计信息
            if result['success']:
                stats = self.analyze_assembly_stats(result['output_dir'])
                result['stats'] = stats
            else:
                result['stats'] = None
            
            results[param_name] = result
        
        # 生成比较报告
        self.generate_comparison_report(results)
        
        # 保存详细结果
        results_file = self.output_dir / 'optimization_results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        print(f"\n✅ 参数优化完成！")
        print(f"📁 结果目录: {self.output_dir}")
        print(f"📊 查看报告: {self.output_dir}/optimization_report.txt")
        
        return results

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='SPAdes基因组组装参数优化')
    parser.add_argument('-1', '--input-r1', required=True, help='输入R1 FASTQ文件')
    parser.add_argument('-2', '--input-r2', required=True, help='输入R2 FASTQ文件')
    parser.add_argument('-o', '--output', default='results/spades_optimization', 
                       help='输出目录 (默认: results/spades_optimization)')
    parser.add_argument('-t', '--threads', type=int, default=8, help='线程数 (默认: 8)')
    parser.add_argument('-m', '--memory', type=int, default=16, help='内存限制GB (默认: 16)')
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.input_r1):
        print(f"❌ 错误: 找不到输入文件 {args.input_r1}")
        sys.exit(1)
    
    if not os.path.exists(args.input_r2):
        print(f"❌ 错误: 找不到输入文件 {args.input_r2}")
        sys.exit(1)
    
    # 创建优化器并运行
    optimizer = AssemblyOptimizer(
        input_r1=args.input_r1,
        input_r2=args.input_r2,
        output_dir=args.output,
        threads=args.threads,
        memory=args.memory
    )
    
    try:
        results = optimizer.run_optimization()
        print("\n🎉 所有组装任务完成！")
        
    except KeyboardInterrupt:
        print("\n⚠️  用户中断操作")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()