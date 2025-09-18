#!/usr/bin/env python3
"""
组装图可视化和统计图表生成脚本
用于第2次课实践操作：基因组组装与质量评估

功能：
1. 生成组装统计对比图表
2. 创建质量评估可视化
3. 绘制组装图结构分析
4. 生成参数优化效果图
5. 创建综合分析仪表板

作者：微生物基因组学课程组
版本：v1.0.0
"""

import os
import sys
import argparse
import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Circle, Rectangle
import matplotlib.patches as mpatches
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

class VisualizationGenerator:
    """组装可视化生成器"""
    
    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置中文字体和图表样式
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei', 'Arial Unicode MS']
        plt.rcParams['axes.unicode_minus'] = False
        plt.style.use('seaborn-v0_8-whitegrid')
        
        # 颜色配置
        self.colors = {
            'primary': '#2E86AB',
            'secondary': '#A23B72',
            'accent': '#F18F01',
            'success': '#C73E1D',
            'warning': '#F4A261',
            'info': '#264653',
            'light': '#E9C46A',
            'dark': '#2A9D8F'
        }
        
        # 图表配置
        self.figure_config = {
            'dpi': 300,
            'bbox_inches': 'tight',
            'facecolor': 'white',
            'edgecolor': 'none'
        }
    
    def load_assembly_stats(self, stats_files):
        """加载组装统计数据"""
        print("📊 加载组装统计数据...")
        
        assembly_data = {}
        
        for name, file_path in stats_files.items():
            if not os.path.exists(file_path):
                print(f"⚠️  文件不存在: {file_path}")
                continue
                
            try:
                if file_path.endswith('.json'):
                    with open(file_path, 'r') as f:
                        data = json.load(f)
                        assembly_data[name] = data
                elif file_path.endswith('.txt'):
                    # 解析QUAST报告格式
                    data = self.parse_quast_report(file_path)
                    assembly_data[name] = data
                else:
                    print(f"⚠️  不支持的文件格式: {file_path}")
                    
            except Exception as e:
                print(f"❌ 加载文件失败 {file_path}: {e}")
        
        print(f"✅ 成功加载{len(assembly_data)}个组装的统计数据")
        return assembly_data
    
    def parse_quast_report(self, report_file):
        """解析QUAST报告文件"""
        stats = {}
        
        with open(report_file, 'r') as f:
            lines = f.readlines()
            
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                    
                if '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        key = parts[0].strip()
                        value = parts[1].strip()
                        
                        # 尝试转换为数字
                        try:
                            if '.' in value:
                                stats[key] = float(value)
                            else:
                                stats[key] = int(value)
                        except ValueError:
                            stats[key] = value
        
        return stats
    
    def create_assembly_comparison_chart(self, assembly_data):
        """创建组装结果对比图表"""
        print("📈 生成组装结果对比图表...")
        
        if not assembly_data:
            print("⚠️  无数据可用于生成对比图表")
            return
        
        # 准备数据
        assemblies = list(assembly_data.keys())
        metrics = ['N50', '# contigs', 'Total length', 'Largest contig']
        
        # 创建子图
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('基因组组装结果对比分析', fontsize=18, fontweight='bold', y=0.98)
        
        # 1. N50对比
        n50_values = []
        for assembly in assemblies:
            n50 = assembly_data[assembly].get('N50', 0)
            if isinstance(n50, str):
                n50 = 0
            n50_values.append(n50)
        
        bars1 = axes[0, 0].bar(assemblies, n50_values, color=self.colors['primary'], alpha=0.8)
        axes[0, 0].set_title('N50 连续性对比', fontsize=14, fontweight='bold')
        axes[0, 0].set_ylabel('N50 长度 (bp)')
        axes[0, 0].tick_params(axis='x', rotation=45)
        
        # 添加数值标签
        for bar, value in zip(bars1, n50_values):
            if value > 0:
                axes[0, 0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(n50_values)*0.01,
                               f'{value:,.0f}', ha='center', va='bottom', fontweight='bold')
        
        # 2. Contig数量对比
        contig_counts = []
        for assembly in assemblies:
            count = assembly_data[assembly].get('# contigs', 0)
            if isinstance(count, str):
                count = 0
            contig_counts.append(count)
        
        bars2 = axes[0, 1].bar(assemblies, contig_counts, color=self.colors['secondary'], alpha=0.8)
        axes[0, 1].set_title('Contig 数量对比', fontsize=14, fontweight='bold')
        axes[0, 1].set_ylabel('Contig 数量')
        axes[0, 1].tick_params(axis='x', rotation=45)
        
        # 添加数值标签
        for bar, value in zip(bars2, contig_counts):
            if value > 0:
                axes[0, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(contig_counts)*0.01,
                               f'{value:,}', ha='center', va='bottom', fontweight='bold')
        
        # 3. 总长度对比
        total_lengths = []
        for assembly in assemblies:
            length = assembly_data[assembly].get('Total length', 0)
            if isinstance(length, str):
                length = 0
            total_lengths.append(length)
        
        bars3 = axes[1, 0].bar(assemblies, total_lengths, color=self.colors['accent'], alpha=0.8)
        axes[1, 0].set_title('基因组总长度对比', fontsize=14, fontweight='bold')
        axes[1, 0].set_ylabel('总长度 (bp)')
        axes[1, 0].tick_params(axis='x', rotation=45)
        
        # 添加数值标签
        for bar, value in zip(bars3, total_lengths):
            if value > 0:
                axes[1, 0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(total_lengths)*0.01,
                               f'{value:,.0f}', ha='center', va='bottom', fontweight='bold')
        
        # 4. 综合质量雷达图
        self.create_quality_radar_chart(axes[1, 1], assembly_data, assemblies)
        
        plt.tight_layout()
        
        # 保存图表
        output_file = self.output_dir / 'assembly_comparison.png'
        plt.savefig(output_file, **self.figure_config)
        plt.close()
        
        print(f"✅ 组装对比图表已保存至: {output_file}")
    
    def create_quality_radar_chart(self, ax, assembly_data, assemblies):
        """创建质量评估雷达图"""
        # 简化为质量分数柱状图
        quality_scores = []
        
        for assembly in assemblies:
            data = assembly_data[assembly]
            
            # 计算综合质量分数 (0-100)
            score = 0
            
            # N50分数 (30%)
            n50 = data.get('N50', 0)
            if isinstance(n50, (int, float)) and n50 > 0:
                n50_score = min(30, n50 / 10000)  # 假设10kb为满分
                score += n50_score
            
            # Contig数量分数 (30%) - 越少越好
            contig_count = data.get('# contigs', 1000)
            if isinstance(contig_count, (int, float)):
                contig_score = max(0, 30 - contig_count / 10)
                score += contig_score
            
            # 总长度分数 (20%)
            total_length = data.get('Total length', 0)
            if isinstance(total_length, (int, float)) and total_length > 0:
                length_score = min(20, total_length / 200000)  # 假设4Mb为满分
                score += length_score
            
            # 最大contig分数 (20%)
            largest_contig = data.get('Largest contig', 0)
            if isinstance(largest_contig, (int, float)) and largest_contig > 0:
                largest_score = min(20, largest_contig / 200000)
                score += largest_score
            
            quality_scores.append(min(100, score))
        
        # 绘制质量分数柱状图
        colors = [self.colors['success'] if score >= 80 else 
                 self.colors['warning'] if score >= 60 else 
                 self.colors['accent'] for score in quality_scores]
        
        bars = ax.bar(assemblies, quality_scores, color=colors, alpha=0.8)
        ax.set_title('综合质量评分', fontsize=14, fontweight='bold')
        ax.set_ylabel('质量分数')
        ax.set_ylim(0, 100)
        ax.tick_params(axis='x', rotation=45)
        
        # 添加分数标签
        for bar, score in zip(bars, quality_scores):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                   f'{score:.1f}', ha='center', va='bottom', fontweight='bold')
        
        # 添加质量等级线
        ax.axhline(y=80, color='green', linestyle='--', alpha=0.7, label='优秀 (80+)')
        ax.axhline(y=60, color='orange', linestyle='--', alpha=0.7, label='良好 (60+)')
        ax.legend(loc='upper right')
    
    def create_parameter_optimization_chart(self, optimization_results):
        """创建参数优化效果图"""
        print("📊 生成参数优化效果图...")
        
        if not optimization_results:
            print("⚠️  无参数优化数据")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('SPAdes参数优化效果分析', fontsize=18, fontweight='bold', y=0.98)
        
        # 准备数据
        param_names = []
        n50_values = []
        contig_counts = []
        runtimes = []
        
        for param_name, result in optimization_results.items():
            if result.get('success') and result.get('stats'):
                param_names.append(param_name)
                n50_values.append(result['stats'].get('n50', 0))
                contig_counts.append(result['stats'].get('total_contigs', 0))
                runtimes.append(result.get('runtime', 0))
        
        if not param_names:
            print("⚠️  无有效的优化结果数据")
            return
        
        # 1. N50优化效果
        bars1 = axes[0, 0].bar(param_names, n50_values, color=self.colors['primary'], alpha=0.8)
        axes[0, 0].set_title('参数对N50的影响', fontsize=14, fontweight='bold')
        axes[0, 0].set_ylabel('N50 (bp)')
        axes[0, 0].tick_params(axis='x', rotation=45)
        
        # 标记最佳结果
        best_n50_idx = n50_values.index(max(n50_values))
        bars1[best_n50_idx].set_color(self.colors['success'])
        
        # 2. Contig数量优化效果
        bars2 = axes[0, 1].bar(param_names, contig_counts, color=self.colors['secondary'], alpha=0.8)
        axes[0, 1].set_title('参数对Contig数量的影响', fontsize=14, fontweight='bold')
        axes[0, 1].set_ylabel('Contig数量')
        axes[0, 1].tick_params(axis='x', rotation=45)
        
        # 标记最佳结果（contig数量越少越好）
        best_contig_idx = contig_counts.index(min(contig_counts))
        bars2[best_contig_idx].set_color(self.colors['success'])
        
        # 3. 运行时间对比
        bars3 = axes[1, 0].bar(param_names, runtimes, color=self.colors['accent'], alpha=0.8)
        axes[1, 0].set_title('参数对运行时间的影响', fontsize=14, fontweight='bold')
        axes[1, 0].set_ylabel('运行时间 (秒)')
        axes[1, 0].tick_params(axis='x', rotation=45)
        
        # 4. 效率vs质量散点图
        # 计算质量分数 (N50/contig_count)
        efficiency_scores = [n50/contig if contig > 0 else 0 for n50, contig in zip(n50_values, contig_counts)]
        
        scatter = axes[1, 1].scatter(runtimes, efficiency_scores, 
                                   c=range(len(param_names)), 
                                   cmap='viridis', s=100, alpha=0.8)
        
        # 添加参数标签
        for i, param in enumerate(param_names):
            axes[1, 1].annotate(param, (runtimes[i], efficiency_scores[i]), 
                              xytext=(5, 5), textcoords='offset points', fontsize=10)
        
        axes[1, 1].set_xlabel('运行时间 (秒)')
        axes[1, 1].set_ylabel('效率分数 (N50/Contig数)')
        axes[1, 1].set_title('运行时间 vs 组装效率', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        
        # 保存图表
        output_file = self.output_dir / 'parameter_optimization.png'
        plt.savefig(output_file, **self.figure_config)
        plt.close()
        
        print(f"✅ 参数优化图表已保存至: {output_file}")
    
    def create_quality_assessment_dashboard(self, busco_results, checkm_results, quast_results):
        """创建质量评估仪表板"""
        print("📊 生成质量评估仪表板...")
        
        fig = plt.figure(figsize=(20, 16))
        gs = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)
        
        fig.suptitle('基因组组装质量评估仪表板', fontsize=20, fontweight='bold', y=0.98)
        
        # 1. BUSCO完整性热图 (占据2x2空间)
        ax1 = fig.add_subplot(gs[0:2, 0:2])
        self.create_busco_heatmap(ax1, busco_results)
        
        # 2. CheckM质量分级 (1x2空间)
        ax2 = fig.add_subplot(gs[0, 2:4])
        self.create_checkm_quality_chart(ax2, checkm_results)
        
        # 3. QUAST统计雷达图 (1x2空间)
        ax3 = fig.add_subplot(gs[1, 2:4])
        self.create_quast_radar_chart(ax3, quast_results)
        
        # 4. 组装连续性分布 (2x2空间)
        ax4 = fig.add_subplot(gs[2:4, 0:2])
        self.create_continuity_distribution(ax4, quast_results)
        
        # 5. 质量等级饼图 (1x1空间)
        ax5 = fig.add_subplot(gs[2, 2])
        self.create_quality_grade_pie(ax5, busco_results, checkm_results)
        
        # 6. 推荐组装排名 (1x1空间)
        ax6 = fig.add_subplot(gs[2, 3])
        self.create_recommendation_ranking(ax6, busco_results, checkm_results, quast_results)
        
        # 7. 质量趋势线 (2x1空间)
        ax7 = fig.add_subplot(gs[3, 2:4])
        self.create_quality_trend_line(ax7, busco_results, checkm_results)
        
        # 保存仪表板
        output_file = self.output_dir / 'quality_assessment_dashboard.png'
        plt.savefig(output_file, **self.figure_config)
        plt.close()
        
        print(f"✅ 质量评估仪表板已保存至: {output_file}")
    
    def create_busco_heatmap(self, ax, busco_results):
        """创建BUSCO结果热图"""
        if not busco_results:
            ax.text(0.5, 0.5, '无BUSCO数据', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('BUSCO完整性热图')
            return
        
        # 准备数据
        assemblies = list(busco_results.keys())
        metrics = ['complete', 'single_copy', 'duplicated', 'fragmented', 'missing']
        
        data_matrix = []
        for assembly in assemblies:
            if busco_results[assembly]:
                row = [busco_results[assembly].get(metric, 0) for metric in metrics]
                data_matrix.append(row)
            else:
                data_matrix.append([0] * len(metrics))
        
        # 创建热图
        im = ax.imshow(data_matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=100)
        
        # 设置标签
        ax.set_xticks(range(len(metrics)))
        ax.set_xticklabels(['完整', '单拷贝', '重复', '片段化', '缺失'], rotation=45)
        ax.set_yticks(range(len(assemblies)))
        ax.set_yticklabels(assemblies)
        
        # 添加数值标签
        for i in range(len(assemblies)):
            for j in range(len(metrics)):
                text = ax.text(j, i, f'{data_matrix[i][j]:.1f}%',
                             ha="center", va="center", color="black", fontweight='bold')
        
        ax.set_title('BUSCO完整性热图', fontsize=14, fontweight='bold')
        
        # 添加颜色条
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('完整性百分比 (%)')
    
    def create_checkm_quality_chart(self, ax, checkm_results):
        """创建CheckM质量图表"""
        if not checkm_results:
            ax.text(0.5, 0.5, '无CheckM数据', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('CheckM质量评估')
            return
        
        assemblies = list(checkm_results.keys())
        completeness = []
        contamination = []
        
        for assembly in assemblies:
            if checkm_results[assembly]:
                completeness.append(checkm_results[assembly].get('completeness', 0))
                contamination.append(checkm_results[assembly].get('contamination', 0))
            else:
                completeness.append(0)
                contamination.append(0)
        
        # 创建散点图
        colors = []
        for comp, cont in zip(completeness, contamination):
            if comp > 90 and cont < 5:
                colors.append(self.colors['success'])  # 高质量
            elif comp > 50 and cont < 10:
                colors.append(self.colors['warning'])  # 中等质量
            else:
                colors.append(self.colors['accent'])   # 低质量
        
        scatter = ax.scatter(contamination, completeness, c=colors, s=100, alpha=0.8)
        
        # 添加质量区域
        ax.axhline(y=90, color='green', linestyle='--', alpha=0.5, label='高质量线')
        ax.axvline(x=5, color='red', linestyle='--', alpha=0.5, label='污染阈值')
        
        # 添加标签
        for i, assembly in enumerate(assemblies):
            ax.annotate(assembly, (contamination[i], completeness[i]), 
                       xytext=(5, 5), textcoords='offset points', fontsize=9)
        
        ax.set_xlabel('污染水平 (%)')
        ax.set_ylabel('完整性 (%)')
        ax.set_title('CheckM质量评估', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def create_quast_radar_chart(self, ax, quast_results):
        """创建QUAST统计雷达图（简化为柱状图）"""
        if not quast_results:
            ax.text(0.5, 0.5, '无QUAST数据', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('QUAST统计分析')
            return
        
        # 选择第一个组装的关键指标
        first_assembly = list(quast_results.keys())[0]
        data = quast_results[first_assembly]
        
        metrics = ['N50', '# contigs', 'Total length', 'Largest contig']
        values = []
        
        for metric in metrics:
            value = data.get(metric, 0)
            if isinstance(value, (int, float)):
                # 标准化值 (0-100)
                if metric == 'N50':
                    normalized = min(100, value / 1000)
                elif metric == '# contigs':
                    normalized = max(0, 100 - value / 10)  # 越少越好
                elif metric == 'Total length':
                    normalized = min(100, value / 50000)
                elif metric == 'Largest contig':
                    normalized = min(100, value / 10000)
                else:
                    normalized = 0
                values.append(normalized)
            else:
                values.append(0)
        
        # 创建柱状图
        bars = ax.bar(metrics, values, color=self.colors['info'], alpha=0.8)
        ax.set_title('QUAST关键指标', fontsize=14, fontweight='bold')
        ax.set_ylabel('标准化分数')
        ax.set_ylim(0, 100)
        ax.tick_params(axis='x', rotation=45)
        
        # 添加数值标签
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                   f'{value:.1f}', ha='center', va='bottom', fontweight='bold')
    
    def create_continuity_distribution(self, ax, quast_results):
        """创建组装连续性分布图"""
        if not quast_results:
            ax.text(0.5, 0.5, '无连续性数据', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('组装连续性分布')
            return
        
        assemblies = list(quast_results.keys())
        n50_values = []
        n90_values = []
        
        for assembly in assemblies:
            data = quast_results[assembly]
            n50_values.append(data.get('N50', 0))
            n90_values.append(data.get('N90', 0))
        
        x = np.arange(len(assemblies))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, n50_values, width, label='N50', 
                      color=self.colors['primary'], alpha=0.8)
        bars2 = ax.bar(x + width/2, n90_values, width, label='N90', 
                      color=self.colors['secondary'], alpha=0.8)
        
        ax.set_xlabel('组装方法')
        ax.set_ylabel('长度 (bp)')
        ax.set_title('组装连续性分布 (N50 vs N90)', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(assemblies, rotation=45)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def create_quality_grade_pie(self, ax, busco_results, checkm_results):
        """创建质量等级饼图"""
        grades = {'高质量': 0, '中等质量': 0, '低质量': 0}
        
        # 基于BUSCO和CheckM结果评估质量等级
        if busco_results and checkm_results:
            for assembly in busco_results.keys():
                busco_data = busco_results.get(assembly, {})
                checkm_data = checkm_results.get(assembly, {})
                
                if busco_data and checkm_data:
                    busco_complete = busco_data.get('complete', 0)
                    checkm_complete = checkm_data.get('completeness', 0)
                    checkm_contam = checkm_data.get('contamination', 100)
                    
                    if busco_complete > 90 and checkm_complete > 90 and checkm_contam < 5:
                        grades['高质量'] += 1
                    elif busco_complete > 70 and checkm_complete > 70 and checkm_contam < 10:
                        grades['中等质量'] += 1
                    else:
                        grades['低质量'] += 1
        
        # 过滤零值
        grades = {k: v for k, v in grades.items() if v > 0}
        
        if grades:
            colors = [self.colors['success'], self.colors['warning'], self.colors['accent']][:len(grades)]
            wedges, texts, autotexts = ax.pie(grades.values(), labels=grades.keys(), 
                                            autopct='%1.1f%%', colors=colors, startangle=90)
            ax.set_title('质量等级分布', fontsize=12, fontweight='bold')
        else:
            ax.text(0.5, 0.5, '无质量数据', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('质量等级分布')
    
    def create_recommendation_ranking(self, ax, busco_results, checkm_results, quast_results):
        """创建推荐组装排名"""
        scores = {}
        
        # 计算综合分数
        if busco_results:
            for assembly in busco_results.keys():
                score = 0
                
                # BUSCO分数
                busco_data = busco_results.get(assembly, {})
                if busco_data:
                    score += busco_data.get('complete', 0) * 0.4
                
                # CheckM分数
                if checkm_results:
                    checkm_data = checkm_results.get(assembly, {})
                    if checkm_data:
                        completeness = checkm_data.get('completeness', 0)
                        contamination = checkm_data.get('contamination', 0)
                        checkm_score = max(0, completeness - 5 * contamination)
                        score += checkm_score * 0.3
                
                # QUAST分数
                if quast_results:
                    quast_data = quast_results.get(assembly, {})
                    if quast_data:
                        n50 = quast_data.get('N50', 0)
                        if isinstance(n50, (int, float)):
                            n50_score = min(30, n50 / 1000)
                            score += n50_score
                
                scores[assembly] = score
        
        if scores:
            # 排序并显示前3名
            sorted_assemblies = sorted(scores.items(), key=lambda x: x[1], reverse=True)[:3]
            
            assemblies = [item[0] for item in sorted_assemblies]
            scores_list = [item[1] for item in sorted_assemblies]
            
            colors = [self.colors['success'], self.colors['warning'], self.colors['accent']][:len(assemblies)]
            bars = ax.barh(assemblies, scores_list, color=colors, alpha=0.8)
            
            ax.set_xlabel('综合分数')
            ax.set_title('推荐排名', fontsize=12, fontweight='bold')
            
            # 添加分数标签
            for bar, score in zip(bars, scores_list):
                ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                       f'{score:.1f}', ha='left', va='center', fontweight='bold')
        else:
            ax.text(0.5, 0.5, '无排名数据', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('推荐排名')
    
    def create_quality_trend_line(self, ax, busco_results, checkm_results):
        """创建质量趋势线"""
        if not busco_results or not checkm_results:
            ax.text(0.5, 0.5, '无趋势数据', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('质量趋势分析')
            return
        
        assemblies = list(busco_results.keys())
        busco_scores = []
        checkm_scores = []
        
        for assembly in assemblies:
            # BUSCO完整性
            busco_data = busco_results.get(assembly, {})
            busco_scores.append(busco_data.get('complete', 0) if busco_data else 0)
            
            # CheckM质量分数
            checkm_data = checkm_results.get(assembly, {})
            if checkm_data:
                completeness = checkm_data.get('completeness', 0)
                contamination = checkm_data.get('contamination', 0)
                checkm_score = max(0, completeness - 5 * contamination)
                checkm_scores.append(checkm_score)
            else:
                checkm_scores.append(0)
        
        x = range(len(assemblies))
        
        ax.plot(x, busco_scores, marker='o', linewidth=2, label='BUSCO完整性', 
               color=self.colors['primary'])
        ax.plot(x, checkm_scores, marker='s', linewidth=2, label='CheckM质量分数', 
               color=self.colors['secondary'])
        
        ax.set_xlabel('组装方法')
        ax.set_ylabel('分数')
        ax.set_title('质量趋势分析', fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(assemblies, rotation=45)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def generate_comprehensive_report(self, assembly_data, optimization_results, 
                                    busco_results, checkm_results, quast_results):
        """生成综合可视化报告"""
        print("🎨 生成综合可视化报告...")
        
        # 1. 组装对比图表
        if assembly_data:
            self.create_assembly_comparison_chart(assembly_data)
        
        # 2. 参数优化图表
        if optimization_results:
            self.create_parameter_optimization_chart(optimization_results)
        
        # 3. 质量评估仪表板
        if busco_results or checkm_results or quast_results:
            self.create_quality_assessment_dashboard(busco_results, checkm_results, quast_results)
        
        # 4. 生成HTML报告索引
        self.create_html_report_index()
        
        print(f"✅ 所有可视化图表已生成完成")
        print(f"📁 输出目录: {self.output_dir}")
    
    def create_html_report_index(self):
        """创建HTML报告索引"""
        html_content = f"""
        <!DOCTYPE html>
        <html lang="zh-CN">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>基因组组装质量评估报告</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }}
                .container {{ max-width: 1200px; margin: 0 auto; background-color: white; padding: 20px; border-radius: 10px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }}
                h1 {{ color: #2E86AB; text-align: center; }}
                h2 {{ color: #A23B72; border-bottom: 2px solid #A23B72; padding-bottom: 5px; }}
                .chart-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; margin: 20px 0; }}
                .chart-item {{ text-align: center; padding: 15px; border: 1px solid #ddd; border-radius: 8px; }}
                .chart-item img {{ max-width: 100%; height: auto; border-radius: 5px; }}
                .timestamp {{ text-align: center; color: #666; font-style: italic; }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>🧬 基因组组装质量评估报告</h1>
                
                <h2>📊 可视化图表</h2>
                <div class="chart-grid">
                    <div class="chart-item">
                        <h3>组装结果对比</h3>
                        <img src="assembly_comparison.png" alt="组装结果对比">
                        <p>展示不同组装方法的N50、contig数量、总长度等关键指标对比</p>
                    </div>
                    
                    <div class="chart-item">
                        <h3>参数优化效果</h3>
                        <img src="parameter_optimization.png" alt="参数优化效果">
                        <p>分析不同SPAdes参数设置对组装质量和运行时间的影响</p>
                    </div>
                    
                    <div class="chart-item">
                        <h3>质量评估仪表板</h3>
                        <img src="quality_assessment_dashboard.png" alt="质量评估仪表板">
                        <p>综合展示BUSCO、CheckM、QUAST等多维度质量评估结果</p>
                    </div>
                    
                    <div class="chart-item">
                        <h3>污染检测分析</h3>
                        <img src="contamination_analysis.png" alt="污染检测分析">
                        <p>GC含量分布、序列异常检测和污染风险评估</p>
                    </div>
                </div>
                
                <h2>📋 报告说明</h2>
                <ul>
                    <li><strong>组装对比</strong>: 比较不同组装策略的效果，帮助选择最佳方法</li>
                    <li><strong>参数优化</strong>: 展示参数调整对组装质量的影响，指导参数选择</li>
                    <li><strong>质量评估</strong>: 多维度评估组装完整性、准确性和污染水平</li>
                    <li><strong>污染检测</strong>: 识别潜在的序列污染和异常区域</li>
                </ul>
                
                <div class="timestamp">
                    <p>报告生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
                </div>
            </div>
        </body>
        </html>
        """
        
        html_file = self.output_dir / 'visualization_report.html'
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"📄 HTML报告已生成: {html_file}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='基因组组装可视化图表生成')
    parser.add_argument('-o', '--output', default='results/visualization',
                       help='输出目录 (默认: results/visualization)')
    parser.add_argument('--assembly-stats', nargs='*', 
                       help='组装统计文件 (格式: name:path)')
    parser.add_argument('--optimization-results', 
                       help='参数优化结果JSON文件')
    parser.add_argument('--busco-results', 
                       help='BUSCO结果JSON文件')
    parser.add_argument('--checkm-results', 
                       help='CheckM结果JSON文件')
    parser.add_argument('--quast-results', 
                       help='QUAST结果JSON文件')
    
    args = parser.parse_args()
    
    # 创建可视化生成器
    generator = VisualizationGenerator(args.output)
    
    # 加载数据
    assembly_data = {}
    optimization_results = {}
    busco_results = {}
    checkm_results = {}
    quast_results = {}
    
    # 加载组装统计数据
    if args.assembly_stats:
        stats_files = {}
        for item in args.assembly_stats:
            if ':' in item:
                name, path = item.split(':', 1)
                stats_files[name] = path
        assembly_data = generator.load_assembly_stats(stats_files)
    
    # 加载其他结果文件
    for result_file, result_dict in [
        (args.optimization_results, optimization_results),
        (args.busco_results, busco_results),
        (args.checkm_results, checkm_results),
        (args.quast_results, quast_results)
    ]:
        if result_file and os.path.exists(result_file):
            try:
                with open(result_file, 'r') as f:
                    result_dict.update(json.load(f))
            except Exception as e:
                print(f"⚠️  加载文件失败 {result_file}: {e}")
    
    # 生成可视化报告
    try:
        generator.generate_comprehensive_report(
            assembly_data, optimization_results, 
            busco_results, checkm_results, quast_results
        )
        print("\n🎉 可视化生成任务完成！")
        
    except Exception as e:
        print(f"\n❌ 发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()