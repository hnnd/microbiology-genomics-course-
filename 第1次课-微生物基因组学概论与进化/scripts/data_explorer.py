#!/usr/bin/env python3
"""
交互式数据探索工具
"""

import pandas as pd
import matplotlib.pyplot as plt

def explore_data():
    """交互式数据探索"""
    try:
        df = pd.read_csv('results/cleaned_genome_data.csv')
        print("=== 微生物基因组数据探索工具 ===\n")
        
        while True:
            print("请选择分析选项:")
            print("1. 查看数据基本信息")
            print("2. 按分类群统计")
            print("3. 查找特定基因组")
            print("4. 极值分析")
            print("5. 退出")
            
            choice = input("\n请输入选项编号 (1-5): ").strip()
            
            if choice == '1':
                print(f"\n数据基本信息:")
                print(f"总记录数: {len(df)}")
                print(f"数据列: {list(df.columns)}")
                print(f"\n基因组大小统计:")
                print(df['Size (Mb)'].describe())
                print(f"\nGC含量统计:")
                print(df['GC%'].describe())
                
            elif choice == '2':
                if 'Group' in df.columns:
                    print(f"\n按分类群统计:")
                    group_counts = df['Group'].value_counts().head(10)
                    print(group_counts)
                else:
                    print("数据中没有分类群信息")
                    
            elif choice == '3':
                search_term = input("请输入要搜索的基因组名称关键词: ").strip()
                if search_term:
                    matches = df[df['#Organism Name'].str.contains(search_term, case=False, na=False)]
                    if len(matches) > 0:
                        print(f"\n找到 {len(matches)} 个匹配结果:")
                        for idx, row in matches.head(10).iterrows():
                            print(f"- {row['#Organism Name']}: {row['Size (Mb)']:.2f} Mb, GC {row['GC%']:.1f}%")
                    else:
                        print("未找到匹配的基因组")
                        
            elif choice == '4':
                print(f"\n极值分析:")
                
                # 最大基因组
                max_size_idx = df['Size (Mb)'].idxmax()
                max_genome = df.loc[max_size_idx]
                print(f"最大基因组: {max_genome['#Organism Name']}")
                print(f"  大小: {max_genome['Size (Mb)']:.2f} Mb")
                print(f"  GC含量: {max_genome['GC%']:.1f}%")
                
                # 最小基因组
                min_size_idx = df['Size (Mb)'].idxmin()
                min_genome = df.loc[min_size_idx]
                print(f"\n最小基因组: {min_genome['#Organism Name']}")
                print(f"  大小: {min_genome['Size (Mb)']:.3f} Mb")
                print(f"  GC含量: {min_genome['GC%']:.1f}%")
                
                # 极端GC含量
                max_gc_idx = df['GC%'].idxmax()
                min_gc_idx = df['GC%'].idxmin()
                
                print(f"\n最高GC含量: {df.loc[max_gc_idx, '#Organism Name']}")
                print(f"  GC含量: {df.loc[max_gc_idx, 'GC%']:.1f}%")
                print(f"  基因组大小: {df.loc[max_gc_idx, 'Size (Mb)']:.2f} Mb")
                
                print(f"\n最低GC含量: {df.loc[min_gc_idx, '#Organism Name']}")
                print(f"  GC含量: {df.loc[min_gc_idx, 'GC%']:.1f}%")
                print(f"  基因组大小: {df.loc[min_gc_idx, 'Size (Mb)']:.2f} Mb")
                
            elif choice == '5':
                print("退出数据探索工具")
                break
                
            else:
                print("无效选项，请重新选择")
            
            print("\n" + "="*50 + "\n")
            
    except FileNotFoundError:
        print("找不到数据文件，请先运行基础分析脚本")
    except Exception as e:
        print(f"发生错误: {e}")

if __name__ == "__main__":
    explore_data()