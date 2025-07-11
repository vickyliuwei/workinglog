#!/usr/bin/env python3
import gzip
import sys
from collections import defaultdict

def load_species_mapping(mapping_file):
    """加载ID到物种的映射"""
    id_to_species = {}
    species_ids = defaultdict(list)
    with open(mapping_file, 'r') as f:
        for line in f:
            id, species = line.strip().split()
            id_to_species[id] = species
            species_ids[species].append(id)
    return id_to_species, species_ids

def calculate_genotype_freq(genotypes):
    """计算基因型频率（忽略缺失数据）"""
    gt_counts = defaultdict(int)
    total = 0
    for gt in genotypes:
        if gt not in ('./.', '.'):
            gt_counts[gt] += 1
            total += 1
    return {gt: count/total for gt, count in gt_counts.items()} if total > 0 else {}

def analyze_vcf(vcf_file, id_to_species, species_ids):
    """分析VCF中的位点共享模式"""
    species_list = sorted(species_ids.keys())
    if len(species_list) != 4:
        raise ValueError("必须包含4个物种")
    
    results = {
        'unique': {s: 0 for s in species_list},
        'shared': defaultdict(int),
        'all_shared': 0
    }
    
    with gzip.open(vcf_file, 'rt') as f:
        # 读取样本头信息
        for line in f:
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                # 创建样本索引
                species_indices = {s: [i for i, id in enumerate(samples) if id in species_ids[s]] 
                                 for s in species_list}
                break
        
        # 处理每个变异位点
        for line in f:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                chrom, pos = fields[0], fields[1]
                variant_id = f"{chrom}:{pos}"
                gts = [fields[i+9].split(':')[0].replace('|', '/') for i in range(len(samples))]
                
                # 计算每个物种的基因型频率
                species_freq = {}
                for s in species_list:
                    s_gts = [gts[i] for i in species_indices[s]]
                    species_freq[s] = calculate_genotype_freq(s_gts)
                
                # 找出有数据的物种
                present_species = [s for s in species_list if species_freq[s]]
                if len(present_species) < 1:
                    continue
                
                # 比较频率差异
                if len(present_species) == 1:
                    results['unique'][present_species[0]] += 1
                else:
                    # 主物种与其他物种比较
                    main_species = present_species[0]
                    other_freq = {}
                    for s in present_species[1:]:
                        for gt, freq in species_freq[s].items():
                            other_freq[gt] = other_freq.get(gt, 0) + freq
                    # 标准化其他物种合并频率
                    other_total = sum(other_freq.values())
                    other_freq = {gt: f/other_total for gt, f in other_freq.items()}
                    
                    # 计算频率差异
                    diff_scores = []
                    for gt in set(species_freq[main_species]) | set(other_freq):
                        f1 = species_freq[main_species].get(gt, 0)
                        f2 = other_freq.get(gt, 0)
                        diff_scores.append(abs(f1 - f2))
                    
                    max_diff = max(diff_scores) if diff_scores else 0
                    
                    if max_diff == 1:  # 主物种特有
                        results['unique'][main_species] += 1
                    elif max_diff == 0:  # 全部共有
                        results['all_shared'] += 1
                    else:  # 部分共有
                        # 两两比较找出具体共享模式
                        shared_pairs = set()
                        for s1, s2 in itertools.combinations(present_species, 2):
                            diff = compare_species_freq(species_freq[s1], species_freq[s2])
                            if diff == 0:
                                shared_pairs.add(frozenset({s1, s2}))
                        
                        if len(shared_pairs) == 1:  # 仅一对物种共享
                            pair = tuple(sorted(next(iter(shared_pairs))))
                            results['shared'][pair] += 1
                        else:  # 更复杂的共享模式
                            results['shared']['complex'] += 1
    
    return results

def compare_species_freq(freq1, freq2):
    """比较两个物种的基因型频率差异"""
    all_gts = set(freq1.keys()) | set(freq2.keys())
    diffs = [abs(freq1.get(gt,0) - freq2.get(gt,0)) for gt in all_gts]
    return max(diffs) if diffs else 1

def generate_report(results, species_list):
    """生成分析报告"""
    print("=== 变异位点分类统计 ===")
    print(f"所有物种共有位点: {results['all_shared']}")
    
    print("\n=== 物种特有变异 ===")
    for s in species_list:
        print(f"{s}特有位点: {results['unique'][s]}")
    
    print("\n=== 两物种共有变异 ===")
    for pair in itertools.combinations(species_list, 2):
        print(f"{pair[0]}-{pair[1]}共有位点: {results['shared'].get(pair, 0)}")
    
    print(f"\n复杂共享模式位点: {results['shared'].get('complex', 0)}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <id_species.map> <input.vcf.gz>")
        sys.exit(1)
    
    mapping_file = sys.argv[1]
    vcf_file = sys.argv[2]
    
    # 加载数据
    id_to_species, species_ids = load_species_mapping(mapping_file)
    species_list = sorted(species_ids.keys())
    print(f"分析物种: {', '.join(species_list)}")
    
    # 分析VCF
    results = analyze_vcf(vcf_file, id_to_species, species_ids)
    
    # 生成报告
    generate_report(results, species_list)

if __name__ == '__main__':
    import itertools
    main()
