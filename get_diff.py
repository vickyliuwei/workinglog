#!/usr/bin/env python3
import gzip
import sys
from collections import defaultdict

def main():
    if len(sys.argv) != 4:
        print("Usage: python genotype_differences.py id.list related_pairs.txt input.vcf.gz")
        sys.exit(1)
    
    # 加载数据
    with open(sys.argv[1]) as f:
        main_ids = [line.strip() for line in f]
    
    related_pairs = defaultdict(list)
    with open(sys.argv[2]) as f:
        for line in f:
            id1, id2 = line.strip().split()
            related_pairs[id1].append(id2)
            related_pairs[id2].append(id1)
    
    # 预加载所有需要的样本
    all_samples = set(main_ids)
    for ids in related_pairs.values():
        all_samples.update(ids)
    
    # 处理VCF
    print("Main_ID\tRelated_ID\tDiff_Count")
    with gzip.open(sys.argv[3], 'rt') as vcf:
        for line in vcf:
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[10:]
                sample_indices = {s:i for i,s in enumerate(samples)}
                break
        
        # 对每个主ID进行处理
        for main_id in main_ids:
            if main_id not in sample_indices:
                continue
                
            related_ids = [rid for rid in related_pairs[main_id] if rid in sample_indices]
            main_idx = sample_indices[main_id]
            
            # 初始化差异计数器
            diff_counts = {rid:0 for rid in related_ids}
            
            # 遍历VCF
            vcf.seek(0)
            for line in vcf:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    main_gt = fields[10+main_idx].split(':')[0]
                    
                    if main_gt in ('./.', '.'):
                        continue
                        
                    for rid in related_ids:
                        rid_idx = sample_indices[rid]
                        rid_gt = fields[10+rid_idx].split(':')[0]
                        
                        if rid_gt not in ('./.', '.') and main_gt != rid_gt:
                            diff_counts[rid] += 1
            
            # 输出结果
            for rid, count in diff_counts.items():
                print(f"{main_id}\t{rid}\t{count}")

if __name__ == '__main__':
    main()
