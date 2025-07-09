import sys
from collections import defaultdict

def parse_attributes(attr_str):
    """解析GFF属性字符串"""
    attrs = {}
    for item in attr_str.split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value
    return attrs

def process_gff(input_path, output_path):
    """处理GFF文件的主函数"""
    # 存储转录本的所有特征位置 {tid: {chrom1, chrom2, ...}}
    transcript_features = defaultdict(set)
    # 存储基因到转录本的映射 {gene_id: [tid1, tid2]}
    gene_transcripts = defaultdict(list)
    # 存储所有行及其相关信息
    lines = []
    
    # 第一次读取：收集信息
    with open(input_path, 'rb') as f:
        for line_bytes in f:
            try:
                line = line_bytes.decode('utf-8')
            except UnicodeDecodeError:
                try:
                    line = line_bytes.decode('latin-1')
                except:
                    line = line_bytes.decode('utf-8', errors='replace')
            
            if line.startswith('#') or not line.strip():
                lines.append(('comment', line))
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                lines.append(('other', line))
                continue
            
            chrom, source, feat_type, start, end, score, strand, phase, attrs = parts
            attrs_dict = parse_attributes(attrs)
            
            if feat_type.lower() == 'gene':
                if 'ID' in attrs_dict:
                    gene_id = attrs_dict['ID']
                    lines.append(('gene', gene_id, line))
                else:
                    lines.append(('other', line))
            
            elif feat_type.lower() in ['mrna', 'transcript']:
                if 'ID' in attrs_dict and 'Parent' in attrs_dict:
                    tid = attrs_dict['ID']
                    gene_id = attrs_dict['Parent']
                    transcript_features[tid].add(chrom)  # 记录转录本自身位置
                    gene_transcripts[gene_id].append(tid)
                    lines.append(('transcript', tid, gene_id, line))
                else:
                    lines.append(('other', line))
            
            elif feat_type.lower() in ['exon', 'cds', 'utr', 'three_prime_utr', 'five_prime_utr']:
                if 'Parent' in attrs_dict:
                    parents = attrs_dict['Parent'].split(',')  # 处理可能有多个Parent的情况
                    for parent in parents:
                        transcript_features[parent].add(chrom)  # 记录特征位置
                    lines.append(('feature', parents, line))
                else:
                    lines.append(('other', line))
            
            else:
                lines.append(('other', line))
    
    # 识别跨区段的转录本
    cross_transcripts = set()
    for tid, chroms in transcript_features.items():
        if len(chroms) > 1:  # 如果转录本或其任何特征位于不同区段
            cross_transcripts.add(tid)
    
    # 第二次处理：写入输出
    with open(output_path, 'w', encoding='utf-8') as out:
        for item in lines:
            if item[0] == 'comment':
                out.write(item[1])
            elif item[0] == 'other':
                out.write(item[1])
            elif item[0] == 'gene':
                gene_id, line = item[1], item[2]
                # 检查该基因是否有至少一个非跨区段转录本
                valid_transcripts = [tid for tid in gene_transcripts[gene_id] 
                                  if tid not in cross_transcripts]
                if valid_transcripts:
                    out.write(line)
            elif item[0] == 'transcript':
                tid, gene_id, line = item[1], item[2], item[3]
                if tid not in cross_transcripts:
                    out.write(line)
            elif item[0] == 'feature':
                parents, line = item[1], item[2]
                # 只有当所有父转录本都不跨区段时才保留
                if not any(parent in cross_transcripts for parent in parents):
                    out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_gff.py input.gff output.gff")
        sys.exit(1)
    
    process_gff(sys.argv[1], sys.argv[2])
    print(f"Processed GFF saved to {sys.argv[2]}")
