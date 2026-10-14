#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
from collections import defaultdict
from io import StringIO

# ==================== 配置 ====================
DIST_FILE = '~/chewbbaca/distance_matrix_all_samples_nie.txt'
ST_FILE = os.path.expanduser('~/chewbbaca/Salmonella_enterica_MLST_result_addMIC_14880.txt')
OUTPUT_FILE = '~/chewbbaca/cgmlst_mlst10_group_assignment.txt'

CGMLST_THRESHOLD = 3
SMALL_ST_MAX_SIZE = 10

# ==================== 并查集 ====================
class UnionFind:
    def __init__(self, elements):
        self.parent = {e: e for e in elements}
        self.rank = {e: 0 for e in elements}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            self.parent[rx] = ry
        elif self.rank[rx] > self.rank[ry]:
            self.parent[ry] = rx
        else:
            self.parent[ry] = rx
            self.rank[rx] += 1

# ==================== 1. 读取距离矩阵 ====================
print("读取距离矩阵...")
df_dist = pd.read_csv(
    DIST_FILE,
    sep='\t',
    index_col=0,
    na_values=['N/A', 'NA', '']
)
df_dist = df_dist.astype(np.float32)
# 确保列顺序与行顺序一致
df_dist = df_dist.loc[df_dist.index, df_dist.index]
samples = df_dist.index.tolist()
n = len(samples)
dist_mat = df_dist.values
print(f"距离矩阵包含 {n} 个样本")

# ==================== 2. 读取 ST 信息 ====================
print("读取 ST 信息...")
df_st = pd.read_csv(ST_FILE, sep='\t', low_memory=False)
df_st.columns = df_st.columns.str.strip()
GENOME_COL = 'Genome_name'
ST_COL = 'MLST'

st_map = {}
for _, row in df_st.iterrows():
    genome = str(row[GENOME_COL]).strip()
    st = row[ST_COL]
    if pd.isna(st) or str(st).strip().upper() in ('N/A', 'NA', 'UNKNOWN', ''):
        st = 'Unknown'
    st_map[genome] = str(st).strip()

# 检查未匹配样本
missing = [s for s in samples if s not in st_map]
if missing:
    print(f"警告：{len(missing)} 个样本在 ST 文件中未找到，将标记为 Unknown。")
    for s in missing[:5]:
        print(f"  示例：{s}")

# ==================== 3. 初始化并查集 ====================
uf = UnionFind(samples)

# ==================== 4. 合并距离 ≤ 阈值的样本 ====================
print(f"合并 cgMLST 距离 ≤ {CGMLST_THRESHOLD} 的样本...")
for i, sample_i in enumerate(samples):
    row = dist_mat[i, :]
    js = np.where((row <= CGMLST_THRESHOLD) & (np.arange(n) > i))[0]
    for j in js:
        uf.union(sample_i, samples[j])

# ==================== 5. 小 ST 整体合并 ====================
print(f"合并样本数 ≤ {SMALL_ST_MAX_SIZE} 的 ST ...")
st_to_samples = defaultdict(list)
for sample in samples:
    st = st_map.get(sample, 'Unknown')
    if st != 'Unknown':
        st_to_samples[st].append(sample)

for st, sample_list in st_to_samples.items():
    if len(sample_list) <= SMALL_ST_MAX_SIZE:
        first = sample_list[0]
        for s in sample_list[1:]:
            uf.union(first, s)

# ==================== 6. 生成最终 group_id ====================
print("生成 group 分配...")
group_root = {sample: uf.find(sample) for sample in samples}
unique_roots = sorted(set(group_root.values()))
root_to_int = {root: idx + 1 for idx, root in enumerate(unique_roots)}
final_groups = {sample: root_to_int[root] for sample, root in group_root.items()}

# ==================== 7. 输出 ====================
df_out = pd.DataFrame({
    'Genome_name': list(final_groups.keys()),
    'group_id': list(final_groups.values())
})
df_out.to_csv(OUTPUT_FILE, sep='\t', index=False)
print(f"完成！结果已保存至 {OUTPUT_FILE}")
print(f"共生成 {len(unique_roots)} 个不可分 group。")
