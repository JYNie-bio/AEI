#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold

# ==================== 配置 ====================
GROUP_FILE = '~/chewbbaca/cgmlst_mlst10_group_assignment.txt'
ST_FILE = os.path.expanduser(
    '~/chewbbaca/Salmonella_enterica_MLST_result_addMIC_14880.txt'
)
OUTPUT_DIR = '~/chewbbaca/split_results'
os.makedirs(OUTPUT_DIR, exist_ok=True)

RANDOM_STATE = 42

# 输出时保留的列
OUTPUT_COLS_BASE = ['Genome_name', 'group_id', 'MLST']


# ==================== 参数解析 ====================
def parse_args():
    parser = argparse.ArgumentParser(
        description='按 group 拆分数据集，支持指定抗生素列。'
    )
    parser.add_argument(
        'columns', nargs='*', default=None,
        help='要处理的抗生素列名，例如 CIP AMP SXT。'
             '不指定则处理默认的 15 个抗生素列。'
    )
    parser.add_argument(
        '--group_file', default=GROUP_FILE,
        help=f'group 分配文件路径（默认 {GROUP_FILE}）'
    )
    parser.add_argument(
        '--st_file', default=ST_FILE,
        help='ST/耐药表型文件路径'
    )
    parser.add_argument(
        '--output_dir', default=OUTPUT_DIR,
        help=f'输出目录（默认 {OUTPUT_DIR}）'
    )
    parser.add_argument(
        '--random_state', type=int, default=RANDOM_STATE,
        help=f'随机种子（默认 {RANDOM_STATE}）'
    )
    return parser.parse_args()


# ==================== 主逻辑 ====================
def main():
    args = parse_args()

    abx_cols = args.columns

    print(f"将处理以下抗生素列: {abx_cols}")

    # 读取 group 分配
    print("读取 group 分配...")
    df_group = pd.read_csv(args.group_file, sep='\t')
    # 列：Genome_name, group_id

    # 读取 ST/耐药表型
    print("读取 ST/耐药表型...")
    df_st = pd.read_csv(args.st_file, sep='\t', low_memory=False)
    df_st.columns = df_st.columns.str.strip()
    print(df_st.columns.tolist())

    # 检查指定列是否存在
    missing_cols = [c for c in abx_cols if c not in df_st.columns]
    if missing_cols:
        print(f"错误：以下列在 ST 文件中不存在: {missing_cols}")
        print(f"ST 文件可用列: {df_st.columns.tolist()}")
        sys.exit(1)

    # 确保 MLST 列存在
    if 'MLST' not in df_st.columns:
        print("错误：ST 文件中未找到 'MLST' 列。")
        sys.exit(1)

    # 合并：只保留需要的列，减少内存
    # 需要 df_group 的 Genome_name, group_id
    # 需要 df_st 的 Genome_name, MLST, 以及所有要处理的 abx_cols
    st_cols_needed = ['Genome_name', 'MLST'] + abx_cols
    df_st_sub = df_st[st_cols_needed].copy()

    df = df_group.merge(df_st_sub, on='Genome_name', how='inner')
    print(f"合并后样本数: {len(df)}")

    # 遍历每个抗生素
    for abx in abx_cols:
        print(f"\n{'='*60}")
        print(f"处理抗生素: {abx}")
        print(f"{'='*60}")

        # 只保留该列标记为 R 或 S 的样本
        df_abx = df[df[abx].isin(['R', 'S'])].copy()
        if len(df_abx) == 0:
            print("  无 R/S 样本，跳过。")
            continue

        # 标签：R=1, S=0
        y = (df_abx[abx] == 'R').astype(int).values
        groups = df_abx['group_id'].values

        n_R = int((y == 1).sum())
        n_S = int((y == 0).sum())
        n_groups = df_abx['group_id'].nunique()
        print(f"  样本数: {len(df_abx)}, R: {n_R}, S: {n_S}, groups: {n_groups}")

        if n_R < 5 or n_S < 5 or n_groups < 5:
            print("  R/S 或 group 数量不足 5，跳过。")
            continue

        # ==================== 第一步：固定 test（20%） ====================
        try:
            sgkf5 = StratifiedGroupKFold(n_splits=5, shuffle=True,
                                         random_state=args.random_state)
            for train_val_idx, test_idx in sgkf5.split(df_abx, y, groups):
                break
        except Exception as e:
            print(f"  5 折划分失败: {e}")
            continue

        df_test = df_abx.iloc[test_idx].reset_index(drop=True)
        df_train_val = df_abx.iloc[train_val_idx].reset_index(drop=True)

        y_test = y[test_idx]
        y_tv = y[train_val_idx]
        groups_tv = groups[train_val_idx]

        test_R = int((y_test == 1).sum())
        test_S = int((y_test == 0).sum())
        print(f"  Test: {len(df_test)}  (R={test_R}, S={test_S})")

        # 只保留需要的列输出
        out_cols = OUTPUT_COLS_BASE + [abx]
        df_test_out = df_test[out_cols]
        test_file = os.path.join(args.output_dir, f'{abx}_test.txt')
        df_test_out.to_csv(test_file, sep='\t', index=False)

        # ==================== 第二步：10 折交叉验证 ====================
        n_R_tv = int((y_tv == 1).sum())
        n_S_tv = int((y_tv == 0).sum())
        max_splits = min(10, n_R_tv, n_S_tv)

        if max_splits < 2:
            print(f"  train_val 中 R={n_R_tv}, S={n_S_tv}，无法交叉验证，跳过。")
            continue

        if max_splits < 10:
            print(f"  警告: train_val 中最小类别数 {max_splits} < 10，"
                  f"实际使用 {max_splits} 折。")
        n_splits = max_splits

        try:
            sgkf10 = StratifiedGroupKFold(n_splits=n_splits, shuffle=True,
                                          random_state=args.random_state)
            splits = list(sgkf10.split(df_train_val, y_tv, groups_tv))
        except Exception as e:
            print(f"  {n_splits} 折划分失败: {e}")
            continue

        for fold, (train_idx, val_idx) in enumerate(splits, 1):
            df_train = df_train_val.iloc[train_idx].reset_index(drop=True)
            df_val = df_train_val.iloc[val_idx].reset_index(drop=True)

            y_train = y_tv[train_idx]
            y_val = y_tv[val_idx]
            train_R = int((y_train == 1).sum())
            train_S = int((y_train == 0).sum())
            val_R = int((y_val == 1).sum())
            val_S = int((y_val == 0).sum())

            print(f"  Fold {fold:2d}: "
                  f"Train {len(df_train):5d} (R={train_R:4d}, S={train_S:4d})  "
                  f"Val {len(df_val):4d} (R={val_R:3d}, S={val_S:3d})")

            # 只保留需要的列输出
            df_train_out = df_train[out_cols]
            df_val_out = df_val[out_cols]

            train_file = os.path.join(args.output_dir, f'{abx}_fold{fold:02d}_train.txt')
            val_file = os.path.join(args.output_dir, f'{abx}_fold{fold:02d}_val.txt')
            df_train_out.to_csv(train_file, sep='\t', index=False)
            df_val_out.to_csv(val_file, sep='\t', index=False)

    print(f"\n全部完成。结果保存在目录: {args.output_dir}")


if __name__ == '__main__':
    main()
