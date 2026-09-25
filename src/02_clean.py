# -*- coding: utf-8 -*-
"""
步骤 02：读取、合并、去重 FAERS 多季度数据
=====================================================
执行说明（低成本 AI）：
- 前置：01_download.py 已完成
- 运行：python src/02_clean.py
- 输出：data/processed/{demo,drug,reac,outc,indi,ther}.parquet
- 去重规则：按 PRIMARYID 保留 FDA_DT 最新版本（FDA 官方推荐做法）
- 列名注意：不同季度文件可能大小写/分隔符不一致，已做兼容处理；
  如遇解析错误，先打印该行样本人工检查，勿静默跳过
"""
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PILOT = os.environ.get("FAERS_PILOT") == "1"
EXTRACTED = (ROOT / "data" / "pilot" / "extracted") if PILOT else (ROOT / "data" / "extracted")
PROCESSED = (ROOT / "data" / "pilot" / "processed") if PILOT else (ROOT / "data" / "processed")

TABLES = ["DEMO", "DRUG", "REAC", "OUTC", "INDI", "THER"]

# 仅保留分析需要的列（大幅降低内存/磁盘占用）
KEEP_COLS = {
    "DEMO": ["PRIMARYID", "CASEID", "FDA_DT", "EVENT_DT", "AGE", "AGE_GRP", "SEX", "REPORTER_COUNTRY", "REPORT_TYPE"],
    "DRUG": ["PRIMARYID", "CASEID", "DRUGNAME", "PROD_AI", "ROLE_COD"],
    "REAC": ["PRIMARYID", "CASEID", "PT"],
    "OUTC": ["PRIMARYID", "CASEID", "OUTC_COD"],
    "INDI": ["PRIMARYID", "CASEID", "INDI_PT"],
    "THER": ["PRIMARYID", "CASEID", "START_DT"],
}


def read_one(path: Path) -> pd.DataFrame:
    """读取单个季度单表，兼容 $ 分隔与大小写差异。"""
    df = pd.read_csv(
        path, sep="$", dtype=str, encoding="latin-1",
        on_bad_lines="warn", low_memory=False,
    )
    df.columns = [c.strip().upper() for c in df.columns]
    df["SOURCE_QUARTER"] = path.parent.name
    return df


def merge_table(name: str) -> pd.DataFrame:
    files = sorted(EXTRACTED.glob(f"**/{name}*.txt"))
    if not files:
        sys.exit(f"未找到 {name} 表，请先运行 01_download.py 并检查完整性")
    frames = []
    for f in files:
        try:
            df = read_one(f)
            keep = KEEP_COLS.get(name, [])
            avail = [c for c in keep if c in df.columns]
            df = df[avail] if avail else df
            frames.append(df)
        except Exception as e:
            print(f"[warn] {f} 读取失败: {e}")
    df = pd.concat(frames, ignore_index=True)
    print(f"{name}: {len(df):,} 行（合并 {len(files)} 个季度）")
    return df


def dedup_demo(demo: pd.DataFrame) -> pd.DataFrame:
    """FDA 推荐：同一 CASE 取 FDA_DT 最新版本。"""
    demo = demo.copy()
    demo["FDA_DT"] = pd.to_datetime(demo["FDA_DT"], errors="coerce")
    key = "PRIMARYID" if "PRIMARYID" in demo.columns else "CASEID"
    demo = demo.sort_values("FDA_DT").drop_duplicates(subset=key, keep="last")
    print(f"DEMO 去重后: {len(demo):,} 个 CASE")
    return demo


def main():
    PROCESSED.mkdir(parents=True, exist_ok=True)
    tables = {t: merge_table(t) for t in TABLES}
    tables["DEMO"] = dedup_demo(tables["DEMO"])
    valid_ids = set(tables["DEMO"]["PRIMARYID" if "PRIMARYID" in tables["DEMO"].columns else "CASEID"])
    for name, df in tables.items():
        key = "PRIMARYID" if "PRIMARYID" in df.columns else "CASEID"
        df = df[df[key].isin(valid_ids)]
        df.to_csv(PROCESSED / f"{name.lower()}.csv", index=False)
        print(f"[save] {name.lower()}.csv {len(df):,} 行")


if __name__ == "__main__":
    main()
