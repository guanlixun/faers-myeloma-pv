# -*- coding: utf-8 -*-
"""
步骤 01：下载并解压 FAERS 季度 ASCII 数据包
=====================================================
执行说明（低成本 AI）：
- 修改 START_YEAR / END_YEAR 后运行：python src/01_download.py
- 数据存放于 data/raw/；解压后的 txt 文件在 data/extracted/
- 网络受限时可手动从 FDA 官网下载 zip 放入 data/raw/ 后运行 --extract-only
"""
import argparse
import io
import re
import zipfile
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
EXTRACTED = ROOT / "data" / "extracted"

START_YEAR, END_YEAR = 2021, 2026   # 覆盖五药上市后完整暴露期；2026 最新可用季度为 Q2（Q3/Q4 尚未发布）
QUARTERS = [1, 2, 3, 4]
# 末端季度截断：END_YEAR 实际可用季度数（1=Q1..4=Q4）。2026-09 探测仅 Q1/Q2 可用。
LAST_YEAR_QUARTERS = 2

# FDA 公开下载页（页面结构可能变动，失败时改为手动下载）
LIST_URL = "https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-latest-quarterly-data-files"


def candidate_urls(year: int, quarter: int):
    """FAERS 季度包候选 URL（FDA 历史上使用过的命名）。"""
    return [
        f"https://fis.fda.gov/content/Exports/faers_ascii_{year}q{quarter}.zip",
        f"https://fis.fda.gov/content/Exports/faers_ascii_{year}Q{quarter}.zip",
    ]


def download_quarter(year: int, quarter: int) -> Path | None:
    """流式下载单个季度 ASCII 包（分块写盘，避免大文件一次性读入内存）。"""
    dest = RAW / f"faers_ascii_{year}q{quarter}.zip"
    if dest.exists() and dest.stat().st_size > 1_000_000:
        print(f"[skip] 已存在 {dest.name} ({dest.stat().st_size/2**20:.1f} MB)")
        return dest
    for url in candidate_urls(year, quarter):
        try:
            with requests.get(url, timeout=(30, 120), stream=True) as r:
                if r.status_code != 200:
                    print(f"[fail] {url} HTTP {r.status_code}")
                    continue
                total = int(r.headers.get("Content-Length", 0))
                tmp = dest.with_suffix(".part")
                size = 0
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        fh.write(chunk)
                        size += len(chunk)
                if size < 1_000_000 or tmp.read_bytes()[:2] != b"PK":
                    print(f"[fail] {url} 内容异常（{size} bytes），跳过")
                    tmp.unlink(missing_ok=True)
                    continue
                tmp.rename(dest)
                print(f"[ok] {url} -> {dest.name} ({size/2**20:.1f} MB)")
                return dest
        except requests.RequestException as e:
            print(f"[fail] {url}: {type(e).__name__}: {e}")
    print(f"[manual] {year}Q{quarter} 自动下载失败，请手动下载后放入 {RAW}")
    return None


def extract_all():
    for z in sorted(RAW.glob("*.zip")):
        target = EXTRACTED / z.stem
        if target.exists():
            continue
        target.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(z) as zf:
            zf.extractall(target)
        print(f"[extract] {z.name} -> {target}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--extract-only", action="store_true")
    args = p.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    EXTRACTED.mkdir(parents=True, exist_ok=True)
    if not args.extract_only:
        for y in range(START_YEAR, END_YEAR + 1):
            qs = QUARTERS if y < END_YEAR else range(1, LAST_YEAR_QUARTERS + 1)
            for q in qs:
                download_quarter(y, q)
    extract_all()
    # 完整性检查：每季度应有 DEMO/DRUG/REAC/OUTC/INDI/THER 等文件
    for d in sorted(EXTRACTED.iterdir()):
        if d.is_dir():
            files = {f.stem[:4].upper() for f in d.glob("*.txt")}
            missing = {"DEMO", "DRUG", "REAC", "OUTC", "INDI", "THER"} - files
            print(f"{d.name}: {'OK' if not missing else '缺 ' + ','.join(missing)}")


if __name__ == "__main__":
    main()
