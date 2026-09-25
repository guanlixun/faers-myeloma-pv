# -*- coding: utf-8 -*-
"""
步骤 03：构建分析数据集（暴露 × 事件 长表）
=====================================================
执行说明（低成本 AI）：
- 前置：02_clean.py 已完成
- 运行：python src/03_build_dataset.py
- 输出：data/processed/analysis_dataset.parquet
  每行 = 一个 CASE × 一个目标药物 × 一个事件类别（命中为 1）
- 关键步骤已内嵌"未匹配药物名探查"报告：
  运行后检查 output/unmatched_drugs.csv，人工确认是否需补充
  config/drug_mapping.py 的 VARIANTS
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.drug_mapping import DRUGS, VARIANT_INDEX, GROUPS
from config.event_definitions import EVENTS, DEATH_OUTCOME_CODE

import os as _os
ROOT = Path(__file__).resolve().parent.parent
_PILOT = _os.environ.get("FAERS_PILOT") == "1"
PROCESSED = (ROOT / "data" / "pilot" / "processed") if _PILOT else (ROOT / "data" / "processed")
OUT = (ROOT / "data" / "pilot" / "output") if _PILOT else (ROOT / "output")

MM_INDICATION_KEYWORDS = ["MYELOMA"]   # INDI_PT 含此关键词视为 MM 适应症


def load(name):
    return pd.read_csv(PROCESSED / f"{name}.csv", dtype=str)


def normalize_name(s: pd.Series) -> pd.Series:
    return s.fillna("").str.upper().str.strip()


def flag_exposures(drug: pd.DataFrame) -> pd.DataFrame:
    """限定 PS（primary suspect）+ 药物名匹配 → CASE × 药物 暴露表。"""
    drug = drug.copy()
    drug["DRUGNAME_N"] = normalize_name(drug["DRUGNAME"])
    if "PROD_AI" in drug.columns:  # 活性成分列可辅助匹配
        drug["PROD_AI_N"] = normalize_name(drug["PROD_AI"])
    else:
        drug["PROD_AI_N"] = ""
    drug["GENERIC"] = drug["DRUGNAME_N"].map(VARIANT_INDEX)
    mask = drug["GENERIC"].isna() & drug["PROD_AI_N"].ne("")
    drug.loc[mask, "GENERIC"] = drug.loc[mask, "PROD_AI_N"].map(VARIANT_INDEX)
    drug = drug[drug["GENERIC"].notna()]
    if "ROLE_COD" in drug.columns:
        drug = drug[drug["ROLE_COD"].eq("PS")]
    key = "PRIMARYID" if "PRIMARYID" in drug.columns else "CASEID"
    exp = drug[[key, "GENERIC"]].drop_duplicates()
    exp["GROUP"] = exp["GENERIC"].map({g: i["group"] for g, i in DRUGS.items()})
    return exp.rename(columns={key: "CASE_KEY"})


def report_unmatched(drug: pd.DataFrame):
    """导出高频但未匹配的药物名，供人工核对映射表。"""
    OUT.mkdir(parents=True, exist_ok=True)
    cand = drug[normalize_name(drug["DRUGNAME"]).str.contains(
        "TECLIST|ELRANAT|TALQUET|IDECAB|CILTACAB|TECVAYLI|ELREXFIO|TALVEY|ABECMA|CARVYKT",
        regex=True, na=False)]
    matched = cand["GENERIC"].notna() if "GENERIC" in cand.columns else pd.Series(False, index=cand.index)
    counts = (normalize_name(cand.loc[~matched, "DRUGNAME"]).value_counts().head(100))
    counts.to_csv(OUT / "unmatched_drugs.csv", header=["count"])
    print(f"[check] 候选药物名频数已导出 output/unmatched_drugs.csv，请人工核对映射表")


def flag_events(reac: pd.DataFrame, outc: pd.DataFrame, indi: pd.DataFrame) -> pd.DataFrame:
    """CASE × 事件类别 命中表。"""
    key = "PRIMARYID" if "PRIMARYID" in reac.columns else "CASEID"
    reac = reac.copy()
    reac["PT_N"] = normalize_name(reac["PT"])
    rows = []
    for ev, info in EVENTS.items():
        pts = {p.upper() for p in info.get("pts", [])}
        hit = reac[reac["PT_N"].isin(pts)][[key]].copy()
        hit["EVENT"] = ev
        rows.append(hit)
        # SOC 层级事件（如全部感染）：FAERS REAC 无 SOC 列，
        # 需 MedDRA 层级文件映射 PT→SOC；若无法获得 MedDRA，
        # 用扩展 PT 清单近似并在论文局限性中声明
        if info.get("soc"):
            print(f"[note] {ev} 定义为 SOC 层级，当前仅用 PT 清单近似，"
                  f"建议人工补充 MedDRA PT→SOC 映射")
    ev_df = pd.concat(rows, ignore_index=True).drop_duplicates()
    # 死亡结局
    okey = "PRIMARYID" if "PRIMARYID" in outc.columns else "CASEID"
    death = outc[outc["OUTC_COD"].eq(DEATH_OUTCOME_CODE)][[okey]].copy()
    death.columns = ["CASE_KEY"]
    death["EVENT"] = "DEATH"
    ev_df = pd.concat([ev_df.rename(columns={key: "CASE_KEY"}), death], ignore_index=True)
    # MM 适应症标记
    ikey = "PRIMARYID" if "PRIMARYID" in indi.columns else "CASEID"
    mm_keys = set(indi.loc[normalize_name(indi["INDI_PT"]).str.contains(
        "|".join(MM_INDICATION_KEYWORDS), na=False), ikey])
    ev_df["MM_INDICATION"] = ev_df["CASE_KEY"].isin(mm_keys)
    return ev_df.drop_duplicates()


def main():
    drug, reac, outc, indi = load("drug"), load("reac"), load("outc"), load("indi")
    exp = flag_exposures(drug)
    report_unmatched(drug)
    ev = flag_events(reac, outc, indi)
    # 暴露集：全部目标药物 CASE×药物（不要求命中目标事件），供 04/06 做分母
    exp.to_csv(PROCESSED / "exposure.csv", index=False)
    print(f"[save] exposure.csv: {len(exp):,} 行, {exp['CASE_KEY'].nunique():,} 个 CASE")
    # 事件命中表全量落盘：CASE×事件，含 MM 标记
    ev.to_csv(PROCESSED / "event_hits.csv", index=False)
    print(f"[save] event_hits.csv: {len(ev):,} 行")
    # 分析数据集：暴露 ∩ 事件（用于直观检查，ROR 统计不再依赖它）
    ds = exp.merge(ev, on="CASE_KEY", how="inner")
    ds.to_csv(PROCESSED / "analysis_dataset.csv", index=False)
    print(f"[save] analysis_dataset.csv: {len(ds):,} 行, "
          f"{ds['CASE_KEY'].nunique():,} 个 CASE")
    print(exp.groupby(["GROUP", "GENERIC"])["CASE_KEY"].nunique())


if __name__ == "__main__":
    main()
