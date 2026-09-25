# -*- coding: utf-8 -*-
"""
步骤 06：敏感性分析（三项，对应 SAP）
=====================================================
执行说明（低成本 AI）：
- 前置：04_ror_analysis.py 已完成
- 运行：python src/06_sensitivity.py
- 输出：output/suppl_table2_sensitivity.csv
- 三项分析：
  S1 仅 MM 适应症记录（MM_INDICATION == True）重复主分析
  S2 排除含合并用药（concomitant）的 CASE（DRUG 表中同 CASE 存在
     ROLE_COD 为 C/SS/I 的记录）后重复主分析
  S3 按报告年份分层（上市早期 vs 近期）比较信号稳定性，
     用于讨论 notoriety bias / Weber effect
- 判读标准（写入论文）：信号方向与主分析一致即视为稳健，
  不要求 CI 完全重叠
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.drug_mapping import DRUGS
from config.event_definitions import EVENTS

import os as _os
ROOT = Path(__file__).resolve().parent.parent
_PILOT = _os.environ.get("FAERS_PILOT") == "1"
PROCESSED = (ROOT / "data" / "pilot" / "processed") if _PILOT else (ROOT / "data" / "processed")
OUT = (ROOT / "data" / "pilot" / "output") if _PILOT else (ROOT / "output")

# 注：ROR 函数与 04_ror_analysis.py 保持一致；如主分析修改公式，此处必须同步
def ror_ci(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    ror = (a / c) / (b / d)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return ror, np.exp(np.log(ror) - 1.96 * se), np.exp(np.log(ror) + 1.96 * se)


def run_ror_subset(exp_sub, ev_sub, all_cases, label, events):
    """对给定暴露子集重复 ROR 主分析（背景=全库排除目标药）。"""
    drug_cases = exp_sub[["CASE_KEY", "GENERIC"]].drop_duplicates()
    target_cases = set(drug_cases["CASE_KEY"])
    bg_hits = ev_sub[ev_sub["CASE_KEY"].isin(all_cases - target_cases)].groupby(
        ["CASE_KEY", "EVENT"]).size().unstack(fill_value=0)
    rows = []
    for generic in DRUGS:
        g_cases = set(drug_cases.loc[drug_cases["GENERIC"] == generic, "CASE_KEY"])
        g_hits = ev_sub[ev_sub["CASE_KEY"].isin(g_cases)].groupby(
            ["CASE_KEY", "EVENT"]).size().unstack(fill_value=0)
        for evn in events:
            a = int(g_hits[evn].gt(0).sum()) if evn in g_hits.columns else 0
            b = len(g_cases) - a
            c = int(bg_hits[evn].gt(0).sum()) if evn in bg_hits.columns else 0
            d = len(all_cases - target_cases) - c
            ror, lo, hi = ror_ci(a, b, c, d)
            rows.append({"analysis": label, "drug": generic, "event": evn,
                         "a": a, "ROR": round(ror, 2),
                         "CI_low": round(lo, 2), "CI_high": round(hi, 2),
                         "signal": bool(lo > 1 and a >= 3)})
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    exp = pd.read_csv(PROCESSED / "exposure.csv", dtype=str)
    ev = pd.read_csv(PROCESSED / "event_hits.csv", dtype=str)
    demo = pd.read_csv(PROCESSED / "demo.csv", dtype=str)
    all_cases = set(demo["PRIMARYID" if "PRIMARYID" in demo.columns else "CASEID"])
    ev["CASE_KEY"] = ev["CASE_KEY"].astype(str)

    # 暴露 CASE 亦需限定在 demo 有效 ID 内
    exp = exp[exp["CASE_KEY"].isin(all_cases)]
    events = list(EVENTS.keys()) + ["DEATH"]

    # MM 标记并入暴露表
    mm_map = ev[["CASE_KEY", "MM_INDICATION"]].drop_duplicates().set_index("CASE_KEY")["MM_INDICATION"]
    exp["MM_INDICATION"] = exp["CASE_KEY"].map(mm_map).fillna(False)

    results = []

    # S1：仅 MM 适应症（MM_INDICATION 为字符串 'True'/'False'）
    ds_mm = exp[exp["MM_INDICATION"].eq("True")]
    results += run_ror_subset(ds_mm, ev, all_cases, "S1_MM_indication_only", events)

    # S2：排除有合并用药的 CASE（DRUG 表中 ROLE_COD 含 C/SS/I 的 CASE）
    drug = pd.read_csv(PROCESSED / "drug.csv", dtype=str)
    dkey = "PRIMARYID" if "PRIMARYID" in drug.columns else "CASEID"
    if "ROLE_COD" in drug.columns:
        concom = set(drug.loc[drug["ROLE_COD"].isin(["C", "SS", "I"]), dkey])
    else:
        concom = set()
    exp_mono = exp[~exp["CASE_KEY"].isin(concom)]
    results += run_ror_subset(exp_mono, ev, all_cases, "S2_no_concomitant", events)

    # S3: 按报告年份分层（FDA_DT 年份中位数分早/近）
    demo["YEAR"] = pd.to_datetime(demo["FDA_DT"], errors="coerce").dt.year
    year_map = demo.set_index(dkey)["YEAR"]
    exp["YEAR"] = exp["CASE_KEY"].map(year_map)
    med_year = int(exp["YEAR"].median())
    for label, sub in [("S3_early", exp[exp["YEAR"] <= med_year]),
                       ("S3_recent", exp[exp["YEAR"] > med_year])]:
        results += run_ror_subset(sub, ev, all_cases, label, events)

    out = pd.DataFrame(results)
    out.to_csv(OUT / "suppl_table2_sensitivity.csv", index=False)
    # 稳健性汇总：主分析为信号的事件，在各敏感性分析中是否仍为信号
    main_sig = pd.read_csv(OUT / "table2_ror_matrix.csv")
    main_sig = main_sig[main_sig["signal"]][["drug", "event"]].assign(main_signal=True)
    chk = out.merge(main_sig, on=["drug", "event"], how="inner")
    summary = (chk.groupby("analysis")["signal"].agg(["sum", "count"])
               .rename(columns={"sum": "信号保持数", "count": "主分析信号数"}))
    print(summary)
    print(f"[save] suppl_table2_sensitivity.csv ({len(out)} 行)")


if __name__ == "__main__":
    main()
