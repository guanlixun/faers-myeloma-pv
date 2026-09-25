# -*- coding: utf-8 -*-
"""
步骤 04：不成比例分析（ROR）主分析 + 头对头比较
=====================================================
执行说明（低成本 AI）：
- 前置：03_build_dataset.py 已完成
- 运行：python src/04_ror_analysis.py
- 输出：output/table2_ror_matrix.csv（主结果表）
        output/table3_headtohead.csv（双抗组 vs CAR-T 组）
        output/figure2_forest.png（森林图）
- 2×2 表定义（以全库其他所有药物为背景）：
    a = 目标药物 & 目标事件报告数
    b = 目标药物 & 其他事件报告数
    c = 其他所有药物 & 目标事件报告数
    d = 其他所有药物 & 其他事件报告数
  ROR = (a/c)/(b/d)，信号阈值：95%CI 下限 > 1 且 a >= 3
- ⚠️ b/d 需要全库背景计数，本脚本基于清洗后的全量 DRUG/REAC 计算
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.drug_mapping import DRUGS, GROUPS
from config.event_definitions import EVENTS

import os as _os
ROOT = Path(__file__).resolve().parent.parent
_PILOT = _os.environ.get("FAERS_PILOT") == "1"
PROCESSED = (ROOT / "data" / "pilot" / "processed") if _PILOT else (ROOT / "data" / "processed")
OUT = (ROOT / "data" / "pilot" / "output") if _PILOT else (ROOT / "output")

MIN_A = 3          # 信号最小报告数
ALPHA = 0.05


def ror_ci(a, b, c, d):
    """ROR 及 95% CI（对数法）。任一格子为 0 时加 0.5 校正。"""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    ror = (a / c) / (b / d)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo, hi = np.exp(np.log(ror) - 1.96 * se), np.exp(np.log(ror) + 1.96 * se)
    return ror, lo, hi


def load(name):
    return pd.read_csv(PROCESSED / f"{name}.csv", dtype=str)


def build_background():
    """全库 CASE × 事件 背景表（所有药物，不止目标药物）。"""
    reac = load("reac")
    key = "PRIMARYID" if "PRIMARYID" in reac.columns else "CASEID"
    reac = reac.copy()
    reac["PT_N"] = reac["PT"].fillna("").str.upper().str.strip()
    frames = []
    for ev, info in EVENTS.items():
        pts = {p.upper() for p in info.get("pts", [])}
        hit = reac[reac["PT_N"].isin(pts)][[key]].copy()
        hit["EVENT"] = ev
        frames.append(hit)
    outc = load("outc")
    okey = "PRIMARYID" if "PRIMARYID" in outc.columns else "CASEID"
    death = outc[outc["OUTC_COD"].eq("DE")][[okey]].copy()
    death.columns = ["CASE_KEY"]
    death["EVENT"] = "DEATH"
    bg = pd.concat(frames, ignore_index=True).rename(columns={key: "CASE_KEY"})
    bg = pd.concat([bg, death], ignore_index=True).drop_duplicates()
    return bg


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    exp = pd.read_csv(PROCESSED / "exposure.csv", dtype=str)
    ev = pd.read_csv(PROCESSED / "event_hits.csv", dtype=str)
    demo = load("demo")
    all_cases = set(demo["PRIMARYID" if "PRIMARYID" in demo.columns else "CASEID"])
    # 目标药物暴露 CASE（分母基础）
    drug_cases = exp[["CASE_KEY", "GENERIC", "GROUP"]].drop_duplicates()
    target_cases = set(drug_cases["CASE_KEY"])
    # CASE × 事件命中矩阵
    ev["CASE_KEY"] = ev["CASE_KEY"].astype(str)
    bg_hits = ev[ev["CASE_KEY"].isin(all_cases - target_cases)].groupby(["CASE_KEY", "EVENT"]).size().unstack(fill_value=0)

    results = []
    events = list(EVENTS.keys()) + ["DEATH"]
    for generic in DRUGS:
        g_cases = set(drug_cases.loc[drug_cases["GENERIC"] == generic, "CASE_KEY"])
        g_hits = ev[ev["CASE_KEY"].isin(g_cases)].groupby(["CASE_KEY", "EVENT"]).size().unstack(fill_value=0)
        for evn in events:
            a = int(g_hits[evn].gt(0).sum()) if evn in g_hits.columns else 0
            b = len(g_cases) - a
            c = int(bg_hits[evn].gt(0).sum()) if evn in bg_hits.columns else 0
            d = len(all_cases - target_cases) - c
            ror, lo, hi = ror_ci(a, b, c, d)
            results.append({
                "drug": generic, "group": DRUGS[generic]["group"], "event": evn,
                "a": a, "b": b, "c": c, "d": d,
                "ROR": round(ror, 2), "CI_low": round(lo, 2), "CI_high": round(hi, 2),
                "signal": bool(lo > 1 and a >= MIN_A),
            })
    t2 = pd.DataFrame(results)
    t2.to_csv(OUT / "table2_ror_matrix.csv", index=False)
    print(f"[save] table2_ror_matrix.csv ({t2['signal'].sum()} 个信号)")

    # ---- Table 1：按药物基线特征（描述性） ----
    demo = demo.copy()
    dkey = "PRIMARYID" if "PRIMARYID" in demo.columns else "CASEID"
    demo[dkey] = demo[dkey].astype(str)
    demo["YEAR"] = pd.to_datetime(demo["FDA_DT"], errors="coerce").dt.year
    demo["AGE_N"] = pd.to_numeric(demo["AGE"], errors="coerce")
    t1_rows = []
    for generic in DRUGS:
        g_cases = set(drug_cases.loc[drug_cases["GENERIC"] == generic, "CASE_KEY"])
        sub = demo[demo[dkey].isin(g_cases)]
        t1_rows.append({
            "drug": generic, "n_cases": len(g_cases),
            "median_age": round(float(sub["AGE_N"].median()), 1) if sub["AGE_N"].notna().any() else None,
            "pct_female": round(100 * sub["SEX"].eq("F").mean(), 1) if "SEX" in sub else None,
            "pct_male": round(100 * sub["SEX"].eq("M").mean(), 1) if "SEX" in sub else None,
            "median_report_year": round(float(sub["YEAR"].median()), 1) if sub["YEAR"].notna().any() else None,
            "top_country": sub["REPORTER_COUNTRY"].value_counts().index[0] if "REPORTER_COUNTRY" in sub and len(sub) else None,
        })
    t1 = pd.DataFrame(t1_rows)
    t1.to_csv(OUT / "table1_baseline.csv", index=False)
    print("[save] table1_baseline.csv")

    # 头对头：双抗组 vs CAR-T 组，比较各事件报告比例（卡方 + 比例差）
    rows = []
    for evn in events:
        for gname, members in GROUPS.items():
            g_cases = set(drug_cases.loc[drug_cases["GROUP"] == gname, "CASE_KEY"])
            g_hits = ev[ev["CASE_KEY"].isin(g_cases)].groupby(["CASE_KEY", "EVENT"]).size().unstack(fill_value=0)
            n_event = int(g_hits[evn].gt(0).sum()) if evn in g_hits.columns else 0
            rows.append({"event": evn, "group": gname, "n_cases": len(g_cases),
                         "n_event": n_event, "pct": round(100 * n_event / max(len(g_cases), 1), 2)})
        b_row, c_row = rows[-2], rows[-1]
        table = np.array([[b_row["n_event"], b_row["n_cases"] - b_row["n_event"]],
                          [c_row["n_event"], c_row["n_cases"] - c_row["n_event"]]])
        if table.min() >= 0 and table.sum() > 0:
            chi2, p, _, _ = stats.chi2_contingency(table, correction=True)
            rows.append({"event": evn, "group": "CHI2_P", "n_cases": None,
                         "n_event": None, "pct": round(p, 4)})
    t3 = pd.DataFrame(rows)
    t3.to_csv(OUT / "table3_headtohead.csv", index=False)
    print(f"[save] table3_headtohead.csv")

    # 森林图（主事件 × 5 药）
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        main_events = ["CRS", "ICANS_NEUROTOXICITY", "INFECTION_ALL",
                       "GPRC5D_OFFTUMOR", "DELAYED_NEUROTOXICITY", "DEATH"]
        fig, axes = plt.subplots(1, len(main_events), figsize=(4 * len(main_events), 5), sharey=True)
        for ax, ev in zip(axes, main_events):
            sub = t2[t2["event"] == ev]
            y = np.arange(len(sub))
            ax.errorbar(sub["ROR"], y,
                        xerr=[sub["ROR"] - sub["CI_low"], sub["CI_high"] - sub["ROR"]],
                        fmt="o", capsize=3)
            ax.axvline(1, ls="--", c="grey")
            ax.set_xscale("log")
            ax.set_title(ev, fontsize=9)
            ax.set_yticks(y, [d[:18] for d in sub["drug"]], fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / "figure2_forest.png", dpi=200, bbox_inches="tight")
        print("[save] figure2_forest.png")
    except Exception as e:
        print(f"[warn] 森林图生成失败: {e}")


if __name__ == "__main__":
    main()
