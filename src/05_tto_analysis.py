# -*- coding: utf-8 -*-
"""
步骤 05：Time-to-Onset（TTO）时序分析 + Weibull 分布检验
=====================================================
执行说明（低成本 AI）：
- 前置：03_build_dataset.py 已完成
- 运行：python src/05_tto_analysis.py
- 输出：output/table4_tto.csv（各药物 × 事件的中位 TTO 与 Weibull 参数）
        output/figure3_tto.png（TTO 分布图）
- TTO = EVENT_DT（DEMO 表）− START_DT（THER 表），单位天
- Weibull 形状参数 β 解读（写入论文）：
    β < 1（95%CI 不含 1）→ 早期失效型（用药后早期高发，随时间下降）
    β ≈ 1 → 随机型（恒定发生）
    β > 1 → 磨损失效型（随暴露时间增加）
- 日期清洗规则：仅保留 0 <= TTO <= 1095（3 年）的记录；
  异常值比例 > 5% 时在论文局限性中声明
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.drug_mapping import DRUGS
from config.event_definitions import EVENTS

import os as _os
ROOT = Path(__file__).resolve().parent.parent
_PILOT = _os.environ.get("FAERS_PILOT") == "1"
PROCESSED = (ROOT / "data" / "pilot" / "processed") if _PILOT else (ROOT / "data" / "processed")
OUT = (ROOT / "data" / "pilot" / "output") if _PILOT else (ROOT / "output")

MAX_TTO_DAYS = 1095


def parse_date(s: pd.Series) -> pd.Series:
    """FAERS 日期为 YYYYMMDD 或部分缺失（YYYYMM/YYYY），仅保留完整 8 位。"""
    s = s.fillna("").str.strip()
    full = s.str.len().eq(8)
    dt = pd.to_datetime(s.where(full), format="%Y%m%d", errors="coerce")
    return dt


def weibull_beta_ci(t: np.ndarray, n_boot: int = 500, seed: int = 42):
    """Weibull 形状参数 β 及 bootstrap 95% CI。"""
    if len(t) < 10:
        return np.nan, np.nan, np.nan
    shape, _, scale = stats.weibull_min.fit(t, floc=0)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        sample = rng.choice(t, size=len(t), replace=True)
        try:
            b, _, _ = stats.weibull_min.fit(sample, floc=0)
            boots.append(b)
        except Exception:
            pass
    if not boots:
        return shape, np.nan, np.nan
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return shape, lo, hi


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    exp = pd.read_csv(PROCESSED / "exposure.csv", dtype=str)
    ev = pd.read_csv(PROCESSED / "event_hits.csv", dtype=str).drop(columns=["MM_INDICATION"], errors="ignore")
    demo = pd.read_csv(PROCESSED / "demo.csv", dtype=str)
    ther = pd.read_csv(PROCESSED / "ther.csv", dtype=str)
    drug = pd.read_csv(PROCESSED / "drug.csv", dtype=str)

    key = "PRIMARYID" if "PRIMARYID" in demo.columns else "CASEID"
    demo = demo[[key, "EVENT_DT"]].copy()
    demo["EVENT_DATE"] = parse_date(demo["EVENT_DT"])

    # THER 表取每个 CASE × 药物 的最早 START_DT
    ther = ther.copy()
    ther["START_DATE"] = parse_date(ther["START_DT"])
    drug_n = drug.copy()
    drug_n["DRUGNAME_N"] = drug_n["DRUGNAME"].fillna("").str.upper().str.strip()
    from config.drug_mapping import VARIANT_INDEX
    drug_n["GENERIC"] = drug_n["DRUGNAME_N"].map(VARIANT_INDEX)
    if "PROD_AI" in drug_n.columns:
        drug_n["PROD_AI_N"] = drug_n["PROD_AI"].fillna("").str.upper().str.strip()
        mask = drug_n["GENERIC"].isna() & drug_n["PROD_AI_N"].ne("")
        drug_n.loc[mask, "GENERIC"] = drug_n.loc[mask, "PROD_AI_N"].map(VARIANT_INDEX)
    drug_n = drug_n[drug_n["GENERIC"].notna()]
    if "ROLE_COD" in drug_n.columns:
        drug_n = drug_n[drug_n["ROLE_COD"].eq("PS")]
    dkey = "PRIMARYID" if "PRIMARYID" in drug_n.columns else "CASEID"
    tkey = "PRIMARYID" if "PRIMARYID" in ther.columns else "CASEID"
    start = (drug_n[[dkey, "GENERIC"]].dropna()
             .merge(ther[[tkey, "START_DATE"]].dropna(), left_on=dkey, right_on=tkey)
             .groupby([dkey, "GENERIC"])["START_DATE"].min().reset_index())

    ev["CASE_KEY"] = ev["CASE_KEY"].astype(str)
    exp["CASE_KEY"] = exp["CASE_KEY"].astype(str)
    # TTO 基础 = 暴露表（CASE×GENERIC）× 事件命中（CASE×EVENT）× 事件日期 × 用药起始日
    tto = (exp.merge(ev[["CASE_KEY", "EVENT"]], on="CASE_KEY", how="inner")
             .merge(demo[[key, "EVENT_DATE"]], left_on="CASE_KEY", right_on=key, how="left")
             .merge(start, left_on=["CASE_KEY", "GENERIC"], right_on=[dkey, "GENERIC"], how="left"))
    tto["TTO"] = (tto["EVENT_DATE"] - tto["START_DATE"]).dt.days
    tto = tto[tto["TTO"].between(0, MAX_TTO_DAYS)]
    print(f"TTO 有效记录: {len(tto):,}")

    rows = []
    for generic in DRUGS:
        for evn in list(EVENTS.keys()) + ["DEATH"]:
            t = tto.loc[(tto["GENERIC"] == generic) & (tto["EVENT"] == evn), "TTO"].dropna().values
            if len(t) < 10:
                continue
            beta, lo, hi = weibull_beta_ci(t)
            rows.append({
                "drug": generic, "event": evn, "n": len(t),
                "median_TTO_days": int(np.median(t)),
                "IQR": f"{int(np.percentile(t,25))}-{int(np.percentile(t,75))}",
                "weibull_beta": round(beta, 2),
                "beta_CI": f"{lo:.2f}-{hi:.2f}",
                "failure_type": ("early" if hi < 1 else "wear-out" if lo > 1 else "random"),
            })
    t4 = pd.DataFrame(rows)
    t4.to_csv(OUT / "table4_tto.csv", index=False)
    print(f"[save] table4_tto.csv ({len(t4)} 行)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        focus = ["CRS", "ICANS_NEUROTOXICITY", "DELAYED_NEUROTOXICITY", "INFECTION_ALL"]
        fig, axes = plt.subplots(1, len(focus), figsize=(4 * len(focus), 4))
        for ax, ev in zip(axes, focus):
            for generic in DRUGS:
                t = tto.loc[(tto["GENERIC"] == generic) & (tto["EVENT"] == ev), "TTO"].dropna()
                if len(t) >= 10:
                    ax.hist(t, bins=30, histtype="step", density=True, label=generic[:12])
            ax.set_title(ev, fontsize=9)
            ax.set_xlabel("days")
        axes[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(OUT / "figure3_tto.png", dpi=200, bbox_inches="tight")
        print("[save] figure3_tto.png")
    except Exception as e:
        print(f"[warn] TTO 图生成失败: {e}")


if __name__ == "__main__":
    main()
