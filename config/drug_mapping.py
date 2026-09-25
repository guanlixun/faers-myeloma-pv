# -*- coding: utf-8 -*-
"""
药物识别映射表 —— FAERS 双抗 vs CAR-T 药物警戒研究
=====================================================
⚠️ 人工核对要点（低成本 AI 执行前必读）：
1. 通用名 / 商品名 / 常见拼写变体必须与 FDA 最新标签逐一核对
2. FAERS 中药物名大小写、连字符、空格混乱，本表统一转为大写后匹配
3. 执行时先用 02_clean.py 中的探查函数导出未匹配的高频药物名，
   人工确认后补充进 VARIANTS 列表
"""

# 五个目标药物：通用名 → (组别, 靶点, 商品名, 拼写变体列表)
DRUGS = {
    "TECLISTAMAB": {
        "group": "BISPECIFIC",      # 双特异性抗体组
        "target": "BCMA",
        "brand": ["TECVAYLI"],
        "variants": ["TECLISTAMAB", "TECLISTAMAB-CQYV", "TECVAYLI"],
    },
    "ELRANATAMAB": {
        "group": "BISPECIFIC",
        "target": "BCMA",
        "brand": ["ELREXFIO"],
        "variants": ["ELRANATAMAB", "ELRANATAMAB-BCMM", "ELREXFIO"],
    },
    "TALQUETAMAB": {
        "group": "BISPECIFIC",
        "target": "GPRC5D",
        "brand": ["TALVEY"],
        "variants": ["TALQUETAMAB", "TALQUETAMAB-TGVS", "TALVEY"],
    },
    "IDECABTAGENE VICLEUCEL": {
        "group": "CAR_T",           # CAR-T 组
        "target": "BCMA",
        "brand": ["ABECMA"],
        "variants": [
            "IDECABTAGENE VICLEUCEL", "IDECABTAGENE", "IDE-CEL",
            "ABECMA", "BB2121",
        ],
    },
    "CILTACABTAGENE AUTOLEUCEL": {
        "group": "CAR_T",
        "target": "BCMA",
        "brand": ["CARVYKTI"],
        "variants": [
            "CILTACABTAGENE AUTOLEUCEL", "CILTACABTAGENE", "CILTA-CEL",
            "CARVYKTI", "CARVYKTI (CILTACABTAGENE AUTOLEUCEL)",
        ],
    },
}

# 反向索引：任一写法 → 标准通用名（02_clean.py 使用）
def build_variant_index():
    idx = {}
    for generic, info in DRUGS.items():
        for v in info["variants"]:
            idx[v.upper().strip()] = generic
    return idx

VARIANT_INDEX = build_variant_index()

# 组间比较定义（04_ror_analysis.py 使用）
GROUPS = {
    "BISPECIFIC": ["TECLISTAMAB", "ELRANATAMAB", "TALQUETAMAB"],
    "CAR_T": ["IDECABTAGENE VICLEUCEL", "CILTACABTAGENE AUTOLEUCEL"],
}
