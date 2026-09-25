# -*- coding: utf-8 -*-
"""
不良事件定义 —— MedDRA Preferred Term (PT) 清单
=====================================================
⚠️ 人工核对要点（低成本 AI 执行前必读）：
1. 本清单为框架初稿，必须由人工对照研究所用 FAERS 季度的
   MedDRA 版本核对每个 PT 的确切拼写（FAERS 使用大写 PT）
2. 核对依据：FDA 各药物标签、已发表药物警戒研究的事件定义
3. 新增/删除 PT 后在论文补充材料 Suppl. Table 1 中同步更新
"""

EVENTS = {
    "CRS": {
        "description": "细胞因子释放综合征",
        "pts": [
            "CYTOKINE RELEASE SYNDROME",
            "CYTOKINE STORM",
        ],
    },
    "ICANS_NEUROTOXICITY": {
        "description": "ICANS 及广义神经毒性",
        "pts": [
            "IMMUNE EFFECTOR CELL-ASSOCIATED NEUROTOXICITY SYNDROME",
            "ENCEPHALOPATHY",
            "CONFUSIONAL STATE",
            "DELIRIUM",
            "APHASIA",
            "SEIZURE",
            "NEUROTOXICITY",
            "CEREBRAL OEDEMA",
        ],
    },
    "INFECTION_ALL": {
        "description": "全部感染（SOC Infections and infestations 的 PT 扩展近似）",
        "pts": [
            "SEPSIS", "SEPTIC SHOCK", "BACTERAEMIA", "BACTERIAL SEPSIS",
            "PNEUMONIA", "PNEUMONIA ASPIRATION", "BRONCHITIS", "BRONCHOPNEUMONIA",
            "LOWER RESPIRATORY TRACT INFECTION", "RESPIRATORY TRACT INFECTION",
            "UPPER RESPIRATORY TRACT INFECTION", "NASOPHARYNGITIS", "SINUSITIS",
            "PHARYNGITIS", "TONSILLITIS", "LARYNGITIS",
            "URINARY TRACT INFECTION", "PYELONEPHRITIS",
            "GASTROENTERITIS", "ENTEROCOLITIS INFECTIOUS", "COLITIS",
            "CLOSTRIDIUM DIFFICILE COLITIS", "DIVERTICULITIS", "PERITONITIS",
            "ABDOMINAL SEPSIS", "CHOLECYSTITIS", "CHOLANGITIS",
            "CELLULITIS", "ERYSIPELAS", "ABSCESS", "ABSCESS LIMB", "ABSCESS SKIN",
            "OSTEOMYELITIS", "SEPTIC ARTHRITIS", "ENDOCARDITIS",
            "MENINGITIS", "ENCEPHALITIS",
            "CANDIDIASIS", "ORAL CANDIDIASIS", "OESOPHAGEAL CANDIDIASIS",
            "SYSTEMIC CANDIDIASIS", "VULVOVAGINAL CANDIDIASIS", "FUNGAL INFECTION",
            "PNEUMOCYSTIS JIROVECII PNEUMONIA", "ASPERGILLOSIS", "MUCORMYCOSIS",
            "HISTOPLASMOSIS", "COCCIDIOIDOMYCOSIS", "CRYPTOCOCCOSIS",
            "CYTOMEGALOVIRUS INFECTION", "CYTOMEGALOVIRUS COLITIS",
            "EPSTEIN-BARR VIRUS INFECTION", "HERPES ZOSTER", "HERPES SIMPLEX",
            "INFLUENZA",
            "COVID-19", "COVID-19 PNEUMONIA", "CORONAVIRUS INFECTION",
            "TUBERCULOSIS", "MYCOBACTERIAL INFECTION", "NOCARDIOSIS",
            "STAPHYLOCOCCAL INFECTION", "STREPTOCOCCAL INFECTION",
            "ENTEROCOCCAL INFECTION", "KLEBSIELLA INFECTION",
            "PSEUDOMONAS INFECTION", "ESCHERICHIASIS",
            "PROTEUS INFECTION", "ACINETOBACTER INFECTION",
            "INFECTION", "POSTOPERATIVE WOUND INFECTION", "WOUND INFECTION",
            "VIRAL INFECTION", "BACTERIAL INFECTION",
        ],
    },
    "SERIOUS_INFECTION": {
        "description": "严重感染（示例子集，需人工补充）",
        "pts": [
            "SEPSIS", "SEPTIC SHOCK", "PNEUMONIA",
            "COVID-19", "COVID-19 PNEUMONIA",
            "CYTOMEGALOVIRUS INFECTION", "PNEUMOCYSTIS JIROVECII PNEUMONIA",
            "PROGRESSIVE MULTIFOCAL LEUKOENCEPHALOPATHY",
        ],
    },
    "GPRC5D_OFFTUMOR": {
        "description": "GPRC5D 相关脱靶毒性（talquetamab 特征谱）",
        "pts": [
            "DYSGEUSIA", "AGEUSIA", "HYPOGEUSIA",
            "DRY MOUTH", "ORAL DISORDER", "STOMATITIS",
            "SKIN EXFOLIATION", "PALMAR-PLANTAR ERYTHRODYSAESTHESIA SYNDROME",
            "RASH MACULO-PAPULAR",
            "NAIL DISORDER", "NAIL DYSTROPHY", "ONYCHOMADESIS",
            "WEIGHT DECREASED", "DECREASED APPETITE",
        ],
    },
    "DELAYED_NEUROTOXICITY": {
        "description": "迟发神经毒性（cilta-cel 标签警示谱）",
        "pts": [
            "PARKINSONISM", "PARKINSON'S DISEASE",
            "CRANIAL NERVE DISORDER", "FACIAL PARALYSIS", "BELL'S PALSY",
            "GUILLAIN-BARRE SYNDROME",
            "MOVEMENT DISORDER", "TREMOR", "ATAXIA", "GAIT DISTURBANCE",
            "PERIPHERAL NEUROPATHY", "POLYNEUROPATHY",
        ],
    },
    "CYTOPENIA": {
        "description": "血液学毒性",
        "pts": [
            "NEUTROPENIA", "FEBRILE NEUTROPENIA",
            "THROMBOCYTOPENIA", "ANAEMIA", "LYMPHOPENIA",
            "PANCYTOPENIA", "BONE MARROW FAILURE",
            "HYPOGAMMAGLOBULINAEMIA",
        ],
    },
    "SECONDARY_MALIGNANCY": {
        "description": "第二原发肿瘤（探索性，含 CAR-T 后 T 细胞肿瘤信号）",
        "pts": [
            "SECOND PRIMARY MALIGNANCY",
            "T-CELL LYMPHOMA", "T-CELL ACUTE LYMPHOBLASTIC LEUKAEMIA",
            "MYELODYSPLASTIC SYNDROME", "ACUTE MYELOID LEUKAEMIA",
        ],
    },
}

# 死亡结局：不使用 PT，用 OUTC 表 OUTC_COD = 'DE' 判定（03_build_dataset.py 处理）
DEATH_OUTCOME_CODE = "DE"
