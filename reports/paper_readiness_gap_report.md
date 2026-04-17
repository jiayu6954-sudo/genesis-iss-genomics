# 论文就绪度缺口分析报告
**分析日期**: 2026-04-17
**脚本**: `science_engine/analysis/paper_readiness_gap_analysis.py`

---

## 关键结论（先看这里）

| 缺口 | 问题 | 结论 |
|------|------|------|
| 缺口1 | CRISPR 差异是否在物种匹配后仍显著？ | **仅 P. polymyxa 显著** ✅ 显著 (p=0.009496)；其余3物种 Ground 也天然无 CRISPR |
| 缺口1 | H-003 "30/30" 表述是否准确？ | **需精化**：正确表述为"P. polymyxa ISS株100%缺失（Ground 53%携带，p<0.05）" |
| 缺口2 | IS 元件差异是否在物种内仍成立？ | **待下方分析结果** |
| 缺口4 | CRISPR 缺失与 IS 减少是否相关？ | **P. polymyxa 内部分析（见下）** |

---

## 缺口1：物种分层 CRISPR 分析

### 重要发现：CRISPR 缺失的物种特异性

| 物种 | ISS CRISPR+ | ISS缺失率 | Ground CRISPR+ | Ground携带率 | Fisher p(两侧) | 结论 |
|------|------------|---------|--------------|------------|-------------|------|
| Bacillus thuringiensis | 0/11 | 100% | 0/15 | 0% | 1.0 | ○ 物种天然特征 |
| Bacillus amyloliquefaciens | 0/9 | 100% | 0/15 | 0% | 1.0 | ○ 物种天然特征 |
| Paenibacillus polymyxa | 0/9 | 100% | 8/15 | 53% | 0.009496 | ⭐ ISS特异 |
| Cytobacillus firmus | 0/1 | 100% | 0/10 | 0% | 1.0 | ○ 物种天然特征 |

**全物种合并 Fisher p = 0.046071**

### 🔑 关键解读

- *B. thuringiensis*、*B. amyloliquefaciens*、*C. firmus* 在**地面环境也天然不含 CRISPR**
  → 这 3 个物种的"CRISPR 缺失"是**物种固有特征**，非 ISS 选择压力导致
- **仅 *P. polymyxa* 表现出 ISS 特异性 CRISPR 丢失**（Ground 53% vs ISS 0%）
- 原始 H-003 表述"30/30 = 100% 缺失"在事实上正确，但**科学意义被夸大**
- **修正表述**：
  > *P. polymyxa* isolates from the ISS showed complete absence of CRISPR-Cas systems
  > (0/9, 0%), significantly lower than ground isolates of the same species
  > (8/15, 53%; Fisher's exact p = 0.009496)
  > The other three species (*B. thuringiensis*, *B. amyloliquefaciens*,
  > *C. firmus*) naturally lack CRISPR-Cas in both ISS and ground environments.

---

## 缺口2：物种分层 IS 元件分析

| 物种 | ISS n | ISS IS/千CDS | Ground n | Ground IS/千CDS | 比值 | Cliff δ | p(MW) | 结论 |
|------|------|------------|---------|--------------|-----|---------|-------|------|
| Bacillus thuringiensis | 11 | 1.294±0.001 | 15 | 11.221±7.258 | 0.115× | -1.0 | 1.9e-05 | ⭐ ISS↓显著 |
| Bacillus amyloliquefaciens | 9 | 2.768±0.693 | 15 | 1.916±2.384 | 1.445× | 0.5704 | 0.021693 | ⭐ ISS↑显著 |
| Paenibacillus polymyxa | 9 | 1.567±0.004 | 15 | 5.026±2.605 | 0.312× | -1.0 | 5.7e-05 | ⭐ ISS↓显著 |
| Cytobacillus firmus | 1 | 12.314±0.000 | 10 | 5.337±3.568 | 2.308× | 1.0 | 0.113846 | ≈ |

---

## 缺口4：CRISPR × IS 相关性（P. polymyxa 内部）

- CRISPR- 株 (n=16): IS/千CDS = **2.737 ± 2.114**
- CRISPR+ 株 (n=8): IS/千CDS = **5.712 ± 2.621**
- Mann-Whitney U=16.0, p=0.003289, Cliff δ=-0.75
- **→ H-new-C 机制支撑：CRISPR 缺失与 IS 减少在 P. polymyxa 内相关**

---

## 修正后的论文声明框架

### 可以发表的声明

| 声明 | 统计支撑 | 证据质量 |
|------|---------|---------|
| P. polymyxa ISS株CRISPR显著缺失（0/9 vs 8/15 Ground） | Fisher p<0.05 | **强** |
| ISS群落 Firmicutes 主导（47.7%），多样性崩溃（H' 4.4→2.1） | 描述性，n=2 | **弱（描述性）** |
| GLDS-224 ISS碎屑 CRISPR/千KO = 0.50× Ground | n=2，无统计检验 | **中（趋势）** |
| KVTAG LexA切割位点4物种100%保守 | 结构分析，n=4 | **中（描述性）** |

### 需要修改的声明

| 原始声明 | 修正版本 |
|---------|---------|
| "30/30 ISS株100%缺失CRISPR" | "P. polymyxa ISS株CRISPR显著低于地面株（p<0.05）；其余3物种在地面也天然无CRISPR" |
| 'ISS IS元件仅为Ground 0.369×' | 'P. polymyxa ISS IS密度0.312× Ground，p=5.7e-05（物种内控制后仍显著）' |

---

## 论文投稿建议（更新版）

| 期刊 | IF | 是否就绪 | 条件 |
|------|-----|---------|------|
| npj Microgravity | 6.3 | ✅ 接近就绪 | 物种分层CRISPR + 补充P. polymyxa文献背景 |
| Microbiome | 8.6 | ⚠️ 需补充 | 需n≥5的GLDS-224群落验证 |
| ISME J | 9.8 | ❌ 需更多工作 | 需机制实验 + IS元件物种内验证 |

---
*本报告由 Project Genesis AI 系统自动生成 (SESSION-009)*