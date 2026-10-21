# Stage-1 standards — verbatim excerpts from the project log

Source: `docs/项目进展与交接.md` of the project (a working log that is appended to over time and is not under version control). Excerpts are copied verbatim, with their line numbers in the log as it stood on 2026-10-21 (file SHA-256 at extraction: 1dd0af36373c818d8fb25f0ae5989a333b9d231c591b3cb066cd5c8854352906). They are reproduced here for transparency; this file has no evidential value for the time at which the text was written — see `../README_stage1_standards.md`.

## 03 配对伪批量预设（stage 1 primary comparison） (log lines 129–134)

```text
- [x] **03 伪批量配对分析（预设，2026-10-08 写定，不得事后更改；已完成，见上文）**：
  1. 主要分析：AT2 谱系（AT2 + AT2/Club-like epithelial + 清理后的增殖上皮）；全部上皮作敏感性分析。
  2. 同一患者同一阶段的多个样本（P4 AAH/AAH1、P21 AIS/AIS1、P23 AIS/AIS1）先合并计数。
  3. AAH/AIS/MIA 合为「前体病变」，与同一患者的 LUAD 配对；P7 无前体病变，不纳入主要分析。
  4. 每个伪批量样本至少 50 个细胞，不足的剔除并列出。
  5. pydeseq2，设计 ~ 患者 + 病变类型，FDR<0.05。CNV 得分只作辅助指标。
```

## 第 3–5 层判定标准 (log lines 92–102)

```text
- **第 3–5 层判定标准（2026-10-08 用户确认，运行前写定；草案全文见 docs/第3-5层判定标准草案.md，以下为确认后的最终版）**。执行顺序：第 3 层 → 第 4 层 → 第 5 层 → 最后对留下的基因跑第 2 层（虚拟敲除，仅作预测）。三态：通过/未通过/无法检验（无法检验不算通过）。第 3–5 层对并集 683 个都算，**入选以预设的 352 为准**，结果标注来源（仅 352 / 两版都有 / 仅 AT2 谱系版）。
  - **下调基因（128 个，按 352 计）**：不进第 3–5 层排序，单列「缺失/恢复清单」，清单里额外报告 DGIdb 激动剂类药物（只探索，不判通过）。
  - **第 3 层（DepMap，方案 A）**：依赖细胞系 Chronos≤−0.5；常见必需=≥90% 全部细胞系依赖；通过=肺腺癌系（98 个模型，53 个有 CRISPR）中≥20% 依赖、非常见必需、在肺腺癌系表达（log2(TPM+1)≥1 的系占≥50%）；否决证据=常见必需，或敲除促进生长（Chronos≥+0.3 的肺腺癌系≥20%）。无存活依赖的基因记「未通过」，但放入「**待补充池**」：后续补 LINCS L1000 的 shRNA/CRISPR 敲低共识特征，检验敲低后 A 上调程序是否下降（下载规模和来源另行报批）。已知局限：DepMap 是 2D 存活依赖，会漏掉驱动浸润但不影响增殖的基因。
  - **第 4 层（I 期队列）**：队列 GSE31210、GSE50081、GSE30219、GSE37745、GSE68465、GSE72094、TCGA-LUAD；肺腺癌、病理 I 期（GSE31210、GSE50081 先筛 I 期；GSE68465 无分期，用 pT1–2 且 pN0 近似；各队列 IA/IB 依据的 TNM 版本不同（GSE31210/GSE50081 为较早版本，TCGA 为 AJCC 7/8 版混合），逐队列注明）；终点优先 RFS/DFS，没有的用 OS（GSE72094、TCGA）；Cox（每 1 SD，校正年龄、性别、分期）+ 随机效应 meta；通过（上调基因）=≥2 个队列风险比>1 且单侧 P<0.05，且 meta 单侧 P 的 BH q<0.05，且无队列显著反向；多探针取平均表达最高的一个；阶段比较（GSE102511、GSE166720、GSE193725）只作辅助。已批准下载平台注释 GPL570、GPL96、GPL15048（下载后核对大小）。
  - **第 5 层（可成药与安全）**：可成药=DGIdb 有抑制类药物（inhibitor、antagonist、blocker、antibody、negative modulator、inverse agonist、suppressor 等；agonist/activator/positive modulator 不算）；A=有已上市药物，B=仅在研（来自 ChEMBL、TTD、GuideToPharmacology、DTC）；**DGIdb 类型缺失但来自 ChEMBL/TTD/GuideToPharmacology/DTC 的相互作用计入，标「类型待核」，最终候选人工核对**。安全：GTEx v8 心（左心室、心房）、肝、肾（皮质、髓质）、全血任一中位 TPM ≥50 为「高表达」否决，10–50 标注「中等」不否决；HPA 25.1 在心肌、肝细胞、肾小管/肾小球、骨髓造血细胞中「High」（可靠性 Approved/Supported/Enhanced）否决（含骨髓），「Medium」标注；HPA 无数据只看 GTEx 并标「蛋白证据缺失」；常见必需基因（第 3 层）标「长期毒性风险」。
  - **注明**：AT2 谱系验证中验证集取的是「AT2 + 气道 Club/secretory」，与发现集的 AT2/Club-like + 增殖上皮不完全对应（验证集没有增殖上皮标签）。
- [x] **第 3–5 层已按确认的标准跑完**（2026-10-08，脚本 08/09/10，汇总 results/11_integrated/integrated_layers.tsv；候选 results/07_candidates/candidates_683.tsv，683=预设 352 ∪ AT2 谱系版 661，上调 434、下调 249；预设 352 内上调 224、下调 128）。平台注释已下载并核对：GPL570（8,471,521 字节）、GPL96（4,522,748 字节）、GPL15048（平台表 897 KB，64,074 行）。
  - **第 3 层（DepMap，肺腺癌 98 个模型，53 个有 CRISPR）**：上调基因通过 53（预设 352 内 29）、未通过 364、无法检验 17；常见必需 70、待补充池 294。
  - **第 4 层（队列）**：GSE31210（168 例，37 事件，RFS）、GSE50081（92，22，DFS）、GSE30219（79，23，DFS，pT1–2N0M0 近似）、GSE37745（34，15，RFS）、GSE68465（261，90，RFS，pT1–2N0 近似，注意死亡未计入 RFS 事件）、GSE72094（249，51，OS）、TCGA（280，68，OS，AJCC 版本分布：7th 156、6th 102、5th 14、3rd 3、4th 1）。其余队列的 TNM 版本 GEO 元数据未注明，**需查原文核对**（未猜测）。上调基因通过 104（预设 352 内 59）、未通过 317、无法检验 13；下调基因 98 个在≥2 队列 HR<1 且单侧 P<0.05（仅报告）。
  - **第 5 层**：上调基因通过 38（预设 352 内 17）；可成药 A 49/B 43；类型待核 30；GTEx 否决 115；HPA High 120；下调基因中有激动剂类药物的 7 个（down_agonists.tsv）。
  - **结果：第 3–5 层全部通过的上调基因只有 1 个——TTK（仅出现在 AT2 谱系版候选，不在预设 352 内）；预设 352 内没有基因同时通过第 3、4、5 层**（通过 2 层的 19 个）。成功标准要求至少 1 个靶点通过全部五层，第 2 层还没跑；按预设的 352 目前没有靶点满足，TTK 来自事后追加的 AT2 谱系验证，须如实标注。主要卡点：第 5 层「无可用抑制类药物」（352 内 47/59 个通过第 4 层的基因）、第 3 层常见必需（PLK1、BIRC5、RRM2 等）、GTEx/HPA 高表达。
```

## 补充分析判定标准（含微环境比较与外部验证规则） (log lines 103–109)

```text
- **补充分析判定标准（2026-10-08 用户确认，运行前写定；草案全文见 docs/补充分析判定标准草案.md，以下为确认后的最终版，含用户的 6 项修改）**。顺序：1 补完第 3 层（LINCS）→ 2 间接路径 → 3 微环境 → 4 探索性（Cox 加增殖评分，不参与判定）；第 2 层（虚拟敲除，仅作预测）最后。入选标准一律不放宽。
  - **批准下载**：Enrichr 库 `LINCS_L1000_CRISPR_KO_Consensus_Sigs`、GEO GSE92742 gene_info（212 KB）、SigCom LINCS API 查询（总计 <100 MB；不下载 5–20 GB 的 Level 5）。另需 GO 细胞周期基因集：用 Enrichr 的 GO_Biological_Process_2023 库（约十 MB，计入 <100 MB 预算），取名称含 cell cycle、mitotic、cell division、DNA replication、chromosome segregation、spindle、cytokinesis、kinetochore 的条目基因（运行前写定，不看结果调整）。
  - **补完第 3 层（LINCS）**：待补充池 294 个（并集 683 内；预设 352 内 148 个）。浸润程序 = 预设 352 内经验证的上调 224 + 下调 128，**只保留 L1000 测到的基因（landmark 或 best-inferred，gene_info 中 pr_is_lm=1 或 pr_is_bing=1），报告剩余基因数**。程序拆为「增殖」（Tirosh 97 个 S/G2M 基因 ∪ GO 细胞周期相关基因）与「非增殖」。每个池内基因的 KO（及可得的 shRNA）共识特征（上调集 U、下调集 D）与程序的重叠用**超几何检验**（背景=L1000 测到的基因，按基因集大小标准化）；得分 S = [−log10 p(D∩程序上调) − log10 p(U∩程序下调)] − [−log10 p(U∩程序上调) − log10 p(D∩程序下调)]；零分布=库内其他基因的特征用同一程序算 S，经验单侧 P，BH 校正；程序基因自身作为 KO 对象的特征不入零分布。**通过=非增殖程序被显著逆转（q<0.05 且 S>0）**；只有增殖部分逆转的标「仅增殖逆转」，不算通过。**CRISPR 敲除不要求目标基因自身下降；shRNA 特征仍要求目标基因在其自身下调集中。** 主判定用跨细胞系共识；肺腺癌系（A549、HCC515 等）单细胞系特征只作辅助。第 3 层通过 = DepMap 通过（方案 A）**或** LINCS 通过；常见必需基因的否决不变。无法检验（基因不在库中）不算通过。
  - **间接路径**：范围 = 预设 352 内（AT2 谱系版并列）上调基因中 DGIdb 无抑制类药物、且在 CollecTRI/DoRothEA 有 ≥5 个靶基因的 TF。调控资源：CollecTRI（decoupler 2.2.0，记录下载日期和版本）+ DoRothEA A/B 级交叉印证；关系计入 = 出现在 CollecTRI 或 DoRothEA A/B，且符号一致。**效应分子必须完整满足第 1 层（含外部验证，以预设 352 版为准，AT2 谱系版并列报告），方向 = TF 方向 × 调控符号，通过来源过滤，并独立通过第 3–5 层。** TF 活性（ULM）只作辅助报告。
  - **微环境**：T/NK、髓系、B/浆、成纤维、血管内皮、平滑肌/周细胞各做 22 对配对伪批量（~ 患者 + 病变类型，FDR<0.05，每侧≥50 个细胞，基因过滤同上皮）；来源过滤反向（群内表达/上皮表达 <0.2 视为疑似上皮来源；固定环境 RNA 名单 SFTPA1/2、SFTPB、SFTPC、SFTPD、NAPSA、SCGB1A1、SCGB3A1、SCGB3A2、SLPI、WFDC2、MUC1、LAMP3、AGER；再加与「SFTPC+SFTPB+SFTPA1 平均表达」Spearman 相关>0.6 的基因）；外部验证在 GSE164789 的同一细胞群按同规则（IAC vs AIS+MIA，单侧 BH q<0.05 且同向）；配体–受体用 LIANA+ 1.10.0（consensus 资源，rank_aggregate），在 GSE164789 同规则验证。第 4 层沿用同一套队列，加该细胞群标志基因平均 z 分数为协变量，调整后仍 q<0.05；第 5 层沿用 DGIdb/HPA/GTEx，**心、肝、肾否决不变；免疫来源靶点的血液/骨髓否决改为：在本研究单核数据中，目标细胞类型表达 ≥ 其他免疫谱系最高者的 3 倍（样本级中位），且 HPA 骨髓造血细胞无「High」（此规则仅用于免疫靶点）**。**微环境靶点第 3 层不适用，按设计无法满足「通过全部五层」，报告为「通过第 1、4、5 层 + 文献证据」，不计入成功标准**；文献证据人工核查（每个候选记录≥1 篇该细胞类型功能研究，不进自动判定）。
  - **探索性**：第 4 层 Cox 加增殖评分（Tirosh 97 个基因，队列内 z 分数平均）作协变量；被检验基因属于该基因集则不调整；报告调整前后 meta HR 与 q；不改变第 4 层判定。
  - **环境**：安装 decoupler 2.2.0 / liana 1.10.0 前先 `pip install --dry-run`，如会改动 scanpy、anndata、numpy、pandas 版本，就另建环境 ggn_lr，不动 ggn_target。
```

## 分期敏感性分析 (log lines 117–117)

```text
- **敏感性分析（按 Peng Table S1 分期；2026-10-08，在看到结果前依据 Table S1 的「Stage」列设定；只报告，不改变任何入选判定）**：Table S1 Tab 1（原表另存 results/01_snRNA/peng_table_s1_tab1.tsv）显示 23 名患者中 II–III 期 5 名：P5（IIA）、P11（IIA）、P20（IIIA）、P24（IIA）、P25（IIB）。排除这 5 名，P7 仍排除（无前体病变），剩 **17 对**（P3、P4、P6、P8、P9、P10、P12、P13、P14、P15、P16、P17、P18、P19、P21、P22、P23）。重跑（脚本 18_sensitivity_stageI.py，设计与过滤同主分析）：a) 上皮 A（AT2 谱系配对伪批量）：显著基因数；与主分析的方向一致率（主分析显著基因中，两个分析都检测到的基因里 log2FC 同号的比例）、log2FC 相关（Pearson 与 Spearman，两者都检测到的基因）；预设 352 个第 1 层基因中仍 FDR<0.05 且同向的数量（含上调/下调）。b) 髓系配对伪批量：SPP1、TREM2、CXCL9 的 log2FC 与 FDR。c) 两侧增殖上皮比例（每个伪批量里 Proliferating epithelial 占 AT2 谱系细胞的比例）的配对比较（双侧 Wilcoxon，均值差，LUAD 更高的患者数），与主分析 22 对并列。**不改变任何入选判定。** 附注：Table S1 的 P24、P25 肿瘤最大径（4.8、5）单位应为 cm；Table S1 无影像资料，病灶是否为 CT 发现的磨玻璃结节未报告（汇总文档的局限已改写）。
```
