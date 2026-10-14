#!/usr/bin/env python3
"""Frozen post-ranking dossiers. No writes to frozen rankings or candidate pool."""
from pathlib import Path
import pandas as pd, numpy as np, json, hashlib, math, shutil
from scipy.stats import nct
from scipy.optimize import brentq
from statsmodels.stats.power import TTestPower
B=Path('/Users/dudu/全部数据库/GGN_target'); R=B/'results/50_screen'
W=Path('/Users/dudu/Documents/Codex/2026-10-09/ggn-target-ggn-target-opt-anaconda3'); O=W/'outputs/step6_dossiers'; O.mkdir(parents=True,exist_ok=True)
inputs={}
def read(p,index=None):
 p=Path(p); inputs[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.read_csv(p,sep='\t',index_col=index,low_memory=False)
def save(d,n):d.to_csv(O/n,sep='\t',index=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def table(d):
 def fmt(v):
  if isinstance(v,(float,np.floating)):return f'{v:.4g}'
  return str(v).replace('|','/').replace('\n',' ')
 return '\n'.join(['| '+' | '.join(map(str,d.columns))+' |','| '+' | '.join(['---']*len(d.columns))+' |']+['| '+' | '.join(fmt(v) for v in row)+' |' for row in d.itertuples(index=False,name=None)])
manifest=json.loads((R/'step4/freeze_manifest.json').read_text()); frozen={f['path']:f['sha256'] for f in manifest['frozen_files']}
assert all(sha(p)==v for p,v in frozen.items())
allr=read(next(iter(frozen))); pool=allr[allr.eligible_topN_pool.eq(True)].copy();assert len(pool)==82
pool=pool.sort_values('rank_eligible_pool_display'); top=pool.head(20).copy()
selected=['PYCR1','BCAT1','CD1A','CD207']; units={'PYCR1':'M1-PYCR1','BCAT1':'M1-BCAT1','CD1A':'M2-CD1A/CD207','CD207':'M2-CD1A/CD207'}
assert set(selected)<=set(top.gene); assert not top[top.gene.isin(selected)].top5_excluded_existing.any()
LIMIT='无亚阶段统计支持，仅方向线索；AAH仅8对'
CAS='本研究自定义的4基因CAS评分可能主要捕捉正常club/分泌细胞，其对癌相关状态的特异性缺乏支持'
KAC='KAC高分细胞在健康肺中罕见。KAC本身为肺泡损伤修复状态，健康肺中罕见不等于对新生过程特异'
DECL='干预方向不可由表达高低推断；本研究不推断抑制该分子是有益还是有害。当前均为探索性、假说生成候选，未证明前驱病变依赖性或人体效用。'
reasons={
'ETV4':'TF遗传操作可行，IHC Approved；无核实的直接药理工具，本轮聚焦代谢与DC两条方向，保留高排序参照。',
'KDELR3':'遗传操作可行，但IHC Uncertain且无核实直接药理工具；ER运输为另一机制，暂不进入本轮。',
'ISG15':'IHC Approved；胞内ISGylation与分泌功能方向未分解，直接抑制方向未知。',
'AMIGO2':'表面分子有定位价值，但IHC Uncertain、直接干预工具及早期方向均待确认。',
'PYCR1':'入选：IHC Approved、酶与遗传扰动可操作，新增原始研究有PYCR1-IN-1工具线索；丙档弱于甲档，必须并列保留。',
'PAFAH1B3':'丙档，IHC Uncertain，当前工具身份/选择性仍待确认；不优于本轮代谢候选。',
'XPR1':'甲档、IHC Supported且遗传扰动可行；磷酸转运不同于氨基酸代谢，列第一备选，不为凑两组强行合并。',
'SRD5A3':'丙档、IHC Uncertain；现有5α还原酶药物不能自动当作SRD5A3选择性工具。',
'CRLF1':'甲档、IHC Approved；细胞因子复合物与受体背景尚未明确，干预方向未知。',
'CD1A':'入选联合单元：甲档、IHC Enhanced，适合与CD207联合验证DC身份/空间分布；不是独立的抑制靶点。',
'CD207':'入选联合单元伴随标志：甲档、IHC Enhanced；与CD1A同一机制，不增加独立候选数。',
'SLC16A7':'IHC Approved；转运底物与代偿未拆分，直接药理工具在冻结证据中未核实。',
'CEACAM5':'IHC Enhanced，但为自定义CAS签名成员，CAS特异性限制直接相关；表达标志不等于依赖。',
'BCAT1':'入选：甲档、遗传扰动和研究性工具线索，支持测试氨基酸代谢；HPA抗体Uncertain且安全未评估，先做抗体/靶向性技术验证。',
'MMP15':'蛋白酶候选，但IHC Uncertain，亚型选择性和基质来源待查。',
'PHLDA3':'存在促癌与保护功能文献冲突，不能确定抑制方向；不选为本轮抑制实验优先项。',
'ENPP5':'IHC Enhanced利于定位，但直接功能方向及工具缺口大。',
'MMP24':'IHC Approved，蛋白酶机制与MMP15相近，缺早期功能和选择性工具确认。',
'ZNF486':'无已核实IHC资源，TF功能与干预工具不足。',
'PTPRN2':'IHC Enhanced可定位，但分泌囊泡相关功能与干预方向待定。'}
groups={'ETV4':'AT2L/转录调控','KDELR3':'AT2L/ER运输','ISG15':'AT2L/干扰素及ISGylation','AMIGO2':'AT2L/黏附','PYCR1':'M1/AT2L氨基酸代谢','PAFAH1B3':'AT2L/脂质水解','XPR1':'AT2L/磷酸转运','SRD5A3':'AT2L/糖基化','CRLF1':'AT2L/细胞因子','CD1A':'M2/Myeloid朗格汉斯/DC','CD207':'M2/Myeloid朗格汉斯/DC','SLC16A7':'AT2L/单羧酸转运','CEACAM5':'AT2L/黏附及CAS标志','BCAT1':'M1/AT2L氨基酸代谢','MMP15':'AT2L/基质蛋白酶','PHLDA3':'AT2L/信号调控','ENPP5':'AT2L/磷酸二酯酶','MMP24':'AT2L/基质蛋白酶','ZNF486':'AT2L/转录调控','PTPRN2':'AT2L/分泌囊泡'}
top['manual_selected']=top.gene.isin(selected); top['candidate_unit']=top.gene.map(units).fillna('未选'); top['mechanism_audit']=top.gene.map(groups);top['manual_reason']=top.gene.map(reasons)
top['selection_source']='82池内主排序前20；同时列消退分层位置';top['clinical_label_limit']=LIMIT;top['CAS_direct_signature_dependency']=top.gene.eq('CEACAM5');top['CAS_limit']=CAS;top['selectivity_final']='尚不能确定'
save(top,'top20_selection_audit.tsv');save(top.head(5),'pure_score_top5.tsv');save(top[top.manual_selected],'selected_gene_dossiers_data.tsv');save(pool[pool.top5_excluded_existing],'pool82_preserved_exclusions.tsv')
# Preserve every existing exclusion in complete candidate table, including records outside pool82.
save(allr[allr.top5_excluded_existing],'all_preserved_exclusions.tsv')
save(allr[allr.pool_before_D7.eq(True)&allr.top5_excluded_existing.eq(True)],'pool199_preserved_exclusions.tsv')
for name in ['main_all_top20.tsv','regression_all_top20.tsv']:
 shutil.copy2(R/'step4'/name,O/name)
# Patient-level power planning: paired continuous assay endpoint, RNA-to-protein standardized effect transfer is an unvalidated assumption.
power=TTestPower(); p_rows=[]; delta_rows=[]
def minimum(d):
 if not np.isfinite(d) or d<=0:return np.inf
 n=max(2,math.ceil(float(power.solve_power(effect_size=d,alpha=.05,power=.8,alternative='two-sided'))))
 while power.power(d,n,alpha=.05)<.8:n+=1
 return n
for g in selected:
 row=top.set_index('gene').loc[g];lin=row.lineage
 c1=read(R/f'step1/main/{lin}_I/pseudobulk_counts.tsv.gz',0);c2=read(R/f'step1/main/{lin}_II/pseudobulk_counts.tsv.gz',0)
 cp1=c1.div(c1.sum(axis=1),axis=0)*1e6;cp2=c2.div(c2.sum(axis=1),axis=0)*1e6
 patients=sorted(set(x.replace('_reference','') for x in cp1.index if x.endswith('_reference'))&set(x.replace('_target','') for x in cp2.index if x.endswith('_target')))
 iac=pd.Series({p:np.log2((cp2.loc[p+'_target',g]+.5)/(cp1.loc[p+'_reference',g]+.5)) for p in patients})
 for stage in ['AAH','AIS','MIA','IAC']:
  v=iac if stage=='IAC' else read(R/f'step1/substage_patient_log2FC_{lin}_{stage}.tsv',0)[g]
  v=v.dropna();n=len(v);mean=float(v.mean());sd=float(v.std(ddof=1));d=mean/sd;t=d*np.sqrt(n)
  lo=brentq(lambda nc:nct.cdf(t,n-1,nc)-.975,-1000,1000)/np.sqrt(n)
  hi=brentq(lambda nc:nct.cdf(t,n-1,nc)-.025,-1000,1000)/np.sqrt(n)
  nn=minimum(d);nc=minimum(lo)
  p_rows.append(dict(gene=g,candidate_unit=units[g],stage=stage,exploratory_n=n,mean_log2CPM_difference=mean,sd_paired_difference=sd,dz=d,dz_ci95_low=lo,dz_ci95_high=hi,min_patient_pairs=nn,conservative_min_pairs=nc,alpha=.05,power=.8,assay='paired continuous IHC H-score or mIF log-density/intensity; transfer assumption, not validated protein effect'))
  delta_rows +=[dict(gene=g,stage=stage,patient=p,paired_log2CPM_difference=x) for p,x in v.items()]
pd_power=pd.DataFrame(p_rows);save(pd_power,'power_stage_effects.tsv');save(pd.DataFrame(delta_rows),'power_patient_differences.tsv')
alloc=[]
for g in selected+['联合四基因面板']:
 pp=pd_power if g=='联合四基因面板' else pd_power[pd_power.gene.eq(g)]
 for typ,col in [('探索性点估计','min_patient_pairs'),('95%CI下界保守','conservative_min_pairs')]:
  a={stage:float(pp[pp.stage.eq(stage)][col].max()) for stage in ['AAH','AIS','MIA','IAC']}; total=sum(a.values());alloc.append(dict(gene_or_panel=g,scenario=typ,**a,paired_adjacent_normal=total,independent_patients=total,total_tissue_specimens=2*total,true_healthy_normal='无法据本研究估计；另设描述性可行性组，不计80%功效'))
allocation=pd.DataFrame(alloc);save(allocation,'power_stage_allocation.tsv')
val=read(R/'step5/pool82_validation_parallel.tsv');save(val[val.gene.isin(selected)],'selected_step5_validation.tsv')
# Contemporary normal reference and corrected selectivity outputs, without changing upstream files.
h=read(R/'step5/HLCA/pool82_expression_donor_summary.tsv');cor=read(R/'step3_correction/corrected_candidates_audit.tsv');save(h[h.gene.isin(selected)],'selected_HLCA_donor_reference.tsv');save(cor[cor.gene.isin(selected)],'selected_corrected_selectivity.tsv')
cal=read(R/'step3_correction/calibration_summary.tsv');save(cal,'selectivity_calibration_reference.tsv')
# Seven downgrade checks are recorded separately for each gene.
down=[]
for g in selected:
 r=top.set_index('gene').loc[g]; metab=g in ['PYCR1','BCAT1']
 vals=[('未触发',f'I LOO={r.I_LOO_status}；既有排除=False'),('不可作为降级依据',f'{r.I_D10}；状态内未检出也不能证明仅组成驱动'),('未触发既有排除',f'D3_I={r.D3_I}; {r.I_support_datasets}；第II步FPKM未达标不能写为独立验证'),('未触发','第I步有证据，非仅IAC出现；但消退仅方向线索'),('风险未解除，降低特异性主张','校正版+AT2L-other+健康HLCA三者合判尚不能确定；不足以触发确定性否决，也不得判安全'),('外部参照受限','HLCA与Peng平台/队列混杂不能分离，禁做合并病例对照；不据此否决内部配对结果'),('触发：方向证据降级为仅晚期模型' if metab else '方向未知：无适用早期扰动支持','已建立肿瘤/细胞系证据不代表前驱病变；不得将未找到证据当作反向证据')]
 for i,(status,evidence) in enumerate(vals,1):down.append(dict(gene=g,criterion=i,status=status,evidence=evidence))
save(pd.DataFrame(down),'seven_downgrade_checks.tsv')
# Resources and experimental proposals: cost/time are rough planning estimates, not quotations.
paths=[]
for g in selected:
 r=top.set_index('gene').loc[g]; immune=g in ['CD1A','CD207']; validated=r.D8_IHC_gene_reliability
 for modality in ['IHC/mIF','细胞系/细胞功能','肺泡类器官','小鼠']:
  if modality=='IHC/mIF':
   feas='可行；抗体先做本地技术验证'+('；BCAT1抗体Uncertain须KO对照/替代抗体，不直接进入大队列' if g=='BCAT1' else '')
   resource=f'HPA记录 {r.D8_IHC_antibody_ids}，可靠性 {validated}；AAH/AIS/MIA/IAC配对FFPE；'+('CD1A/CD207/HLA-DR/CD45及上皮标志联合定位' if immune else '候选蛋白+上皮标志SFTPC/EPCAM，KAC与修复状态多标志辅助；病理医师标定ROI')
   endpoints='患者级CD1A+CD207+ DC密度、空间距离；单阳性单独记录' if immune else '患者级上皮内H-score/强度；阳性细胞比例为次要终点，区分细胞数与细胞内强度'
   controls='配对邻近正常、独立健康正常描述性组、损伤修复肺；阳性组织、无一抗/单染、批内参考、盲法；BCAT1/PYCR1加敲低细胞块验证抗体'
   interpretation='仅证实定位/蛋白差异，不能决定抑制或细胞清除；不得把ROI数当样本数'
   budget='技术小试2–5万元，4–8周；正式队列费用另按功效标本数和当地报价';cheap='是（技术小试；不含合作取材成本）'
  elif modality=='细胞系/细胞功能':
   feas='代谢基因遗传扰动可行，癌细胞系外推有限' if not immune else '常规肺癌细胞系不适用；需人源DC功能系统，条件性可行'
   resource=('至少两条sgRNA/siRNA、抗性cDNA救援；靶点表达阳性肺癌细胞并含正常上皮；'+('PYCR1-IN-1仅研究性线索' if g=='PYCR1' else 'BCAT1研究性WQQ-345工具线索；不可用gabapentin默认替代选择性抑制')) if not immune else '人外周血来源DC/组织DC；CD1A/CD207双标、抗原呈递和T细胞共培养；需先确认细胞身份'
   endpoints='靶蛋白抑制、存活/克隆、凋亡及代谢读出；正常上皮损伤并列' if not immune else '抗原呈递、T细胞反应、细胞生存；先观察再决定干预方向'
   controls='非靶向、两独立序列、救援、载体；药物需剂量反应与遗传一致性' if not immune else '无抗原、同型、单培养、正常DC；不能以CD207+细胞消除等同CD207分子抑制'
   interpretation='只回答所用模型的功能依赖；药理/遗传不一致需排查脱靶；未覆盖AAH/AIS'
   budget='约3–10万元，2–4个月（DC系统约5–15万元）';cheap='否'
  elif modality=='肺泡类器官':
   feas='有条件可行，需合作取得前驱病灶及正常供者' if not immune else '纯上皮类器官不适用，需人DC-上皮共培养'
   resource='患者来源前驱病灶/正常AT2类器官、培养基及ECM、CRISPRi/siRNA递送；'+('同供者DC与免疫共培养' if immune else '研究性化合物及靶向性读出')
   endpoints='已建立类器官之后干预；纵向体积/存活、上皮分化与损伤修复；按供者汇总'
   controls='供者配对正常、损伤修复刺激、非靶向和救援、空白载体；多个独立供者，孔数不替代供者数'
   interpretation='检测前驱相关依赖与修复毒性；阳性仍不等于人体结节消退；需保留真正前驱组织身份'
   budget='约10–30万元，4–9个月；免疫共培养可能更高';cheap='否'
  else:
   feas='条件性可行，现成品系/药物跨物种活性未核实' if g!='CD1A' else '普通小鼠缺乏直接CD1A同源模型；需人源化系统，不列常规可行路径'
   resource=f'冻结同源记录={r.D8_mouse_ortholog_record}；现成敲除品系未核实；'+('CD207/Langerin系谱/细胞消除模型须另核实；不是CD1A/CD207分子药理替代' if immune else '肺上皮条件性遗传扰动、适合前驱阶段的模型及微CT/病理随访')
   endpoints='病灶已建立后随机开始干预，测原病灶负荷/消退及进展；并测正常肺修复、体重/器官毒性'
   controls='时间匹配车辆/同窝对照、盲法病理、基线病灶负荷、药效靶向；患者方向不能由鼠瘤缩小反推'
   interpretation='区分减少发生与既有病灶消退；小鼠为单位；CD1A人源化差异与免疫重建需额外验证'
   budget='约20–60万元以上，9–18个月以上；人源化/新建品系另计';cheap='否'
  paths.append(dict(gene=g,candidate_unit=units[g],path=modality,feasibility=feas,resources=resource,endpoints=endpoints,controls=controls,interpretation=interpretation,time_cost=budget,lowest_cost=cheap))
pathdf=pd.DataFrame(paths);save(pathdf,'four_validation_paths.tsv')
# Stable prose sources: extra literature updates are annotations, never a rewrite of frozen D6.
sources='''- [PYCR1研究工具与肺癌模型原始研究](https://www.nature.com/articles/s12276-025-01577-z)：报告PYCR1-IN-1及遗传扰动；不是前驱病灶临床疗效证据。
- [BCAT1肺癌遗传扰动原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC9483610/)；[BCAT1与EGFR-TKI耐药研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC11324870/)：后者提供WQQ-345研究工具线索，耐药癌模型外推受限。
- [Gabapentin与BCAT1非特异性问题原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC6528808/)：不把gabapentin作为已验证的选择性BCAT1工具。
- [Choe原文](https://pubmed.ncbi.nlm.nih.gov/40824728/)：以聚类与标志基因界定CAS，无本研究自定义四基因评分方法。
- [正常人远端气道SCGB3A2分泌细胞研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297319/)：正常分泌细胞表达使“癌特异”解释不足；不意味着所有SCGB3A2+细胞都是经典club细胞。
- HPA抗体证据来自冻结D8记录，可靠性不等于本实验室FFPE验证；未采购、未做湿实验。'''
requests=[
 dict(dataset='Fu 2025',accession='NODE OEP002580; GSA HRA002624/HRA002983',url='https://www.nature.com/articles/s41422-025-01200-w',route='由机构PI通过NODE/GSA受控访问；优先RNA HRA002983，WGS按研究必要性申请',request='原始/计数RNA、患者-样本映射、AAH/AIS/MIA/IAC病理、邻近/健康正常定义、批次、吸烟及治疗、影像有无、随访时点/结局和删失；核对954患者/1008样本各模态实际可用数',benefit='扩大分期复现/随访关联；AAH极少不能补足AAH证据',priority='1'),
 dict(dataset='Chen 2025',accession='EGAD50000000637; DAC EGAC50000000385',url='https://ega-archive.org/dacs/EGAC50000000385',route='机构PI向EGA所列DAC申请并签署DUA；获批后受控下载',request='RNA/WES及患者级映射、配对正常和病理复核、批次、临床协变量；明确94前驱(25 AAH/69 AIS)与73 MIA的各模态样本覆盖',benefit='AAH/AIS/MIA阶段锚点；无IAC则不能独立检验前体到浸润',priority='1'),
 dict(dataset='HTAN Lung PCA/PCAPP',accession='phs002371; GSE319666 processed',url='https://data.humantumoratlas.org/data-access',route='先核对HTAN项目/病例/文件manifest，公开level3/4先用；原始数据经dbGaP审批后走GC/SB-CGC；不是所有HTAN均受控',request='RNA与LCM区域、患者/病灶配对、AAH/AIS/浸润病理原始标签、批次、临床及空间图像manifest；明确HTA2 PCAPP与HTA3 Lung PCA不同项目',benefit='独立区域与患者层面验证；GSE319666与HTAN重复病例不算两独立队列',priority='2'),
 dict(dataset='GSE164789作者健康肺锚点',accession='GSE164789; preprint 10.21203/rs.3.rs-1084110/v1',url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164789',route='准备给通讯作者的数据请求，未发送；确认预印本健康供者样本未上传GEO及共享许可',request='真正无癌健康供者正常肺原始/计数矩阵、细胞注释、供者ID与年龄性别吸烟、解剖区域、处理/测序批次、质控及伦理共享范围；不得用62例A侧邻近组织替代',benefit='同平台同流程健康参照；仍需核对批次/取材，不能预先保证完全可比或声称是全领域唯一锚点',priority='1')]
save(pd.DataFrame(requests),'controlled_data_request_list.tsv')
request_md='# 受控数据与健康肺申请清单（草案，未发送/未提交）\n\n目标：复现冻结候选与阶段关联，分离正常修复表达和前驱病变关联；不预设抑制有益。所有申请由有资格机构/PI完成。\n\n共同准备：研究摘要、固定候选/分析方案、申请人机构与PI、伦理批件或豁免、DUA、获授权成员、加密本地存储与访问日志、保存/销毁计划、不重识别承诺。只申请必要字段，不申请直接身份标识。\n'
for req in requests:
 request_md+='\n## '+req['dataset']+'\n\n'+'\n\n'.join(f'**{k}**：{v}' for k,v in req.items() if k!='dataset')+'\n'
request_md+='\n## 作者申请正文草案\n\n我们拟检验冻结的前驱病变表达候选，并评估癌相关状态评分与正常修复状态的区分能力。贵组预印本提及健康供者正常肺，而公开GSE164789所用A侧为肿瘤邻近组织。希望确认健康供者样本能否共享去标识的计数矩阵、细胞注释与供者/批次/取材区域信息，以及必要的原始测序数据访问途径。我们将遵守伦理、使用范围及署名/引用要求，并在分析前核对样本身份、流程和批次。若不能共享原始数据，受限统计输出或作者端执行冻结分析也有帮助。\n'
(O/'data_requests.md').write_text(request_md)
report=['# 第6步：冻结后候选档案与实验验证方案','\n本轮选择4个基因、3个候选单元、2条机制方向：PYCR1、BCAT1、CD1A/CD207。这是D6+D8与适切性人工判断，不是统计学证明这几项优于未选项。候选池、分数、权重、门槛及排除项保持冻结。',
'## 证据界限（档案与论文必须保留）',DECL,
'第I步：82条相应发现主张中22条方向一致，其中16条counts与FPKM非参数两法一致。第II步：7条只可表述为“counts方向一致，FPKM非参数核对未达正向标准”，不得写成“第II步获得独立验证”。两步强度不对称。',
'ETV4/KDELR3/ISG15/PYCR1/PAFAH1B3/SRD5A3/TREM2的MW q依次为0.063/0.221/0.083/0.180/0.113/0.074/0.085，全部无法判定。此处按用户已核验数值保留，未改写上游验证表。',
CAS+'。三条依据：健康HLCA的CAS高分16.30%，KAC为0.91%；第二阶段GSE164789为+0.94个百分点、5/14、P=0.69，未复现（该项引用用户已核验结果，本步未读取旧阶段结果或重新计算）；SCGB3A2为正常气道club/分泌细胞相关分泌蛋白。Choe原文用聚类与标志基因，无本研究评分方法。',
'HLCA全细胞与Peng探针法单核平台不同，阈值各自内部确定，百分比不得直接相减。这不是证伪CAS，而是特异性未获支持。SCGB3A2也存在于正常远端气道分泌细胞，不能把所有高分细胞直接定为经典club细胞。',KAC+'。',
'## 挑选口径与纯分数对照',
'“前20”按82池内冻结主排序前20解释，沿用第4–5步eligible_pool_top20；所有82条都有第I步资格，因此此池内消退名单顺序相同。另附全集主排序及消退分层各前20以避免混淆。BCAT1、CD1A/CD207不在全集前20，未声称属于全集前20。本轮全部取自主排序的池内前20，消退名单只作并列信息。',
table(top[['gene','rank_eligible_pool_display','rank_all_display','rank_regression_display','score_main','cohort_suitability_tier','manual_selected','manual_reason']]),
'纯池内分数前5为ETV4、KDELR3、ISG15、AMIGO2、PYCR1。人工选择未替代该基准。XPR1是第一备选；若后续改选须留新的决定记录，不修改冻结名次。',
'既有11条排除记录全部保留，其中82池中的C12orf49、CRTAP、ZDHHC16、TREM2仍排除。第5步范围校正后的“反向=0”不撤销其他队列显著反向或LOO排除。',
'## 两条机制的归组与最弱环节',
'**M1：AT2L氨基酸代谢，PYCR1+BCAT1。** 前者关联脯氨酸合成，后者关联支链氨基酸转氨；依据同一细胞大类与氨基酸代谢功能归组，二者是平行节点，未证明存在单一线性通路或共同依赖。共同可测试环节是改变氨基酸代谢后前驱上皮的存活与修复代价。最弱环节：缺少真正前驱模型，且蛋白定位、代谢通量、选择性均待验证。',
'**M2：朗格汉斯/DC状态，CD1A+CD207。** 以联合定位和功能验证作为一个单元。CD1A单阳性并非朗格汉斯细胞专属，CD207也不能在所有背景下等同同一细胞；先核实双阳性、HLA-DR/CD45和组织空间。共同可研究环节是抗原呈递及上皮–DC关系，当前无证据决定增强还是抑制。最弱环节：细胞积聚可能是保护反应或伴随现象，不能由标志基因上调推出清除细胞有益。',
'其余20条逐项归属见top20_selection_audit.tsv。MMP15/MMP24为同类基质蛋白酶；ETV4/ZNF486同为转录调控大类但不是同一已证通路；其余AT2L候选共享来源不等于共享机制。不将所有上皮基因强行归成一条机制。']
for g in selected:
 r=top.set_index('gene').loc[g]; vr=val[val.gene.eq(g)].iloc[0]; immune=g in ['CD1A','CD207']; cr=cor[cor.gene.eq(g)]; hr=h[h.gene.eq(g)&h.cell_type.isin(['AT2','Club (non-nasal)','DC1','DC2','Migratory DCs'])]
 report += ['## '+g+'：'+units[g],reasons[g],f'来源：82池主排序第{int(r.rank_eligible_pool_display)}；全集主排序第{int(r.rank_all_display)}；消退分层第{int(r.rank_regression_display)}。适切性{r.cohort_suitability_tier}档（甲含GSE164789的AIS/MIA配对支持；丙仅bulk组合且依赖GSE288479部分实性浸润癌的GG区域，不能当成真正前驱病灶独立支持）。临床标签：'+('两者（消退方向线索+阻断方向候选）' if g=='PYCR1' else '消退方向线索')+'；'+LIMIT+'。',
 f'第I步log2FC={r.I_log2FoldChange:.3f}，95%CI [{r.I_log2FC_ci95_low:.3f},{r.I_log2FC_ci95_high:.3f}]，BH q={r.I_padj:.3g}，{int(r.I_n_pairs)}对，方向一致率{r.I_direction_consistency:.3f}；第II步log2FC={r.II_log2FoldChange:.3f}，BH q={r.II_padj:.3g}。第II步上游合格={r.II_upstream_eligible}，未合格者的该数值仅描述，不提升阻断标签。',
 f'冻结后GSE282617第I步counts={vr.GSE282617_I_direction_verdict}，q={vr.GSE282617_I_q_greater_lineage:.4g}；FPKM非参数={vr.GSE282617_I_MW_direction_verdict}，MW q={vr.GSE282617_I_MW_q_greater_lineage:.4g}。第II步对应主张={vr.GSE282617_II_claim_scope}；counts={vr.GSE282617_II_direction_verdict}；MW q={vr.GSE282617_II_MW_q_greater_lineage:.4g}；PYCR1只能写counts方向一致，FPKM非参数核对未达正向标准。空间核对GSE307534：{vr.GSE307534_inference_verdict}；GSE189487：{vr.GSE189487_inference_verdict}，不等于细胞来源/因果验证。',
 f'D3_I={r.D3_I}，独立支持：{r.I_support_datasets}；D4来源比I={r.I_source_ratio:.3g}，D5增殖196命中={r.I_prolif196}；D1′={r.D1_prime:.4f}、D2′={r.D2_prime:.4f}、D3′={r.D3_prime_used:.4f}，冻结总分={r.score_main:.6f}。',
 f'D6人工核对后直接抑制证据={r.D6_inhibitory_evidence_after_manual_review}。新增研究性工具仅作第6步注释，不改变冻结D6，不能视为临床可用药。D8 IHC={r.D8_IHC_gene_reliability}。D9：{r.D9_direction_annotation}；现有癌细胞扰动仅晚期模型，免疫候选方向未知。',
 f'D10：{r.I_D10}；支持状态={r.I_D10_supported_states}；LOO={r.I_LOO_status}，最大P={r.I_LOO_max_p:.3g}，缺失P数={r.I_LOO_missing_p_count}；没有由状态分析证明组成贡献已排除。',
 f'D7器官安全档：{r.D7_safety_assessment_tier}；{r.D7_safety_interpretation}。未触发否决不代表长期干预安全。',
 '**选择性：尚不能确定。** 2.6.2前体校正仅21次置换、未达200次，未检出膨胀不代表消除病灶分布混杂；浸润候选子集校准失败（6.77%，95%CI 5.14–8.39），不作为选择性推断。'+('此候选为Myeloid；AT2斜率及AT2L-other是跨谱系背景，不能当成髓系状态特异检验。' if immune else '需与AT2L-other及健康肺正常AT2共同看待，显著关联不等于差异的差异。'),
 table(cr[['comparison','n_valid_pairs','median_centered_slope_difference','q_BH_within_all_detected_genes']])+'\n（校正版数字列示审计；LUAD项不可作推断，Myeloid项非本谱系选择性证据。）',
 f'AT2L-other前体对正常：log2FC={r.D7_other_Precursor_log2FoldChange:.3f}，95%CI [{r.D7_other_Precursor_effect_ci95_low:.3f},{r.D7_other_Precursor_effect_ci95_high:.3f}]，q={r.D7_other_Precursor_padj:.3g}；IAC对正常log2FC={r.D7_other_LUAD_log2FoldChange:.3f}，q={r.D7_other_LUAD_padj:.3g}。该对照不是健康肺，也不构成候选×状态交互作用检验。',
 table(hr[['cell_type','donors_observed','donor_median_detection','donor_median_mean_log1pCP10k']])+'\nHLCA为供者级描述性参照，检测率不是跨平台风险差，未据此判定安全。',
 'CAS依赖审计：不是自定义四基因签名成员，入选的D1/D3不依赖CAS身份；如后续用CAS高分定位作为支持，该证据受以下限制：'+CAS+'。平台与阈值不可比，不能由百分比差异证伪CAS。'+KAC+'。',
 '四路径如下；预算为人民币量级估计，不是已询价费用，均不含机构管理、样本获取及失败重做。']
 for _,pp in pathdf[pathdf.gene.eq(g)].iterrows():
  report += ['### '+pp['path']+'；最低成本路径：'+pp.lowest_cost]+[f'**{label}**：{pp[col]}' for col,label in [('feasibility','可行性'),('resources','资源'),('endpoints','终点'),('controls','对照'),('interpretation','判读'),('time_cost','周期与成本')]]
 report += ['### 七条降级检查',table(pd.DataFrame(down).query('gene == @g')[['criterion','status','evidence']]),'需验证阶段才能判定：药理效应是否由遗传扰动支持、是否仅减少发生而不能使既有病灶消退。',DECL]
report += ['## IHC/mIF功效与分期标本需求',
'**实际蛋白检测所需最小样本量目前不可直接确定。** 以下是条件性合作规划：假设未来预先指定的患者级连续蛋白终点（IHC H-score、mIF强度或对数密度）具有与RNA配对差相同的标准化效应dz。该假设尚未验证，不能把表中数字称为已验证的IHC功效。效应来自本研究探索性结果，有筛选偏倚及高估风险。',
'方法：AAH/AIS/MIA采用本研究第1步患者配对log2(CPM+0.5)差；IAC仅为本步规划，使用第1步原始伪bulk按各样本全基因总量转CPM，与同患者正常做相同差。dz=均值/配对差SD；95%CI反演非中心t，包含效应和方差的不确定性；双侧α=0.05、power=0.8，配对t功效求解后向上取整并检查前一个整数。按CI下界重算保守数，下界≤0记inf（不能得到有限保守样本量），绝不取绝对值或人为补一个数字。',
'每一期为不同患者，一位患者贡献一块该期病灶和一块配对邻近正常；每名患者仅计一次，不同病灶/视野/细胞不是独立n。正常列为邻近正常组织，不是无癌健康供者。真正健康肺缺少配对蛋白效应，无法给80%功效的最小数；须在技术小试后单独估计。IAC对正常用于描绘阶段谱，不能替代直接IAC对前体的阻断检验。',
 table(pd_power[['gene','stage','exploratory_n','dz','dz_ci95_low','dz_ci95_high','min_patient_pairs','conservative_min_pairs']]),
 table(allocation),
'MIA探索性仅4对，低于发现检验最低8对；功效计算只作不稳定的规划估计，不能新增阶段显著性主张。CD1A/CD207表中的功效针对各单标志的连续表达代理；双阳性DC密度与空间距离是新终点，没有已知效应，需小试重估，不能套用表中n。联合面板每期取四基因需要数的最大值，所有基因用同一患者队列，不能把四套患者数相加。这里只保证各单项名义α=.05下的边际功效，不保证四基因×四期全家族错误率或所有终点同时成功；正式方案须指定唯一主要比较或重新做多重性调整。表内为可分析数，无脱片/失败加成；预计损耗r时招募ceil(n/(1-r))。巨大n或inf表示当前证据不足，不代表临床可招募。先做抗体/方差技术小试，再用盲态蛋白方差重估正式样本量，不改变候选排序。',
'## 明确不能得出的结论',
'1. 不能声称发现病灶发生/侵袭的必要或充分因果事件。\n2. 不能从上调推出应抑制；CD1A/CD207清除可能损害保护性免疫。\n3. 不能声称已建立AAH/AIS/MIA会消退，或已证明阻断其浸润。\n4. 不能将晚期/耐药细胞系、移植瘤或研究性化合物等同早期人体获益。\n5. 不能称CAS四基因评分具有癌状态特异性；也不能说已证伪原文CAS群体。\n6. 不能把健康KAC罕见等同肿瘤特异，或把D7未否决等同安全。\n7. 不能把HLCA与Peng百分比直接相减、将平台队列混杂当疾病差异。\n8. 不能把第II步counts同向写作独立验证；FPKM非参数七项均未达正向标准。\n9. 不能由状态内未显著判定仅细胞组成驱动，亦不能说组成贡献已分解。\n10. 不能把组织学分期等同影像GGN，或将实性区域自动当浸润。\n11. 不能把不同公开数据的重复患者算独立复现。\n12. 不能把RNA规划n当真实蛋白样本量；AAH八对没有支持临床消退结论。\n13. 不能把CD1A与CD207算两个独立机制，也不能把两种代谢酶称为已证统一轴。\n14. 不能由药理效应确定靶向性、由预防肿瘤发生推出治疗已存在病灶。',
'## 数据申请与下一步',
'完整清单及作者请求草案见data_requests.md和controlled_data_request_list.tsv，尚未发送或提交。优先Chen的早期分期与GSE164789健康锚点；Fu扩大总体/随访但不能凭总体n推断AAH充分；HTAN先核对项目与病例重叠，再按所需原始文件申请。',
'## 主要外部依据（原始研究/官方入口）',sources,
'[Fu原始研究与数据登记](https://www.nature.com/articles/s41422-025-01200-w)；[Chen EGA DAC](https://ega-archive.org/dacs/EGAC50000000385)；[HTAN官方访问规则](https://data.humantumoratlas.org/data-access)。',
'## 审计与复现',
'56_dossiers.py只读取本研究50_screen既有结果并生成本步文件。未调用旧results/00–40；第二阶段CAS数值仅显式引用用户核验。候选和排除规则未更改，见QA.json、input_sha256.tsv和output_sha256.tsv。']
full_report='\n\n'.join(report)+'\n'
(O/'step6_dossiers_report.md').write_text(full_report)
for gene in selected:
 start=full_report.index('## '+gene+'：')
 end=full_report.find('\n\n## ',start+3)
 if end<0:end=len(full_report)
 (O/f'dossier_{gene}.md').write_text('# '+gene+' 候选档案\n\n'+full_report[full_report.index('## 证据界限'):full_report.index('## 挑选口径')]+full_report[start:end]+'\n\n'+sources+'\n\n功效、完整限制与申请清单见同目录step6_dossiers_report.md。\n')
# Basic arithmetic and scope checks, without tests that recapitulate prose.
for rr in p_rows:
 for dkey,nkey in [('dz','min_patient_pairs'),('dz_ci95_low','conservative_min_pairs')]:
  d,n=rr[dkey],rr[nkey]
  if np.isfinite(n):
   assert power.power(d,n,alpha=.05)>=.8
   assert n==2 or power.power(d,n-1,alpha=.05)<.8
assert len(pool[pool.top5_excluded_existing])==4
assert set(pool[pool.top5_excluded_existing].gene)=={'C12orf49','CRTAP','ZDHHC16','TREM2'}
assert len(down)==28 and len(paths)==16
assert len(allr[allr.pool_before_D7.eq(True)&allr.top5_excluded_existing.eq(True)])==11
assert all(sha(p)==v for p,v in frozen.items())
qa=dict(frozen_hashes_unchanged=True,frozen_sha256=frozen,pool_n=82,selected_genes=selected,independent_candidate_units=3,mechanism_groups=2,existing_exclusions_in_pool=4,existing_exclusions_in_pool199=11,all_existing_exclusion_records=int(allr.top5_excluded_existing.sum()),four_path_rows=16,downgrade_rows=28,power_minimal_integer_checks_passed=True,protein_effect_transfer_validated=False,CAS_stage2_numbers_source='user-verified statement; not recomputed',selection_scope='pool82 top20, not all5207 top20')
(O/'QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2));save(pd.DataFrame([dict(path=k,sha256=v) for k,v in inputs.items()]),'input_sha256.tsv')
shutil.copy2(__file__,O/'56_dossiers.py')
prose=W/'work/step6/paper_limits.txt'
if prose.exists():shutil.copy2(prose,O/'paper_limitations_ready.md')
save(pd.DataFrame([dict(file=p.name,sha256=sha(p)) for p in sorted(O.iterdir()) if p.is_file() and p.name!='output_sha256.tsv']),'output_sha256.tsv')
print(json.dumps(qa,ensure_ascii=False,indent=2));print(allocation.to_string(index=False))
