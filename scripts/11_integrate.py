#!/usr/bin/env python
"""汇总第 3–5 层（上调基因）：results/11_integrated/integrated_layers.tsv、pass_layers_3_4_5.tsv。第 3 层 = DepMap 通过（方案 A）；LINCS 无通过者，不改变结果。
用法：python scripts/11_integrate.py"""
import pandas as pd
P='results'
c=pd.read_csv(f'{P}/07_candidates/candidates_683.tsv',sep='\t',index_col=0)
l3=pd.read_csv(f'{P}/08_layer3_depmap/layer3_depmap.tsv',sep='\t',index_col=0)
l4=pd.read_csv(f'{P}/09_layer4_clinical/layer4_clinical.tsv',sep='\t',index_col=0)
l5=pd.read_csv(f'{P}/10_layer5_safety_druggable/layer5.tsv',sep='\t',index_col=0)
t=c[['log2FoldChange','tier','set','in_main352','direction']].copy()
t['L3']=l3['status']; t['L4']=l4['status']; t['L5']=l5['status']
t['L3_frac_dep_luad']=l3.get('frac_dep_luad'); t['L3_common_essential']=l3.get('common_essential')
t['L4_meta_HR']=l4['meta_HR']; t['L4_meta_q']=l4['meta_q']; t['L4_n_support']=l4['n_support']; t['L4_n_cohorts']=l4['n_cohorts']
for k in ['drug_grade','approved_drugs','other_drugs','type_to_verify','gtex_max_tpm','gtex_max_tissue','gtex_level','hpa_level','hpa_high_tissues','reason']:
    t['L5_'+k]=l5[k]
up=t[t.direction=='up'].copy()
for L in ['L3','L4','L5']: up[L+'_pass']=up[L]=='通过'
up['n_pass_3_5']=up[['L3_pass','L4_pass','L5_pass']].sum(axis=1)
t.loc[up.index,'n_pass_3_5']=up['n_pass_3_5']
t.to_csv(f'{P}/11_integrated/integrated_layers.tsv',sep='\t')
allp=up[up.n_pass_3_5==3].sort_values(['in_main352','tier','L4_meta_q'],ascending=[False,True,True])
allp.to_csv(f'{P}/11_integrated/pass_layers_3_4_5.tsv',sep='\t')
for lab,sub in [('预设352(上调224)',up[up.in_main352]),('并集683(上调434)',up)]:
    print(lab,'通过层数分布',sub['n_pass_3_5'].value_counts().sort_index().to_dict(),'；L3&L4',int((sub.L3_pass&sub.L4_pass).sum()),'；L3&L5',int((sub.L3_pass&sub.L5_pass).sum()),'；L4&L5',int((sub.L4_pass&sub.L5_pass).sum()))
print('第3-5层全部通过：',len(allp),'；其中预设352：',int(allp.in_main352.sum()),'；',list(allp.index))
