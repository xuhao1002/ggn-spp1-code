import numpy as np, pandas as pd
from statsmodels.duration.hazard_regression import PHReg
from scipy import stats
R='/Users/dudu/全部数据库/GGN_target/raw/'
def stage3(s):
    s=str(s)
    for k,v in [('IV','III-IV'),('III','III-IV'),('II','II'),('I','I')]:
        if f'Stage {k}' in s: return v
    return np.nan
def cox(df,time,event,covs,tag,one_sided=True):
    d=df.dropna(subset=[time,event,'stas']+covs).copy(); d=d[d[time]>0]
    X=pd.DataFrame({'stas':d.stas.astype(float)})
    for c in covs:
        if d[c].dtype==object or str(d[c].dtype)=='category':
            X=X.join(pd.get_dummies(d[c],prefix=c,drop_first=True).astype(float))
        else: X[c]=d[c].astype(float)
    m=PHReg(d[time].values,X.values,status=d[event].values.astype(int)).fit(disp=0)
    b,se=m.params[0],m.bse[0]; z=b/se; p1=stats.norm.sf(z)
    print(f'{tag}: n={len(d)} 事件={int(d[event].sum())} STAS+={int(d.stas.sum())} HR={np.exp(b):.2f} (95%CI {np.exp(b-1.96*se):.2f}-{np.exp(b+1.96*se):.2f}) 单侧P={p1:.3f} 双侧P={2*p1 if p1<.5 else 2*(1-p1):.3f}')
    return b,se
# TCGA
lt=pd.read_excel(R+'STAS_labels/STAS_labels_TCGA_from_PanLR_20261011.xlsx'); y=lt.groupby('患者编号').label.agg(['max','size'])
cl=pd.read_csv(R+'TCGA_LUAD/TCGA-LUAD.clinical.tsv.gz',sep='\t'); cl=cl[cl['sample_type.samples']=='Primary Tumor'].drop_duplicates('submitter_id').set_index('submitter_id')
sv=pd.read_csv(R+'TCGA_LUAD/TCGA-LUAD.survival.tsv.gz',sep='\t').drop_duplicates('_PATIENT').set_index('_PATIENT')
t=pd.DataFrame({'stas':y['max'],'ns':y['size']}).join(sv[['OS.time','OS']]).join(cl[['age_at_index.demographic','gender.demographic','ajcc_pathologic_stage.diagnoses']])
t['stage']=t['ajcc_pathologic_stage.diagnoses'].map(stage3); t['age']=t['age_at_index.demographic']; t['sex']=(t['gender.demographic']=='male').astype(int)
print('== TCGA OS')
a=cox(t,'OS.time','OS',['stage','age','sex'],'TCGA 全部分期 校正分期/年龄/性别')
tI=t[t.stage=='I']; b1=cox(tI,'OS.time','OS',['age','sex'],'TCGA 仅I期 校正年龄/性别')
cox(t,'OS.time','OS',[],'TCGA 未校正')
cox(t[t.ns==1],'OS.time','OS',['stage','age','sex'],'TCGA 仅单切片 校正')
# CPTAC
lc=pd.read_excel(R+'STAS_labels/STAS_labels_CPTAC_from_PanLR_20261011.xlsx'); yc=lc.groupby('患者编号').STAS.agg(['max','size'])
c=pd.read_csv(R+'STAS_clinical/CPTAC-3.clinical.tsv.gz',sep='\t',dtype=str).drop_duplicates('submitter_id').set_index('submitter_id')
u=pd.DataFrame({'stas':yc['max'],'ns':np.minimum(yc['size'],4)}).join(c[['vital_status.demographic','days_to_death.demographic','days_to_last_follow_up.diagnoses','ajcc_pathologic_stage.diagnoses','gender.demographic','progression_or_recurrence.diagnoses','days_to_recurrence.diagnoses']])
u['dead']=(u['vital_status.demographic']=='Dead').astype(float)
u.loc[~u['vital_status.demographic'].isin(['Dead','Alive']),'dead']=np.nan
u['time']=pd.to_numeric(np.where(u.dead==1,u['days_to_death.demographic'],u['days_to_last_follow_up.diagnoses']),errors='coerce')
u['stage']=u['ajcc_pathologic_stage.diagnoses'].map(stage3); u['sex']=(u['gender.demographic']=='male').astype(int)
print('== CPTAC OS'); print('可用',u.dropna(subset=['time','dead','stage']).shape[0],'随访中位(天)',u.time.median())
e=cox(u,'time','dead',['stage','sex','ns'],'CPTAC 校正分期/性别/切片数')
cox(u,'time','dead',[],'CPTAC 未校正')
cox(u[u.ns==1],'time','dead',['stage','sex'],'CPTAC 仅单切片 校正')
# 复发
u['rec']=(u['progression_or_recurrence.diagnoses']=='yes').astype(float); u.loc[u['progression_or_recurrence.diagnoses']=='not reported','rec']=np.nan
u['rtime']=pd.to_numeric(np.where(u.rec==1,u['days_to_recurrence.diagnoses'],u['days_to_last_follow_up.diagnoses']),errors='coerce')
print('== CPTAC 复发（已报告者）'); cox(u,'rtime','rec',['stage','sex','ns'],'CPTAC 复发 校正')
# 合并 OS
(b_t,s_t),(b_c,s_c)=a,e; w=np.array([1/s_t**2,1/s_c**2]); bp=(w*np.array([b_t,b_c])).sum()/w.sum(); sp=w.sum()**-.5
print(f'== 固定效应合并 OS: HR={np.exp(bp):.2f} (95%CI {np.exp(bp-1.96*sp):.2f}-{np.exp(bp+1.96*sp):.2f}) 单侧P={stats.norm.sf(bp/sp):.3f}')
