#!/opt/anaconda3/envs/ggn_target/bin/python
"""Step4 only: fixed D1/D2/D3 ranks, two simultaneous freezes, no validation/dossiers.
Requires an explicit scoring clarification record before computing ranks.
"""
from pathlib import Path
import argparse,json,hashlib,datetime,os,math
import numpy as np,pandas as pd
B=Path('/Users/dudu/全部数据库/GGN_target');R=B/'results/50_screen';O=R/'step4';O.mkdir(exist_ok=True)
SCRIPT=Path(__file__)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(df,path):df.to_csv(path,sep='\t',index=False,float_format='%.17g')
def order(df,score='score_main'):
 return df.sort_values([score,'gene','lineage'],ascending=[False,True,True],kind='mergesort').copy()
def percentiles(values,eligible,mode):
 # Main: each discovery-step eligible set. Sensitivity: all5207 raw then mask.
 values=pd.to_numeric(values,errors='coerce')
 if mode=='within_step_eligible':
  assert values[eligible].notna().all()
  out=pd.Series(0.,index=values.index);out.loc[eligible]=values[eligible].rank(method='average')/int(eligible.sum());return out
 ref=values if mode=='raw_all_then_mask' else values.where(eligible,0.)
 ref=ref.fillna(-np.inf)
 pct=ref.rank(method='average',pct=True)
 return pct.where(eligible & values.notna(),0.)
def rank_sensitivity(c,score,scope):
 s=order(c[scope],score);return s[['gene','lineage',score]].head(20)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--freeze',action='store_true');args=ap.parse_args()
 decision=O/'scoring_clarification.json'
 if not decision.exists():raise SystemExit('No ranking calculated: scoring clarification is pending. Complete the explicit decision record first.')
 dec=json.loads(decision.read_text())
 assert dec['axis_rule']=='I_for_I_eligible_else_II' and dec['percentile_rule']=='within_step_eligible'
 assert dec.get('explicit_user_confirmation') is True
 if not args.freeze:print(json.dumps(dec,ensure_ascii=False,indent=2));return
 if (O/'freeze_manifest.json').exists():
  m=json.loads((O/'freeze_manifest.json').read_text())
  assert all(sha(x['path'])==x['sha256'] for x in m['frozen_files'])
  print('Existing freeze verified; no new ranking/freeze created.');return
 if list(R.glob('ranking_frozen_*.tsv')):
  raise SystemExit('Existing frozen ranking files detected without this batch manifest; inspect them before creating any additional freeze.')
 src=R/'step3/candidates_step3.tsv';c=pd.read_csv(src,sep='\t',low_memory=False)
 assert len(c)==5207 and not c.duplicated(['gene','lineage']).any()
 assert c.I_upstream_eligible.sum()==672 and c.pool_before_D7.sum()==199 and c.pool_after_D7.sum()==82
 assert (c.pool_after_D7&c.D7_safety_assessment_tier.eq('安全未评估')).sum()==36
 assert not any(x.startswith('D7_KAC_Precursor_') or x.startswith('D7_KAC_LUAD_') for x in c)
 # Pre-score implementation lock. No scores are inspected before this is saved.
 stamp=datetime.datetime.now().astimezone().strftime('%Y%m%dT%H%M%S%z')
 methods={'time':stamp,'protocol_sha256':sha(B/'docs/新研究_靶点优先级筛选_预设方案.md'),'input_sha256':sha(src),'script_sha256':sha(SCRIPT),'decision':dec,'N_percentile_reference':{'D1':672,'D2':4757,'D3_each_axis':5207},'D1_D2':'Original log2FC times positive-patient fraction. Average rank divided by step-specific eligible N:672/4757. Ineligible step zero, no imputation. All5207 raw-then-mask reference is sensitivity only.','D3':'I and histology-II raw counts and percentiles independently retained over5207. Zero support gives percentile0. Select I for I-upstream-eligible records; otherwise II. Never add axes, regional excluded. Current II=0: equivalent to max ONLY after masking axes whose discovery step is ineligible; NOT max of raw axes (2109 II-only records have positive raw I). Future II data never changes the discovery-step rule.','weights_main':[.33,.33,.34],'weights_sensitivity':[[.35,.25,.40],[.50,.20,.30]],'regression':'Filter I-upstream-eligible672 AFTER computing scores; same values/weights/denominator.','eligible_pool':'Existing axis-matched199 and original organ-gate82 only; HPAmissing36 remain eligible but distinctly unassessed. No repair-result/D6/D8/D9 score/gate.','ties':'Average ranks represented as integer half-ranks; denominators672/4757/5207 scaled to LCM; exact weighted integer numerator divided by200*LCM equals the fixed formula and prevents machine-epsilon tie breaking. Equal scores retain equal minimum rank; display order gene then lineage, no clinical tiebreak.','sensitivity_broad_prolif':'Exclude row from eligible pool if it is broad1311 in every axis that presently qualifies it; do not re-estimate percentiles or main scores.','stop':'No frozen-cohort read, no step5/6, no dossiers, no manually chosen<=5.'}
 (O/'score_methods_before_calculation.json').write_text(json.dumps(methods,ensure_ascii=False,indent=2))
 with (R/'run_log.txt').open('a') as f:f.write('\n'+stamp+' STEP4 METHODS LOCK before score computation SHA256='+sha(O/'score_methods_before_calculation.json')+'; no frozen validation.\n')
 for tag,out in [('I','D1_prime'),('II','D2_prime')]:
  prod=c[tag+'_log2FoldChange']*c[tag+'_direction_consistency']
  raw=c[tag+'_D_raw'];ok=prod.notna()&raw.notna()
  assert np.allclose(prod[ok],raw[ok],atol=1e-12,rtol=1e-10)
  c[out]=percentiles(raw,c[tag+'_upstream_eligible'],dec['percentile_rule'])
 for col,out in [('D3_I','D3_I_prime'),('D3_II_histology','D3_II_prime')]:
  c[out]=c[col].rank(method='average',pct=True).where(c[col]>0,0.)
 c['D3_score_axis']=np.where(c.I_upstream_eligible,'I','II_histology')
 c['D3_prime_used']=np.where(c.I_upstream_eligible,c.D3_I_prime,c.D3_II_prime)
 c['D3_raw_used']=np.where(c.I_upstream_eligible,c.D3_I,c.D3_II_histology)
 # Average ranks are integer or half-integer. Compute the same weighted formula
 # with integer numerators to avoid machine-epsilon tie-breaking artefacts.
 counts=np.array([672,4757,5207],dtype=np.int64)
 common=math.lcm(*map(int,counts))
 half_ranks=np.rint(c[['D1_prime','D2_prime','D3_prime_used']].to_numpy()*counts*2).astype(np.int64)
 assert np.allclose(half_ranks/(counts*2),c[['D1_prime','D2_prime','D3_prime_used']].to_numpy(),atol=1e-14,rtol=0)
 scaled=half_ranks*(common//counts)
 c['score_exact_numerator']=scaled@np.array([33,33,34])
 c['score_exact_denominator']=200*common
 c['score_main']=c.score_exact_numerator/c.score_exact_denominator
 assert np.allclose(c.score_main,.33*c.D1_prime+.33*c.D2_prime+.34*c.D3_prime_used,atol=1e-14,rtol=0)
 for name,w in [('score_sensitivity_035_025_040',(35,25,40)),('score_sensitivity_050_020_030',(50,20,30))]:
  c[name]=(scaled@np.array(w))/(200*common)
 c['D1_prime_all_reference']=percentiles(c.I_D_raw,c.I_upstream_eligible,'raw_all_then_mask')
 c['D2_prime_all_reference']=percentiles(c.II_D_raw,c.II_upstream_eligible,'raw_all_then_mask')
 all_half=np.rint(c[['D1_prime_all_reference','D2_prime_all_reference','D3_prime_used']].to_numpy()*5207*2).astype(np.int64)
 c['score_sensitivity_all_reference']=(all_half@np.array([33,33,34]))/(200*5207)
 assert np.isfinite(c.score_main).all() and c.score_main.between(0,1).all()
 assert c.loc[~c.I_upstream_eligible,'D1_prime'].eq(0).all() and c.loc[~c.II_upstream_eligible,'D2_prime'].eq(0).all()
 assert c.D3_II_prime.eq(0).all()
 c['eligible_topN_pool']=c.pool_after_D7
 c['score_rank_min_all']=c.score_main.rank(method='min',ascending=False).astype(int)
 c=order(c);c['rank_all_display']=np.arange(1,len(c)+1)
 ep=order(c[c.eligible_topN_pool]);ep['rank_eligible_pool_display']=np.arange(1,len(ep)+1)
 rm={(r.gene,r.lineage):r.rank_eligible_pool_display for r in ep.itertuples()};c['rank_eligible_pool_display']=[rm.get((r.gene,r.lineage),np.nan) for r in c.itertuples()]
 regression=order(c[c.I_upstream_eligible]);regression['rank_regression_display']=np.arange(1,len(regression)+1)
 rmap={(r.gene,r.lineage):r.rank_regression_display for r in regression.itertuples()};c['rank_regression_display']=[rmap.get((r.gene,r.lineage),np.nan) for r in c.itertuples()]
 regression=order(c[c.I_upstream_eligible]).copy()
 c['in_main_all_top20']=c.rank_all_display.le(20)
 c['in_regression_all_top20']=c.rank_regression_display.le(20)
 c['in_eligible_pool_top20']=c.rank_eligible_pool_display.le(20)
 c['rank_status']='第4步冻结；后续验证不得改变排序'
 regression['rank_status']=c['rank_status'].iloc[0]
 c['freeze_id']=stamp;regression['freeze_id']=stamp
 # All currently eligible rows belong to I; filtering to the regression list cannot change their relative order.
 assert list(zip(ep.gene,ep.lineage))==list(zip(regression[regression.eligible_topN_pool].gene,regression[regression.eligible_topN_pool].lineage))
 assert np.array_equal(c.set_index(['gene','lineage']).loc[regression.set_index(['gene','lineage']).index,'score_main'].values,regression.score_main.values)
 def stops(r):
  why=[]
  if not (r.I_upstream_eligible or r.II_upstream_eligible):why.append('无上游合格步')
  if not r.pool_before_D7:why.append('匹配发现步的D3<1')
  if r.D7_safety_veto==True:why.append('器官安全否决')
  elif r.D7_data_incomplete==True:why.append('关键定量安全数据缺失')
  if r.top5_excluded_opposite:why.append('显著反向：不得top5')
  if r.top5_excluded_I_LOO or r.top5_excluded_II_LOO:why.append('LOO单患者驱动：不得top5')
  return '；'.join(why) if why else '未触发既有门槛/排除；非已完成候选适宜性验证'
 c['gate_or_top5_notes']=c.apply(stops,axis=1)
 regression=c.set_index(['gene','lineage']).loc[regression.set_index(['gene','lineage']).index].reset_index()
 paths=[R/f'ranking_frozen_{stamp}.tsv',R/f'ranking_frozen_regression_{stamp}.tsv']
 # Write both to temporary paths, verify both, then publish as one logged freeze batch.
 for df,p in zip([c,regression],paths):save(df,p.with_suffix('.tmp'))
 for df,p in zip([c,regression],paths):
  chk=pd.read_csv(p.with_suffix('.tmp'),sep='\t',low_memory=False);assert len(chk)==len(df) and chk[['gene','lineage']].equals(df[['gene','lineage']].reset_index(drop=True))
 for p in paths:os.replace(p.with_suffix('.tmp'),p)
 records=[{'path':str(p),'rows':len(df),'sha256':sha(p)} for p,df in zip(paths,[c,regression])]
 manifest={'freeze_id':stamp,'frozen_files':records,'script_sha256':sha(SCRIPT),'methods_sha256':sha(O/'score_methods_before_calculation.json'),'candidate_pool':82,'unassessed_in_pool':36,'both_tables_same_score':True,'no_step5_6':True}
 (O/'freeze_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
 with (R/'freeze_log.txt').open('a') as f:
  for r in records:f.write(stamp+'\t'+r['path']+'\t'+str(r['rows'])+'\tSHA256='+r['sha256']+'\n')
 save(c,R/'candidates_full.tsv')
 save(c.head(20),O/'main_all_top20.tsv');save(regression.head(20),O/'regression_all_top20.tsv')
 save(c[c.eligible_topN_pool].head(20),O/'eligible_pool_top20.tsv');save(c[c.eligible_topN_pool].head(5),O/'eligible_pool_pure_score_top5.tsv')
 save(c.head(5),O/'main_all_pure_score_top5.tsv');save(regression.head(5),O/'regression_all_pure_score_top5.tsv')
 # Sensitivity scopes remain distinct: full rankings vs permitted eligible pool. Never replace main ranks.
 scopes={'main_all':pd.Series(True,index=c.index),'regression_all':c.I_upstream_eligible,'eligible82':c.eligible_topN_pool};sr=[]
 for sc,mask in scopes.items():
  ref=set(map(tuple,order(c[mask]).head(20)[['gene','lineage']].values))
  for name in ['score_sensitivity_035_025_040','score_sensitivity_050_020_030','score_sensitivity_all_reference']:
   alt=order(c[mask],name).head(20);keys=set(map(tuple,alt[['gene','lineage']].values));sr.append(dict(scope=sc,analysis=name,reference_n=len(ref),alternative_n=len(keys),intersection=len(keys&ref),overlap_fraction=len(keys&ref)/len(ref),jaccard=len(keys&ref)/len(keys|ref)))
   save(alt,O/f'{sc}_{name}_top20.tsv')
 broad_axis=(c.I_axis_matched_D3_pass&c.I_prolif1311.eq(False))|(c.II_axis_matched_D3_pass&c.II_prolif1311.eq(False))
 broad_pool=c.eligible_topN_pool&broad_axis;bp=order(c[broad_pool]);save(bp.head(20),O/'broad1311_eligible_top20.tsv');ref=set(map(tuple,order(c[c.eligible_topN_pool]).head(20)[['gene','lineage']].values));alt=set(map(tuple,bp.head(20)[['gene','lineage']].values));sr.append(dict(scope='eligible82',analysis='broad_prolif1311_exclusion_only',reference_n=len(ref),alternative_n=len(alt),intersection=len(ref&alt),overlap_fraction=len(ref&alt)/len(ref),jaccard=len(ref&alt)/len(ref|alt)))
 broad_upstream=(c.I_upstream_eligible&c.I_prolif1311.eq(False))|(c.II_upstream_eligible&c.II_prolif1311.eq(False))
 for sc,mask in [('main_all',broad_upstream),('regression_all',c.I_upstream_eligible&c.I_prolif1311.eq(False))]:
  broad=order(c[mask]).head(20);save(broad,O/f'broad1311_{sc}_top20.tsv')
  original=c if sc=='main_all' else c[c.I_upstream_eligible]
  ref=set(map(tuple,order(original).head(20)[['gene','lineage']].values));alt=set(map(tuple,broad[['gene','lineage']].values))
  sr.append(dict(scope=sc,analysis='broad_prolif1311_exclusion_only',reference_n=len(ref),alternative_n=len(alt),intersection=len(ref&alt),overlap_fraction=len(ref&alt)/len(ref),jaccard=len(ref&alt)/len(ref|alt)))
 save(pd.DataFrame(sr),O/'ranking_sensitivity.tsv')
 save(order(c,'score_sensitivity_all_reference'),O/'all_reference_sensitivity_full.tsv')
 save(c[c.eligible_topN_pool&~broad_pool][['gene','lineage','I_prolif1311','II_prolif1311']],O/'broad1311_removed_from82.tsv')
 summary={'N_all':len(c),'N_regression':len(regression),'N_axis_pool':int(c.pool_before_D7.sum()),'N_organ_pool':82,'N_HPA_unassessed':36,'N_broad_pool':len(bp),'main_top20_eligible':int(c.head(20).eligible_topN_pool.sum()),'regression_top20_eligible':int(regression.head(20).eligible_topN_pool.sum()),'main_top20_top5_exclusion':int(c.head(20).top5_excluded_existing.sum()),'regression_top20_top5_exclusion':int(regression.head(20).top5_excluded_existing.sum()),'eligible_top20_top5_exclusion':int(c[c.eligible_topN_pool].head(20).top5_excluded_existing.sum()),'score_changed_in_regression':0,'freeze_id':stamp,'manually_selected_candidates':0,'dossiers_created':0,'no_frozen_cohort_read':True}
 (O/'step4_QA.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
 with (R/'run_log.txt').open('a') as f:f.write('\n'+stamp+' STEP4 FREEZE COMPLETE all5207/regression672, two SHA256 in freeze_log; pool82 preserved incl36 unassessed; sensitivities report only; no manual top5/dossiers/frozen-cohort analysis.\n')
 print(json.dumps(manifest,ensure_ascii=False,indent=2));print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
