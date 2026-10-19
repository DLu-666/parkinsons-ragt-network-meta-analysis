from pathlib import Path
from collections import defaultdict,deque,Counter
import csv,json,math,itertools,hashlib
import numpy as np

W=Path(__file__).resolve().parent;D=W
NAMES={'motor':'Motor impairment','balance':'Balance','mobility':'Mobility','fog':'Freezing of gait','qol':'Quality of life','10MWT':'10MWT','6MWT':'6MWT'}
LEVELS=['No concerns','Some concerns','Major concerns'];GRADES=['High','Moderate','Low','Very low']
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def write(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
ledger=json.loads((D/'extracted_ledger.json').read_text(encoding='utf8'))
studies={str(s['id']):s for s in ledger}
reason={
'1':'H&Y stages 3-4 were more advanced than in several comparator links, a potential modifier of training response.',
'5':'H&Y below stage 3 differed from the more advanced cohorts in other links.',
'8':'Auditory and visual cues accompanied robot training and previous rehabilitation could continue; cointerventions may modify relative effects.',
'11':'The older cohort and ten sessions spread over three months differed from the predominant four-week regimens.',
'12':'The contrast added an exoskeleton to progressive multimodal exercise over eight weeks; it was not an isolated four-week gait-practice contrast.'}
ind={s:(2 if s in reason else 1) for s in studies if s!='6'}
indrows=[{'study':s,'report':studies[s]['label'],'indirectness_score':ind[s],
          'study_judgment':'Moderate' if ind[s]==2 else 'Low',
          'reason':reason.get(s,'Population, gait-training comparison and end-of-intervention outcomes addressed the broad review question; no additional study-level relevance concern was assigned.'),
          'source':studies[s]['source'],'severity':studies[s]['HY'],'device':studies[s]['device'],'dose':studies[s]['dose']} for s in ind]
write(D/'certainty_study_relevance.csv',indrows)
edges=read(W/'study_conductance.csv');risk={(r['study'],r['domain']):r['Overall'] for r in read(D/'ROBUST_RCT_by_outcome.csv')}
results=[r for o in NAMES for r in read(D/'comparative_effects.csv') if r['analysis']=='main' and r['domain']==o]
prediction={(r['domain'],r['treatment1'],r['treatment2']):r for r in read(D/'predictive_intervals.csv')}
checks={(r['domain'],r['treatment1'],r['treatment2']):r for r in read(D/'node_splitting.csv')}
loop_outcomes={r['domain'] for r in checks.values() if r['status']=='Node split fitted'}

def contributions(domain,a,b):
    ee=[e for e in edges if e['domain']==domain]
    nodes=sorted({t for e in ee for t in (e['treatment1'],e['treatment2'])})
    group=defaultdict(list)
    for e in ee:group[tuple(sorted((e['treatment1'],e['treatment2'])))].append(e)
    weights={k:sum(float(e['conductance']) for e in v) for k,v in group.items()}
    L=np.zeros((len(nodes),len(nodes)))
    for (s,t),w in weights.items():
        i,j=nodes.index(s),nodes.index(t);L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
    c=np.zeros(len(nodes));c[nodes.index(a)]=1;c[nodes.index(b)]=-1
    potential=np.linalg.pinv(L)@c;assert np.max(np.abs(L@potential-c))<1e-8
    flow={}
    for (s,t),w in weights.items():
        ff=w*(potential[nodes.index(s)]-potential[nodes.index(t)])
        if abs(ff)>1e-10:flow[(s,t) if ff>0 else (t,s)]=abs(ff)
    ec=defaultdict(float)
    def shortest():
        queue=deque([(a,[])]);seen={a}
        while queue:
            u,path=queue.popleft()
            if u==b:return path
            for v in nodes:
                if v not in seen and flow.get((u,v),0)>1e-9:
                    seen.add(v);queue.append((v,path+[(u,v)]))
    total=0
    while True:
        pp=shortest()
        if pp is None:break
        amount=min(flow[e] for e in pp);total+=amount
        for e in pp:flow[e]-=amount;ec[tuple(sorted(e))]+=amount/len(pp)
    assert abs(total-1)<1e-7
    sc=defaultdict(float)
    for e,cc in ec.items():
        for x in group[e]:sc[x['study']]+=cc*float(x['conductance'])/weights[e]
    assert abs(sum(sc.values())-1)<1e-7
    return sc

def round_level(score):return min(2,max(0,int(math.floor(score+0.5))-1))
def crossings(lo,hi,mu):
    boundary=-.2 if mu>=0 else .2
    return int(lo<0<hi)+int(lo<boundary<hi)

ratings=[];primary=[];contrib=[]
for r in results:
    o,a,b=r['domain'],r['treatment1'],r['treatment2'];key=(o,a,b)
    sc=contributions(*key);rev=contributions(o,b,a)
    assert all(abs(sc[s]-rev[s])<1e-8 for s in set(sc)|set(rev))
    high=sum(v for s,v in sc.items() if risk[(s,o)]=='High Risk')
    risk_score=1+2*high;within=round_level(risk_score)
    indirect_score=sum(ind[s]*v for s,v in sc.items());indirect=round_level(indirect_score)
    for s,cc in sc.items():contrib.append({'domain':o,'treatment1':a,'treatment2':b,'study':s,'report':studies[s]['label'],
                                        'contribution_percent':100*cc,'risk':risk[(s,o)],'indirectness_score':ind[s]})
    mu,lo,hi=(float(r[k]) for k in ['mean','lower','upper']);pr=prediction[key];ns=checks[key]
    plo,phi=float(pr['predictive_lower']),float(pr['predictive_upper'])
    ci_cross=crossings(lo,hi,mu);pi_cross=crossings(plo,phi,mu)
    imp=0 if lo>=-.2 and hi<=.2 else ci_cross
    het=max(0,pi_cross-ci_cross)
    assert pi_cross>=ci_cross,(key,ci_cross,pi_cross)
    if ns['status']=='Node split fitted':
        p=float(ns['tail_probability']);assert p>=.1
        incoh=0;itext=f'Independent CET-TT split: posterior tail probability {p:.3f}; no discrepancy detected, with low power.'
    elif o not in loop_outcomes:
        incoh=2;itext='No independently testable loop; major concern reflects inability to assess incoherence, not a detected contradiction.'
    else:
        incoh=1;itext='No valid comparison-specific split. Another comparison in this network was tested; no global test was available. Some concern retained.'
    reporting=1
    parts=[within,reporting,max(indirect,incoh),max(imp,het)]
    grade=GRADES[min(3,sum(parts))]
    relevant=[studies[s]['label'] for s,v in sc.items() if ind[s]==2 and v>1e-8]
    row={'domain':o,'outcome':NAMES[o],'treatment1':a,'treatment2':b,'direct_studies':int(pr['direct_studies']),
         'mean':mu,'lower':lo,'upper':hi,'predictive_lower':plo,'predictive_upper':phi,
         'high_risk_contribution_percent':100*high,'within_study_score':risk_score,'indirectness_score':indirect_score,
         'within_study_bias':LEVELS[within],'reporting_bias':LEVELS[reporting],'indirectness':LEVELS[indirect],
         'imprecision':LEVELS[imp],'heterogeneity':LEVELS[het],'incoherence':LEVELS[incoh],
         'applied_within':within,'applied_reporting':reporting,'applied_indirectness_incoherence':max(indirect,incoh),
         'applied_imprecision_heterogeneity':max(imp,het),'total_downgrade':sum(parts),'certainty':grade,
         'within_reason':f'High-risk contribution {100*high:.2f}%; CINeMA Average score {risk_score:.3f}.',
         'reporting_reason':'Reporting bias suspected: published English-language evidence; unpublished and registry results have not been comprehensively reconciled. One downgrade; no statistical proof claimed.',
         'indirectness_reason':f'Contribution-weighted relevance score {indirect_score:.3f}. Moderate-relevance contributors: '+(', '.join(relevant) if relevant else 'none')+'. Indirect-only status alone did not cause a downgrade.',
         'imprecision_reason':f'CrI {lo:.3f} to {hi:.3f}; {ci_cross} crossings of zero and the boundary opposite the mean. Reference SMD +/-0.2 was author-approved during revision, not a validated MCID.',
         'heterogeneity_reason':f'Prediction interval {plo:.3f} to {phi:.3f}; {pi_cross} crossings versus {ci_cross} for the CrI. Prediction is conditional on the consistent model.',
         'incoherence_reason':itext,
         'overall_reason':'Randomized evidence began High. Within-study and reporting penalties were added to the greater of indirectness/incoherence and the greater of imprecision/heterogeneity; floor Very low.'}
    ratings.append(row)
    primary.append({'Indicator':NAMES[o],'Comparison':a+' versus '+b,'Number of studies':row['direct_studies'],
                    'Within-study bias':row['within_study_bias'],'Reporting bias':row['reporting_bias'],
                    'Indirectness':row['indirectness'],'Imprecision':row['imprecision'],
                    'Heterogeneity':row['heterogeneity'],'Incoherence':row['incoherence'],'Confidence rating':grade})
assert len(ratings)==47
write(D/'certainty_assessment.csv',ratings);write(D/'CINeMA app input.csv',primary)
write(D/'study_contributions.csv',contrib)
write(D/'study_conductance.csv',edges)
summary={'comparisons':len(ratings),'certainty_counts':dict(Counter(r['certainty'] for r in ratings)),
         'RAGT_comparisons':sum(r['treatment1']=='RAGT' for r in ratings),
         'domain_counts':{key:dict(Counter(r[key] for r in ratings)) for key in ['within_study_bias','reporting_bias','indirectness','imprecision','heterogeneity','incoherence']},
         'untestable_networks':[o for o in NAMES if o not in loop_outcomes],
         'contribution_checks':'Each row sums to 100%; flow conservation and reverse-contrast invariance passed.',
         'rules_sha256':hashlib.sha256((W/'Certainty assessment methods.md').read_bytes()).hexdigest()}
(W/'assessment_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
print(json.dumps(summary,indent=2))
