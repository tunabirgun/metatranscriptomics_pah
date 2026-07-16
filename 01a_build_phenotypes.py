# Build per-cohort phenotype tables (GSM -> group) from GEO series matrices.
# Applies locked eligibility rules; verifies counts against the registered inventory.
import gzip, glob, os, csv, re
RAW="02_acquisition/raw"; OUT="03_preprocessing/pheno"; os.makedirs(OUT, exist_ok=True)
def sm(gse):
    p=glob.glob(f"{RAW}/*/{gse}/{gse}_series_matrix.txt.gz")[0]
    d={}
    with gzip.open(p,'rt',encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('!Sample_'):
                k=line.split('\t',1)[0]; v=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
                d.setdefault(k,[]).append(v)
    return d
def col(d,key,idx=0): return d[key][idx]
def write(gse, gsms, titles, groups, subtypes, reasons):
    with open(f"{OUT}/{gse}_pheno.tsv","w",newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(["gsm","title","group","subtype","reason"])
        for a in zip(gsms,titles,groups,subtypes,reasons): w.writerow(a)
    keep=[g for g in groups if g in ("PAH","control")]
    from collections import Counter
    c=Counter(groups)
    print(f"{gse}: PAH={c['PAH']} control={c['control']} drop={sum(v for k,v in c.items() if k not in ('PAH','control'))}  -> {dict(c)}")

# GSE113439
d=sm("GSE113439"); ds=col(d,'!Sample_characteristics_ch1',1)
g=[];sub=[];rs=[]
for x in ds:
    x=x.replace('disease state:','').strip()
    if 'normal control' in x: g.append('control');sub.append('control');rs.append('')
    elif 'CTEPH' in x: g.append('drop_CTEPH');sub.append('CTEPH');rs.append('Group4 CTEPH')
    elif 'idiopathic PAH' in x: g.append('PAH');sub.append('IPAH');rs.append('')
    elif 'CHD' in x: g.append('PAH');sub.append('CHD-PAH');rs.append('')
    elif 'CTD' in x: g.append('PAH');sub.append('CTD-PAH');rs.append('')
    else: g.append('drop_other');sub.append(x);rs.append('unmapped')
write("GSE113439",col(d,'!Sample_geo_accession'),col(d,'!Sample_title'),g,sub,rs)

# GSE117261
d=sm("GSE117261"); cg=col(d,'!Sample_characteristics_ch1',1); st=col(d,'!Sample_characteristics_ch1',2)
g=[];sub=[];rs=[]
for c,s in zip(cg,st):
    s=s.replace('pah_subtype:','').strip(); c=c.replace('clinical_group:','').strip()
    if c=='FD': g.append('control');sub.append('FD');rs.append('')
    elif s in ('IPAH','APAH','FPAH'): g.append('PAH');sub.append(s);rs.append('')
    else: g.append('drop_nonGrp1');sub.append(s);rs.append('Other/WHO4 non-Grp1')
write("GSE117261",col(d,'!Sample_geo_accession'),col(d,'!Sample_title'),g,sub,rs)

# GSE53408
d=sm("GSE53408"); sg=col(d,'!Sample_characteristics_ch1',0)
g=[];sub=[];rs=[]
for x in sg:
    if 'normal control' in x: g.append('control');sub.append('control');rs.append('')
    else: g.append('PAH');sub.append('PAH');rs.append('')
write("GSE53408",col(d,'!Sample_geo_accession'),col(d,'!Sample_title'),g,sub,rs)

# GSE15197
d=sm("GSE15197"); ts=col(d,'!Sample_characteristics_ch1',0)
g=[];sub=[];rs=[]
for x in ts:
    if 'Normal lung' in x: g.append('control');sub.append('control');rs.append('')
    elif 'idiopathic pulmonary fibrosis' in x: g.append('drop_Grp3');sub.append('PH-IPF');rs.append('Group3 PH-IPF')
    elif 'pulmonary arterial hypertension' in x: g.append('PAH');sub.append('PAH');rs.append('')
    else: g.append('drop_other');sub.append(x[:20]);rs.append('unmapped')
write("GSE15197",col(d,'!Sample_geo_accession'),col(d,'!Sample_title'),g,sub,rs)

# GSE48149 (source_name_ch1)
d=sm("GSE48149"); src=col(d,'!Sample_source_name_ch1')
g=[];sub=[];rs=[]
for x in src:
    x=x.strip()
    if x=='NL': g.append('control');sub.append('control');rs.append('')
    elif x in ('IPAH (PPH)','SSc-PAH'): g.append('PAH');sub.append('IPAH' if 'PPH' in x else 'SSc-PAH');rs.append('')
    else: g.append('drop_fibrosis');sub.append(x);rs.append('IPF/SSc-PF non-Grp1')
write("GSE48149",col(d,'!Sample_geo_accession'),col(d,'!Sample_title'),g,sub,rs)

# GSE254617
d=sm("GSE254617"); cond=col(d,'!Sample_characteristics_ch1',2); st=col(d,'!Sample_characteristics_ch1',3)
g=[];sub=[];rs=[]
for c,s in zip(cond,st):
    c=c.replace('condition:','').strip(); s=s.replace('ph subtype:','').strip()
    if c=='FD': g.append('control');sub.append('FD');rs.append('')
    elif s in ('IPAH','APAH','HPAH'): g.append('PAH');sub.append(s);rs.append('')
    elif s=='PVOD': g.append('PAH');sub.append('PVOD');rs.append('Group1prime flagged')
    else: g.append('drop_nonGrp1');sub.append(s);rs.append('Unknown/WHOGRP4')
write("GSE254617",col(d,'!Sample_geo_accession'),col(d,'!Sample_title'),g,sub,rs)

# GSE272776
d=sm("GSE272776"); dis=col(d,'!Sample_characteristics_ch1',3)
g=[];sub=[];rs=[]
for x in dis:
    x=x.replace('disease:','').strip()
    g.append('PAH' if x=='PAH' else 'control');sub.append(x);rs.append('')
write("GSE272776",col(d,'!Sample_geo_accession'),col(d,'!Sample_title'),g,sub,rs)

# GSE208592 (VALIDATION)
d=sm("GSE208592"); cond=col(d,'!Sample_characteristics_ch1',1)
g=[];sub=[];rs=[]
for x in cond:
    x=x.replace('condition:','').strip()
    g.append('PAH' if x=='PAH' else 'control');sub.append(x);rs.append('')
write("GSE208592",col(d,'!Sample_geo_accession'),col(d,'!Sample_title'),g,sub,rs)
print("\nInventory targets: 113439=14/11 117261=54/25 53408=12/11 15197=18/13 48149=18/9 254617=88/52 272776=8/8 | val 208592=15/18")
