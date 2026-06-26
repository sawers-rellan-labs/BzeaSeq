#!/usr/bin/env python
"""27M-wideseq Jaccard donor search for ONE NIL (whole genome), missing-aware.

Windows from calls_taxa_r5 (state>0); NIL genotypes from the 27.6M-position allelic counts
(alt read >=1 -> alt-present=1, ref-only -> 0, uncovered dropped); panel = wideseq_ref
(217 teosinte accessions, ~60% missing) extracted per window via bcftools -r. Per block:
missing-aware symmetric Jaccard of non-B73 presence vs each accession -> nearest; aggregate
genome-wide (weighted by block Mb) -> taxa composition.

Usage: ~/anaconda3/bin/python agent/jaccard_27m_one.py PN14_SID1284
"""
import sys, subprocess, time, csv
import numpy as np
from collections import defaultdict

name = sys.argv[1]
B = "/Volumes/BZea/bzeaseq"
COUNTS = f"{B}/allelic_counts/{name}.allelicCounts.tsv"
WREF = B + "/wideseq_ref/wideseq_chr{}.vcf.gz"
CALLS = "data/rtiger_50K/calls_taxa_r5.csv"
META = "data/donor_id/refpanel_metadata.csv"
MIN_SNP = 50
t0 = time.time()

# panel samples (VCF column order) + taxa
samps = subprocess.run(["bcftools","query","-l",WREF.format(1)],capture_output=True,text=True).stdout.split()
meta = {r["sample"]: r for r in csv.DictReader(open(META))}
taxon = np.array([meta.get(s,{}).get("maizegdb_prefix","NA") for s in samps])
print(f"panel: {len(samps)} accessions")

# windows for this NIL (state>0)
import pandas as pd
calls = pd.read_csv(CALLS)
W = calls[(calls.name==name)&(calls.state>0)][["chr","start_bp","end_bp"]].values.tolist()
known = None
# infer known donor accession/taxon from counts dir is N/A here; read from calls 'donor' col
kd = calls[(calls.name==name)]["donor"].iloc[0]
print(f"{name}: {len(W)} donor blocks | known donor taxon = {kd}")

# --- NIL genotypes in the windows: one awk pass over the 27.6M-line counts ---
win_lines = "".join(f"chr{int(c)}\t{int(s)}\t{int(e)}\n" for c,s,e in W)
awk = r'''
NR==FNR{c[FNR]=$1;s[FNR]=$2;e[FNR]=$3;nw=FNR;next}
/^@/||$1=="CONTIG"{next}
{for(i=1;i<=nw;i++) if($1==c[i]&&$2>=s[i]&&$2<=e[i]){tot=$3+$4; if(tot>0) print $1"\t"$2"\t"($4>=1?1:0); break}}
'''
p = subprocess.run(["awk","-F","\t",awk,"/dev/stdin",COUNTS],input=win_lines,capture_output=True,text=True)
q = {}
for line in p.stdout.splitlines():
    c,pos,a = line.split("\t"); q[(c,int(pos))] = int(a)
print(f"NIL covered SNPs across windows: {len(q):,}  ({time.time()-t0:.0f}s)")

# --- per-block panel extract + missing-aware Jaccard ---
by_tax, by_acc = defaultdict(float), defaultdict(float)
def parse_gt(g):
    g=g.split(':')[0]
    if '.' in g: return -1
    return 1 if '1' in g else 0
for c,s,e in W:
    cs=f"chr{int(c)}"; mb=(e-s)/1e6
    rec = subprocess.run(["bcftools","view","-H","-r",f"{cs}:{int(s)}-{int(e)}",WREF.format(int(c))],
                         capture_output=True,text=True).stdout.splitlines()
    Pq=[]; Prow=[]
    for ln in rec:
        f=ln.split("\t"); pos=int(f[1]); key=(cs,pos)
        if key not in q: continue
        Pq.append(q[key]); Prow.append([parse_gt(g) for g in f[9:]])
    if len(Pq) < MIN_SNP:
        continue
    qv=np.array(Pq); P=np.array(Prow)            # (L,), (L x 217)
    q1=(qv==1); a1=(P==1); valid=(P!=-1)
    inter=(q1[:,None]&a1&valid).sum(0); union=((q1[:,None]|a1)&valid).sum(0)
    J=np.where(union>0,inter/union,0.0); best=int(J.argmax())
    by_tax[taxon[best]]+=mb; by_acc[samps[best]]+=mb
    print(f"  chr{int(c)}:{int(s/1e6)}-{int(e/1e6)}Mb ({mb:.0f}Mb, {len(Pq)} SNPs) -> {samps[best]} [{taxon[best]}] J={J[best]:.3f}")

tot=sum(by_tax.values())
print(f"\n=== {name} (known Zd, λ high; 50K wrongly called Zl) — 27M result ===")
if tot==0: print("no callable blocks"); sys.exit()
comp={k:round(v/tot,3) for k,v in sorted(by_tax.items(),key=lambda x:-x[1])}
dom=max(by_tax,key=by_tax.get)
print(f"taxa composition (Mb-weighted): {comp}")
print(f"dominant taxon: {dom}  | known: {kd}  | MATCH: {dom==kd}")
print(f"top individual: {max(by_acc,key=by_acc.get)}  | {time.time()-t0:.0f}s total")
