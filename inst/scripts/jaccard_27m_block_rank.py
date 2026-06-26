#!/usr/bin/env python
"""Rank ALL 217 panel accessions for the largest donor block of one NIL (27M wideseq),
to distinguish mislabel vs panel-gap. Prints top hits + where the known-taxon accessions fall.
Usage: ~/anaconda3/bin/python agent/jaccard_27m_block_rank.py PN14_SID1284 Zd
"""
import sys, subprocess, csv
import numpy as np, pandas as pd

name, known = sys.argv[1], sys.argv[2]
B="/Volumes/BZea/bzeaseq"
COUNTS=f"{B}/allelic_counts/{name}.allelicCounts.tsv"
WREF=B+"/wideseq_ref/wideseq_chr{}.vcf.gz"
samps=subprocess.run(["bcftools","query","-l",WREF.format(1)],capture_output=True,text=True).stdout.split()
meta={r["sample"]:r for r in csv.DictReader(open("data/donor_id/refpanel_metadata.csv"))}
taxon=np.array([meta.get(s,{}).get("maizegdb_prefix","NA") for s in samps])

calls=pd.read_csv("data/rtiger_50K/calls_taxa_r5.csv")
blk=calls[(calls.name==name)&(calls.state>0)].assign(mb=lambda d:(d.end_bp-d.start_bp)/1e6).sort_values("mb",ascending=False).iloc[0]
c,s,e=int(blk.chr),int(blk.start_bp),int(blk.end_bp)
print(f"{name}: largest block chr{c}:{s/1e6:.0f}-{e/1e6:.0f}Mb ({blk.mb:.0f}Mb), known donor taxon {known}")

# NIL covered SNPs in the block (early-exit awk on this chr)
awk=r'''/^@/||$1=="CONTIG"{next} $1==CH && $2>=S && $2<=E{tot=$3+$4; if(tot>0)print $2"\t"($4>=1?1:0)} $1==CH && $2>E{exit}'''
p=subprocess.run(["awk","-F","\t","-v",f"CH=chr{c}","-v",f"S={s}","-v",f"E={e}",awk,COUNTS],capture_output=True,text=True)
q={int(l.split("\t")[0]):int(l.split("\t")[1]) for l in p.stdout.splitlines()}
print(f"NIL covered SNPs in block: {len(q):,}")

# panel for the region
rec=subprocess.run(["bcftools","view","-H","-r",f"chr{c}:{s}-{e}",WREF.format(c)],capture_output=True,text=True).stdout.splitlines()
def gt(g):
    g=g.split(':')[0]; return -1 if '.' in g else (1 if '1' in g else 0)
Pq,Prow=[],[]
for ln in rec:
    f=ln.split("\t"); pos=int(f[1])
    if pos not in q: continue
    Pq.append(q[pos]); Prow.append([gt(x) for x in f[9:]])
qv=np.array(Pq); P=np.array(Prow)
q1=(qv==1); a1=(P==1); valid=(P!=-1)
inter=(q1[:,None]&a1&valid).sum(0); union=((q1[:,None]|a1)&valid).sum(0)
ncov=valid.sum(0)
J=np.where(union>0,inter/union,0.0)
df=pd.DataFrame({"accession":samps,"taxon":taxon,"J":J.round(3),"n_covalid":ncov}).sort_values("J",ascending=False).reset_index(drop=True)
print(f"co-observed SNPs (block): {q1.sum()} NIL non-B73 of {len(qv)} covered\n")
print("=== TOP 15 accessions ===")
print(df.head(15).to_string())
print(f"\n=== all known-taxon ({known}) accessions, with rank ===")
kd=df[df.taxon==known].copy(); kd["rank"]=kd.index+1
print(kd.to_string())
print(f"\nbest {known} J = {kd.J.max():.3f} (rank {int(kd['rank'].min())}) | top overall {df.taxon[0]} {df.accession[0]} J={df.J[0]:.3f}")
