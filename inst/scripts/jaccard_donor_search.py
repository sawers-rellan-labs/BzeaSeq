#!/usr/bin/env python
"""Jaccard donor search on SNP50K (Dilkes WideSeq method) — with confidence filtering.

Per NIL: ancestry calls (calls_taxa_r5) define introgression windows; within each, the NIL's
per-SNP genotype (alt observed=1, ref-only=0, uncovered=NA) is symmetric-Jaccard-matched vs the
dense 217-accession teosinte 50K panel. Per block we record the best accession, its taxon, the
best J, and the MARGIN (best J minus the best J of any *other* taxon). A block is "confident" if
best J >= FLOOR and margin >= MARGIN; a NIL is "resolved" if it has >=1 confident block, and its
dominant taxon is taken over confident blocks (Mb-weighted). Validate vs the pedigree donor taxon.

Usage: ~/anaconda3/bin/python agent/jaccard_donor_search.py
"""
import numpy as np, pandas as pd, glob, os, csv
from collections import defaultdict

PANEL="data/donor_id/refpanel_gt.tsv"; SAMP="data/donor_id/refpanel_samples.txt"
META="data/donor_id/refpanel_metadata.csv"; CALLS="data/rtiger_50K/calls_taxa_r5.csv"
COUNTS="data/rtiger_50K/counts"; OUT="results/donor_id"; MAF, MIN_SNP = 0.05, 20
FLOOR, MARGIN = 0.40, 0.05            # default confidence thresholds (see grid below)

samps=[l.strip() for l in open(SAMP)]
meta={r["sample"]:r for r in csv.DictReader(open(META))}
taxon_of=np.array([meta.get(s,{}).get("maizegdb_prefix","NA") for s in samps])
chrom,pos,rows=[],[],[]
for line in open(PANEL):
    c=line.rstrip("\n").split("\t"); chrom.append(c[0]); pos.append(int(c[1]))
    rows.append([1 if "1" in g else 0 for g in c[2:]])
A=np.array(rows,dtype=np.int8); chrom=np.array(chrom); pos=np.array(pos)
posidx={(chrom[i],pos[i]):i for i in range(len(pos))}
afreq=A.mean(1); informative=(afreq>=MAF)&(afreq<=1-MAF)
print(f"panel {A.shape[0]} SNPs x {A.shape[1]} accessions | informative {informative.sum()}")
name2={os.path.basename(p)[:-4]:(os.path.basename(os.path.dirname(p)),p) for p in glob.glob(f"{COUNTS}/*/*.tsv")}
calls=pd.read_csv(CALLS)

def load_query(path):
    q=np.full(A.shape[0],-1,dtype=np.int8)
    for line in open(path):
        c=line.rstrip("\n").split("\t"); i=posidx.get((c[0],int(c[1])))
        if i is None: continue
        ref,alt=int(c[3]),int(c[5])
        if alt>=1: q[i]=1
        elif ref>=1: q[i]=0
    return q

def search_nil(name):
    donor_acc,path=name2[name]; known=donor_acc[:2]; q=load_query(path)
    blks=[]
    for _,b in calls[(calls.name==name)&(calls.state>0)].iterrows():
        cstr=f"chr{int(b.chr)}"
        sel=(chrom==cstr)&(pos>=b.start_bp)&(pos<=b.end_bp)&informative&(q!=-1)
        L=int(sel.sum())
        if L<MIN_SNP: continue
        q1=(q[sel]==1); Asub=(A[sel]==1)
        inter=(q1[:,None]&Asub).sum(0); union=(q1[:,None]|Asub).sum(0)
        with np.errstate(invalid="ignore"): J=np.where(union>0,inter/union,0.0)
        best=int(J.argmax()); bt=taxon_of[best]; bj=float(J[best])
        other=J.copy(); other[taxon_of==bt]=-1; margin=bj-float(other.max())
        blks.append(dict(taxon=bt, acc=samps[best], J=bj, margin=margin,
                         mb=(b.end_bp-b.start_bp)/1e6, chr=int(b.chr),
                         start=int(b.start_bp), end=int(b.end_bp), nsnp=L))
    return known, blks

def dominant(blks, conf=False):
    by=defaultdict(float)
    for x in blks:
        if conf and not (x["J"]>=FLOOR and x["margin"]>=MARGIN): continue
        by[x["taxon"]]+=x["mb"]
    if not by: return None,0.0,0.0
    tot=sum(by.values()); dom=max(by,key=by.get); return dom, by[dom]/tot, tot

os.makedirs(OUT,exist_ok=True)
names=[n for n in sorted(name2) if ((calls.name==n)&(calls.state>0)).any()]
print(f"full run: {len(names)} NILs", flush=True)
allblk=[]; summ=[]
for i,n in enumerate(names):
    known,blks=search_nil(n)
    for x in blks: allblk.append({**x,"name":n,"known":known})
    udom,ufrac,_=dominant(blks,False)
    cdom,cfrac,cmb=dominant(blks,True)
    summ.append(dict(name=n,known_taxon=known,n_blocks=len(blks),
                     dom_taxon=udom,dom_frac=round(ufrac,3),
                     uncond_match=(udom==known) if udom else None,
                     resolved=(cdom is not None),conf_dom_taxon=cdom,
                     conf_frac=round(cfrac,3),conf_match=(cdom==known) if cdom else None,
                     mean_bestJ=round(np.mean([x["J"] for x in blks]),3) if blks else 0))
    if (i+1)%300==0: print(f"  {i+1}/{len(names)}",flush=True)

bdf=pd.DataFrame(allblk); sdf=pd.DataFrame(summ)
bdf.to_csv(f"{OUT}/nil_block_matches.csv",index=False); sdf.to_csv(f"{OUT}/nil_donor_summary.csv",index=False)

# block-level J/margin distribution split by correctness (calibration)
bdf["correct"]=bdf.taxon==bdf.known
print("\nblock best-J by correctness (quantiles):")
for lab,g in [("correct",bdf[bdf.correct]),("wrong",bdf[~bdf.correct])]:
    print(f"  {lab:7} J  p10/50/90 = {g.J.quantile(.1):.2f}/{g.J.quantile(.5):.2f}/{g.J.quantile(.9):.2f}"
          f"  margin p50 = {g.margin.quantile(.5):.3f}  n={len(g)}")

# unconditional vs confidence-filtered match rate
u=sdf[sdf.uncond_match.notna()]
print(f"\nUNCONDITIONAL: dom==known {int(u.uncond_match.sum())}/{len(u)} ({u.uncond_match.mean():.0%})")
print("\nCONFIDENCE grid (resolved NILs ; match% among resolved):")
print(f"{'FLOOR':>6}{'MARGIN':>7}{'resolved':>10}{'match%':>8}")
for fl in (0.3,0.4,0.5):
    for mg in (0.02,0.05,0.10):
        m=bdf.copy();
        conf=m[(m.J>=fl)&(m.margin>=mg)]
        # per-NIL dominant over confident blocks
        d=conf.groupby(["name","known"]).apply(lambda g: g.groupby("taxon").mb.sum().idxmax(), include_groups=False)
        if len(d)==0: print(f"{fl:>6}{mg:>7}{0:>10}{'-':>8}"); continue
        dd=d.reset_index(name="cdom"); dd["match"]=dd.cdom==dd.known
        print(f"{fl:>6}{mg:>7}{len(dd):>10}{dd.match.mean():>7.0%}")

r=sdf[sdf.resolved]
print(f"\n@ default FLOOR={FLOOR} MARGIN={MARGIN}: resolved {len(r)}/{len(sdf)} ({len(r)/len(sdf):.0%}) | "
      f"match among resolved {r.conf_match.mean():.0%}")
print("by known taxon (resolved n ; match% among resolved):")
for t,g in r.groupby("known_taxon"):
    print(f"  {t}: resolved {len(g):4d}  match {g.conf_match.mean():.0%}")
print(f"\nwrote {OUT}/nil_donor_summary.csv + nil_block_matches.csv")
