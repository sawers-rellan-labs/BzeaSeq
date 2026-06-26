# Donor identification by Jaccard similarity to a teosinte panel

*Lab notebook — BZea NIL donor-of-origin without founder sequences. (2026-06-25)*

> **Provenance.** This analysis was run from the zealtiger NIL-genotyping workspace (ancestry
> calls + allelic counts staged there and on the `/Volumes/BZea/` mount); it lives here in
> **`sawers-rellan-labs/BzeaSeq`** because the wideseq / ancestry-calling lineage is the donor-ID
> home. Scripts: `inst/scripts/jaccard_donor_search.py`, `jaccard_27m_one.py`,
> `jaccard_27m_block_rank.py`. Outputs: `results/donor_id/`. Input data paths in **Setup** below
> are as run (zealtiger `data/…` + the mount); they are not vendored into this repo.

## Precedent

Brian Dilkes's **WideSeq** preprint (Khangura, Kaur, San Miguel, Dilkes; DOI
**10.64898/2026.02.20.707111**; local `data/2026.02.20.707111.source.xml`). For an introgression
of *unknown origin*, they binarize the non-reference SNP calls in the window (non-B73 = 1, ref = 0)
and compute the **Jaccard index** against every accession in a **public panel** (maize HapMap3,
1,210 accessions); the highest-Jaccard accession is the donor haplotype. Worked examples:
b94 introgression → **LH52 98.2%** (164 SNPs); Z004E0189 → CML247 98.7%. Two flavors: Fig 10
"non-B73 calls only" (containment) and Fig 11 "both positions." Crucially it needs **no founder
sequence or panel of one's own** — only a public set of *known genotypes* to match against.

## Motivation

We want the **teosinte donor per introgression** for the BZea NILs, but we have **no donor founder
sequences and no founder panel of our own**. The GOOGA composite-pseudo-reference approach
(PLoS Genet pgen.1011072) needs per-founder assemblies → not available. Jaccard-to-a-public-panel
is the route that fits our resources: match called introgression SNPs to the teosinte reference
panel we already have.

## Setup

**Query (per NIL):** ancestry calls (`data/rtiger_50K/calls_taxa_r5.csv`, RTIGER per-taxon r=5)
define the introgression windows (state>0). Within each window the NIL's **per-SNP genotype** is
the query — alt allele observed = 1, ref-only = 0, uncovered = excluded. (Not bin allele-frequency:
Jaccard is a set operation on per-SNP presence.)

- 50K query: `data/rtiger_50K/counts/<donor>/<name>.tsv` (ref/alt counts at ~49K array SNPs).
- 27M query: `/Volumes/BZea/bzeaseq/allelic_counts/<name>.allelicCounts.tsv` (GATK
  CollectAllelicCounts at the 27.6M teosinte-vs-B73 wideseq positions; **no BAMs needed**).

**Reference panel:** 217 teosinte accessions, taxa from `bzea_50K_cohort_ref_metadata.csv`
(`maizegdb_prefix` Zx/Zv/Zd/Zl/Zh + Zn/teo_mix).

- 50K panel: `bzea_50K_ref_panel.vcf.gz` — **dense (0.0% missing**, verified: only 0/0,1/1,0/1,1/0).
- 27M panel: `wideseq_ref/wideseq_chr*.vcf.gz` — same 217 accessions, tabix-indexed, **~60% missing**.

**Method:** per block, missing-aware Jaccard of non-B73 presence (query vs each accession over the
co-observed intersection); nearest accession = donor for that block; aggregate genome-wide weighted
by block Mb → per-NIL **donor-individual** table + **species/taxa** rollup. Validate the dominant
taxon against the pedigree donor taxon.

**Panel taxa coverage (the key constraint):**

| taxon | panel accessions |
|---|---|
| Zv parviglumis | 88 |
| Zx mexicana | 79 |
| Zn nicaraguensis | 14 |
| teo_mix | 12 |
| Zl luxurians | 9 |
| Zd diploperennis | 5 |
| Zh huehuetenangensis | ~1 |

Scripts: `inst/scripts/jaccard_donor_search.py` (50K, all NILs), `inst/scripts/jaccard_27m_one.py`
(27M one genome), `inst/scripts/jaccard_27m_block_rank.py` (rank all accessions for one block).

## Expected results

1. Jaccard recovers the donor **taxon** where the panel represents it (**Zx, Zv**) and degrades
   where it doesn't (**Zd, Zl, Zh**). Match rate should track panel accession count.
2. The 27M panel (~1000× more SNPs/window) **would help only if the limit is discrimination**
   (SNP density). If the limit is **panel representation** (no close relative of the donor exists in
   the panel), more SNPs will not help.

## Obtained results

### 1. 50K full run — 1,385 NILs

Dominant taxon vs pedigree donor taxon. Two Jaccard normalizations tried (≈ a wash):

| taxon | n | symmetric Jaccard | Fig-10 containment |
|---|---|---|---|
| Zv | ~360 | **92%** | 84% |
| Zx | ~565 | 76% | **80%** |
| Zl | ~115 | 57% | 57% |
| Zd | ~240 | 45% | 42% |
| Zh | 61 | 2% | **11%** |
| **overall** | | **70%** (966/1379) | 69% (919/1324) |

Match rate **tracks panel representation almost linearly** (Zv 88 accs→92% … Zh 1 acc→2%),
exactly as expected. The normalization choice is secondary (symmetric penalizes alt-rich
accessions; containment rewards them — opposite biases, net wash).

### 2. 27M test on the highest-coverage *failed* Zd — `PN14_SID1284` (λ=0.87)

50K wrongly called it **Zl** (0.58). Whole-genome 27M run (6 blocks, **1,030,400 covered SNPs**,
~3 min):

```
chr1 22K SNPs -> RIL2 [Zl] J=0.569 ; chr2 40K -> RIL2 [Zl] 0.304 ; chr3 277K -> RIL2 [Zl] 0.201
chr3 646K -> RIL2 [Zl] 0.194 ; chr5 20K -> RIL2 [Zl] 0.265 ; chr7 26K -> Ame21884 [Zd] 0.278
=> Zl 0.97 / Zd 0.03 ; MATCH False
```

**27M did NOT rescue it** — it reproduced the 50K answer *more* confidently. More SNPs did not help.

### 3. Diagnostic — rank all 217 accessions on the chr3:67–150 Mb block (646K SNPs, 42,661 NIL non-B73)

```
RIL2     Zl 0.194   RIL003 Zl 0.194   RIMH001 Zh 0.192
5D3      Zd 0.190 (rank 4)  Ame2317 Zd 0.189   Ame21884 Zd 0.184   5D7 Zd 0.182
... Zh cluster 0.15–0.18 ...   best Zd = 0.190 (rank 4) vs top Zl 0.194 (Δ 0.004)
```

**Verdict: panel-gap, NOT mislabel.** (a) The top 7 — Zl, Zh *and* Zd — sit within **0.012**; best
Zd trails the "winner" by **0.004** (noise). (b) The absolute best J is **~0.19** (vs ~0.98 for a
real Dilkes match) — *no* accession is a true relative of this donor. The diploperennis donor simply
isn't in the panel; the whole perennial Zd/Zl/Zh clade collapses to a flat, unresolved plateau.

**Conclusion:** the limit for the distal taxa is **panel representation, not SNP density** — 646K
SNPs reproduce the 50K plateau exactly. This NIL is a correctly-labeled Zd that is *unresolvable*,
not *wrong*. The headline "Zd 45%" therefore **over-counts failures**: many are low-J plateau ties.
This motivates a **confidence filter** (next section) to separate confident calls from unresolved.

## Confidence-filtered re-score (50K)

*Hypothesis:* most Zd/Zl/Zh "misses" are low-J / low-margin plateau ties, not confident
mis-assignments. A filter — require the block's best J above a floor **and** beating the
runner-up taxon by a margin — should leave the confident calls (Zx/Zv) intact while reclassifying
the plateau ties as **"unresolved"** rather than counting them wrong.

**Result: the confidence filter does NOT work — hypothesis rejected.**

Block best-J by correctness (using the pedigree taxon as truth) is essentially the *same* for
correct and wrong calls:

| blocks | J p10 / p50 / p90 | margin p50 | n |
|---|---|---|---|
| **correct** (taxon==pedigree) | 0.29 / **0.52** / 0.79 | 0.054 | 7,665 |
| **wrong** | 0.28 / **0.50** / 0.78 | 0.031 | 5,405 |

J (and margin) barely discriminate. The whole grid confirms it:

| FLOOR | MARGIN | resolved NILs | match% among resolved |
|---|---|---|---|
| 0.3 | 0.02 | 1354 | 69% |
| 0.3 | 0.10 | 912 | 72% |
| 0.4 | 0.05 | 1226 | 68% |
| 0.5 | 0.10 | 862 | 73% |

Even the strictest filter (0.5/0.10) only reaches 73% while discarding ~40% of NILs. At the default
**FLOOR=0.4, MARGIN=0.05**: resolved 1226/1385 (89%), **match among resolved 68% — no better than
the 70% unconditional.** By taxon among resolved: **Zv 90%, Zx 76%, Zl 48%, Zd 30%, Zh 0%.**

**Interpretation — corrected against the Chen 2022 *Zea* phylogeny.** Confidence filtering at the
*species* level fails, but the reason is **not** "confident-but-wrong sister matches" as I first
wrote (and I had the sisters wrong). Re-scoring at the **section** level changes the whole picture.

## Correction & section-level re-analysis (Chen 2022)

Chen et al. 2022, *Genome sequencing reveals evidence of adaptive variation in the genus Zea*
(Nat. Genet., **doi:10.1038/s41588-022-01184-y**; bioRxiv 2022.06.03.494450). Relevant topology:

- **Section Luxuriantes** (perennial/distal clade): *Z. diploperennis* (Zd), *Z. perennis* (Zp),
  *Z. luxurians* (Zl), *Z. nicaraguensis* (Zn).
- **Section Zea / Z. mays**: *Z. parviglumis* (Zv), *Z. mexicana* (Zx), **Z. huehuetenangensis (Zh)**,
  *Z. mays* (Zm).
- **Z. mays, luxurians, and diploperennis diverged nearly *contemporaneously* (~120 kya)** — a
  near-polytomy; **Z. perennis split from Z. diploperennis ~48 kya** (Zd's true sister is **Zp**,
  not Zl); **Z. nicaraguensis is ~a subspecies of Z. luxurians** (Zl↔Zn).

So my earlier "sisters Zx↔Zv, Zd↔Zl" was wrong: Zd↔Zp and Zl↔Zn; and **Zh is in *Z. mays*** (section
Zea), not the perennial clade. Two consequences I had missed:

1. **Panel metadata bug:** the 8 *Z. perennis* accessions are mislabeled `maizegdb_prefix = Zh`
   (`taxa_label = Z. huehuetenangensis`); only the `reference_group` column is correct. So my
   `Zh` category was mostly *perennis*, and there is exactly **1 true Zh accession** (RIMH001). The
   per-taxon scoring used the buggy prefix as truth → mis-grounded for Zh/Zd.
2. **The species-level "failures" are within-section near-polytomy + sparse/mislabeled panel — not
   the method picking a wrong lineage.**

**Section-level accuracy** (matched accession's *correct* species → section, vs the NIL's donor
section; using `reference_group`):

| donor taxon | species-level | section-level |
|---|---|---|
| Zd | 45% | **91%** |
| Zh | 2% | **93%** |
| Zl | 57% | **96%** |
| Zv | 92% | **99%** |
| Zx | 76% | **100%** |
| **overall** | 70% | **97%** |

For **Zd** NILs the matched-accession species (Mb-weighted): *diploperennis* **0.40** (true taxon,
plurality) / *perennis* 0.22 (true sister) / *luxurians* 0.19 → **78% in the correct section, true
species on top.** The Zh "2%" was an artifact: huehuetenangensis is *Z. mays*, and Zh NILs correctly
match Zv/Zx (section Zea) **93%** of the time — "wrong" only because there's 1 true Zh accession.

## Conclusions (revised)

- **The method recovers the donor's *section* with 97% accuracy** (Zd 91, Zh 93, Zl 96, Zv 99,
  Zx 100) — it reliably places every NIL in the right half of the *Zea* tree.
- **Species-level resolution is real only within section Zea (Zv 92%, Zx 76%).** Within section
  **Luxuriantes (Zd/Zl/Zn/Zp) the species are a near-polytomy** (Chen 2022: ~simultaneous
  divergence) and **cannot be resolved** at the species level from short introgressions — confirmed
  by 27M (reproduces the plateau) and by confidence filtering (J doesn't separate correct/wrong
  *because the species genuinely overlap*). This is **biology, not a method or confidence defect.**
- **Fix the panel metadata** (re-label the 8 *Z. perennis* accessions from `Zh`→`Zp`; use
  `reference_group`, not `maizegdb_prefix`) before any per-taxon reporting.
- **Practical use:** report the donor's **section** for all NILs (reliable); report **species** for
  **Zv/Zx**; for Luxuriantes report "section Luxuriantes (diploperennis/luxurians/perennis,
  species-unresolved)". The per-NIL composition in `results/donor_id/nil_donor_summary.csv` is usable
  as a section call genome-wide and a species call for Zv/Zx.

## Verdict (F. Rodríguez): the SNP50K skim data is NOT diagnostic for donor identification

Bottom line from this side-quest: **this data cannot point to a donor in the reference panel, and
within the perennial clade it cannot even point to the species. It is not diagnostic for donor ID.**

- **No donor individual is identifiable.** The best Jaccard for any panel accession is ~0.19–0.5
  (a real donor match, à la Brian, is ~0.98). The "top" accession is just the marginal peak of a
  noisy distribution — not a donor.
- **No species within section Luxuriantes** (Zd/Zl/Zn/Zp): near-polytomy → 45–57% species, no
  better with 27M SNPs or any confidence threshold.
- The **only** signals that survive are coarse and do **not** constitute donor ID: the **section**
  (Z. mays vs perennial, 97%) and the **Zv↔Zx subspecies** split (76–92%, the one place subspecies
  resolves). Useful as broad ancestry context; useless for naming a donor accession or a Luxuriantes
  species.

## Why is Brian's WideSeq diagnostic (~98%) and ours isn't? — plausible explanations

Ranked, most decisive first:

1. **Donor IN the panel vs NOT in the panel (the crux).** Brian's introgressions came from *known
   maize inbred lines* (b94→LH52, →Mo17; Z004E0189→CML247) that **are members of HapMap3**. The
   true donor is literally in the reference panel → near-identity (98%). Our teosinte NIL donors are
   specific wild accessions almost certainly **absent** from the 217-accession diversity panel →
   the best possible hit is a *relative*, which for teosinte is distant (J~0.19). You can't match to
   something that isn't there.
2. **Inbred (homozygous, ~clonal) vs outbred teosinte (heterozygous, diverse).** A maize inbred is a
   fixed near-homozygous genotype that matches its panel copy almost perfectly. A teosinte
   "accession" is an **outcrossing, heterozygous population**, not one genotype — so even the *same*
   accession sampled twice differs, and within-accession/within-taxon diversity dilutes the Jaccard.
   A clean high-J match to a single teosinte donor may not exist *even in principle*.
3. **Panel size & representation.** HapMap3 ≈ 1,210 maize lines densely tile the elite/founder maize
   gene pool → almost any maize donor has a close relative. Our 217 teosinte, with **1–14 accessions
   per distal taxon**, sparsely sample teosinte diversity (esp. Luxuriantes).
4. **Phylogenetic structure / divergence.** Maize inbreds = one shallow gene pool, resolvable line-
   to-line. Teosinte section Luxuriantes = **near-polytomy (~120 kya)** → little diagnostic structure
   to assign even a species, let alone an accession.
5. **Marker ascertainment.** wideseq SNPs mark **teosinte-vs-B73** divergence — many are shared
   across *all* teosinte (fixed vs maize) → uninformative for *which* teosinte. HapMap3 SNPs segregate
   *among maize lines* → directly diagnostic for Brian's question.
6. **Query + panel sparsity (secondary).** NIL skim ~0.4–0.9× *and* the 27M teosinte panel ~60%
   missing → thin, noisy co-observed overlap lowers J. (Secondary: section-level still worked, so this
   is not the main driver.)

**The first two are the heart of it:** Brian matched introgressions to their *actual inbred donors,
present in a huge panel*; we are matching *outbred teosinte* introgressions to a *small panel that
does not contain the donors*. Same method, fundamentally different data regime → his is diagnostic,
ours is not.

Outputs: `results/donor_id/nil_donor_summary.csv`, `results/donor_id/nil_block_matches.csv`.