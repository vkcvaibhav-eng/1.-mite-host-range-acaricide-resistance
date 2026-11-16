#!/usr/bin/env python3
"""Stage 5e: write the BEAST 2 input file for one concatenated matrix.

Model (protocol: relaxed clock, relative dating)
- partitions as merged by ModelFinder in the pilot run; GTR + gamma (4 categories) with estimated base
  frequencies in each; relative partition rates estimated, weighted mean fixed at 1
- one uncorrelated lognormal relaxed clock over all partitions, mean rate fixed at 1, so branch lengths are
  in expected substitutions per site and node ages are relative; trees are rescaled to root age 1 afterwards
- birth-death tree prior
- clades listed in the constraints file are kept monophyletic (backbone of Matsuda et al. 2018 and outgroup)

Usage: python scripts/05e_beast_xml.py tetranychoidea [chain_length] [constraints.csv]
Output: analysis/beast/<name>.xml (one file; chains differ only by the seed given on the command line)
"""
import csv
import re
import sys
from pathlib import Path
from xml.sax.saxutils import quoteattr

ROOT = Path(__file__).resolve().parents[1]
name = sys.argv[1] if len(sys.argv) > 1 else "tetranychoidea"
chain = int(sys.argv[2]) if len(sys.argv) > 2 else 100_000_000
cons_file = Path(sys.argv[3]) if len(sys.argv) > 3 else ROOT / "data" / "clean" / f"beast_constraints_{name}.csv"
aln_file = ROOT / "data" / "clean" / "alignments" / f"concat_{name}.fasta"
scheme = ROOT / "results" / "trees" / "pilot" / f"pilot_{name}.best_scheme.nex"
out = ROOT / "analysis" / "beast"
out.mkdir(parents=True, exist_ok=True)

# ---- alignment
seqs, cur = {}, None
for line in open(aln_file):
    line = line.strip()
    if line.startswith(">"):
        cur = line[1:]
        seqs[cur] = []
    elif line:
        seqs[cur].append(line)
seqs = {k: "".join(v).upper() for k, v in seqs.items()}
taxa = list(seqs)
nsite = len(next(iter(seqs.values())))

# ---- partitions: charset name -> BEAST filter string and number of sites
parts = []
for m in re.finditer(r"charset\s+(\S+)\s*=\s*([^;]+);", scheme.read_text()):
    pid, ranges = m.group(1), m.group(2).split()
    n = 0
    for r in ranges:
        a, b, step = re.match(r"(\d+)-(\d+)(?:\\(\d+))?", r).groups()
        n += len(range(int(a), int(b) + 1, int(step or 1)))
    parts.append((pid, ",".join(ranges), n))
assert sum(n for _, _, n in parts) == nsite, "partitions do not cover the alignment"

# ---- monophyly constraints
clades = []
if cons_file.exists():
    for r in csv.DictReader(open(cons_file)):
        members = r["taxa"].split()
        assert set(members) <= set(taxa), f"unknown taxon in clade {r['clade']}"
        clades.append((r["clade"], members))

X = []
w = X.append
NS = ":".join(["beast.pkgmgmt", "beast.base.core", "beast.base.inference", "beast.base.evolution.alignment",
               "beast.base.evolution.tree", "beast.base.evolution.tree.coalescent", "beast.base.inference.util",
               "beast.base.evolution.operator", "beast.base.inference.operator", "beast.base.evolution.sitemodel",
               "beast.base.evolution.substitutionmodel", "beast.base.evolution.likelihood",
               "beast.base.evolution.branchratemodel", "beast.base.evolution.speciation",
               "beast.base.inference.distribution", "beast.base.inference.parameter"])
RP = "beast.base.inference.parameter.RealParameter"
w('<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
w(f'<!-- written by scripts/05e_beast_xml.py; {len(taxa)} taxa, {nsite} sites, {len(parts)} partitions, '
  f'{len(clades)} monophyly constraints -->')
w(f'<beast version="2.7" namespace="{NS}">')
w('<data id="aln" spec="Alignment" dataType="nucleotide">')
for t in taxa:
    w(f'  <sequence id="seq_{t}" spec="Sequence" taxon="{t}" totalcount="4" value="{seqs[t]}"/>')
w('</data>')
for pid, filt, _ in parts:
    w(f'<data id="{pid}" spec="FilteredAlignment" data="@aln" filter={quoteattr(filt)}/>')
w('<taxonset id="taxa" spec="TaxonSet">')
for t in taxa:
    w(f'  <taxon id="{t}" spec="Taxon"/>')
w('</taxonset>')

w(f'<run id="mcmc" spec="MCMC" chainLength="{chain}">')
store = max(chain // 2000, 1000)
w(f'<state id="state" spec="State" storeEvery="{store}">')
w('  <tree id="tree" spec="beast.base.evolution.tree.Tree" name="stateNode" taxonset="@taxa"/>')
w(f'  <parameter id="BDBirthRate" spec="{RP}" name="stateNode" lower="0.0" upper="1000.0">10.0</parameter>')
w(f'  <parameter id="BDDeathRate" spec="{RP}" name="stateNode" lower="0.0" upper="1.0">0.5</parameter>')
w(f'  <parameter id="ucldStdev" spec="{RP}" name="stateNode" lower="0.0" upper="5.0">0.3</parameter>')
w('  <stateNode id="rateCategories" spec="beast.base.inference.parameter.IntegerParameter" '
  f'dimension="{2 * len(taxa) - 2}">1</stateNode>')
for pid, _, _ in parts:
    for r in ["AC", "AG", "AT", "CG", "GT"]:
        w(f'  <parameter id="rate{r}.{pid}" spec="{RP}" name="stateNode" lower="0.0">1.0</parameter>')
    w(f'  <parameter id="gammaShape.{pid}" spec="{RP}" name="stateNode" lower="0.05" upper="50.0">0.5</parameter>')
    w(f'  <parameter id="freq.{pid}" spec="{RP}" name="stateNode" dimension="4" lower="0.0" upper="1.0">0.25</parameter>')
    w(f'  <parameter id="mutationRate.{pid}" spec="{RP}" name="stateNode" lower="0.0">1.0</parameter>')
w('</state>')

w('<init id="startTree" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@tree" '
  'taxa="@aln" rootHeight="0.5">')
w(f'  <populationModel spec="ConstantPopulation"><parameter name="popSize" spec="{RP}">1.0</parameter></populationModel>')
for cid, _ in clades:
    w(f'  <constraint idref="mono.{cid}"/>')
w('</init>')

w('<distribution id="posterior" spec="CompoundDistribution">')
w('<distribution id="prior" spec="CompoundDistribution">')
w('  <distribution id="BirthDeath" spec="beast.base.evolution.speciation.BirthDeathGernhard08Model" '
  'birthDiffRate="@BDBirthRate" relativeDeathRate="@BDDeathRate" tree="@tree"/>')
U = "beast.base.inference.distribution.Uniform"
w(f'  <distribution id="prior.BDBirthRate" spec="Prior" x="@BDBirthRate"><distr spec="{U}" lower="0.0" upper="1000.0"/></distribution>')
w(f'  <distribution id="prior.BDDeathRate" spec="Prior" x="@BDDeathRate"><distr spec="{U}" lower="0.0" upper="1.0"/></distribution>')
w('  <distribution id="prior.ucldStdev" spec="Prior" x="@ucldStdev"><distr spec="Gamma" alpha="0.5396" beta="0.3819"/></distribution>')
for pid, _, _ in parts:
    for r in ["AC", "AG", "AT", "CG", "GT"]:
        beta = "20.0" if r == "AG" else "10.0"
        w(f'  <distribution id="prior.rate{r}.{pid}" spec="Prior" x="@rate{r}.{pid}"><distr spec="Gamma" alpha="0.05" beta="{beta}"/></distribution>')
    w(f'  <distribution id="prior.gammaShape.{pid}" spec="Prior" x="@gammaShape.{pid}"><distr spec="Exponential" mean="1.0"/></distribution>')
for cid, members in clades:
    w(f'  <distribution id="mono.{cid}" spec="beast.base.evolution.tree.MRCAPrior" monophyletic="true" tree="@tree">')
    w(f'    <taxonset id="clade.{cid}" spec="TaxonSet">' + "".join(f'<taxon idref="{t}"/>' for t in members) + '</taxonset>')
    w('  </distribution>')
w('</distribution>')

w('<distribution id="likelihood" spec="CompoundDistribution" useThreads="true">')
for i, (pid, _, _) in enumerate(parts):
    w(f'  <distribution id="treeLikelihood.{pid}" spec="TreeLikelihood" data="@{pid}" tree="@tree">')
    w(f'    <siteModel id="siteModel.{pid}" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape.{pid}" mutationRate="@mutationRate.{pid}">')
    w(f'      <substModel id="gtr.{pid}" spec="GTR" rateAC="@rateAC.{pid}" rateAG="@rateAG.{pid}" rateAT="@rateAT.{pid}" rateCG="@rateCG.{pid}" rateGT="@rateGT.{pid}">')
    w(f'        <parameter id="rateCT.{pid}" spec="{RP}" estimate="false" name="rateCT" lower="0.0">1.0</parameter>')
    w(f'        <frequencies id="freqs.{pid}" spec="Frequencies" frequencies="@freq.{pid}"/>')
    w('      </substModel>')
    w('    </siteModel>')
    if i == 0:
        w('    <branchRateModel id="relaxedClock" spec="beast.base.evolution.branchratemodel.UCRelaxedClockModel" rateCategories="@rateCategories" tree="@tree">')
        w('      <distr id="clockLogNormal" spec="beast.base.inference.distribution.LogNormalDistributionModel" S="@ucldStdev" meanInRealSpace="true">')
        w(f'        <parameter id="clockLogNormal.M" spec="{RP}" estimate="false" name="M">1.0</parameter>')
        w('      </distr>')
        w(f'      <parameter id="ucldMean" spec="{RP}" estimate="false" name="clock.rate">1.0</parameter>')
        w('    </branchRateModel>')
    else:
        w('    <branchRateModel idref="relaxedClock"/>')
    w('  </distribution>')
w('</distribution>')
w('</distribution>')

SC = "ScaleOperator"
w(f'<operator id="op.BDBirthRate" spec="{SC}" parameter="@BDBirthRate" scaleFactor="0.75" weight="3.0"/>')
w(f'<operator id="op.BDDeathRate" spec="{SC}" parameter="@BDDeathRate" scaleFactor="0.75" weight="3.0"/>')
w(f'<operator id="op.treeScaler" spec="{SC}" tree="@tree" scaleFactor="0.9" weight="3.0"/>')
w(f'<operator id="op.rootScaler" spec="{SC}" tree="@tree" rootOnly="true" scaleFactor="0.8" weight="3.0"/>')
w('<operator id="op.uniform" spec="beast.base.evolution.operator.Uniform" tree="@tree" weight="30.0"/>')
w('<operator id="op.subtreeSlide" spec="SubtreeSlide" tree="@tree" size="0.02" weight="15.0"/>')
w('<operator id="op.narrow" spec="Exchange" tree="@tree" weight="15.0"/>')
w('<operator id="op.wide" spec="Exchange" tree="@tree" isNarrow="false" weight="3.0"/>')
w('<operator id="op.wilsonBalding" spec="WilsonBalding" tree="@tree" weight="3.0"/>')
w(f'<operator id="op.ucldStdev" spec="{SC}" parameter="@ucldStdev" scaleFactor="0.5" weight="3.0"/>')
w('<operator id="op.catRandomWalk" spec="IntRandomWalkOperator" parameter="@rateCategories" windowSize="1" weight="10.0"/>')
w('<operator id="op.catSwap" spec="SwapOperator" intparameter="@rateCategories" weight="10.0"/>')
w('<operator id="op.catUniform" spec="UniformOperator" parameter="@rateCategories" weight="10.0"/>')
for pid, _, _ in parts:
    for r in ["AC", "AG", "AT", "CG", "GT"]:
        w(f'<operator id="op.rate{r}.{pid}" spec="{SC}" parameter="@rate{r}.{pid}" scaleFactor="0.5" weight="0.2"/>')
    w(f'<operator id="op.gammaShape.{pid}" spec="{SC}" parameter="@gammaShape.{pid}" scaleFactor="0.5" weight="0.2"/>')
    w(f'<operator id="op.freq.{pid}" spec="DeltaExchangeOperator" parameter="@freq.{pid}" delta="0.01" weight="0.2"/>')
w('<operator id="op.mutationRates" spec="DeltaExchangeOperator" delta="0.75" weight="2.0">')
for pid, _, _ in parts:
    w(f'  <parameter idref="mutationRate.{pid}"/>')
w('  <weightvector spec="beast.base.inference.parameter.IntegerParameter" dimension="'
  f'{len(parts)}" estimate="false" lower="0" upper="0">' + " ".join(str(n) for _, _, n in parts) + '</weightvector>')
w('</operator>')

every = max(chain // 10000, 1000)
w(f'<logger id="tracelog" spec="Logger" fileName="{name}.$(seed).log" logEvery="{every}" model="@posterior" sanitiseHeaders="true" sort="smart">')
for ref in ["posterior", "likelihood", "prior", "BirthDeath", "BDBirthRate", "BDDeathRate", "ucldStdev"]:
    w(f'  <log idref="{ref}"/>')
w('  <log id="treeStat" spec="beast.base.evolution.tree.TreeStatLogger" tree="@tree"/>')
w('  <log id="rateStat" spec="beast.base.evolution.RateStatistic" branchratemodel="@relaxedClock" tree="@tree"/>')
for pid, _, _ in parts:
    w(f'  <log idref="treeLikelihood.{pid}"/>')
    w(f'  <log idref="mutationRate.{pid}"/>')
    w(f'  <log idref="gammaShape.{pid}"/>')
    for r in ["AC", "AG", "AT", "CG", "GT"]:
        w(f'  <log idref="rate{r}.{pid}"/>')
    w(f'  <log idref="freq.{pid}"/>')
w('</logger>')
w(f'<logger id="treelog" spec="Logger" fileName="{name}.$(seed).trees" logEvery="{every}" mode="tree">')
w('  <log id="treeWithRates" spec="beast.base.evolution.TreeWithMetaDataLogger" branchratemodel="@relaxedClock" tree="@tree"/>')
w('</logger>')
w(f'<logger id="screenlog" spec="Logger" logEvery="{every * 10}">')
for ref in ["posterior", "likelihood", "prior"]:
    w(f'  <log idref="{ref}"/>')
w('</logger>')
w('</run>')
w('</beast>')
path = out / f"{name}.xml"
path.write_text("\n".join(X) + "\n")
print(f"{path.relative_to(ROOT)}: {len(taxa)} taxa, {nsite} sites, partitions "
      + ", ".join(f"{p} ({n})" for p, _, n in parts) + f"; {len(clades)} constraints; chain {chain:,}")
