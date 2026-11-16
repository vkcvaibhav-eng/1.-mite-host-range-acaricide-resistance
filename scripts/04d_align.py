"""Stage 4d: align the chosen sequences and build the concatenated matrices.

Input:   data/clean/accessions.csv
Outputs (one set per superfamily tree, as planned for stage 5):
         data/clean/alignments/concat_<group>.fasta
         data/clean/alignments/partitions_<group>.nex    COI codon positions, 18S, 28S
         data/clean/alignments/<marker>_<group>.fasta    single-marker alignments
         results/tables/stage4_alignment_summary.csv
         logs/versions_stage4.txt

COI: translated with the invertebrate mitochondrial code; near-complete proteins (>= 400 aa) are aligned with
MAFFT L-INS-i, shorter ones are added with --addfragments, and the nucleotides are threaded back on the
protein alignment. Where a taxon has two COI segments they are merged (segment 1 has priority).
18S and 28S: MAFFT --auto with --adjustdirection, then trimAl -gappyout.
"""
import subprocess
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/stage4"
OUT = ROOT / "data/clean/alignments"
WORK.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
TRIMAL = "trimal"
GROUPS = {"tetranychoidea": ["Tetranychidae", "Tenuipalpidae", "Cheyletidae", "Trombiculidae"],
          "eriophyoidea": ["Eriophyidae", "Phytoptidae"]}


def run(cmd, stdin=None):
    return subprocess.run(cmd, capture_output=True, text=True, check=True, input=stdin).stdout


def read_fasta(text):
    seqs, name = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0].replace("_R_", "", 1) if line.startswith(">_R_") else line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, seqs):
    Path(path).write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))


def align_coi(rows, tag):
    nuc = {f"{r.taxon.replace(' ', '_')}|{r.segment}": r.sequence for r in rows.itertuples()}
    prot = {k: str(Seq(v).translate(table=5)) for k, v in nuc.items()}
    long_ = {k: v for k, v in prot.items() if len(v) >= 400}
    short = {k: v for k, v in prot.items() if len(v) < 400}
    write_fasta(WORK / f"coi_{tag}_long.faa", long_)
    aln = run(["mafft", "--quiet", "--thread", "4", "--localpair", "--maxiterate", "1000",
               str(WORK / f"coi_{tag}_long.faa")])
    if short:
        (WORK / f"coi_{tag}_long.aln").write_text(aln)
        write_fasta(WORK / f"coi_{tag}_short.faa", short)
        aln = run(["mafft", "--quiet", "--thread", "4", "--addfragments", str(WORK / f"coi_{tag}_short.faa"),
                   str(WORK / f"coi_{tag}_long.aln")])
    paln = read_fasta(aln)
    # thread nucleotides on the protein alignment
    codon = {}
    for k, p in paln.items():
        s, i, out = nuc[k], 0, []
        for aa in p:
            if aa == "-":
                out.append("---")
            else:
                out.append(s[i:i + 3])
                i += 3
        assert i == len(s), (k, i, len(s))
        codon[k] = "".join(out)
    # merge segments per taxon
    merged = {}
    for k in sorted(codon):
        taxon, seg = k.split("|")
        if taxon not in merged:
            merged[taxon] = list(codon[k])
        else:
            base = merged[taxon]
            for j in range(0, len(base), 3):
                if base[j:j + 3] == ["-", "-", "-"]:
                    base[j:j + 3] = codon[k][j:j + 3]
    merged = {k: "".join(v) for k, v in merged.items()}
    # drop codon columns that are empty or held by a single taxon
    n = len(next(iter(merged.values())))
    keep = [j for j in range(0, n, 3) if sum(s[j:j + 3] != "---" for s in merged.values()) >= 2]
    return {k: "".join(s[j:j + 3] for j in keep) for k, s in merged.items()}


def align_rrna(rows, marker, tag):
    seqs = {r.taxon.replace(" ", "_"): r.sequence for r in rows.itertuples()}
    raw = WORK / f"{marker}_{tag}.fasta"
    write_fasta(raw, seqs)
    aln = run(["mafft", "--quiet", "--thread", "4", "--adjustdirection", "--auto", str(raw)])
    a = read_fasta(aln)
    alnf = WORK / f"{marker}_{tag}.aln"
    write_fasta(alnf, {k: v.upper() for k, v in a.items()})
    trimmed = WORK / f"{marker}_{tag}.trim"
    run([TRIMAL, "-in", str(alnf), "-out", str(trimmed), "-gappyout", "-fasta"])
    t = read_fasta(trimmed.read_text())
    return {k: v.upper() for k, v in t.items()}, len(next(iter(a.values())))


def main():
    acc = pd.read_csv(ROOT / "data/clean/accessions.csv", keep_default_na=False)
    summary = []
    for tag, fams in GROUPS.items():
        sub = acc[acc["family"].isin(fams)]
        blocks, raw_len = {}, {}
        blocks["COI"] = align_coi(sub[sub["marker"] == "COI"], tag)
        for m in ("18S", "28S"):
            rows = sub[sub["marker"] == m]
            if len(rows) >= 3:
                blocks[m], raw_len[m] = align_rrna(rows, m, tag)
        taxa = sorted(set().union(*[set(b) for b in blocks.values()]))
        parts, concat, pos = [], {t: [] for t in taxa}, 1
        for m, b in blocks.items():
            L = len(next(iter(b.values())))
            write_fasta(OUT / f"{m}_{tag}.fasta", b)
            for t in taxa:
                concat[t].append(b.get(t, "-" * L))
            if m == "COI":
                for c in range(3):
                    parts.append(f"  charset COI_pos{c + 1} = {pos + c}-{pos + L - 1}\\3;")
            else:
                parts.append(f"  charset r{m} = {pos}-{pos + L - 1};")
            present = sum(set(s) - set("-N?") != set() for s in b.values())
            summary.append(dict(tree=tag, marker=m, taxa_with_marker=present, taxa_in_matrix=len(taxa),
                                aligned_length=L, length_before_trimming=raw_len.get(m, L),
                                missing_cells_pct=round(100 * (1 - sum(len(s) - s.count("-") - s.count("N")
                                                                         for s in b.values()) / (L * len(taxa))), 1)))
            pos += L
        write_fasta(OUT / f"concat_{tag}.fasta", {t: "".join(v) for t, v in concat.items()})
        (OUT / f"partitions_{tag}.nex").write_text("#nexus\nbegin sets;\n" + "\n".join(parts) + "\nend;\n")
        summary.append(dict(tree=tag, marker="all", taxa_with_marker=len(taxa), taxa_in_matrix=len(taxa),
                            aligned_length=pos - 1))
    s = pd.DataFrame(summary)
    s.to_csv(ROOT / "results/tables/stage4_alignment_summary.csv", index=False)
    print(s.to_string(index=False))
    mafft_v = subprocess.run(["mafft", "--version"], capture_output=True, text=True).stderr.strip()
    trimal_v = run([TRIMAL, "--version"]).strip()
    import Bio
    (ROOT / "logs/versions_stage4.txt").write_text(
        f"MAFFT {mafft_v}\n{trimal_v}\nBiopython {Bio.__version__}\npandas {pd.__version__}\n")


if __name__ == "__main__":
    main()
