# Tetranychinae backbone constraint

`constraint_matsuda2018.nwk` is the topology of Figure 1 (maximum-likelihood tree of 652 genes, nucleotides) of
Matsuda, Kozaki, Ishii and Gotoh (2018), PLoS ONE 13: e0203136, read from the published figure
(`data/raw/matsuda2018/`; the authors deposited reads only, DDBJ DRA007145, and no tree file).

Rules used when reading the figure
- Nodes shown with less than 100% bootstrap support are collapsed (inside Panonychus and inside the
  Oligonychus clade a3-1; both genera are left as unresolved clades).
- The three places where the nucleotide and amino-acid trees disagree (text of the Results) are collapsed:
  Eotetranychus uchidai and Schizotetranychus gilvus + S. shii; the position of S. recki in clade a1-4;
  the inside of Panonychus.
- Inside the Tetranychus urticae group the figure is too dense to read safely; T. urticae, T. turkestani and
  T. pueraricola are left unresolved next to the clade T. ezoensis + T. parakanzawai + T. kanzawai.
- Names follow this project's taxon names: T. okinawanus = T. gloveri, Eotetranychus asiaticus =
  E. sexmaculatus, Tetranychina harti = Petrobia harti; the red and green forms of T. urticae are one taxon.
- Species of the figure that are not in the alignment are pruned by `scripts/05b_constrained_iqtree.sh`
  before use; species of the alignment that are not in the figure are placed freely by the data.
