# One command rebuilds every table and figure: make all
# Targets are filled in as each stage is completed.

.PHONY: all feasibility resistance hosts sequences phylogeny clean-results

all: feasibility resistance hosts sequences phylogeny

feasibility:
	python scripts/01_feasibility.py

clean-results:
	rm -rf results/tables/* results/figures/*

resistance:
	python scripts/02_resistance.py
	python scripts/02b_introduction_years.py
	python scripts/02d_verification_summary.py

# stage 3: needs the WFO backbone classification.csv (doi:10.5281/zenodo.20782718, _DwC_backbone_R.zip);
# set WFO_CLASSIFICATION to its path
hosts:
	python scripts/03a_candidate_species.py
	python scripts/03b_host_records.py
	python -W ignore scripts/03c_host_breadth.py

# stage 4: needs MAFFT and trimAl on the PATH; 04b downloads from GenBank (about 20 minutes)
sequences:
	python scripts/04a_taxa.py
	python scripts/04b_fetch_sequences.py 2026-11-16
	python -W ignore scripts/04c_select_sequences.py
	python scripts/04d_align.py

# stage 5: needs iqtree2 (2.3.6) on the PATH; 05a about 15 minutes, 05b about 90 minutes on 2 cores.
# 05e writes the BEAST 2 input; the BEAST runs themselves are not part of make (see logs/decisions.md, entry 71)
phylogeny:
	./scripts/05a_pilot_iqtree.sh
	./scripts/05b_constrained_iqtree.sh
	./scripts/05c_eriophyoidea_iqtree.sh
	python -W ignore scripts/05d_ml_summary.py
	python -W ignore scripts/05f_beast_constraints.py
	python scripts/05e_beast_xml.py tetranychoidea
