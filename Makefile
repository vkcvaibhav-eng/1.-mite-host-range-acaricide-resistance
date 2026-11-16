# One command rebuilds every table and figure: make all
# Targets are filled in as each stage is completed.

.PHONY: all feasibility clean-results

all: feasibility resistance hosts

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
