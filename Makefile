# One command rebuilds every table and figure: make all
# Targets are filled in as each stage is completed.

.PHONY: all feasibility clean-results

all: feasibility resistance

feasibility:
	python scripts/01_feasibility.py

clean-results:
	rm -rf results/tables/* results/figures/*

resistance:
	python scripts/02_resistance.py
	python scripts/02b_introduction_years.py

resistance_old:
	python scripts/02_resistance.py
