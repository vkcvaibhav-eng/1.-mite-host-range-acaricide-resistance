# One command rebuilds every table and figure: make all
# Targets are filled in as each stage is completed.

.PHONY: all feasibility clean-results

all: feasibility

feasibility:
	python scripts/01_feasibility.py

clean-results:
	rm -rf results/tables/* results/figures/*
