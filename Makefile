# Flint: reproduce the released statistics, reports, and figures from the
# retained participant tables. Requires Python 3.11 or newer.

PYTHON ?= python3.11
VENV := .venv
PY := $(VENV)/bin/python

.PHONY: all venv test verify cv report atlas figures clean

all: verify cv report atlas figures   ## Run every stage in dependency order

venv: $(PY)                            ## Create the virtual environment
$(PY):
	$(PYTHON) -m venv $(VENV)
	$(PY) -m pip install -r requirements.txt

test: $(PY)                            ## Unit and published-table regression tests
	$(PY) -m unittest discover -s tests

verify: $(PY)                          ## Recompute all 880 association rows and the bootstrap
	$(PY) scripts/verify.py

cv: $(PY)                              ## Leave-one-out and repeated five-fold cross-validation
	$(PY) scripts/cross_validate_readouts.py

report: $(PY)                          ## Render CROSS_VALIDATION.md from results/cross_validation/
	$(PY) scripts/write_cv_report.py

atlas: $(PY)                           ## Render PHENOTYPE_ATLAS.md from results/context/phenotypes/
	$(PY) scripts/build_phenotype_catalog.py

figures: $(PY)                         ## Regenerate every figure in figures/
	$(PY) scripts/plot_results.py

clean:
	rm -rf $(VENV) scripts/__pycache__ tests/__pycache__
