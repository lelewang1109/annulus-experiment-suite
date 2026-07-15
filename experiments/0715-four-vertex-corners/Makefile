PYTHON ?= python3

.PHONY: setup run quick check

setup:
	$(PYTHON) -m pip install -r requirements.txt

run:
	$(PYTHON) scripts/run_all.py

quick:
	$(PYTHON) scripts/run_all.py --skip-slow

check:
	$(PYTHON) scripts/validate_repository.py
