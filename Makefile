PY ?= .venv/bin/python
PIP ?= .venv/bin/pip

.PHONY: setup data clinical expression cohorts analysis simulate readme site validate api test lint docker clean

setup:
	python3 -m venv .venv
	$(PIP) install --upgrade pip
	$(PIP) install -e ".[dev]"

clinical:
	$(PY) scripts/fetch_clinical.py

expression:
	$(PY) scripts/fetch_expression.py

data: clinical expression

cohorts:
	PYTHONPATH=src $(PY) -c "from replica.config import load_config; from replica.data import build_cohorts; build_cohorts(load_config())"

analysis: cohorts
	PYTHONPATH=src $(PY) scripts/run_analysis.py

simulate:
	PYTHONPATH=src $(PY) scripts/run_simulations.py

readme:
	PYTHONPATH=src $(PY) scripts/build_readme.py

site:
	PYTHONPATH=src $(PY) scripts/build_site.py

validate:
	PYTHONPATH=src $(PY) scripts/validate_artifacts.py

pipeline: data analysis site validate
	@echo "Replica pipeline complete: artifacts/results.json and site/index.html"

api:
	PYTHONPATH=. $(PY) -m uvicorn app.api:app --host 0.0.0.0 --port 8000

test:
	PYTHONPATH=src $(PY) -m pytest -q

lint:
	$(PY) -m ruff check src scripts tests app

docker:
	docker build -t replica:1.0.0 .

clean:
	rm -rf .pytest_cache .ruff_cache
