.PHONY: help setup format lint test check clean ui ui-test

.DEFAULT_GOAL := help

help:  ## Show available make commands
	@echo "=== Makefile ==="
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

setup:  ## Initial project setup — sync dependencies + install package in editable mode (run once)
	uv sync
	uv pip install -e . --reinstall

format:  ## Auto-format code with ruff
	uv run ruff format src/

lint:    ## Run ruff linting + auto-fix issues
	uv run ruff check --fix src/

test:    ## Run tests with pytest (quiet mode)
	uv run pytest

check:   ## Run full quality check (NON-destructive) — recommended before committing
	uv run ruff format --check src/
	uv run ruff check src/
	uv run pytest

ui:      ## Start the web UI for running the pipeline (LAN, no login)
	bash web/start.sh

ui-test: ## Self-test the web UI (add FULL=1 to also run a real pipeline into a temp dir)
	python3 web/selftest.py $(if $(FULL),--full,)

clean:   ## Remove build artifacts, caches, __pycache__, and any .egg-info folders (anywhere)
	rm -rf build/ dist/ *.egg-info/
	rm -rf .pytest_cache/ .ruff_cache/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true