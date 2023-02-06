# Gratitude to IIASA Climate Assessment, which this is based upon
FILES_TO_FORMAT_PYTHON=scripts/*.py

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([0-9a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

.PHONY: help
help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

.PHONY: checks
checks:   ## run all the checks
	@echo "\n\n=== black ==="; black --check scripts/*.py || echo "--- black failed ---" >&2; \
		echo "\n\n=== isort ==="; isort --check-only --quiet scripts/*.py || echo "--- isort failed ---" >&2; \
		echo "\n\n=== flake8 ==="; flake8 scripts/*.py || echo "--- flake8 failed ---" >&2; \
		echo

.PHONY: format
format:  ## re-format files
	make isort
	make black

.PHONY: black
black:   ## use black to autoformat code
	black --target-version py37 $(FILES_TO_FORMAT_PYTHON)

isort:   ## format the code
	isort $(FILES_TO_FORMAT_PYTHON)
