Repository Guidelines
=====================

Project Structure & Module Organization
--------------------------------------
- Core Python package lives in `antibody_abtigen/` (`cli.py`, `pipeline.py`, `mapping.py`, `sabdab.py`, `structure.py`, `yaml_converter.py`, `pymol_env.py`); `run.py` is a thin entrypoint for `uv run python run.py`.
- `scripts/` contains runnable helpers for SLURM (`run_full_pipeline.sh`) and quick sanity checks (`run_test.sh`).
- Data/cache/output defaults: `data/` (downloads, cache), `output/` (CSV/JSON + DP_* folders), `output_yamls/` (converted YAMLs). Keep large generated files out of commits.
- `docs/dev/` holds research notes and conversion guides; update or add brief summaries when behavior changes.

Build, Test, and Development Commands
-------------------------------------
- Install deps: `uv sync` (optionally `uv sync --extra pymol` if PyMOL is available).
- Quick dry-run (CI parity, no structure downloads): `uv run python run.py --limit 10 --dry-run --no-pymol`.
- Small functional run: `uv run antibody-abtigen build --output ./output --data-dir ./data --limit 5 --no-pymol`.
- YAML conversion only: `uv run antibody-abtigen to-yaml --input ./output --output ./output_yamls`.
- HPC/SLURM pipeline: `bash scripts/run_full_pipeline.sh` (expects `.venv` + optional `pymol-env` conda env).
- Linting: `uv run ruff check antibody_abtigen/ run.py --ignore E501`.

Coding Style & Naming Conventions
---------------------------------
- Python 3.10+; prefer type hints and module-level docstrings. Use snake_case for functions/variables, UpperCamelCase for classes, and lowercase-kebab for CLI flags.
- Keep lines reasonably short; `ruff` ignores E501 but wrap when sensible. Favor clear, explicit variable names for PDB/UniProt identifiers.
- CLI additions should expose both short and long flags following existing click patterns.

Testing Guidelines
------------------
- No unit tests yet; smoke-test with the dry-run command above before opening a PR. For structure downloads, cap `--limit` and default to `--no-pymol` in CI.
- If adding tests, place under `tests/` and use `pytest`; prefer fixtures that mock network calls to SAbDab/UniProt instead of real downloads.
- Validate outputs by checking `output/dataset_summary.csv` and `output/processing_log.json` are produced and non-empty.

Commit & Pull Request Guidelines
--------------------------------
- Commits: short imperative subject (`Add dry-run smoke test`), include rationale in the body when not obvious. Group related changes; avoid committing generated data/cache outputs.
- PRs: provide a concise description, linked issue (if any), commands run with results (e.g., snippets from dry-run/ruff), and note any data regeneration steps. Include screenshots only if UI/visual docs are affected.

Security & Configuration Tips
-----------------------------
- Do not hardcode API tokens; current pipeline uses public SAbDab/PDB endpoints. Keep private credentials in your environment and out of Git.
- When enabling PyMOL, prefer an isolated `pymol-env` conda environment to avoid system-wide installs.
- Large downloads live in `data/` and `output/`; ensure `.gitignore` stays effective before committing.***
