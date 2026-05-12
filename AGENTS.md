# AGENTS.md

Agent guidance for the St. Jude Cloud bioinformatics workflows repository.

## What this repo is

A collection of WDL 1.1 pipelines and tool wrappers for genomic analysis. Primary artifacts are `.wdl` files; there is no traditional build step. Python and R exist only as scripts embedded in Docker images under `scripts/`.

## Layout

```
workflows/        # End-to-end pipelines (rnaseq/, dnaseq/, chipseq/, methylation/, etc.)
tools/            # Individual WDL task wrappers (one tool per file)
data_structures/  # WDL struct definitions
docker/           # Dockerfiles + package.json version metadata per image
scripts/          # Python and R scripts (not standalone; embedded in Docker images)
tests/            # pytest-workflow YAML test definitions
  tools/          # Per-tool test YAMLs
  workflows/      # Per-workflow test YAMLs
  input_json/     # Input JSON fixtures
  input/          # Binary fixtures (BAMs, FASTQs, etc.) — tracked via git-lfs
template/         # task-examples.wdl (canonical patterns) + common-parameter-meta.txt
developer_scripts/
  run_sprocket_or_miniwdl.sh   # Unified test runner
  update_container_tags.sh     # Rewrites container tags for branch testing (do not commit)
```

## Developer commands

### Setup

```bash
uv sync          # Python deps (canonical; use this, not pip)
git lfs pull     # Required after clone to get test fixture files
```

### WDL lint / format / check (Sprocket)

```bash
sprocket lint .                    # Lint all WDL
sprocket format tools/bwa.wdl      # Format a file
sprocket check .                   # Thorough check (validates inputs, etc.)
sprocket validate <wdl> <inputs>   # Validate inputs against a workflow
```

`sprocket.toml` disables the `ContainerUri` rule and sets `deny_notes = true` (notes in WDL cause check failure).

### Python (scripts/ only)

```bash
uv run ruff format --check --diff scripts/
uv run ruff format scripts/
uv run ruff check scripts/
```

### R (scripts/ only)

```r
styler::style_dir()   # format
lintr::lint_dir()     # lint
```

### Tests

```bash
# All tests
uv run pytest --kwdof --wt $(nproc)

# Single test by tag
pytest --tag bwa

# Single test by name
pytest -k bwa_aln

# Choose runner (defaults to sprocket)
RUNNER=miniwdl pytest --tag bwa
```

`pytest.ini` always injects `--git-aware --symlink`. Use `--kwdof` to keep outputs on failure.

Every test calls `./developer_scripts/run_sprocket_or_miniwdl.sh` internally — do not call runners directly in test commands.

### Run a WDL task directly (outside pytest)

```bash
# sprocket
sprocket run --output-dir output --target bwa_aln tools/bwa.wdl <input_file>

# miniwdl
miniwdl run --task bwa_aln --verbose --dir output/. -i tests/tools/input_json/bwa_aln.json tools/bwa.wdl
```

## CI pipeline

All jobs trigger on push. Key facts:
- `reference` and `slow` tags are **excluded** from the CI test matrix.
- CI deletes `slow`-tagged tests from YAML files in-place before running — do not replicate this destructively locally.
- CI builds Docker images and runs `update_container_tags.sh` to rewrite tags before pytest — do not commit the rewritten WDL files.
- Every test runs against both `sprocket` and `miniwdl` runners as a matrix.
- `sprocket` is built from source (`cargo install sprocket --locked`) in CI.
- `requirements-ci.txt` pins older versions than `pyproject.toml`; use `uv` locally.

## WDL conventions (non-obvious)

- All WDL files must be `version 1.1`.
- All tasks must include `set -euo pipefail` when using pipes or multiple commands.
- Multi-core tasks: accept `use_all_cores: Boolean = false` and `ncpu: Int = 2`; use `$(nproc)` when `use_all_cores` is true. `use_all_cores` must be the last Boolean in the input block; `ncpu` precedes memory/disk inputs.
- Resource inputs (`ncpu`, `modify_memory_gb`, `modify_disk_size_gb`) must be overridable.
- Single output → `outfile_name`; multiple outputs or extension matters → `prefix`.
- BAM/BAI companion pairs must be localized via `ln -s` to CWD and cleaned up with `rm` at task end.
- Scripts go in `scripts/`; never embed Python/R directly in WDL `command` blocks.
- Imports within the repo must use **relative paths**. External imports must pin a **tagged release** (not `main`/`master`) — enforced by `pull-check.yaml`.
- Deprecated tasks: add `deprecated: true` and `warning: "**[DEPRECATED]**..."` to `meta`. Never add `deprecated: false`.
- Update the relevant `CHANGELOG.md` when modifying any WDL under a subdirectory.

## Parameter meta conventions

Prefer names from `template/common-parameter-meta.txt`. Key ones:
- `bam` (not `input_bam`, `in_bam`)
- `bam_index` (not `bai`)
- `read_one_fastq_gz` / `read_two_fastq_gz` (not `read1`/`read2`)
- `paired_end` (not `paired`)

## Docker image versioning

Each image lives in `docker/<tool>/` with `Dockerfile` and `package.json`. `version` = underlying tool version; `revision` starts at `0`, increments for image-only changes, resets to `0` on tool version upgrades. Prefer BioContainers images over creating new custom images.

## Testing quirks

- Binary fixtures in `tests/input/` are in **git-lfs** — run `git lfs pull` before testing.
- Test files prefixed with `_` (e.g., `_test_methylation-preprocess.yaml`) are disabled.
- Sprocket outputs land in `output/runs/*/_latest/outputs.json`, then get copied to `output/outputs.json`. miniwdl outputs go directly to `output/`. Test YAMLs reference `output/outputs.json`.
- Fixture inputs are downsampled to small chromosomes (chrY/chrM, chr9/chr22) for speed.

## Existing instruction sources

- `.github/instructions/wdl.instructions.md` — applies to `**/*.wdl`; references CONTRIBUTING.md, best-practices.md, template/task-examples.wdl, and template/common-parameter-meta.txt.
- `template/task-examples.wdl` — canonical WDL task patterns; read before writing new tasks.
- `template/common-parameter-meta.txt` — required/banned parameter_meta strings.
- `CONTRIBUTING.md` — general coding style.
- `best-practices.md` — WDL-specific best practices.
