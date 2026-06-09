# Packaging plan — `structbio` (lightweight, git-installable)

> Status: draft / proposed · Date: 2026-06-09
> Turn the `structural-bioinfo` Claude plugin's scripts into a single
> `pip install git+…`-able Python package with a unified `structbio` CLI, **without**
> breaking the plugin or the existing plugin→repo sync workflow.

## Context (today)

- The toolkit ships as a Claude Code plugin: `plugins/structural-bioinfo/skills/<skill>/`,
  each skill = `SKILL.md` + `scripts/*.py` (+ `references/`, sometimes `assets/`, `examples/`).
- **~8,300 LOC** of Python across 10 skills / 12 scripts. **Every script already has both
  `def main()` and an `if __name__ == "__main__"` guard** — so exposing CLIs is near-zero refactor.
- Dependencies are *tiered*: about half the tools are stdlib-only; others need numpy, rdkit,
  mdtraj, the `mkdssp` binary, or torch+esm.
- Duplication is significant (the main "future development" pain): manual PDB column parsing in
  **6** skills, raw RCSB/HTTP in **6**, a self-contained HTML report shell in **5**, mmCIF parsing
  in 4, Kabsch/RMSD in 2, plus a conda re-exec/env-discovery hack in 3.
- Source-of-truth is `~/.claude/plugins/structural-bioinfo/`; `sync-claude-plugin` mirrors it
  one-way into `plugins/structural-bioinfo/` here and regenerates the top-level README.

## Goal & non-goals

**Goal:** `pip install git+…`, a unified `structbio <subcommand>` CLI, and a single source of
truth — while the plugin keeps working with zero install.

**Non-goals (out of scope for "lightweight"):** PyPI / bioconda publishing, release CI, a semver
release process, and the full DRY refactor (that is the optional Phase 2). No speculative abstraction.

## Key design decisions

### 1. The package lives *inside* the plugin dir
Put `structbio/` and `pyproject.toml` **under `plugins/structural-bioinfo/`**, next to `skills/`.

| Concern | Why inside-plugin wins |
|---|---|
| Sync workflow | Already mirrors `~/.claude/plugins/structural-bioinfo/` → repo; the package rides along |
| Zero-install plugin users | Package ships *with* the plugin; scripts bootstrap it via `sys.path` |
| pip users | `pip install "git+…#subdirectory=plugins/structural-bioinfo"` |
| Single source of truth | Logic lives once in `structbio/`; plugin `scripts/*.py` become ~6-line shims |

The plugin loader only reads `.claude-plugin/`, `skills/`, etc., so `structbio/` and
`pyproject.toml` are ignored — no conflict. (Repo-root would be the conventional Python layout, but
it would not ship with a marketplace plugin install and is not covered by the existing sync.)

### 2. A single `structbio` dispatcher (not per-tool commands)
One console entry point, `structbio`, dispatches to subcommands. Dispatch is by lazy
`importlib.import_module`, so a subcommand's heavy deps (rdkit/mdtraj/torch) load **only when that
subcommand runs** — nothing heavy loads for `structbio holo-search`.

```python
# structbio/cli.py
import sys, importlib

COMMANDS = {
    "holo-search":             "structbio.holo_search",
    "pdb-extract":             "structbio.pdb_extract",
    "epitope-host":            "structbio.epitope_host",
    "epitope-ss":              "structbio.epitope_ss",
    "af-ss":                   "structbio.af_ss",
    "pmhc-search":             "structbio.pmhc.search",
    "tcr-split":               "structbio.pmhc.split",
    "esm-featurize":           "structbio.esm_featurize",
    "pharmacophore":           "structbio.pharmacophore.analyze",
    "pharmacophore-report":    "structbio.pharmacophore.report",
    "pharmacophore-dashboard": "structbio.pharmacophore.dashboard",
}

def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        _print_help(COMMANDS); return 0
    cmd, rest = argv[0], argv[1:]
    if cmd not in COMMANDS:
        sys.stderr.write(f"structbio: unknown subcommand '{cmd}'\n"); _print_help(COMMANDS); return 2
    mod = importlib.import_module(COMMANDS[cmd])   # lazy → heavy deps load per subcommand
    sys.argv = [f"structbio {cmd}"] + rest         # each tool's own argparse takes over
    return mod.main()
```

`structbio <sub> --help` forwards to the tool's own argparse help (with `prog="structbio <sub>"`),
because the dispatcher rewrites `sys.argv` and calls the tool's existing `main()` unchanged.

| Subcommand | Module | From skill |
|---|---|---|
| `structbio holo-search` | `structbio.holo_search` | pdb-holostructure-search (`search_holostructures.py`) |
| `structbio pdb-extract` | `structbio.pdb_extract` | pdb-extractor (`download_and_extract.py`) |
| `structbio epitope-host` | `structbio.epitope_host` | search-epitope-host (`search_epitope_host.py`) |
| `structbio epitope-ss` | `structbio.epitope_ss` | epitope-secondary-structure (`epitope_ss.py`) |
| `structbio af-ss` | `structbio.af_ss` | af-secondary-structure (`af_ss.py`) |
| `structbio pmhc-search` | `structbio.pmhc.search` | pdb-pmhc-tcr-search (`search_pmhc_tcr.py`) |
| `structbio tcr-split` | `structbio.pmhc.split` | pdb-pmhc-tcr-search (`split_tcr_pmhc.py`) |
| `structbio esm-featurize` | `structbio.esm_featurize` | esm-featurize (`featurize.py`) |
| `structbio pharmacophore` | `structbio.pharmacophore.analyze` | pharmacophore-analyzer (`pharmacophore_analysis.py`) |
| `structbio pharmacophore-report` | `structbio.pharmacophore.report` | pharmacophore-report-generator (`extract_report_data.py`) |
| `structbio pharmacophore-dashboard` | `structbio.pharmacophore.dashboard` | pharmacophore-dashboard (`generate_report.py`) |

### 3. Lean core + optional extras
Heavy/conflicting deps are extras; imports stay lazy (which is what the conda re-exec hack was
working around — this deletes it).

## Target layout

```
plugins/structural-bioinfo/
├── pyproject.toml                      # NEW — metadata, single entry point, extras
├── structbio/                          # NEW — the engine (single source of truth)
│   ├── __init__.py   cli.py            # dispatcher
│   ├── holo_search.py  pdb_extract.py  epitope_host.py
│   ├── epitope_ss.py   af_ss.py        esm_featurize.py
│   ├── pmhc/           __init__.py  search.py  split.py
│   └── pharmacophore/  __init__.py  analyze.py  report.py  dashboard.py
│       └── assets/     styles.css   app.js      # moved here, bundled as package-data
├── skills/<skill>/
│   ├── SKILL.md                        # unchanged content (+ a pip-install note)
│   └── scripts/<x>.py                  # becomes a thin bootstrap shim
└── .claude-plugin/  README.md  LICENSE
```

## Standalone-vs-pip, resolved (the shim)

Each `skills/<s>/scripts/<x>.py` becomes a shim that prefers the installed package and falls back to
the in-tree copy, so `python scripts/x.py …` keeps working with no install:

```python
#!/usr/bin/env python3
import os, sys
try:
    from structbio.pharmacophore.analyze import main      # pip-installed
except ModuleNotFoundError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))  # in-tree
    from structbio.pharmacophore.analyze import main
if __name__ == "__main__":
    main()
```

## `pyproject.toml` sketch

```toml
[build-system]
requires = ["setuptools>=61"]
build-backend = "setuptools.build_meta"

[project]
name = "structbio"            # placeholder name — rename if preferred
version = "0.1.0"
requires-python = ">=3.9"
dependencies = ["numpy"]      # light; used by 3 tools → keeps install one-step

[project.optional-dependencies]
pharmacophore = ["rdkit"]              # pharmacophore                (conda recommended)
md            = ["mdtraj"]             # af-ss                        (conda recommended)
esm           = ["esm", "torch"]       # esm-featurize (heavy; verify the PyPI name)
all           = ["rdkit", "mdtraj", "esm", "torch"]
# external, not pip-installable: mkdssp → `conda install -c conda-forge dssp`  (epitope-ss)

[project.scripts]
structbio = "structbio.cli:main"

[tool.setuptools]
packages = ["structbio", "structbio.pmhc", "structbio.pharmacophore"]
[tool.setuptools.package-data]
"structbio.pharmacophore" = ["assets/*.css", "assets/*.js"]
```

**Tool → install matrix:** core (`pip install structbio`) covers holo-search, pdb-extract,
epitope-host, pmhc-search/tcr-split, plus the numpy-using report & dashboard; `[pharmacophore]`→rdkit,
`[md]`→mdtraj, `[esm]`→torch+esm. `epitope-ss` additionally needs the `mkdssp` binary (conda-forge).

## Phase 1 — package-ize (the lightweight deliverable)

Done in the plugin **source** (`~/.claude/plugins/structural-bioinfo/`), then synced. Verify each
step before the next.

1. **Skeleton** — `structbio/` + `__init__.py`s + `cli.py` + `pyproject.toml`.
2. **Move logic in** — relocate each `scripts/<x>.py` → its `structbio/…` module (mechanical; every
   file already has `main()`). Fix the dashboard's 3 sibling imports to package-relative, move
   `assets/` into `structbio/pharmacophore/assets/`, and load them via `importlib.resources` instead
   of the current `__file__`-relative path.
   → verify: `python -m structbio.cli <sub> --help` for all 11 subcommands.
3. **Shims** — replace each `scripts/<x>.py` with the bootstrap shim above.
   → verify: `python skills/<s>/scripts/<x>.py --help` still works with **no** pip install.
4. **`pip install -e` smoke test** — `structbio --help` lists subcommands; `structbio <sub> --help`
   forwards correctly; re-run the **ATR dashboard end-to-end** and diff the HTML vs. current output.
5. **Docs + hygiene** — README install section gets the `git+…#subdirectory=…` line + extras matrix
   and `structbio <sub>` examples (the `python scripts/…` examples stay valid via shims). Add
   `*.egg-info/`, `build/`, `dist/` to the sync rsync excludes so editable-install artifacts aren't
   committed.

After Phase 1: one copy of the logic; `pip install git+…` works with extras; plugin still works with
zero install; sync workflow unchanged.

## Phase 2 — optional, incremental DRY (the "future development" payoff)

Do these only as you touch each area; each is independently shippable and verified against current
outputs (golden-file diff on the ATR/CXCR4 reports). Ordered by payoff:

| Extract → shared module | Removes duplication in |
|---|---|
| `structbio/io/pdb.py` (ATOM/HETATM parsing) | **6 skills** |
| `structbio/io/rcsb.py` + `http.py` | **6 skills** |
| `structbio/report/html.py` (HTML shell) | **5 skills** |
| `structbio/io/mmcif.py` | 4 skills |
| `structbio/align.py` (Kabsch / RMSD) | 2 skills |

Add a thin `tests/` (pytest) as modules are extracted — especially PDB parsing and Kabsch, where
bugs are silent.

## Usage / distribution (after Phase 1)

```bash
# core tools
pip install "git+https://github.com/skblnw/structural-bioinfo-skills.git#subdirectory=plugins/structural-bioinfo"
# with extras
pip install "structbio[pharmacophore] @ git+https://github.com/skblnw/structural-bioinfo-skills.git#subdirectory=plugins/structural-bioinfo"

structbio holo-search Q13535 -o atr_holo.csv
structbio pdb-extract --csv atr_holo.csv --uniprot Q13535 --output-dir raw/
structbio pharmacophore raw/9L45/9L45_clean.pdb -d pharmacophores/
structbio pharmacophore-dashboard pharmacophores/ --raw-dir raw/ --out report.html
```

## Risks / tradeoffs

- **Zero-install claim softens:** preserved via the shim, but the cleanest path becomes a one-time
  `pip install`. The heavy skills already required conda installs, so this mostly affects the 4
  stdlib tools.
- **Two run surfaces** (shim + dispatcher) — mitigated because both import the same module; the
  `--help` matrix in steps 3–4 catches drift.
- **`esm` PyPI name** needs confirming (EvolutionaryScale `esm` vs `fair-esm`) before pinning `[esm]`.
- **Rollback is trivial:** Phase 1 is additive + shims; restoring the original `scripts/*.py` reverts.

## Effort

Phase 1 is genuinely small given every script already has `main()` — mostly moving files + 12 tiny
shims + one `pyproject.toml` + `cli.py` + the dashboard asset-loading tweak. Phase 2 is the larger,
optional work, taken one module at a time.

## Open decisions

1. **Package name** — `structbio` (placeholder). Prefer `sbio` / `struct_bioinfo` / other?
2. **numpy in core vs. a `[report]` extra** — currently in core for one-step installs; flip it to
   keep a truly stdlib-only core.
3. ~~Command style~~ — **resolved: single `structbio` dispatcher with subcommands.**
4. **Scope now** — Phase 1 only, or Phase 1 + start Phase 2's `io/pdb.py`?
