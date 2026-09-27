# Reproducing V52-R1

## Reference documents and scope

The scientific references are the author's article V52 (one `main.tex`, Appendices A–E including C.3) and the standalone mathematical supplement V52 (one `mathematical_catalogue.tex`, plus bibliography). Their original SHA-256 values are in [provenance/input_source_hashes.json](provenance/input_source_hashes.json). Both were freshly compiled to recover numbering; [validation/fresh_numbering.json](validation/fresh_numbering.json) records the comparison with the original label map. No catalogue editor or historical manuscript generator is run.

`publication/PUBLICATION_INDEX.json` lists all labels, source lines, document namespaces, numbers, pages and coverage classifications. `publication/SCRIPT_REFERENCE_MAP.json` and `.csv` add program/function/test, inputs, outputs, verification method, tolerances and actual run status. Unmapped analytical claims are not silently represented as computationally proved.

## Tested environment; installation is separate

Python 3.13.5; Linux x86_64; SymPy 1.14.0, mpmath 1.3.0, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8. The precise interpreter/platform and transitive distribution versions are recorded with the runs. Only that environment was executed. Windows/WSL instructions below are examples for the author, not a cross-platform compatibility claim.

Create a virtual environment and install `requirements.txt` using an internet connection, then run the calculations without one. On Linux:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

For Windows Python, the activation command is `.venv\Scripts\activate` in `cmd.exe`. A user already working in WSL can use the Linux commands instead. No API key, token or environment credential is required.

**What was actually done here.** A virtual environment with `include-system-site-packages = false` was created. DNS prevented a fresh pip download. The installed distributions were copied by their package metadata into that environment; no packages are bundled in this release. The sandbox's installed Matplotlib `__init__.py` contained a platform-added serializer after the upstream module. That appended hook was removed in the isolated copy to prevent importing `caas_jupyter_tools`; the scientific programs were not modified for it. This is an explicitly adjusted local dependency, not an unmodified upstream wheel installation. [validation/environment_preparation.json](validation/environment_preparation.json) documents the procedure. A fresh internet-enabled installation remains an author-side confirmation step.

`run_all.py` starts child processes with the active interpreter, a local `sitecustomize.py` socket-connection guard and one BLAS/OpenMP thread per process. Dependency installation is not attempted by the runner. No remote reads or writes occur during calculations. The offline guard and isolated import paths were checked separately.

## Commands and exit codes

```bash
python run_all.py --mode quick --output results/quick
python run_all.py --mode full --workers 4 --output results/full --article-source "/path/to/main.tex"
python run_all.py --mode full --workers 4 --output results/full --article-source "/path/to/main.tex" --resume
```

Quick mode is a deliberately smaller scope, not a shortened claim of full reproduction. Full mode includes the separately supplied article transcription test. `PASS` means the current command returned 0; `FAIL` means a calculation/assertion failed; `BLOCKED` means a mandatory input or dependency was unavailable; `NOT RUN` means outside the selected quick scope. Overall exit codes are 0, 1 and 2 respectively for success, a failed test, and blocked mandatory work. A quick run may return 0 while explicitly listing full-only stages as `NOT RUN`.

Use a new output directory for an independent regeneration. An existing directory requires `--resume` and matching algorithm/input hashes. Heavy non-COM, radial and angular records carry their own checkpoint fingerprints. Resume validates data and reruns algebraic audits; it never obtains a new PASS by reading an old `success=true`. Fast tests and derived comparisons run again. A source edit or changed input normally requires a new output directory. Cosmetic Python docstrings/comments are excluded from the runner's semantic hash, but complete source hashes remain in the release manifest.

All child stdout/stderr, command, return code and duration are saved under the selected `logs/`. Long steps update checkpoint files and elapsed-time messages approximately every minute. The non-COM generator saves one channel at a time; the radial generator saves one pair/bare component at a time; the angular generator saves one matrix position at a time. No reference result is overwritten. Partial directories are not complete outputs.

## Calculation graph

The short branch is: reduced identities → local scalar vertices → printed-catalogue checks; background diagnostics → explicit article transcription; Figures 1–3; canonical mixed normalization; preserved finite-time/nonlinear regressions.

The heavy radiation branch is: exact ADM generator (81 bare components and 54 pair sources) → 81 assembled contacts → permutation and radial controls → exact printed-vs-regenerated comparison. The 54/81 records also feed helicity/transfer/Wigner projection → 74 positions → separate homogeneous exchange completion → angular controls. The same radial data feed the joint soft-limit checks. Historical manuscript-building scripts are not dependencies.

The independent fixed-kinematics branch is: exact non-COM engine → 81 channel checkpoints → recombined coefficients and nonzero registry → independent finite-window quadrature witnesses. Neither the published number 57 nor an absent old JSON table is used to generate amplitudes.

## Normalizations, regulators and error criteria

The radiation calculations use U=0, no independent matter, g_R=k=1 for the tabular dimensionless coefficients and r=k'/k>0. COM momenta and half-angle t=tan(theta/2)>0 are explicit in the generator. Linear external tensor polarizations have squared norm 2. The internal crossed tensor basis has norm 2q²; the relative orientation sign is retained in contractions. The isolated homogeneous scalar has a negative kinetic sign and is contracted separately from the five homogeneous tensor directions. Gaussian covariance and volume are not silently fixed to numbers.

The non-COM engine uses the different, **unit-norm** deterministic polarization frame documented in `parameters/noncom.json`: four momenta of magnitude 2, scalar speed 1/sqrt(3), tensor speed 1, signs ++--. Its 57 count is obtained by exact nonzero grouped Laurent coefficients. Cross-comparing the two polarization conventions without their conversion factors would be incorrect.

Time endpoints remain 0<x_i<x_f<infinity. The independent quadrature uses the three intervals in `parameters/noncom_windows.json`, orders 48 and 80, convergence tolerance 1e-8, and a nonzero witness threshold max(1e-9, 100 times the two-order difference). The two-order difference is a convergence diagnostic, **not a rigorous error bound**. Common nonzero external-wave and spatial normalization factors are stripped in the witness coefficient, as in its driver.

Exact symbolic tests require zero in the specified rational/algebraic field. Sampled substitutions are marked as samples, not all-variable identities. Printed decimal comparisons use 0.51 times the last displayed decimal unit. Numerical DOP853 tolerances and initial conditions are documented in [parameters/embedded_parameters.md](parameters/embedded_parameters.md); the figure configurations are read at runtime. No random seed is needed: there is no stochastic calculation.

The soft check uses r=1+mu*q before q tends to zero. Finite cutoffs, finite wave-packet bands, state choice and separate zero modes remain part of the construction. Neither passing coefficients nor cancellation of leading boundary terms establishes a regulator-free unitary limit.

## Figures

Original article PDFs are preserved in `data/reference_figures/`. Regenerated curves, CSV data, PDF and PNG files are under each run. Figure 1 uses the same c=1.2 and initial data, solves the unreduced system and checks the constraint, determinant and round-trip integration. Figure 2 reconstructs the analytic flow/nullclines and a regular early-time trajectory; its explicit new numerical seed/mesh and extended trajectory are documented. The original displays only an early segment. Scientific content was compared visually; byte equality of rendered PDFs is not required. Figure 3 keeps the supplied plotting formula and constants. No article figure was automatically replaced.

## Standalone supplement build

The published comparison source is never regenerated. Copy its source to an empty build folder before compiling:

```bash
python tools/build_supplement.py --output results/supplement_build
```

This optional document step requires pdfLaTeX; BibTeX is used when available, otherwise the supplied `.bbl` is retained. It is not a dependency of the scientific calculations. No fonts, styles or new macros are introduced into the documents. Fresh compilation logs and numbering checks are retained. Check `KNOWN_LIMITATIONS.md` before claiming a complete independent verification of the article.

## Using the supplied heavy checkpoints explicitly

The bundled final full validation used real heavy records generated in preceding runs of this preparation. One initial run was deliberately rescheduled for parallel angular work. A second independent portable run completed all 54 pair-source and 81 bare-contact generations; those 135 records match the first generation. It was deliberately stopped before a corrected path-dependent control. Neither interrupted invocation is labelled a full PASS. The final frozen-code invocation used `--resume`, rechecked each generator fingerprint and actually reran the audits and all downstream controls. See `validation/preparation_run_history.json` and the final run summary.

To repeat that **resumed** validation without mislabelling it as regeneration from zero:

```bash
python tools/seed_checkpoints.py --from-results reference_results/v52_r1 --output results/resumed
python run_all.py --mode full --resume --workers 4 --output results/resumed --article-source "/path/to/V52/article/source/main.tex"
```

The staging helper copies records, clears derived homogeneous-exchange fields for explicit recalculation, and writes `NOT RUN`; it cannot certify them. The source and staged hashes and cleared fields are recorded. Use the earlier command with an empty new output directory and **without** `--resume` for fresh generation. Original computation logs, final tests and copied reference outputs are distinguished in the validation records. All new source/result hashes are retained; no inherited PASS substitutes for an executed test.
