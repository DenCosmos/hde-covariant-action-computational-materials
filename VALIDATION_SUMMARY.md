# V52-R1 validation summary

**Status:** local unpublished computational candidate. The final full runner returned **PASS / 0** for its documented computational scope. This does not assert an independent proof of every analytical statement of V52.

The frozen portable full invocation took 421.092 seconds and executed 30 stage commands. It explicitly resumed heavy records newly generated in this preparation, not inherited published answers. Short checks and downstream comparisons were rerun. The final invocation must not be advertised as a from-zero generation. Two earlier fresh radial computations gave identical 54 source and 81 bare-contact records.

| Calculation | Actually obtained | Scope |
|---|---|---|
| Fixed non-COM, C.3 (356) | 81 components; 57 exact nonzero; 366 nonzero frequency groups | Grouped exact coefficients, specified momenta and normalization |
| Finite-window witnesses | 57 nonzero on each of 3 intervals; maximum 48/80 order difference 1.154e-12 | Pointwise witnesses, not all-window or zero-regulator proof |
| Radial reconstruction | 54 pair sources; 81 bare and assembled contacts; 41 nonzero full contacts | Radiation COM configuration, distinct from non-COM 57 |
| Contact permutations | 243 full identities; 972 including 729 part identities | Exact rational/algebraic comparisons |
| Printed supplement comparison | 2296 PASS; 0 FAIL | Unchanged V52 coefficient source versus regenerated expressions |
| Angular and homogeneous exchange | 74 positions; 12 nonhomogeneous + 12 homogeneous ordered terms each; 222 structural checks | Coefficient dictionaries at j=0,2,3,4, not scattering eigenvalues |
| Angular / radial controls | 77 / 18 PASS | Wigner controls include an explicitly chosen angle |
| Joint soft limit | 209 PASS | r=1+mu*q, nine hard pairs and 81 contacts; retained boundary terms |
| Preserved engine / mixed canonical | 338 / 11 PASS | Finite-time/nonlinear tests; separate canonical 2x2 normalization identities |
| Seven original working programs | All rerun successfully | 53 analytic, 62 local-vertex, 15 entropy, 186 printed-structure and 35 background symbolic checks, plus numerical/source/figure checks |
| Figures 1–3 | Generated with numerical data | Figures 1–2 are declared new implementations; originals not replaced |

## Final command records

| Stage | Status | Return code | Seconds |
|---|---|---:|---:|
| analytic | PASS | 0 | 4.001 |
| local_vertices | PASS | 0 | 2.001 |
| entropy | PASS | 0 | 1.001 |
| printed_catalogue | PASS | 0 | 6.001 |
| background | PASS | 0 | 6.002 |
| figure01 | PASS | 0 | 3.004 |
| figure02 | PASS | 0 | 2.001 |
| figure03 | PASS | 0 | 1.001 |
| mixed_normalization | PASS | 0 | 1.001 |
| finite_time_identities | PASS | 0 | 54.009 |
| source_transcription | PASS | 0 | 1.001 |
| radial_pair_9_18 | PASS | 0 | 1.002 |
| radial_pair_0_9 | PASS | 0 | 1.002 |
| radial_pair_18_27 | PASS | 0 | 1.002 |
| radial_pair_27_36 | PASS | 0 | 1.001 |
| radial_pair_45_54 | PASS | 0 | 1.001 |
| radial_pair_36_45 | PASS | 0 | 1.002 |
| radial_bare_0_27 | PASS | 0 | 1.002 |
| radial_bare_27_54 | PASS | 0 | 1.001 |
| radial_bare_54_81 | PASS | 0 | 1.004 |
| noncom | PASS | 0 | 5.002 |
| contacts | PASS | 0 | 27.004 |
| permutations | PASS | 0 | 17.003 |
| radial_controls | PASS | 0 | 1.001 |
| printed_vs_recomputed | PASS | 0 | 139.021 |
| angular | PASS | 0 | 52.008 |
| global_exchange | PASS | 0 | 14.002 |
| angular_controls | PASS | 0 | 1.001 |
| soft | PASS | 0 | 69.012 |
| noncom_windows | PASS | 0 | 15.003 |

## Evidence and remaining limits

[Actual results](validation/ACTUAL_RESULTS.json), [full run](reference_results/v52_r1/run_summary.json), [run history](validation/preparation_run_history.json), [two fresh radial generations](validation/two_fresh_radial_generations.json), [maps](SCRIPT_REFERENCE_MAP.md). The maps contain 8695 relation rows and cover 1018 labelled publication elements with test, computation or context links; 1669 labels are indexed in total. Unmapped elements are not silently marked PASS.

The input archive and all 247 extracted originals remained unchanged. The main article and its Appendix C.3 remain separate from the supplement and are not included in the public repository. The supplement compiled autonomously without changing its TeX/bibliography. The exact source-to-printed-number mapping has no differences from V52.

Fresh internet pip installation was BLOCKED by network resolution; the executed virtual environment was isolated and populated from local distributions, with its sandbox-specific Matplotlib change disclosed. The local CFF check passed, but official cffconvert/schema validation was NOT RUN. The later licensing update records the author-selected restrictive LICENSE. The actual repository URL, release date and public-access confirmation are still pending. Zenodo and its DOI are deferred, not requirements of the current GitHub-only route. This update did not rerun the scientific calculations; see `validation/license_update_checks.json`. None of these limitations is presented as a successful public release.

The full arbitrary-background independent-matter gravitational scattering correction/eigenvalues are not reproduced. This package restores the specified radiation scalar/tensor exchanges and verifies canonical mixed-field identities. It does not establish a regulator-free unitary limit. See [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md).
