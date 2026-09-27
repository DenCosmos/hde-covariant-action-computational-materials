# Script–publication reference map

Scientific text **V52**, computational candidate **V52-R1**. Fresh numbering: 1669 labelled elements in two separate namespaces; 1018 have a computation, test or input/context relation. 8695 explicit relation rows are recorded. A relation is not a promise of an independent proof of every referenced formula.

Use [REFERENCE_MAP.csv](publication/REFERENCE_MAP.csv) or [REFERENCE_MAP.json](publication/REFERENCE_MAP.json) to search by printed number **and document**, LaTeX label, test identifier or function. The [complete publication index](publication/PUBLICATION_INDEX.json) also records elements without a separate test; the [reverse program index](publication/PROGRAM_INDEX.json) lists every related element for each program. All printed numbers and pages were obtained from [fresh compilation](validation/fresh_numbering.json), not historical filename arithmetic.

Results below are the actually executed [final full run](reference_results/v52_r1/run_summary.json). Heavy records were newly computed earlier in the same preparation and validated on resume; checks and derived comparisons ran again. This is not a claim that the final invocation recalculated every heavy record from zero. Output paths are relative to the selected results directory.

## 1. `article_appsDE_verify_reduced_identities.py`

**Role:** symbolic identities. **Status:** PASS.
**Method:** Exact auxiliary elimination, differentiation, and polynomial identities
**Input:** Embedded symbolic actions; no catalogue input
**Outputs:** `analytic_v40/checks.csv`; `analytic_v40/run_summary.json`.
**Origin:** `GitHub_bundle/scripts/verify_analytic.py`. **Limit:** Algebraic identities, not independent existence/compactness proofs or regulator-free unitarity.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `quadratic_and_crossing` in `scripts/article_appsDE_verify_reduced_identities.py` | article §D.1, `eq:scalar-finite-k-complete-action` (441)<br>article §D.1, `eq:scalar-finite-k-complete-coefficients` (442)<br>article §D.1, `eq:scalar-finite-k-complete-equation` (443)<br>article §7.1, `eq:crossing-linear-leading-coefficients` (210) |
| `canonical_and_late` in `scripts/article_appsDE_verify_reduced_identities.py` | article §D.2, `eq:power-finite-k-linear-equation` (448)<br>article §D.2, `eq:power-exact-Q-coefficients` (449)<br>article §D.2, `eq:power-canonical-mass-rational` (450)<br>article §D.2, `eq:power-canonical-mass-limits` (451) |
| `background_global` in `scripts/article_appsDE_verify_reduced_identities.py` | article §7.7, `eq:power-global-dimensionless-system` (256)<br>article §7.7, `eq:v33-lyapunov` (254)<br>article §7.7, `eq:v33-second-desitter-instability` (259) |
| `gauge_generator` in `scripts/article_appsDE_verify_reduced_identities.py` | article §5, `sec:hamiltonian` (5) |
| `scalar_polynomials_and_sources` in `scripts/article_appsDE_verify_reduced_identities.py` | article §E, `app:quartic-explicit-sources` (E) |
| `nonresonant_majorant` in `scripts/article_appsDE_verify_reduced_identities.py` | article §E.3, `eq:local-majorant-tree-polynomial` (476)<br>article §E.3, `eq:local-majorant-explicit-Ckappa` (477) |

## 2. `supplement_sec1_6_verify_local_vertices.py`

**Role:** symbolic generation and checks. **Status:** PASS.
**Method:** Finite Lagrange inversion versus repeated implicit differentiation
**Input:** Implicit equation s L+L^5 Uprime+4 Phi=0; embedded generic derivatives; radiation U=0 comparison
**Outputs:** `local_vertices_v40/derivatives.txt`; `local_vertices_v40/results.json`; `local_vertices_v40/checks.csv`.
**Origin:** `GitHub_bundle/scripts/verify_local_vertices.py`. **Limit:** 27 derivatives and five scalar local vertices; no fifth/sixth-order metric vertices.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `finite_derivative` in `scripts/supplement_sec1_6_verify_local_vertices.py` | article §C.9, `eq:local-sixth-all-derivatives` (419)<br>article §C.9, `eq:local-sixth-lagrange-proof` (421)<br>supplement §1.6, `eq:cat-local-W-1-0` (551)<br>supplement §1.6, `eq:cat-local-W-0-1` (552)<br>supplement §1.6, `eq:cat-local-W-2-0` (553)<br>supplement §1.6, `eq:cat-local-W-1-1` (554)<br>supplement §1.6, `eq:cat-local-W-0-2` (555)<br>supplement §1.6, `eq:cat-local-W-3-0` (556)<br>supplement §1.6, `eq:cat-local-W-2-1` (557)<br>supplement §1.6, `eq:cat-local-W-1-2` (558)<br>supplement §1.6, `eq:cat-local-W-0-3` (559)<br>supplement §1.6, `eq:cat-local-W-4-0` (560)<br>supplement §1.6, `eq:cat-local-W-3-1` (561)<br>supplement §1.6, `eq:cat-local-W-2-2` (562)<br>supplement §1.6, `eq:cat-local-W-1-3` (563)<br>supplement §1.6, `eq:cat-local-W-0-4` (564)<br>supplement §1.6, `eq:cat-local-W-5-0` (565)<br>supplement §1.6, `eq:cat-local-W-4-1` (566)<br>supplement §1.6, `eq:cat-local-W-3-2` (567)<br>supplement §1.6, `eq:cat-local-W-2-3` (568)<br>supplement §1.6, `eq:cat-local-W-1-4` (569)<br>supplement §1.6, `eq:cat-local-W-0-5` (570)<br>supplement §1.6, `eq:cat-local-W-6-0` (571)<br>supplement §1.6, `eq:cat-local-W-5-1` (572)<br>supplement §1.6, `eq:cat-local-W-4-2` (573)<br>supplement §1.6, `eq:cat-local-W-3-3` (574)<br>supplement §1.6, `eq:cat-local-W-2-4` (575)<br>supplement §1.6, `eq:cat-local-W-1-5` (576)<br>supplement §1.6, `eq:cat-local-W-0-6` (577) |
| `main` in `scripts/supplement_sec1_6_verify_local_vertices.py` | article §C.9, `eq:local-sixth-complete-vertices` (423)<br>supplement §1.6, `eq:cat-local-L-2` (578)<br>supplement §1.6, `eq:cat-local-L-3` (579)<br>supplement §1.6, `eq:cat-local-L-4` (580)<br>supplement §1.6, `eq:cat-local-L-5` (581)<br>supplement §1.6, `eq:cat-local-L-6` (582) |

## 3. `article_sec7_3_verify_entropy.py`

**Role:** symbolic identities. **Status:** PASS.
**Method:** Symbolic entropy derivatives and dimensional checks
**Input:** Embedded entropy parameters r,nu,gamma,m,enthalpy positive
**Outputs:** `logs/entropy.log`.
**Origin:** `GitHub_bundle/scripts/verify_entropy.py`. **Limit:** Not a general nonequilibrium thermodynamic derivation; dimensional checks are consistency checks.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/article_sec7_3_verify_entropy.py` | article §7.3, `eq:generalized-apparent-entropy` (218)<br>article §7.3, `eq:generalized-entropy-background` (219)<br>article §7.3, `eq:generalized-entropy-matching` (220)<br>article §7.3, `eq:generalized-entropy-no-ghost` (221)<br>article §7.3, `eq:generalized-entropy-critical-radius` (222)<br>article §7.3, `eq:generalized-entropy-time-derivative` (223) |

## 4. `supplement_sec1_verify_printed_catalogue.py`

**Role:** transcription/structural checks. **Status:** PASS.
**Method:** Parse 331 rational functions and reconstruction/vertex tables; exact identities in printed coefficients
**Input:** supplement/source; newly calculated local_vertices_v40
**Outputs:** `catalogue_checks.json`.
**Origin:** `GitHub_bundle/scripts/verify_catalogue.py`. **Limit:** Printed-coefficient checks, not independent ADM derivation; the separate strict comparator supplies that comparison.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/supplement_sec1_verify_printed_catalogue.py` | supplement §1.2, `sec:cat-contacts` (1.2)<br>supplement §1.2.14, `sec:cat-component-table` (1.2.14)<br>supplement §1.3, `sec:cat-pairs` (1.3)<br>supplement §1.6, `sec:cat-local-vertices` (1.6) |

## 5. `article_appD4_compute_background_diagnostics.py`

**Role:** numerical diagnostics and symbolic checks. **Status:** PASS.
**Method:** DOP853 integration at stored tolerances; auxiliary linear algebra and asymptotic slopes
**Input:** Embedded parameter grid; see parameters/embedded_parameters.md; runtime output argument only
**Outputs:** `n02_background_diagnostics/numerical_results.json`; `n02_background_diagnostics/saddle_rates.csv`; `n02_background_diagnostics/crossing_convergence.csv`; `n02_background_diagnostics/late_rates.csv`.
**Origin:** `GitHub_bundle/scripts/n02_background_diagnostics.py`. **Limit:** Numerical examples and tolerances, not proofs of global dynamics.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `symbolic_reduction` in `scripts/article_appD4_compute_background_diagnostics.py` | article §D.1, `eq:scalar-finite-k-complete-coefficients` (442)<br>article §D.2, `eq:power-canonical-mass-rational` (450) |
| `saddle_numerics` in `scripts/article_appD4_compute_background_diagnostics.py` | article §7.4, `eq:power-law-full-background` (231)<br>article §7.4, `eq:power-law-full-eigenvalues` (232)<br>article §7.4, `tab:v31-numerical-dynamics` (3) |
| `crossing_numerics` in `scripts/article_appD4_compute_background_diagnostics.py` | article §D.4, `eq:numerical-crossing-point-v40` (460)<br>article §D.4, `eq:numerical-crossing-initial-background-v40` (461)<br>article §D.4, `tab:crossing-reproduction-v40` (5) |
| `late_symbolic` in `scripts/article_appD4_compute_background_diagnostics.py` | article §D.2, `eq:power-canonical-mass-limits` (451) |
| `late_numerics` in `scripts/article_appD4_compute_background_diagnostics.py` | article §D.4, `eq:numerical-late-initial-data-v40` (462)<br>article §D.4, `eq:numerical-curvature-envelope-v40` (463)<br>article §D.4, `tab:late-rate-convergence-v40` (6) |

## 6. `article_V52_verify_source_transcription.py`

**Role:** editorial/structural and numerical transcription. **Status:** PASS.
**Method:** Resolve references; parse 18 table values; compare to displayed-decimal tolerances
**Input:** Explicit --article-source main.tex or folder, references.bib, three figures; supplementary source; newly computed numerical CSVs
**Outputs:** `release_checks.json`.
**Origin:** `GitHub_bundle/scripts/verify_release.py`. **Limit:** No article redistributed. Saddle maxima are expected test constants; other 18 values parsed. No .aux needed.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/article_V52_verify_source_transcription.py` | article §D.4, `app:numerical-protocols` (D.4)<br>article §D.4, `tab:crossing-reproduction-v40` (5)<br>article §D.4, `tab:late-rate-convergence-v40` (6)<br>article §7.4, `tab:v31-numerical-dynamics` (3) |

## 7. `article_fig01_plot_pole_crossing.py`

**Role:** new numerical figure reconstruction. **Status:** PASS.
**Method:** DOP853 integration of unreduced Euler equation plus differentiated constraint; determinant/residual/reversibility checks
**Input:** parameters/figure01.json read at runtime
**Outputs:** `figure01/pole_crossing_c12.csv`; `figure01/pole_crossing_c12.pdf`; `figure01/checks.json`.
**Origin:** `NEW implementation from V52 equations; original figure retained`. **Limit:** One specified regular trajectory. Scientific curves, not byte identity to original PDF.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `coefficients` in `scripts/article_fig01_plot_pole_crossing.py` | article §6.4, `eq:standard-pole-ratio` (139) |
| `rhs` in `scripts/article_fig01_plot_pole_crossing.py` | article §6.4, `eq:unreduced-pole-action` (136)<br>article §6.4, `eq:pole-A-explicit` (137)<br>article §6.4, `eq:pole-unreduced-constraint` (140)<br>article §6.4, `eq:pole-constraint-coefficients` (141)<br>article §6.4, `eq:pole-cauchy-determinant` (145) |
| `main` in `scripts/article_fig01_plot_pole_crossing.py` | article §6.4, `fig:pole-crossing-c12` (1) |

## 8. `article_fig02_plot_phase_portrait.py`

**Role:** new numerical figure reconstruction. **Status:** PASS.
**Method:** DOP853 phase trajectory, fixed points, eigenvalues, nullclines, half-seed sensitivity test
**Input:** parameters/figure02.json read at runtime; seed and plot mesh are explicit new reconstruction choices
**Outputs:** `figure02/phase_portrait_c12.pdf`; `figure02/phase_portrait_trajectory.csv`; `figure02/phase_portrait_field.csv`; `figure02/phase_portrait_nullclines.csv`; `figure02/checks.json`.
**Origin:** `NEW implementation from V52 phase system; original figure retained`. **Limit:** Background phase portrait, not perturbation stability. Extended trajectory; original only displays an early segment.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `vector_field` in `scripts/article_fig02_plot_phase_portrait.py` | article §7.6, `eq:compact-phase-system` (242) |
| `main` in `scripts/article_fig02_plot_phase_portrait.py` | article §7.6, `fig:phase-portrait` (2)<br>article §7.6, `eq:compact-fixed-points` (243)<br>article §7.6, `eq:compact-eigenvalues` (244)<br>article §7.6, `eq:early-regular-manifold` (245) |

## 9. `article_fig03_plot_strong_coupling.py`

**Role:** figure generation. **Status:** PASS.
**Method:** Evaluate the two cubic operator estimates and their minimum
**Input:** Embedded c grid and dimensionless normalization from original plotting script
**Outputs:** `figures/strong_coupling_scale_standard.csv`; `figures/strong_coupling_scale_standard.pdf`.
**Origin:** `GitHub_bundle/scripts/plot_strong_coupling.py`. **Limit:** Operator estimates are not an independently determined physical unitarity threshold.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/article_fig03_plot_strong_coupling.py` | article §A.11, `fig:strong-coupling-standard` (3)<br>article §A.11, `eq:main-standard-cubic-scale` (315)<br>article §A.11, `eq:main-standard-cubic-numerical` (317) |

## 10. `article_appC3_compute_noncom_matrix_elements.py`

**Role:** restored exact computation. **Status:** PASS.
**Method:** Exact multilinear ADM expansion, auxiliary and Legendre contact, all 18 ordered exchanges per external component
**Input:** parameters/noncom.json read at runtime; lib/article_appC_noncom_engine.py
**Outputs:** `noncom/coefficients.json`; `noncom/component_registry.csv`; `noncom/verification.json`.
**Origin:** `NEW checkpointed driver of preserved original_scripts/noncom_original.py engine; replaces missing-data-only verifier`. **Limit:** All 81 separate non-COM combinations; 57 tested only after computation. No infinite-time or regulator-removal claim.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/article_appC3_compute_noncom_matrix_elements.py` | article §C.3, `eq:radiation-noncom-check-momenta` (356)<br>article §C.2, `eq:v34r-exact-scalar-density` (342)<br>article §C.2, `eq:v34r-exact-gravity-density` (343)<br>article §6.10, `eq:v33-legendre-contact` (191) |
| `audit` in `scripts/article_appC3_compute_noncom_matrix_elements.py` | article §C.3, `eq:radiation-noncom-check-momenta` (356)<br>article §C.3, `eq:radiation-contact-decomposition` (351)<br>article §C.3, `eq:radiation-ordered-pair-reconstruction` (354) |

## 11. `article_appC3_verify_noncom_finite_windows.py`

**Role:** new numerical witness checks. **Status:** PASS.
**Method:** Independent 48/80-order nested quadrature over three positive finite intervals
**Input:** Newly computed noncom/coefficients.json; parameters/noncom_windows.json read at runtime
**Outputs:** `noncom/finite_windows.csv`; `noncom/finite_window_checks.json`.
**Origin:** `NEW independent nested Gauss-Legendre evaluation of the restored V52 finite-window kernel`. **Limit:** Nonzero witnesses and order-convergence comparison, not a rigorous quadrature bound or all-window identity proof.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/article_appC3_verify_noncom_finite_windows.py` | article §C.3, `eq:radiation-noncom-check-momenta` (356)<br>article §6.10, `eq:v33-radial-kernel` (192)<br>article §C.4, `eq:v34-elementary-single` (359)<br>article §C.4, `eq:v34-elementary-ordered` (360) |

## 12. `supplement_sec1_3_generate_radial_sources.py`

**Role:** restored exact ADM generation. **Status:** PASS.
**Method:** Exact rational field ADM expansion; 81 bare components, 54 unordered external-pair species combinations
**Input:** Embedded symbolic r>0,t=tan(theta/2)>0; --task/--start/--stop/--resume control generation
**Outputs:** `radial_generic/bare/*.json`; `radial_generic/pair/*.json`.
**Origin:** `GitHub_bundle/original_scripts/generate_radial_v40.py; preserved radiation_adm_v40.py`. **Limit:** TT norm 2 externally; internal cross norm 2q^2; homogeneous scalar negative kinetic sign preserved. Odd-reflection bare zeros are exact symmetry exclusions.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `real_basis` in `scripts/supplement_sec1_3_generate_radial_sources.py` | article §C.5, `eq:v34-linear-polarizations` (364)<br>article §C.3, `eq:radiation-vertex-coefficient-basis` (353) |
| `pair_data` in `scripts/supplement_sec1_3_generate_radial_sources.py` | article §C.3, `eq:radiation-cubic-source-definition` (352)<br>article §C.2, `eq:v34r-nonzero-constraints` (345)<br>article §C.2, `eq:v34r-global-reduction` (348)<br>supplement §1.3, `eq:cat-pair-polynomials` (56) |
| `main` in `scripts/supplement_sec1_3_generate_radial_sources.py` | article §C.3, `eq:radiation-catalogue-external-momenta` (350)<br>supplement §1.1, `eq:cat-momenta` (1) |

## 13. `supplement_sec1_2_assemble_contacts.py`

**Role:** restored exact contact computation. **Status:** PASS.
**Method:** Contract three auxiliary partitions and signed Legendre products, H4=-L4+aux+Legendre
**Input:** 54 new pair records and 81 new bare records
**Outputs:** `radial_generic/contacts/*.json`; `radial_generic/assembly_summary.json`.
**Origin:** `GitHub_bundle/original_scripts/assemble_radial_v40.py`. **Limit:** Radiation COM frame; not the separate fixed non-COM configuration.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/supplement_sec1_2_assemble_contacts.py` | article §C.3, `eq:radiation-contact-decomposition` (351)<br>article §6.10, `eq:v33-legendre-contact` (191)<br>supplement §1.5.2, `eq:cat-partition-contact-sum` (498)<br>supplement §1.5.2, `eq:cat-partition-velocity-nonzero` (499)<br>supplement §1.5.2, `eq:cat-partition-velocity-homogeneous` (500) |

## 14. `article_appC3_verify_contact_permutations.py`

**Role:** restored exact symmetry checks. **Status:** PASS.
**Method:** Within-pair swaps and pair transposition with r-dependent scaling, angle sign and complex conjugation
**Input:** 81 newly assembled contact dictionaries
**Outputs:** `radial_generic/permutation_checks.json`; `radial_generic/symmetry_check.json`; `radial_generic/reconstruction.csv`.
**Origin:** `GitHub_bundle/original_scripts/check_catalogue_v40.py`. **Limit:** 243 full-contact identities plus 729 component identities = 972. Counts are distinct and not substituted for each other.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/article_appC3_verify_contact_permutations.py` | article §C.3, `eq:radiation-contact-permutations` (355)<br>supplement §1.2, `eq:cat-contact-symmetries` (3)<br>supplement §1.2.14, `sec:cat-component-table` (1.2.14) |

## 15. `article_appC3_verify_radial_controls.py`

**Role:** restored exact/sample control checks. **Status:** PASS.
**Method:** 18 scalar/tensor source and homogeneous checks, including six sampled crossing controls
**Input:** New pair records; explicit generic control formulas
**Outputs:** `radial_generic/independent_controls.csv`.
**Origin:** `GitHub_bundle/original_scripts/verify_radial_controls_v40.py`. **Limit:** Some controls are exact substitutions at selected points, not identities at every angle.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/article_appC3_verify_radial_controls.py` | article §C.3, `eq:radiation-cubic-source-definition` (352)<br>article §C.7, `eq:v34-global-source-polynomials` (401) |

## 16. `article_appC5_compute_angular_coefficients.py`

**Role:** restored mixed angular computation. **Status:** PASS.
**Method:** Exact helicity conversion, Wigner polynomial, q-transfer Laurent coefficients; all scalar and tensor internal exchanges
**Input:** New 54 pair and 81 contact records; default j=0,2,3,4
**Outputs:** `angular_v40/j*.json`; `angular_v40/block_registry_*.csv`.
**Origin:** `GitHub_bundle/original_scripts/angular_coefficients_v40.py`. **Limit:** 74 coefficient positions, not 74 finite energy scattering eigenamplitudes. Finite cutoffs and state parameters remain explicit.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `wigner` in `scripts/article_appC5_compute_angular_coefficients.py` | article §C.5, `eq:v34-wigner-polynomial` (373)<br>supplement §1.5, `eq:cat-angular-explicit-wigner-coefficients` (492) |
| `transfer` in `scripts/article_appC5_compute_angular_coefficients.py` | article §C.5, `eq:v34-transfer-variable` (370)<br>supplement §1.5, `eq:cat-angular-source-transfer` (490) |
| `main` in `scripts/article_appC5_compute_angular_coefficients.py` | article §C.5, `eq:v34-helicity-transform` (367)<br>article §C.5, `eq:angular-block-complete-order` (379)<br>article §C.5, `eq:angular-exchange-finite-reconstruction` (377)<br>article §C.5, `eq:angular-contact-finite-reconstruction` (378)<br>article §C.6, `eq:boundary-finite-coefficient-reconstruction` (392)<br>supplement §1.5, `eq:cat-angular-exchange-convolution` (494)<br>supplement §1.5, `eq:cat-angular-contact-convolution` (495)<br>supplement §1.5, `eq:cat-boundary-explicit-basis` (496) |

## 17. `article_appC7_complete_global_exchange.py`

**Role:** restored homogeneous mixed exchange computation. **Status:** PASS.
**Method:** Project AA,AB,BA,BB products for scalar and five homogeneous tensor directions; both time orderings
**Input:** New radial source and 74 angular dictionaries; symbolic Gaussian covariance
**Outputs:** `angular_v40/complete_summary.json`; `angular_v40/block_registry.csv`; `angular_v40/j*_*.json`.
**Origin:** `GitHub_bundle/original_scripts/complete_angular_global_v40.py; pure to_z extracted from editorial generator without executing it`. **Limit:** 12 nonhomogeneous plus 12 homogeneous ordered exchange blocks per position; no covariance or volume silently selected.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `angular_moment` in `scripts/article_appC7_complete_global_exchange.py` | article §C.5, `eq:v34-angular-moments` (371) |
| `main` in `scripts/article_appC7_complete_global_exchange.py` | article §C.7, `eq:v34r-global-free-flow` (400)<br>article §C.7, `eq:v34-global-source-polynomials` (401)<br>article §C.7, `eq:v34-global-contraction` (402)<br>article §C.7, `eq:v34-isolated-zero-distribution` (403)<br>article §C.6, `eq:v34r-boundary-table-mixed-dictionary` (390)<br>article §C.6, `eq:v34r-boundary-table-square-dictionary` (391) |

## 18. `article_appC5_verify_angular_controls.py`

**Role:** restored angular controls. **Status:** PASS.
**Method:** Wigner d comparison at theta=pi/3 for 74 positions plus three scalar contact moments
**Input:** 74 new angular records and symbolic Wigner formula
**Outputs:** `angular_v40/independent_controls.csv`.
**Origin:** `GitHub_bundle/original_scripts/verify_angular_controls_v40.py`. **Limit:** 77 checks. Wigner point tests alone are not an all-angle identity proof; coefficients were separately generated symbolically.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/article_appC5_verify_angular_controls.py` | article §C.5, `eq:v34-wigner-polynomial` (373)<br>article §C.5, `eq:v34-angular-contact-series` (375) |

## 19. `article_appC6_verify_soft_limit.py`

**Role:** restored exact soft-limit checks. **Status:** PASS.
**Method:** Exact Laurent soft coefficients and five boundary leading groups
**Input:** New radial pairs and contact records; r=1+mu*q before q->0
**Outputs:** `soft_v40/summary.json`; `soft_v40/checks.csv`; `soft_v40/five_groups.txt`.
**Origin:** `GitHub_bundle/original_scripts/verify_soft_v40.py`. **Limit:** Does not prove existence of an unregulated unitary limit. Equal radii are not substituted before the joint soft limit.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/article_appC6_verify_soft_limit.py` | article §6.10, `eq:v34r-soft-scalar-vertex` (195)<br>article §6.10, `eq:v33-soft-singularity` (196)<br>article §C.6, `eq:boundary-five-groups-leading` (393)<br>article §C.6, `eq:boundary-scalar-five-coefficients` (394) |

## 20. `supplement_secs1_2_1_5_compare_recomputed_coefficients.py`

**Role:** new independent printed-coefficient comparison. **Status:** PASS.
**Method:** Exact restored r,z functions compared with printed R maps, reconstruction, 46 crossed auxiliary terms and homogeneous lapse polynomials
**Input:** Unchanged V52 supplementary TeX; independently generated 54 pairs and 81 contacts
**Outputs:** `printed_vs_recomputed/checks.json`.
**Origin:** `NEW parser/comparator; preserved algebra and pure half-angle transform reused`. **Limit:** Catalogue is never regenerated; tests carry exact per-coefficient source labels. No fictitious intermediate data.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `read_printed` in `scripts/supplement_secs1_2_1_5_compare_recomputed_coefficients.py` | supplement §1.2, `sec:cat-contacts` (1.2)<br>supplement §1.3, `sec:cat-pairs` (1.3) |
| `main` in `scripts/supplement_secs1_2_1_5_compare_recomputed_coefficients.py` | supplement §1.2.14, `sec:cat-component-table` (1.2.14)<br>supplement §1.5.2, `sec:cat-contact-partitions` (1.5.2)<br>supplement §1.5.2, `eq:cat-partition-homogeneous-polynomials` (501)<br>supplement §1.5.2, `eq:cat-partition-auxiliary-minus` (503) |

## 21. `article_appC1_verify_mixed_normalization.py`

**Role:** new exact normalization checks. **Status:** PASS.
**Method:** Independent canonical action boundary check, Euler differentiation, symplectic identity, determinant and simple-pole polynomial remainders
**Input:** Generic real symmetric W and real antisymmetric Omega; positive kinetic matrix; simple positive root, positive nonzero norm
**Outputs:** `mixed_normalization/checks.json`.
**Origin:** `NEW implementation directly from unambiguous V52 equations; not found historical code`. **Limit:** 11 identities; no full independent-matter mixed scattering spectrum or all-background gravitational correction is computed.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `main` in `scripts/article_appC1_verify_mixed_normalization.py` | article §C.1, `eq:v34-canonical-matrices` (326)<br>article §C.1, `eq:v34-antisymmetric-matrix` (327)<br>article §C.1, `eq:v34-canonical-frequency` (328)<br>article §C.1, `eq:v34-canonical-action` (329)<br>article §C.1, `eq:v34-canonical-evolution` (330)<br>article §C.1, `eq:v34-symplectic-norm` (331)<br>article §C.1, `eq:v34-frequency-determinant` (332)<br>article §C.1, `eq:v34-frequency-roots` (333)<br>article §C.1, `eq:v34-spectral-residue` (334)<br>article §C.1, `eq:v34-group-velocity` (335) |

## 22. `article_appC_verify_finite_time_and_nonlinear_identities.py`

**Role:** restored exact regression driver. **Status:** PASS.
**Method:** 338 exact tests, including direct ADM points, finite-time/scalar angular/homogeneous and nonlinear identities
**Input:** Embedded exact polynomials and special ADM configurations in lib/article_appC_noncom_engine.py
**Outputs:** `engine_regressions/checks.json`; `engine_regressions/results.json`.
**Origin:** `NEW explicit driver of 13 preserved noncom_original.py test functions`. **Limit:** Per-function scope only. Scalar j=0..8 moments are separate from the 74 mixed angular positions. No regulator-free S matrix.

| Function (file) | Publication elements: document, subsection, label and printed number |
|---|---|
| `test_radiation_expansion` in `lib/article_appC_noncom_engine.py` | article §C.2, `eq:v34r-exact-scalar-density` (342)<br>article §C.2, `eq:v34r-exact-gravity-density` (343)<br>article §C.2, `eq:v34r-extrinsic-density` (344)<br>article §C.2, `eq:v34r-nonzero-constraints` (345) |
| `test_general_Legendre_transform` in `lib/article_appC_noncom_engine.py` | article §6.10, `eq:v33-legendre-contact` (191) |
| `test_complete_example` in `lib/article_appC_noncom_engine.py` | article §C.3, `eq:radiation-noncom-check-momenta` (356) |
| `test_time_integrals` in `lib/article_appC_noncom_engine.py` | article §C.4, `eq:v34-elementary-single` (359)<br>article §C.4, `eq:v34-elementary-ordered` (360)<br>article §C.4, `eq:v34-single-time-series` (361)<br>article §C.4, `eq:v34-ordered-time-series` (362)<br>article §C.4, `eq:v34-time-order-identity` (363) |
| `test_signed_Legendre_transform` in `lib/article_appC_noncom_engine.py` | article §6.10, `eq:v33-legendre-contact` (191)<br>article §C.2, `eq:v34r-global-reduction` (348) |
| `test_homogeneous_and_state` in `lib/article_appC_noncom_engine.py` | article §C.7, `eq:v34r-global-free-flow` (400)<br>article §C.7, `eq:v34-global-contraction` (402) |
| `test_generic_scalar_contact` in `lib/article_appC_noncom_engine.py` | article §C.4, `eq:v34-scalar-leading-contact` (357)<br>article §C.4, `eq:v34-scalar-equal-radii-contact` (358) |
| `test_direct_com` in `lib/article_appC_noncom_engine.py` | article §C.3, `eq:radiation-contact-decomposition` (351)<br>article §C.4, `eq:v34-scalar-equal-radii-contact` (358) |
| `test_angular_projection` in `lib/article_appC_noncom_engine.py` | article §C.5, `eq:v34-angular-exchange-series` (374)<br>article §C.5, `eq:v34-angular-contact-series` (375) |
| `test_radial_kernel` in `lib/article_appC_noncom_engine.py` | article §C.3, `eq:radiation-cubic-source-definition` (352) |
| `test_channel_count_and_leakage` in `lib/article_appC_noncom_engine.py` | article §6.10, `eq:v33-channel-counts` (194)<br>article §C.8, `eq:v34-finite-time-leakage` (411) |
| `test_nonlinear` in `lib/article_appC_noncom_engine.py` | article §C.10, `eq:v34-crossing-series-H` (428)<br>article §C.10, `eq:v34-crossing-series-L` (429)<br>article §C.10, `eq:v34r-Bianchi-monotonicity` (437) |
| `test_simple_wave` in `lib/article_appC_noncom_engine.py` | article §C.10, `eq:v34r-simple-wave-equation` (431)<br>article §C.10, `eq:v34r-simple-wave-characteristics` (432) |

## Reading the scientific scope

The non-COM row anchored at article Eq. (356) concerns its fixed kinematics and the associated 57-count paragraph. The 243 full-contact identities are distinguished from 729 additional component-part identities. The 74 angular outputs are symbolic coefficient dictionaries with mixed scalar/tensor nonhomogeneous and separately completed homogeneous exchange. They are not 74 regulator-free scattering eigenvalues or a numerical S-matrix. The canonical mixed-normalization identities are not a reproduction of the full independent-matter scattering problem. All unassigned index entries retain a noncomputational designation.
