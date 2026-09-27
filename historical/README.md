# Historical scripts and adaptation provenance

The supplied V52 bundle contains exactly **7 working scripts and 23 historical Python files**. The 23 history copies below have `.py.txt` extensions and are **not** part of the scientific run command. Original bytes are retained only in the private work archive. Any sanitized paths are listed in `provenance/history_redactions.json`; scientific expressions are not replaced with expected answers.

[Machine-readable origins](../provenance/FILE_ORIGINS.json), [CSV](../provenance/FILE_ORIGINS.csv), [actual imports](../provenance/DEPENDENCIES.json).

| Historical file | Role | Supported adaptation | Execution of original |
|---|---|---|---|
| `angular_coefficients_v40.py` | scientific calculation, regression or reusable algebra | `scripts/article_appC5_compute_angular_coefficients.py` | NOT RUN |
| `apply_replacements.py` | historical editorial/build orchestration; not a scientific entry point | `none; editorial excluded` | NOT RUN |
| `assemble_radial_v40.py` | scientific calculation, regression or reusable algebra | `scripts/supplement_sec1_2_assemble_contacts.py` | NOT RUN |
| `audit_release_v40.py` | historical editorial/build orchestration; not a scientific entry point | `none; editorial excluded` | NOT RUN |
| `build_article_v40.py` | historical editorial/build orchestration; not a scientific entry point | `none; editorial excluded` | NOT RUN |
| `build_local_catalogue_v40.py` | historical editorial/build orchestration; not a scientific entry point | `none; editorial excluded` | NOT RUN |
| `build_mathematical_catalogue_v40.py` | historical editorial/build orchestration; not a scientific entry point | `none; editorial excluded` | NOT RUN |
| `check_catalogue_v40.py` | scientific calculation, regression or reusable algebra | `scripts/article_appC3_verify_contact_permutations.py` | NOT RUN |
| `complete_angular_global_v40.py` | scientific calculation, regression or reusable algebra | `scripts/article_appC7_complete_global_exchange.py` | NOT RUN |
| `generate_radial_v40.py` | scientific calculation, regression or reusable algebra | `scripts/supplement_sec1_3_generate_radial_sources.py` | NOT RUN |
| `n02_background_diagnostics.py` | scientific calculation, regression or reusable algebra | `scripts/article_appD4_compute_background_diagnostics.py` | NOT RUN |
| `noncom_original.py` | scientific calculation, regression or reusable algebra | `lib/article_appC_noncom_engine.py` | NOT RUN |
| `radial_coeff_algebra.py` | scientific calculation, regression or reusable algebra | `lib/radial_coeff_algebra.py` | NOT RUN |
| `radiation_adm_v40.py` | scientific calculation, regression or reusable algebra | `lib/radiation_adm_v40.py` | NOT RUN |
| `run_verifications_v40.py` | historical editorial/build orchestration; not a scientific entry point | `none; editorial excluded` | NOT RUN |
| `verify_analytic_v40.py` | scientific calculation, regression or reusable algebra | `scripts/article_appsDE_verify_reduced_identities.py` | NOT RUN |
| `verify_angular_controls_v40.py` | scientific calculation, regression or reusable algebra | `scripts/article_appC5_verify_angular_controls.py` | NOT RUN |
| `verify_editorial_changes_v40.py` | historical editorial/build orchestration; not a scientific entry point | `none; editorial excluded` | NOT RUN |
| `verify_entropy_v44.py` | scientific calculation, regression or reusable algebra | `scripts/article_sec7_3_verify_entropy.py` | NOT RUN |
| `verify_local_vertices_v40.py` | scientific calculation, regression or reusable algebra | `scripts/supplement_sec1_6_verify_local_vertices.py` | NOT RUN |
| `verify_noncom_v40.py` | historical coefficient-file reader, superseded by fresh derivation and audit | `scripts/article_appC3_compute_noncom_matrix_elements.py` | NOT RUN |
| `verify_radial_controls_v40.py` | scientific calculation, regression or reusable algebra | `scripts/article_appC3_verify_radial_controls.py` | NOT RUN |
| `verify_soft_v40.py` | scientific calculation, regression or reusable algebra | `scripts/article_appC6_verify_soft_limit.py` | NOT RUN |

The historical catalogue generators are not used to overwrite the published comparison object. All parameterized coefficient generators start from existing ADM/source algorithms. The missing non-COM table is regenerated, not supplied as an expected answer. Reconstructed figures 1–2 and the mixed-canonical/printed-coefficient comparison routines are identified as new code, not rediscovered historical files.
