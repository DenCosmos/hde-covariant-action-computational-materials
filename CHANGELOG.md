# Changelog

## V52-R1 — author-selected licensing update, 27 September 2026

Added the author-approved restrictive LICENSE (all rights reserved with limited permission for unmodified, non-commercial scholarly verification), using the author name already present in the source metadata. Added LICENSE_SCOPE.md; updated README, CFF/BibTeX descriptions and the local licensing decision record. No invented SPDX identifier or public URL is used. Removed the superseded LICENSE_PENDING.md.

The current publication guide now covers GitHub without a Zenodo prerequisite. The metadata helper respects `publication_target: github` and does not require a DOI in that route; it preserves the selected licence scope and only inserts an actual optional licence URL. The archive helper has an explicit `--github-only` option. These are publication-tool changes, not scientific algorithm changes.

Scientific scripts, modules, parameters, results, historical source programs, the supplement PDF/TeX and publication-to-code maps are unchanged. Prior scientific validation records remain historical evidence and were not represented as newly executed. Checksums and static packaging/metadata checks are regenerated for this archive. GitHub's independent service terms, including forking and AI-related grants, are disclosed in LICENSE_SCOPE.md. No remote action was performed. The earlier Zenodo archive is not updated by this step.

## V52-R1 — local preparation, 27 September 2026

The computational version is new; the scientific article and mathematical supplement remain V52. No publication date is asserted.

Restored the supplied exact ADM dependencies and radial/angle pipeline, replacing broken historical paths and missing-data-only checks with actual calculations. Added an explicitly new fixed non-COM driver with exact frequency grouping, checkpointing and independent finite-window numerical witnesses. Extracted only the pure half-angle transform from the historical catalogue editor. Historical manuscript generation is not executed.

Adapted the seven existing working checks to themed filenames and explicit inputs. Restored the eight radial/contact/angular/soft drivers, added strict independent comparison with the unchanged published coefficient functions, and added canonical mixed-normalization tests. A new driver reruns the preserved finite-time, scalar-angular, homogeneous and nonlinear checks.

Reconstructed Figures 1 and 2 from the stated equations with explicit configuration/data and diagnostics; retained all original figures and the supplied Figure 3 algorithm. Added a quick/full runner, offline calculation guard, per-command logs, fingerprinted heavy checkpoints, document namespaces, fresh compilation-derived numbering and bidirectional maps.

No scientific formula was altered to match code. Changes to algorithms/interfaces and before/after comparisons are in `provenance/ADAPTATION_REGISTRY.json` and validation records. In particular, the global-completion summary now counts actual positions rather than hard-coding 74 for a custom j list; angular checkpoint compatibility is now validated. No arithmetic coefficients were changed. The Matplotlib adjustment is to the isolated execution environment, not the scientific code.

The inaccurate catalogue attribution following article Eq. (356) is addressed only by a proposed local replacement in `publication/PROPOSED_ARTICLE_REFERENCE_CORRECTION.md`. It was not silently applied to the article. Publication metadata and licensing remain author-controlled.
