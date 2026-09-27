# Known limits of V52-R1

## Scientific limits retained from V52

The finite-time radiation calculation is for U=0 without independent matter. It includes scalar and tensor external species and their internal exchanges, but it is not the full arbitrary-background mixed matter–gravity scattering problem. The 74 angular positions are coefficient dictionaries with two radial arguments, not a finite on-shell scattering matrix or its eigenvalues. Homogeneous Gaussian covariance, separate spatial-volume normalization and finite regulators remain explicit. The existence of a unitary limit as regulators are removed has **not been established**.

The canonical mixed-field program verifies the stated action, symplectic, determinant and simple-pole identities under their hypotheses. It does not calculate the independent-matter mixed scattering eigenvalues or the full background correction of order E²/M_Pl². The article itself distinguishes the conserved-current gravitational exchange from that full correction. No supplied program independently reproduces that complete arbitrary-background correction; the package does not claim otherwise.

The old 243 count concerns full-contact permutation identities. The restored checker additionally checks bare, auxiliary and Legendre components, giving 972 tests in total. These numbers must not be interchanged. Wigner controls at a selected angle are point tests; the actual coefficient generator uses exact polynomial projection. Numerical examples do not prove existence, compactness, global stability or all asymptotic assertions. Local scalar vertices through sixth order are not metric vertices of the same orders.

## Specific restoration qualifications

Figure 1 and Figure 2 programs are explicitly **new implementations of already stated equations**, not rediscovered original plotting files. Their sampling choices and numerical checks are saved; original PDFs remain unchanged. Figure 2 extends the regular trajectory beyond the short early segment shown in V52. The non-COM driver is new but calls the preserved exact engine; its old missing coefficient file has not been fabricated. The comparison source for the printed catalogue remains unchanged. The pure half-angle transformation was extracted from a historical editorial generator without executing the generator.

The source-transcription checker requires the separately supplied article; it cannot check unpublished article files absent from a reader's copy. Its two saddle maxima are expected constants, while the 18 crossing/late values are parsed from the supplied text. The programme does not prove all equations in the 93-page article. The full label index explicitly records unassigned analytical/context elements.

## Environment and publication limits

A fresh internet-enabled pip installation could not be performed in this execution environment. An isolated virtual environment was populated from installed packages, and a sandbox-only Matplotlib hook was removed there. This is disclosed, not presented as a stock-wheel installation. Other Python versions and operating systems were not executed.

No remote repository, release, Zenodo record or DOI was created. The author has selected the custom restrictive terms in [LICENSE](LICENSE) for original code, text and data; see [LICENSE_SCOPE.md](LICENSE_SCOPE.md) for their scope and independent GitHub platform rights. Zenodo is deferred. Public visibility and author metadata remain subject to the author's final confirmation. The final checked commit hash belongs in an external release record, release description, Zenodo metadata and article citation, not in the commit it identifies.

CITATION.cff validation level is recorded in `validation/citation_validation.json`. Do not call a local YAML/field check an official CFF schema validation. The author should run the official validator after filling actual metadata when that level is not recorded as PASS.

The public historical copies are sanitized provenance, not executable supported entry points. Original files and failed/diagnostic initial logs remain in the local work archive. A successfully prepared archive is not evidence that a scientific stage with FAIL, BLOCKED or NOT RUN has passed. Consult `VALIDATION_SUMMARY.md` and the saved run summaries for actual execution status.
