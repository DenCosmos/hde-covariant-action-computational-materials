# Proposed local reference correction — NOT applied

Document: scientific article V52, Appendix C.3, immediately after Eq. (356), label `eq:radiation-noncom-check-momenta`.

The existing sentence attributes the 57 nonzero non-COM elements to coefficients listed in the printed catalogue. That catalogue is a different COM/radial configuration. The recovered non-COM calculation is now a separate output of V52-R1, using the explicitly stated Eq. (356) momenta.

Minimal replacement, contingent on the saved non-COM verification and actual release metadata:

```latex
The separate coefficient calculation for the kinematics of
Eq.~\eqref{eq:radiation-noncom-check-momenta}, provided in the
computational materials~\cite{YerokhinV52R1Computational}, gives
57 nonzero full matrix elements after combining terms with identical
frequency-dependent exponentials.
```

Use the finalized `CITATION.bib` entry `YerokhinV52R1Computational`, after the author reserves and confirms the actual version DOI and publishes that object. The DOI of the article is a different identifier and must not be copied into this record.

After metadata confirmation, also update the material-availability paragraph to name computational release V52-R1, tag v52-r1, the exact version DOI and the checked commit. Keep the existing mathematical-supplement citation consistent with its unchanged scientific V52 source; the new computational citation can be added separately. Do not claim that a reserved DOI is already published. Keep each document in its original one-main-TeX-file structure. No catalogue, code or full local work archive is added to the article or arXiv source archive.

The exact printed equations, their signs, conditions and caveats remain unchanged. The source files in this delivery have **not** received this proposed replacement. After any later approved editorial edit, regenerate the article's numbering map and update the source fingerprint rather than calling the changed text the verified original.
