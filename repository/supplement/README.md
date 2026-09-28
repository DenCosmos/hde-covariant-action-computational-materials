# Standalone mathematical supplement

[`V52_mathematical_supplement.pdf`](V52_mathematical_supplement.pdf) is the corrected English mathematical catalogue. The existing filename is retained to preserve repository links; the document itself has no version designation. [`source/mathematical_catalogue.tex`](source/mathematical_catalogue.tex) is the complete standalone source, including the text, equations, three extended tables and bibliography.

From `source/`, run:

```bash
pdflatex -interaction=nonstopmode -halt-on-error mathematical_catalogue.tex
pdflatex -interaction=nonstopmode -halt-on-error mathematical_catalogue.tex
pdflatex -interaction=nonstopmode -halt-on-error mathematical_catalogue.tex
```

BibTeX is not required. `source/references.bib` and `source/mathematical_catalogue.bbl` are synchronized optional bibliography exports; the source no longer reads them. The article source, article PDF, article `.aux`, Python results and historical generators are not needed for compilation. The explicitly named TeX packages must be installed.

The catalogue accompanies Danylo Yerokhin's manuscript *A minimal local covariant action for holographic dark energy: constraints, perturbations, and nonlinear dynamics*. References explicitly identified as belonging to the main article use its numbering; all other section and equation numbers belong to this supplement. The local equation numbers remain 1--1145.

The English source contains the reviewed translation and editorial corrections synchronized with the corrected Russian catalogue. These clarify indexing, angular reconstruction, normalization and terminology without changing the coefficient data. The computational programs and their reference results are unchanged. Do not regenerate this reviewed source with historical catalogue generators.

See the repository [LICENSE](../LICENSE), [verification scope](../VALIDATION_SUMMARY.md), [known limitations](../KNOWN_LIMITATIONS.md) and [change log](../CHANGELOG.md). This document update does not assert a new computational release, arXiv identifier, DOI or journal reference.
