# Standalone mathematical supplement — V52

`V52_mathematical_supplement.pdf` is the compiled document. `source/mathematical_catalogue.tex` contains the entire text, equations and three extended tables in one file. `source/references.bib` and `source/mathematical_catalogue.bbl` are the only other source dependencies. The article source, article PDF, article `.aux`, Python results and historical generators are not needed for compilation.

From `source/`, run `pdflatex mathematical_catalogue.tex` three times. To update bibliography metadata, run `pdflatex mathematical_catalogue.tex`, `bibtex mathematical_catalogue`, and `pdflatex mathematical_catalogue.tex` three times. No network is needed once TeX and the explicitly named packages are installed.

The document accompanies Danylo Yerokhin's main manuscript V52, dated 27 September 2026. It is a locally prepared supplement, not a separately published journal paper. Main-article external numbers are fixed to V52 and were verified from its independent build. Local equation numbers run from 1 to 1145. Do not regenerate this source from `original_scripts/`.

The supplement source has no macros, class files or style files added for the split. The article citation is intentionally an unpublished-manuscript record without an invented arXiv identifier. Complete release metadata and choose a licence before public distribution. See the package-level README for verification scope and the pre-existing non-centre-of-mass provenance issue.
