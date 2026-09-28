# V52-R1 computational materials

**A minimal local covariant action for holographic dark energy: constraints, perturbations, and nonlinear dynamics**  
Danylo Yerokhin — computational package **V52-R1**; the mathematical supplement contains the reviewed English editorial corrections.

This repository provides the author-controlled computational materials and mathematical supplement for the associated article. The English catalogue includes the reviewed corrections described in [CHANGELOG.md](CHANGELOG.md). No new computational release or DOI is asserted by this catalogue update. The author-selected restrictive licence is unchanged. See the [GitHub publication guide in Russian](PUBLICATION_GUIDE_RU.md).

## Copyright and permitted use

**Copyright © 2026 Danylo Yerokhin. All rights reserved.** See [LICENSE](LICENSE) for the limited permission to download, retain, inspect and execute unmodified materials for personal, academic and non-commercial verification of the associated published results. Other uses requiring the holder's permission are not authorized by that grant. This is not an open-source software licence. Citation does not grant reproduction or redistribution rights.

Third-party rights and GitHub's independently applicable platform terms remain in effect. Those platform terms include public-repository forking and certain AI-training rights; this LICENSE does not override them. See [licence scope and platform limitations](LICENSE_SCOPE.md).

## First computation

Install dependencies before going offline. The tested interpreter was Python 3.13.5 on Linux; other Python versions and operating systems were not tested.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_all.py --mode quick --output results/quick
```

The dependency download is an installation step, **not** part of the calculations. In the preparation environment, network installation was unavailable; an isolated virtual environment was populated from installed distributions. This is **not a verified clean PyPI installation**. The exact procedure and environment modification are disclosed in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

Quick mode reruns the short symbolic checks, background diagnostics, all three figures and the finite-time regression identities. It explicitly reports the heavy computations as `NOT RUN`. To regenerate the documented complete computational scope from scratch:

```bash
python run_all.py --mode full --workers 4 --output results/full \
  --article-source "/path/to/V52/article/source/main.tex"
```

Supply the original article source directory, including `references.bib` and its three figures. The article is intentionally **not distributed here**. A missing article source produces `BLOCKED` and exit code 2 in full mode; it does not prevent independent calculations. The calculations do not require LaTeX. Use a new output directory for a fresh computation. `--resume` reuses only matching heavy checkpoints and reruns the short tests.

## What is in this repository

The [mathematical supplement](supplement/V52_mathematical_supplement.pdf) and its [standalone source](supplement/source/mathematical_catalogue.tex) contain the reviewed English corrections. The document has no version designation; its existing PDF filename is retained for link compatibility. The coefficient data and equation numbering are unchanged. The bibliography is embedded in the source; see the [compilation instructions](supplement/README.md). [scripts](scripts/) contains the executable themed programs; [lib](lib/) contains the shared exact algebra. [parameters](parameters/) distinguishes live JSON configurations from descriptions of embedded constants.

[reference_results](reference_results/) separates inherited V51/V52 outputs from newly produced V52-R1 results. [validation](validation/) contains the actual commands, return codes and preparation checks. [historical archive](historical/) is a **non-executable, sanitized historical archive**, not another working pipeline. Original unredacted historical files are retained only in the author's local work archive.

The reader's entry point is [SCRIPT_REFERENCE_MAP.md](SCRIPT_REFERENCE_MAP.md). It links each program/function/test to V52 labels and freshly compiled printed numbers. Machine-readable maps and the complete label index are in [publication](publication/). The namespaces `article:` and `supplement:` are always distinct.

## Scope and limits

The fixed non-center-of-mass computation after article Eq. (356) is separate from the printed center-of-mass catalogue. Its 81 components, grouped frequency structures and finite-window witnesses are saved, with the nonzero count determined after computation. The radial reconstruction separately produces 81 bare contacts, 54 pair sources and the mixed scalar/tensor angular coefficient dictionaries at j=0,2,3,4. See [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) for the actual achieved scope and [VALIDATION_SUMMARY.md](VALIDATION_SUMMARY.md) for recorded results.

These are **finite-time, state- and regulator-dependent coefficients**, not a proof of a regulator-free unitary S matrix. Canonical mixed-field normalization identities are not a calculation of arbitrary-background mixed scattering eigenvalues. No numerical example substitutes for the analytical existence or boundedness arguments in the article.

[Reproduction details](REPRODUCIBILITY.md) · [Changes and provenance](CHANGELOG.md) · [Author publication procedure](PUBLICATION_GUIDE_RU.md) · [Citation metadata](CITATION.cff)
