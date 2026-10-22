# Manuscript Changelog

## Final manuscript revisions — 22 October 2026

- Updated the manuscript’s associated software-release DOI to the v1.2.0 archive: https://doi.org/10.5281/zenodo.23566896.
- Updated the manuscript package’s release-related metadata to identify software version v1.2.0.
- Enlarged Tables 2 and 3 in the Word manuscript to improve readability and make better use of the available text width.
- Regenerated the PDF from the revised Word manuscript and checked the relevant table pages and title page.
- Retained the manuscript’s existing scientific results, numerical comparisons, references, and figure content unless otherwise specified in the revision history below.

## Layout and reference-style revision — 9 October 2026

- Reflowed text, figures, and tables in a clean, single-column academic layout to reduce avoidable page-bottom whitespace.
- Standardized references to unnumbered, alphabetized author–year entries, consistent with in-text author–year citations.
- Confirmed that Harfst et al. appears alphabetically between Hairer and Harris.
- Preserved the requested author block, clickable ORCID, and PDF author metadata.

## Pre-release corrections — 9 October 2026

- Regenerated the previously blank Figure 8 from the saved long-duration eccentric-orbit result and inserted the plot into the manuscript. Visually checked the two panels and the dip near orbit 724.
- Corrected Section 5.2: Yoshida4 has the smallest energy error; Verlet and Yoshida4 are both at machine precision for angular momentum, consistent with Table 2.
- Updated the SciPy bibliography entry to Virtanen, Gommers, Oliphant et al. (2020), preserving the in-text Virtanen et al. citation. Added Hairer, Nørsett, and Wanner (1993) at the first DOP853 mention.
- Shortened the Figure 3 and Figure 7 captions to keep them with their figures, moved the Figure 7 explanatory paragraph into the main text, and reduced the size of Figure 1.
- Visually checked equations and scientific-notation exponents on pages 2–4 and in Tables 2–3.
- Removed “Reproducible” from the title until a public archive and DOI were available.
- At this stage, public repository publication and DOI assignment were still pending.

## Version 1.0 — 10 October 2026

- Bolded the manuscript title and added a version/date line beneath the author block.
- Updated the abstract to report the Table 3 phase-error ranking across all three tested resolutions.
- Clarified the substantial velocity-Verlet position/phase separation in the 1,000-orbit experiment.
- Stated that the RK4 phase-error sweep does not establish clean fourth-order asymptotic scaling.
- Changed the in-text IAU resolution citation to (IAU, 2012).
- Replaced the archive to-do sentence with a statement that the code was not publicly released when this manuscript version was prepared.

## Scientific and reproducibility revisions

- Replaced the mixed one-year/five-year comparison in Section 5.2 with Table 2, comparing all four methods over the same nominal five-year interval using a two-day step: 913 steps over 1,826 days.
- Reran the eccentric Kepler problem with \(e = 0.9\) for 1,000 periods at 5,000, 10,000, and 20,000 steps per orbit.
- Used 200 diagnostic samples per orbit in the long-duration experiment.
- Added Table 3 with the nine method/resolution combinations and regenerated Figure 8.
- Explained the sharp energy-error dip near orbit 724 as a near-zero crossing in the sampled, phase-dependent energy-error curve, rather than an integration blow-up.
- Added a common five-year comparison script, a long-duration resolution-sweep script, updated result metadata, a DOI audit note, and `CITATION.cff` to the reproducibility materials.

## Author information and reference corrections

- Replaced the title-page author block with the requested information:
  - Daniel Durodola
  - Department of Physics, University of Lagos, Lagos, Nigeria
  - ORCID: 0009-0009-7171-4901 · danieldurodola02@gmail.com
- Made the ORCID identifier clickable, linking to https://orcid.org/0009-0009-7171-4901.
- Set the PDF author metadata to `Daniel Durodola` and saved the manuscript under the filename `Durodola_NBody_Benchmark.pdf`.
- Corrected reference [8] (Harfst et al.) to the special-purpose-supercomputers paper, including its article title, journal, pages, and DOI: https://doi.org/10.1016/j.newast.2006.11.003.
- Linked the corrected DOI in the PDF.
- Audited the other DOI-bearing references; no further mismatch was identified. The IAU Resolution B2 citation remains an official institutional URL rather than a DOI.

## Software release and citation metadata — v1.2.0

- The software has been archived on Zenodo as version v1.2.0.
- Version-specific software DOI: https://doi.org/10.5281/zenodo.23566896.
- Earlier archived software releases remain available:
  - v1.1.0: https://doi.org/10.5281/zenodo.23480336
  - v1.0.0: https://doi.org/10.5281/zenodo.23286569
- Updated the release badge and DOI badge references for the v1.2.0 release.
- Updated `CITATION.cff` to identify version 1.2.0 and include the Zenodo DOI.
- The software DOI identifies the archived software release; it does not, by itself, establish that the associated manuscript has been peer-reviewed or published.

## Document layout

- Used available space on pages 4 and 6 for Tables 2 and 3 to avoid unnecessary page breaks.
- Kept the manuscript at 12 pages in the previously checked layout.
- Enlarged Tables 2 and 3 in the Word document and regenerated the PDF to reflect the revised table sizing.


