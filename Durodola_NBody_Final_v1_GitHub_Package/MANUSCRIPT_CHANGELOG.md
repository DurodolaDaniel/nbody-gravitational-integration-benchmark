# Final manuscript change log — 2026-10-09

1. **Author block (title page):** replaced with exactly the three requested lines:
   - Daniel Durodola
   - Department of Physics, University of Lagos, Lagos, Nigeria
   - ORCID: 0009-0009-7171-4901 · danieldurodola02@gmail.com
2. **ORCID link:** made `0009-0009-7171-4901` clickable to `https://orcid.org/0009-0009-7171-4901`.
3. **PDF metadata and filename:** set PDF author metadata to `Daniel Durodola`; saved final manuscript as `Durodola_NBody_Benchmark.pdf`. The previous author-name ordering did not occur elsewhere in the manuscript text, running headers, or footers.
4. **Reference [8] (Harfst et al.):** corrected the article title, journal, pages, and DOI to the special-purpose-supercomputers paper, DOI `10.1016/j.newast.2006.11.003`; linked the DOI in the PDF. The previous DOI `10.1016/j.parco.2007.01.001` belonged to a different paper.
5. **DOI audit:** checked the other DOI-bearing references; no further mismatch was identified. The IAU Resolution B2 citation remains an official institutional URL rather than a DOI.
6. **Section 5.2:** removed the mixed one-year/five-year comparison and added Table 2 with all four methods evaluated over the same nominal five-year interval at a two-day step (913 steps; 1,826 days).
7. **Section 5.4 and Figure 8:** reran the e = 0.9 Kepler problem for 1,000 periods at 5,000, 10,000, and 20,000 steps per orbit. Diagnostics use 200 samples per orbit. Added Table 3 with the nine method/resolution combinations, regenerated Figure 8, and explained the sharp energy-error dip near 724 periods as a near-zero crossing in the sampled, phase-dependent energy-error curve—not an integration blow-up.
8. **Whitespace/layout:** used existing blank space on pages 4 and 6 for the two new tables, keeping the manuscript at 12 pages and avoiding unnecessary page breaks.
9. **Repository reproducibility materials:** added a script for the common five-year comparison, a long-duration resolution-sweep script, updated result metadata, a DOI audit note, and `CITATION.cff` with the requested ORCID.

**External release status:** this package has not been uploaded to GitHub and has no Zenodo DOI. No public URL or DOI has been invented. An authenticated publishing step is still required before describing the repository as publicly archived.


## Layout and reference-style revision — 9 October 2026

- Reflowed text, figures, and tables in a clean single-column academic layout to reduce avoidable page-bottom whitespace.
- Standardized references to unnumbered, alphabetized author–year entries, consistent with in-text author–year citations.
- Confirmed Harfst et al. appears alphabetically between Hairer and Harris.
- Preserved the requested author block, clickable ORCID, and PDF author metadata.


## Pre-release corrections — 9 October 2026

- Regenerated the previously blank Figure 8 from the saved long-duration eccentric-orbit result and inserted the plot into the manuscript; visually verified the two panels and the dip near orbit 724.
- Corrected §5.2: Yoshida4 has the smallest energy error; Verlet and Yoshida4 are both at machine precision for angular momentum, consistent with Table 2.
- Updated the SciPy bibliography entry to Virtanen, Gommers, Oliphant et al. (2020), preserving the in-text Virtanen et al. citation. Added Hairer, Nørsett, and Wanner (1993) at the first DOP853 mention.
- Shortened Figure 3 and Figure 7 captions to keep them with their figures; moved the Figure 7 explanatory paragraph into the body. Reduced Figure 1 size.
- Visually checked equations and scientific-notation exponents on pages 2–4 and in Tables 2–3. Removed “Reproducible” from the title until a public archive and DOI actually exist.
- **Release blocker:** public GitHub upload and Zenodo DOI are still pending because authenticated account publishing is not available in this environment. No public URL or DOI has been fabricated.

## Version 1.0 — 10 October 2026

- Bolded the manuscript title and added a version/date line under the author block.
- Updated the abstract to report the Table 3 phase-error ranking across all three tested resolutions.
- Clarified the severe velocity-Verlet position/phase separation in the 1,000-orbit experiment.
- Stated that the RK4 phase-error sweep does not establish clean fourth-order asymptotic scaling.
- Changed the in-text IAU resolution citation to (IAU, 2012).
- Replaced the archive to-do sentence with a plain statement that the code was not publicly released when this manuscript version was prepared.
