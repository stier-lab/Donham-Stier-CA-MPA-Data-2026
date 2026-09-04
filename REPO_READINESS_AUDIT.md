# Data Repo Readiness Audit

Date: 2026-09-03

Scope: public data-processing repository for PISCO, KFM/MBON, LTER, Landsat,
and MPA source-data harmonization.

## Status

- Public GitHub repo verified: `stier-lab/Donham-Stier-CA-MPA-Data-2026`.
- GitHub repo is archived/read-only after the 2026-09-03 cleanup.
- Current role: source-data provenance and harmonized CSV production.
- Companion public analysis repo:
  `stier-lab/sbc-2026-donham-kelp-mpa-cascade`.

## Cleanup Decisions

- Updated README language from the earlier Conservation Letters framing to the
  current Journal of Applied Ecology manuscript framing.
- Removed tracked `drive-recovery/` files from current HEAD because they were
  duplicate recovered raw/local materials that contradicted the repo boundary
  that raw monitoring data stay out of git. Local copies are preserved in
  ignored `local_archive/2026-09-03_data_repo_harsh_cleanup/drive-recovery/`.
- Removed tracked historical pipeline logs from current HEAD. Local copies are
  preserved in ignored `local_archive/2026-09-03_data_repo_harsh_cleanup/logs/`;
  future `logs/*.log` outputs remain ignored.
- Removed `firebase-debug.log` from the tracked public tree; a local archived
  copy is retained under `local_archive/2026-09-03_public_repo_cleanup/`.
- Moved ignored `.DS_Store` scratch state into the same local archive and added
  ignore rules for local archives and future Firebase debug logs.
- Replaced placeholder Dryad DOI and private Google Drive setup language with a
  provider/project-data-steward input setup note. Also softened generated-output
  documentation so it does not claim a final Dryad archive before submission.

## Boundaries

- Do not track raw monitoring data, local symlinks, manuscript drafts,
  citation-source PDFs, cover letters, or generated manuscript exports here.
- Keep source-data and harmonization documentation here.
- Keep statistical analysis and public result summaries in the analysis repo.
