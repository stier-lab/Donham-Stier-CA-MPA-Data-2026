# Data Repo Readiness Audit

Date: 2026-09-03

Scope: public data-processing repository for PISCO, KFM/MBON, LTER, Landsat,
and MPA source-data harmonization.

## Status

- Public GitHub repo verified: `stier-lab/Donham-Stier-CA-MPA-Data-2026`.
- GitHub repo is archived/read-only as of 2026-09-03; this local cleanup commit
  can be pushed only after unarchiving the repo or choosing a new remote.
- Current role: source-data provenance and harmonized CSV production.
- Companion public analysis repo:
  `stier-lab/sbc-2026-donham-kelp-mpa-cascade`.

## Cleanup Decisions

- Updated README language from the earlier Conservation Letters framing to the
  current Journal of Applied Ecology manuscript framing.
- Removed `firebase-debug.log` from the tracked public tree; a local archived
  copy is retained under `local_archive/2026-09-03_public_repo_cleanup/`.
- Moved ignored `.DS_Store` scratch state into the same local archive and added
  ignore rules for local archives and future Firebase debug logs.

## Boundaries

- Do not track raw monitoring data, local symlinks, manuscript drafts,
  citation-source PDFs, cover letters, or generated manuscript exports here.
- Keep source-data and harmonization documentation here.
- Keep statistical analysis and public result summaries in the analysis repo.
