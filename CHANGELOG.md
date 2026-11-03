# Changelog

Weekly releases are tagged `vYYYY.MM.DD` by `.github/workflows/weekly-release.yml`.

## Unreleased

- New name and site: Supremum. The landing page is `docs/index.html`; the generated entry
  table moves to `docs/browse.html` (`tools/export.py` now writes `browse.html`).
- priority_claims: 8 verdicts re-checked against primary sources (2026-11-03). Third-party
  pre-emption is now 9 of 231 claims (3.9%); equal-or-better before the claim, 43 of 231.

## v2026.10.13

- Provenance, fidelity and feasibility fields on every entry (SCHEMA.md).
- Ingest the Erdős AI-contributions wiki and the DS1 Ramsey tables.
- Ingest the Packing Center tables and Cohn's kissing-number and sphere-packing tables, with
  cross-links to existing entries.
- Validator-verified leaderboard records mirrored as current (CHANGES_20261013.md).
- Weekly self-versioning workflow, CI on push and pull request (Python 3.12, 3.13).
- JOSS readiness: LICENSE, pyproject.toml, install and support sections, issue templates,
  CITATION.cff licence split.

## v0.1.0 - 2026-10-12

- Initial public release: audited entries, validator, exporter, badges, proposal tool, tests.
