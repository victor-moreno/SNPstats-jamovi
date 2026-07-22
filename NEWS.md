History:
- 260420 v0.2.0 First public release
- 260501 v0.3.0 Added PGS submodule
- 260508 v0.4.0 Combined submenus for SNPstats and categorical response
- 260717 v0.5.0 Refactored to eliminate table refresh

Issues:
- snpPGS gets argument caseLevel without default. Function fails in tests unless snpPGS.h.R is patched and package reinstalled. No problem for jamovi

Plan:
- parallel speed-up
- remove genetics denepdency because package is marked obsolete
