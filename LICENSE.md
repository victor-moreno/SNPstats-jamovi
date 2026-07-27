# License

SNPstats is free software: you can redistribute it and/or modify it under the
terms of the **GNU General Public License version 3** as published by the Free
Software Foundation. The full text is in [`COPYING`](COPYING) and at
<https://www.gnu.org/licenses/gpl-3.0.html>.

Copyright © 2026 Victor Moreno (Catalan Institute of Oncology)

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Why GPL and not a permissive license

SNPstats is a jamovi module, and that decides the question. It is built on
`jmvcore` (GPL >= 2) — every analysis class inherits from `jmvcore::Analysis`,
so there is no jamovi module that does not link to it — and on `haplo.stats`
(GPL >= 2), which provides the haplotype EM and the haplotype GLM. A work that
requires GPL libraries to function is distributed under the GPL.

Removing the `genetics` dependency (v1.0.0) does not change this. That was done
because the package is marked obsolete upstream and a module whose dependency
leaves the pinned jamovi snapshot fails to install rather than degrading — an
availability fix, not a licensing one. The replacements in `R/snp_genetics.R`
were written from the published statistical definitions, not ported from
`genetics`, and are covered by this license along with the rest of the module.
