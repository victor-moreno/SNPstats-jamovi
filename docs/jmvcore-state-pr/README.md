# jmvcore `State` fix — ready to submit

The fix is committed on branch **`fix/state-cannot-be-constructed`** inside the
`jamovi-src/` clone (`git@github.com:jamovi/jamovi.git`). It is *not* pushed —
this sandbox has no SSH access.

## Contents

| file | what it is |
|---|---|
| `PR.md` | the pull-request description, ready to paste |
| `reproduce.R` | standalone reproduction; needs only jmvcore, no module or data |
| `0001-jmvcore-fix-State.patch` | the commit as a patch, if you would rather apply it elsewhere |

## To submit

```bash
cd jamovi-src
git log --oneline -1          # fccf50c3 jmvcore: fix State so it can be constructed
git show --stat HEAD          # jmvcore/R/state.R + jmvcore/tests/testthat/test-state.R

# push to your fork and open the PR against jamovi/jamovi:main
git remote add fork git@github.com:<you>/jamovi.git
git push fork fix/state-cannot-be-constructed
```

Then open the PR with the body from `PR.md`.

**Before it can be merged, jamovi requires a Contributor License Agreement.**
See `jamovi-src/CONTRIBUTING.md` — there is an individual form and an entity form
(the latter if you are contributing on behalf of ICO).

## Verifying it yourself

```bash
# build the patched jmvcore into a scratch library
R_LIBS=".Rlib-arm" R CMD INSTALL --library=<scratch> jamovi-src/jmvcore
R_LIBS="<scratch>:.Rlib-arm" Rscript docs/jmvcore-state-pr/reproduce.R
R_LIBS="<scratch>:.Rlib-arm" Rscript -e \
  'library(testthat); library(jmvcore); test_file("jamovi-src/jmvcore/tests/testthat/test-state.R")'
```

`git -C jamovi-src stash push -- jmvcore/R/state.R` reverts to the buggy version
if you want to see the failure first; `git -C jamovi-src stash pop` restores it.

## Note on the working tree

`jamovi-src/docker/jamovi-Dockerfile` was already modified before this work and
is deliberately **not** part of the commit.
