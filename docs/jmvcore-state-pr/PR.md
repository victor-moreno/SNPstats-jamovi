# PR: jmvcore — fix `State` so it can be constructed

**Repository:** https://github.com/jamovi/jamovi
**Branch:** `fix/state-cannot-be-constructed`
**Base:** `main` (verified against `e56b7091`; `origin/main` has no later change
to `jmvcore/R/state.R` — the file's last touch is `06b955d1 Moved jmvcore and
compiler into tree`)
**Files:** `jmvcore/R/state.R`, `jmvcore/tests/testthat/test-state.R` (new)

---

## Summary

`jmvcore::State$new()` throws unconditionally:

```
Error in super$initialize(options, name, "", FALSE, clearWith) :
  argument "refs" is missing, with no default
```

`State$initialize` passes **five positional arguments** to
`ResultsElement$initialize`, which takes **six** and gives none of them a
default.

This makes **`type: State` unusable in any module**. The compiler accepts it —
`State` is in the results type enum (`compiler.js:287`) — and emits a call the
class cannot satisfy, so the analysis fails to construct, before `.init()` runs.
The module is dead rather than degraded, and nothing about the error points at
the real cause.

## Reproduction

Needs nothing but jmvcore (`reproduce.R` in this directory runs it):

```r
library(jmvcore)
State$new(options = Options$new(), name = "myState", clearWith = list("someOption"))
#> Error in super$initialize(options, name, "", FALSE, clearWith) :
#>   argument "refs" is missing, with no default
```

That is exactly the call `jamovi-compiler` generates for

```yaml
- name: myState
  type: State
  clearWith:
      - someOption
```

`compiler.js` `resultsify()` writes `options=options` and then only the
properties the `.r.yaml` declares (`compiler.js:407-417`), so `title`, `visible`
and `refs` are never supplied.

## Root cause

`jmvcore/R/state.R`:

```r
initialize=function(options, name, title, visible, clearWith) {
    super$initialize(options, name, '', FALSE, clearWith)
}
```

`jmvcore/R/results.R`:

```r
initialize=function(options, name, title, visible, clearWith, refs) { ... }
```

Three problems in the one method, the second masking the third:

1. **`refs` is neither accepted nor forwarded.** Every sibling takes
   `refs=character()` and passes `refs=refs`.
2. **The parent is called positionally.** Every sibling calls it by name — which
   is what has kept them working as the parent signature grew.
3. **No parameter has a default**, so the compiler's three-argument call would
   fail on `title`/`visible` even with (1) and (2) fixed. Currently masked: the
   body ignores both (passing the literals `''` and `FALSE`), so lazy evaluation
   never forces them and `refs` is the only error seen.

`State` is the only `ResultsElement` subclass with either problem:

| class | `super$initialize` call |
|---|---|
| Array, Group, Html, Image, Notice, Preformatted, Table | by name, all six arguments |
| Action, Output | by name, own signature |
| **State** | **`super$initialize(options, name, '', FALSE, clearWith)`** |

## The change

Brought into line with its siblings — defaults on every parameter, `refs`
accepted and forwarded, parent called by name. `title` and `visible` stay
hardcoded in the call, since a `State` carries saved analysis data and is never
shown:

```r
initialize=function(
    options,
    name='',
    title='',
    visible=FALSE,
    clearWith='*',
    refs=character()) {

    super$initialize(
        options=options,
        name=name,
        title='',
        visible=FALSE,
        clearWith=clearWith,
        refs=refs)
}
```

## Tests

`jmvcore/tests/testthat/test-state.R` (new) covers the compiler's call shape,
`refs` round-tripping, the untitled/not-visible guarantee, and
`setState()`/`$state`.

All five assertions fail before this change:

```
── State can be constructed the way the compiler constructs it ──
   argument "refs" is missing, with no default
── State accepts and forwards refs, like every other results element ──
   unused argument (refs = "someref")
```

The rest of the suite is unaffected. Running `test_dir()` on
`jmvcore/tests/testthat` before and after:

| | failed | errors | passed |
|---|---|---|---|
| before | 2 | 8 | 174 |
| after | **0** | **5** | **181** |

The 5 remaining errors are pre-existing and unrelated — `test-parseaddress.R`,
`test-table.R` and `test-utils.R` call unexported internals, which resolve under
`test_check('jmvcore')` but not under a bare `test_dir()`. They are identical
before and after.

## Why it matters

Found while looking for a way to avoid recomputing an expensive fit on every
option click, in a module whose haplotype × covariate interaction calls
`haplo.stats::haplo.glm` — about 1.6 s for a four-SNP model on 2,838 subjects.
Its result rows cannot be pre-created in `.init()` (the haplotype set is not
known until the EM has run), so `rowCount` is always 0 after a restore and the
fit re-runs on *every* option click, including ones that cannot affect it.

Caching the fitted rows — roughly 5 KB of text — in a `State` element whose
`clearWith` lists the options the fit depends on is exactly the right shape for
this: the payload persists across clicks and is invalidated automatically by the
same declaration mechanism the tables already use. With the fix applied locally,
that works precisely as documented; state survives the save/restore cycle and
clears when `clearWith` fires.

For completeness, the alternatives were measured:

| carrier | state survives restore? |
|---|---|
| `Image` | yes |
| `Html` | inconsistent — survived on a top-level element, not on a nested one |
| `Table` | no |
| `State` | cannot be constructed |

Only `Image` is reliable today, and using a hidden image purely as a state box is
a hack we would rather not ship.

## Notes

- No API break: every existing call to `State$new()` is a strict subset of the
  new signature, and the class's observable behaviour (untitled, not visible) is
  unchanged.
- Happy to adjust the test style or split the fix and tests if you would prefer.
