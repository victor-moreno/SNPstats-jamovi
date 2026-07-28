# jmvcore — fix `State` so it can be constructed

Problem found while looking for a way to avoid recomputing an expensive fit on every
option click. Haplotypes are estimated via EM algorithm and its result rows cannot be pre-created in `.init()` (the haplotype set is not known until the EM has run), so `rowCount` is always 0 after a restore and the
fit re-runs on *every* option click, including ones that cannot affect it.

Caching the fitted rows — roughly 5 KB of text — in a `State` element whose
`clearWith` lists the options the fit depends on is exactly the right shape for
this: the payload persists across clicks and is invalidated automatically by the
same declaration mechanism the tables already use. With the fix applied locally,
that works precisely as documented; state survives the save/restore cycle and
clears when `clearWith` fires.

The problem is that **no call to `jmvcore::State$new()` succeeds.**
`State$initialize` passes **five positional arguments** to
`ResultsElement$initialize`, which takes **six** and gives none of them a
default — and `State` itself gives none of its own parameters a default either,
so the failure is a cascade:

```r
State$new()                                        # argument "options" is missing
State$new(options=o)                               # argument "name" is missing
State$new(options=o, name="s")                     # argument "clearWith" is missing
State$new(options=o, name="s", clearWith=list())   # argument "refs" is missing   <-- the compiler's call
State$new(options=o, name="s", refs="r")           # unused argument (refs = "r")
```

The fourth line is the one that matters in practice. This makes **`type: State`
unusable in any module**: the compiler accepts it — `State` is in the results
type enum (`compiler.js:287`) — and `resultsify()` emits `options=options`
followed by only the properties the `.r.yaml` declares (`compiler.js:407-417`),
so `title`, `visible` and `refs` are never supplied. The analysis then fails to
construct, before `.init()` runs.

There are 3 problems:

1. **`refs` is neither accepted nor forwarded.** Every sibling takes
   `refs=character()` and passes `refs=refs`.
2. **The parent is called positionally.** Every sibling calls it by name — which
   is what has kept them working as the parent signature grew.
3. **No parameter has a default**, which is what produces the cascade above.
   `title` and `visible` are the exception: the body ignores them, passing the
   literals `''` and `FALSE`, so R's lazy evaluation never forces them and they
   never raise — they simply cannot be supplied usefully either.

## Fix

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

## Is `options` still required after the fix?

Yes, deliberately. `options` is the one genuinely required argument — the parent
stores it and everything downstream, including `clearWith` evaluation and
translation, depends on it — so it keeps no default and `State$new()` still
raises `argument "options" is missing`. That is **unchanged** by this patch: the
unpatched class raised exactly the same error for the same call.

It also matches every sibling. `Html`, `Preformatted`, `Image`, `Notice`,
`Group`, `Array`, `Action` and `Output` all leave `options` without a default and
all raise on a bare `$new()`. The single exception is `Table`, which explicitly
opts in with

```r
if (missing(options))
    options <- Options$new()
```

so that `Table$new()` works standalone in its own tests. If you would like
`State` to do the same for symmetry with `Table`, that is a one-line addition —
but it seemed out of scope for a bug fix, and the compiler always supplies
`options` anyway.

## Tests

`jmvcore/tests/testthat/test-state.R` (new) covers the compiler's call shape,
`refs` round-tripping, the untitled/not-visible guarantee, and
`setState()`/`$state`. All five assertions fail before this change.

The rest of the suite is unaffected — `test_dir()` on `jmvcore/tests/testthat`
before and after:

| | failed | errors | passed |
|---|---|---|---|
| before | 2 | 8 | 174 |
| after | **0** | **5** | **181** |

The 5 remaining errors are pre-existing and unrelated: `test-parseaddress.R`,
`test-table.R` and `test-utils.R` call unexported internals, which resolve under
`test_check('jmvcore')` but not under a bare `test_dir()`. Identical before and
after.