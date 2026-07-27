# Bug report: `jmvcore::State` cannot be constructed — `type: State` is unusable

**Package:** jmvcore 2.7.38
**R:** 4.5.3 (2026-03-11), macOS arm64
**Severity:** `type: State` is unusable in any module; the analysis fails to
construct, so the whole analysis is dead rather than degraded.
**Reported by:** Victor Moreno (Catalan Institute of Oncology) — SNPstats module

---

## Summary

`jmvcore::State$new()` throws unconditionally:

```
Error in super$initialize(options, name, "", FALSE, clearWith) :
  argument "refs" is missing, with no default
```

`State$initialize` calls the parent with **five positional arguments**, but
`ResultsElement$initialize` takes **six** and none of them has a default. It is
the only `ResultsElement` subclass in jmvcore that calls the parent positionally,
and the only one that does not accept and forward `refs`.

The practical consequence is that declaring `type: State` in a `.r.yaml`
compiles cleanly — `State` is in the compiler's accepted type enum — and then
the module fails at runtime, on construction, before `.init()` runs.

---

## Reproduction

No module needed:

```r
library(jmvcore)
State$new(options = Options$new(), name = "myState", clearWith = list("someOption"))
#> Error in super$initialize(options, name, "", FALSE, clearWith) :
#>   argument "refs" is missing, with no default
```

That call is exactly what the jamovi UI compiler generates for

```yaml
- name: myState
  type: State
  clearWith:
    - someOption
```

namely, in the generated `*.h.R`:

```r
self$add(jmvcore::State$new(
    options=options,
    name="myState",
    clearWith=list("someOption")))
```

so any module with a `type: State` element reproduces it. In our case the
analysis object could not be created at all:

```
snpStatsClass$new(options = ..., data = ...)
#> Error: argument "refs" is missing, with no default
```

---

## Root cause

`jmvcore::State$initialize` (2.7.38):

```r
function (options, name, title, visible, clearWith)
{
    super$initialize(options, name, "", FALSE, clearWith)
}
```

`jmvcore::ResultsElement$initialize`:

```r
function (options, name, title, visible, clearWith, refs)
```

with **no defaults on any of the six**. The call supplies five positionally, so
`refs` is missing and the error is raised as soon as `ResultsElement$initialize`
touches it.

There are three separate problems in that one method:

1. **`refs` is neither accepted nor forwarded.** Every sibling class takes
   `refs = character()` and passes `refs = refs`.
2. **The parent is called positionally.** Every sibling calls it by name, which
   is what has kept them working as the parent signature grew.
3. **No parameter has a default**, so the compiler-generated three-argument call
   (`options`, `name`, `clearWith`) would fail on `title`/`visible` even if (1)
   and (2) were fixed. This one is currently masked: the body ignores `title`
   and `visible` — passing the literals `""` and `FALSE` — so R's lazy
   evaluation never forces them, and `refs` is the only error you see.

### Every `ResultsElement` subclass, for contrast

Surveyed programmatically in 2.7.38:

| class | `super$initialize` call |
|---|---|
| Array | `super$initialize(options = options, name = name, title = title, visible = visible, clearWith = clearWith, refs = refs)` |
| Group | *(same, by name, 6 args)* |
| Html | *(same)* |
| Image | *(same)* |
| Notice | *(same)* |
| Preformatted | *(same)* |
| Table | *(same)* |
| Action | by name, own signature |
| Output | by name, own signature |
| **State** | **`super$initialize(options, name, "", FALSE, clearWith)`** |

`State` is the outlier on both counts.

For direct comparison, `Html$initialize`:

```r
function (options, name = "", title = "", visible = TRUE, clearWith = "*",
          refs = character(), content = "")
{
    super$initialize(options = options, name = name, title = title,
                     visible = visible, clearWith = clearWith, refs = refs)
    ...
}
```

---

## Proposed fix

Bring `State` into line with its siblings — accept `refs`, give the parameters
defaults, and call the parent by name. `State` is meant to be invisible and
untitled, so keep those hardcoded in the call rather than in the signature:

```r
initialize = function(options, name = "", title = "", visible = FALSE,
                      clearWith = "*", refs = character()) {
    super$initialize(options = options, name = name, title = "",
                     visible = FALSE, clearWith = clearWith, refs = refs)
}
```

A minimal one-line alternative, if you prefer to change as little as possible,
is to add `refs` and forward it:

```r
initialize = function(options, name, title, visible, clearWith, refs = character()) {
    super$initialize(options, name, "", FALSE, clearWith, refs)
}
```

This also works — `title` and `visible` are never forced — but it leaves the
class relying on lazy evaluation to tolerate the arguments the compiler does not
supply, so we would suggest the first form.

---

## Verification

Both the failure and the fix were checked end to end against a real module
(SNPstats), by replacing the method in the loaded namespace:

```r
St <- get("State", envir = asNamespace("jmvcore"))
St$set("public", "initialize",
       function(options, name = "", title = "", visible = FALSE,
                clearWith = "*", refs = character()) {
         super$initialize(options = options, name = name, title = "",
                          visible = FALSE, clearWith = clearWith, refs = refs)
       }, overwrite = TRUE)
```

Result:

```
=== 1. UNPATCHED jmvcore 2.7.38 ===
  analysis with a `type: State` element -> ERROR: argument "refs" is missing, with no default

=== 2. APPLYING THE PROPOSED FIX ===
  State$initialize replaced

=== 3. PATCHED ===
  analysis constructed OK
  state set ( 16 rows ) and saved
  after unrelated click + restore: state survived = TRUE  rows = 16
  after a clearWith-triggering change: state cleared = TRUE
```

So with the fix, `State` behaves exactly as documented: the payload survives the
save/restore cycle, and it is invalidated automatically when the element's
`clearWith` fires.

---

## Why this matters to us

We hit this looking for a way to avoid recomputing an expensive fit on every
option click. Our haplotype × covariate interaction analysis calls
`haplo.stats::haplo.glm`, which for a four-SNP model on 2,838 subjects costs
about **1.6 s**. Its result rows cannot be pre-created in `.init()` — the
haplotype set is not known until the EM has run — so `rowCount` is always 0
after a restore and the fit re-runs on *every* option click, including ones that
cannot affect it.

The natural fix is to cache the fitted rows (about **5 KB** as text) in a
`State` element whose `clearWith` lists the options the fit actually depends on:
the payload then persists across clicks and is invalidated automatically by the
same declaration mechanism the tables already use. That is precisely what
`type: State` appears designed for, and it is currently unreachable.

We checked the alternatives, for completeness:

| carrier | state survives restore? |
|---|---|
| `Image` | yes |
| `Html` | inconsistent — survived on a top-level element, not on a nested one |
| `Table` | **no** |
| `State` | **cannot be constructed** |

Only `Image` is reliable, and using a hidden image purely as a state box is a
hack we would rather not ship.

Happy to test a patched build if that is useful.

---

## Attachment

`docs/jmvcore-State-bug-verify.R` in this repository is the script that produced
the output above. It needs the SNPstats module installed, but the two-line
reproduction at the top of this report needs nothing but jmvcore.
