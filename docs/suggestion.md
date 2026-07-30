# Suggestion for jmvcore: let a table restore its own shape

**Status:** proposal, not implemented here. Nothing in SNPstats depends on it.
**Context:** jmvcore 2.7.38, R 4.5.3, module SNPstats 1.0.0.
**Related:** `docs/jmvcore-addrow-pr/` (the quadratic `addRow`) — complementary,
not the same issue.

---

## Summary

A results table whose **shape** (row count, row keys, dynamically-added columns)
is only knowable *after* the analysis has computed can never be restored from
the saved state, so it is rebuilt on every option click — including clicks on
options that are not in its `clearWith`.

The information needed to restore it is already in the state file. jmvcore reads
it, then discards it. The fix is small and can be made opt-in.

## The symptom, from a module author's side

Some tables cannot predict their shape in `.init()`. In SNPstats the haplotype
tables are the clean example: their rows are the haplotypes found by
`haplo.stats::haplo.em`, and the interaction tables additionally add one column
per level of the interaction covariate. Neither is derivable from the options —
only from the data, after the EM has run.

The consequence is that `Table$isFilled()` is never TRUE after a restore for
these tables, so the module must refill them on every click. To keep that
affordable, the module caches the *computed products* of the fits in
`Table$state` (which **is** restored) and rewrites the rows from cache.

Measured on `data/CRCgenet-SNPs.tsv`, 8 SNPs + 1 covariate, all five haplotype
tables enabled (29 rows each):

```
first run:        34.7 s
unrelated click:   0.293 s total   (haplo.em calls: 0, haplo.glm calls: 0)
  snp_prepare      0.052 s
  everything else  0.241 s   <- the rebuild this proposal would remove
```

So the workaround does work, and the residual cost is small *here*. What it
costs is not primarily time:

- every such table needs a bespoke caching layer, and the cached state must stay
  under jamovi's 500 kB limit or it silently stops persisting (in SNPstats,
  `vcov(haplo.glm)` had to be trimmed to its coefficient block: 250 kB → 28 kB);
- dynamically-added columns must be re-added on every run as well;
- tables that *can* predict their shape must do so with a parallel
  implementation in `.init` (`.pre_rows`, row-count predictors) that has to stay
  in step with what `.run` actually produces;
- and the rebuild is quadratic in row count, because of `addRow` (see the other
  write-up).

## Root cause

`Table$fromProtoBuf` restores by iterating **the `.init` skeleton** and looking
each row up in the saved protobuf:

```r
for (i in seq_along(private$.rowNames)) {          # the skeleton's rows
    rowName <- private$.rowNames[[i]]
    fromRowIndex <- indexOf(rowName, tablePB$rowNames)
    if (!is.na(fromRowIndex)) {
        for (j in seq_along(private$.columns)) {   # the skeleton's columns
            toCol <- private$.columns[[j]]
            fromColIndex <- columnPBIndicesByName[[toCol$name]]
            if (!is.null(fromColIndex))
                toCol$getCell(i)$fromProtoBuf(cells[[fromColIndex]][[fromRowIndex]])
        }
    }
}
```

Restore is a **left join from the skeleton onto the saved results**. Saved rows
with no skeleton counterpart are dropped; so are saved columns. Yet
`tablePB$rowNames` and every saved column with its cells — the complete shape —
are in hand at that moment.

Two things make this worth changing rather than working around:

1. **The guards above that loop have already decided the shape is still valid.**
   `fromProtoBuf` returns early if `clearWith` fired:

   ```r
   if (someChanges && identical("*", private$.clearWith)) return()
   if (any(oChanges %in% private$.clearWith))            return()
   for (clearName in private$.clearWith)
       if (any(vChanges %in% private$.options$option(clearName)$vars)) return()
   ```

   Reaching the copy loop *means* "nothing this table depends on changed" — and
   the shape is discarded anyway.

2. **`Table$state` is already restored under exactly those guards** (it is
   `super$fromProtoBuf`, called after them). So a module may persist its
   computed intermediates across a click but not its displayed output. That
   inconsistency is the whole defect, and it is what forces the caching
   workaround: authors reach for `state` precisely because it is the only thing
   that survives.

There is no way for the module to close this gap itself. `.init` runs *before*
the restore, so it cannot consult the saved shape to decide what to pre-create:

```
after init():    rows: 0   state: NULL      <- .init sees nothing
after .load():   rows: 0   state: present   <- saved state appears only here
after run():     rows: 29
```

`.init` also cannot know whether the work is needed at all — only whether the
option is ticked, which stays true for every subsequent click.

## Proposal

`.r.yaml` already accepts `rows: 0` (a literal count) and `rows: (snps)` (an
option expression). Add a third form for "not derivable from options at all":

```yaml
- name: haploFreqTable
  type: Table
  rows: dynamic
  clearWith: [snps, response, haploFreqMin, haploMinCount, completeCases]
```

and in `Table$fromProtoBuf`, once the `clearWith` guards have passed, adopt the
saved shape instead of intersecting with it:

```r
# after the guards, before the existing copy loop
if (private$.rowsDynamic) {
    for (nm in tablePB$rowNames)
        if (is.na(indexOf(nm, private$.rowNames)))
            self$addRow(rowKey = <key for nm>)
    for (colPB in columnsPB)
        if (! colPB$name %in% names(private$.columns))
            self$addColumn(name  = colPB$name, title = colPB$title,
                           type  = colPB$type,  format = colPB$format)
}
```

The existing copy loop then fills them, unchanged.

### Why this layer

- **Opt-in**: no behaviour change for any existing module or table.
- **No new state**: it uses what the state file already carries; nothing extra is
  computed, stored or versioned.
- **It makes `clearWith` mean what it already claims**, and removes the
  `state`-vs-shape inconsistency described above.
- **It removes the workaround class, not just its cost**: no caching of computed
  products, no 500 kB budget, no row-count prediction in `.init` that can drift
  from what `.run` produces.

### The same flag should cover `Array`

`Array$fromProtoBuf` has the identical left join — it iterates
`private$.itemNames` and looks each up in `arrayPB$elements`, dropping saved
items the skeleton lacks. So arrays built with `addItem` (per-SNP result groups,
per-stratum LD groups) have the same defect and the same fix.

### The one implementation snag

`Table$asProtoBuf` stores only

```r
table$rowNames <- private$.rowNames
```

i.e. the JSON-serialised keys, not the keys themselves. Two options:

- **jmvcore-only**: parse them back. `.rowNames` is `sapply(.rowKeys, toJSON)`,
  so this round-trips for scalar keys and is enough for the common case.
- **Cleaner, but cross-repo**: add a `rowKeys` field to the protobuf schema.
  That is a coordinated jamovi + jmvcore change.

## Two smaller, independent asks

- **`deleteRows(rowNo)`.** Today `deleteRows()` takes no arguments — all or
  nothing. A module that over-predicts a row count therefore cannot trim the
  surplus, only rebuild wholesale. This also closes off the obvious
  "over-allocate rows in `.init`" workaround, because `Table$isFilled()` returns
  FALSE if *any* visible cell is `NULL`/`NA`, so a single blank trailing row
  re-triggers the very refill it was meant to avoid.
- **The quadratic `addRow`** (`docs/jmvcore-addrow-pr/`). Complementary: shape
  restore avoids the rebuild, the `addRow` fix makes the unavoidable *first*
  build cheap. 2016 rows: 41.3 s → 0.14 s.

## What does not work (checked, so it need not be re-derived)

- **A post-restore hook.** Adding rows after `.load` is too late by
  construction — restore is the thing that fills them. `postInit` is nearly that
  hook and still would not help; in 2.7.38 `Analysis$init()` sets the status to
  `"inited"`, so `Analysis$run()` never calls `postInit()` at all.
- **Running `.load` before `.init`.** `.init` creates the objects
  `fromProtoBuf` merges into, so the order cannot simply be flipped.
- **Exposing "is this a restore?" or the previous row count to `.init`.** It
  would work, but it inverts the responsibility: every module would reimplement
  shape bookkeeping that the framework already holds in the state file.
