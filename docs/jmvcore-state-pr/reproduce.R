# Minimal reproduction of the jmvcore::State constructor bug.
# Needs nothing but jmvcore — no module, no jamovi, no data.
#
#   Rscript reproduce.R
#
# Before the fix both calls fail; after it, both succeed.

library(jmvcore)

cat("jmvcore:", as.character(packageVersion("jmvcore")), "\n")
cat("R:      ", R.version.string, "\n\n")

# ── 1. Exactly what jamovi-compiler emits for a `type: State` element ────────
#
#     - name: myState
#       type: State
#       clearWith:
#           - someOption
#
# compiler.js resultsify() writes `options=options` followed by only the
# properties the .r.yaml declares, so title / visible / refs are never supplied.
cat("1. call shapes\n   BEFORE the fix every one of these fails, in a cascade.\n   AFTER it, only the bare call does - `options` is required by design.\n")
o <- Options$new()
shapes <- list(
  list('State$new()',                                       function() State$new()),
  list('State$new(options=o)',                              function() State$new(options = o)),
  list('State$new(options=o, name="s")',                    function() State$new(options = o, name = "s")),
  list('State$new(options=o, name="s", clearWith=list())',  function() State$new(options = o, name = "s", clearWith = list())),
  list('State$new(options=o, name="s", refs="r")',          function() State$new(options = o, name = "s", refs = "r")))
for (sh in shapes) {
  res <- tryCatch({ sh[[2]](); "OK" }, error = function(e) paste("FAILS:", conditionMessage(e)))
  cat(sprintf("   %-52s %s\n", sh[[1]], res))
}
cat("\n   The 4th shape is what jamovi-compiler emits (options= plus only the\n")
cat("   properties the .r.yaml declares), so it is the one modules hit.\n")
cat("   `options` having no default is deliberate and matches every sibling\n")
cat("   except Table; this fix is about name / clearWith / refs.\n\n")

# ── 2. refs, which every other results element accepts ──────────────────────
cat("2. passing refs, as every other results element allows\n")
r2 <- tryCatch({
  State$new(options = Options$new(), name = "myState", refs = "someref")
  "OK"
}, error = function(e) paste("FAILS:", conditionMessage(e)))
cat("   ", r2, "\n\n")

# ── 3. the same call against every sibling class, for contrast ──────────────
cat("3. siblings, same call shape\n")
ns <- asNamespace("jmvcore")
for (nm in c("Html", "Preformatted", "Table", "Image", "Notice", "Group", "State")) {
  gen <- get(nm, envir = ns)
  res <- tryCatch({
    gen$new(options = Options$new(), name = "x", title = "", visible = FALSE,
            clearWith = list(), refs = character())
    "constructs"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("   %-13s %s\n", nm, res))
}

cat("\nRoot cause: State$initialize calls super$initialize() with five positional\n")
cat("arguments; ResultsElement$initialize takes six and none has a default.\n")
