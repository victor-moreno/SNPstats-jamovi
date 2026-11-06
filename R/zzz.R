# Loaded last (R sources load in C-locale order), so this runs after snpPGS.h.R.
# A `type: Level` option cannot have a yaml default, so the compiler emits a bare
# `caseLevel,` formal and snpPGS() fails from R when caseLevel is omitted. Set here,
# the default survives every build route; patch_h.sh only keeps the committed .h.R tidy.
formals(snpPGS)["caseLevel"] <- list(NULL)
