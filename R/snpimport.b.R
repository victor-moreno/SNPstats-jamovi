
# snpImport -- import genotypes from PLINK and VCF files for a list of SNPs.
#
# The reading and decoding live in plink_read.R and payload.R, which know
# nothing about jamovi. This file is orchestration only: pull the options
# apart, call those, fill the tables, and hand the built columns to
# .performOpen, which is the only thing that produces data -- this analysis
# never writes into the sheet it runs in (see .performOpen and the r.yaml
# comment on openNewState for why).
#
# setRowNums is not relevant here for the same reason: there is no Output
# growing the current dataset to reset rows on.

snpImportClass <- R6::R6Class(
  "snpImportClass",
  inherit = snpImportBase,
  private = list(

    .cache     = NULL,
    .cacheDone = FALSE,
    .cacheErr  = NULL,

    # The rows of both tables are created here, not in .run(). jamovi rebuilds
    # the analysis object on every option click, and Table$fromProtoBuf only
    # restores cells into rows that already exist by the end of .init() -- rows
    # added in .run() are never restored. That is why the tables blanked and
    # refilled on every click: the restore had nothing to put cells into.
    # (Same finding as ../SNPstats/R/snpStats.b.R:29-34.)
    #
    # The engine builds that object fresh per request and only then restores,
    # so every table here starts empty and fromProtoBuf never adds a row --
    # which is why a clearWith shorter than the option panel cannot duplicate
    # rows. .syncRows below does not rely on that; a 2026-09 review read the
    # unconditional adds as a duplication bug, and a guard that makes the
    # invariant local is cheaper than an argument.
    #
    # This also means .prepared() decodes the payload once per click rather
    # than once per session: the cache is on the object, and the object is new.
    # It cannot be deferred out of .init() -- which SNPs get a summary row
    # depends on MAF, HWE and missingness, i.e. on the decode -- and the work
    # is bounded by what the browser is allowed to send (PAYLOAD_SAFE, 3.6 MiB).
    #
    # Only the .bim text is parsed here, never the .bed: this runs on every
    # click, so it has to stay cheap.
    .init = function() {

      private$.setInstructions()

      # Row created here, unconditionally, so it exists to restore into on the
      # next rebuild whatever .performOpen wrote to it (see the r.yaml comment
      # on openNewState) -- same reason every other table's rows are created in
      # .init rather than .run.
      if (self$results$openNewState$rowCount == 0)
        self$results$openNewState$addRow(rowKey = 1, values = list(fired = ""))

      if (!private$.haveData()) {
        # An empty report and an empty summary table next to the instructions
        # are just noise.
        self$results$provenance$setVisible(FALSE)
        self$results$summary$setVisible(FALSE)
        return()
      }
      self$results$provenance$setVisible(self$options$showProvenance)
      self$results$summary$setVisible(self$options$showSummary)

      d <- private$.prepared()
      if (is.null(d)) return()

      if (self$options$showProvenance)
        private$.syncRows(self$results$provenance,
                          as.list(private$.provKeys()),
                          list(item = private$.provKeys()))

      # Only the SNPs that survive filtering get a row: a summary listing SNPs
      # that are not in the spreadsheet reads as a bug, whatever the status
      # column says. What was dropped, and why, goes in the import report.
      if (self$options$showSummary) {
        j <- which(d$keep)
        private$.syncRows(self$results$summary, as.list(j), list(
          snp     = d$ids[j],
          chr     = d$bim$chr[j],
          bp      = as.numeric(d$bim$bp[j]),
          alleles = paste0(d$bim$a1[j], "/", d$bim$a2[j])))
      }
    },

    # Create exactly the rows named by `keys`, whatever the table already holds.
    #
    # In the engine that is always an empty table (a fresh analysis object per
    # request, restored only after .init() returns), so this reduces to
    # .addRows. It is written this way so the invariant does not depend on
    # reading enginer.cpp: a table that somehow arrived with rows is rebuilt
    # rather than appended to, and one that already has the right keys is left
    # alone so restored cells survive.
    .syncRows = function(tbl, keys, values = NULL) {
      if (tbl$rowCount > 0L) {
        if (identical(tbl$rowKeys, keys)) return(invisible(FALSE))
        tbl$deleteRows()
      }
      private$.addRows(tbl, keys, values)
      invisible(TRUE)
    },

    # Add many rows without jmvcore's quadratic cost.
    #
    # Table$addRow recomputes `.rowNames <- sapply(.rowKeys, toJSON)` on *every*
    # call (jmvcore/R/table.R:371), so building n rows is O(n^2) serialisations:
    # measured 5.4 s for 500 rows, paid on every option click because rows have
    # to be created in .init(). This does the same work with the serialisation
    # hoisted out of the loop. Falls back to addRow if jmvcore's internals ever
    # move, so a version bump degrades to slow rather than broken.
    .addRows = function(tbl, keys, values = NULL) {
      pr <- tbl$.__enclos_env__$private
      # Checking that the fields exist is not the same as checking that they
      # still mean what this assumes. A jmvcore that kept the names and changed
      # the semantics would corrupt every table silently, so the fast path is
      # proved against jmvcore's own addRow() once per session before it is
      # used, and a disagreement demotes this to the slow path rather than
      # producing wrong rows. See .addRowsUsable().
      ok <- !is.null(pr$.rowKeys) && !is.null(pr$.columns) &&
            is.numeric(pr$.rowCount) && .addRowsUsable()
      if (!ok) {
        for (i in seq_along(keys))
          tbl$addRow(rowKey = keys[[i]],
                     values = lapply(values, function(v) v[i]))
        return(invisible(FALSE))
      }
      .add_rows_fast(tbl, keys, values)
      invisible(TRUE)
    },

    # What counts as "loaded" depends on the format: a VCF names its own
    # samples and carries its own variant fields, a .tped needs its .tfam, and
    # a .bed needs both companions.
    .haveData = function() {
      if (!nzchar(self$options$genoContent)) return(FALSE)
      switch(self$options$sourceFormat,
             vcf  = TRUE,
             ped  = nzchar(self$options$variantContent),
             tped = nzchar(self$options$sampleContent),
             nzchar(self$options$variantContent) &&
               nzchar(self$options$sampleContent))
    },

    # Dispatch on the source format. Each reader returns the same shape -- fam,
    # a .bim-like variant table and an allele-2 dosage matrix -- so nothing
    # downstream knows or cares which format the data came from.
    #
    # The allele *labels* differ legitimately between formats: a VCF names REF
    # first, a .bim picks by frequency. The emitted genotype strings are
    # identical either way, which is what test-formats.R checks.
    .readSource = function() {
      fmt <- self$options$sourceFormat
      geno <- self$options$genoContent

      if (identical(fmt, "vcf")) {
        v <- read_vcf(payload_lines(geno, "VCF"))
        return(list(fam = v$fam, bim = v$bim, dose = v$dose,
                    raw_bytes = nchar(geno) / 4 * 3, skipped = v$skipped,
                    unreadable = v$unreadable, format = "VCF"))
      }

      if (identical(fmt, "ped")) {
        # A .ped carries its own sample fields, so only the .map comes with it.
        pd <- read_ped(payload_lines(geno, ".ped"),
                       payload_lines(self$options$variantContent, ".map"))
        return(list(fam = pd$fam, bim = pd$bim, dose = pd$dose,
                    raw_bytes = nchar(geno) / 4 * 3, skipped = 0L,
                    unreadable = pd$unreadable,
                    format = "PLINK text (.ped/.map)"))
      }

      if (identical(fmt, "tped")) {
        t <- read_tped(payload_lines(geno, ".tped"),
                       payload_lines(self$options$sampleContent, ".tfam"))
        return(list(fam = t$fam, bim = t$bim, dose = t$dose,
                    raw_bytes = nchar(geno) / 4 * 3, skipped = 0L,
                    unreadable = t$unreadable,
                    format = "PLINK text (.tped/.tfam)"))
      }

      fam <- read_fam(payload_lines(self$options$sampleContent, ".fam"))
      bim <- read_bim(payload_lines(self$options$variantContent, ".bim"))
      raw <- payload_bytes(geno, "genotype data")

      # Whether the three files were the same three files. The length check
      # below cannot answer that: a variant's offset is computed from the .fam's
      # sample count and the .bim's line number, and the expectation here is
      # re-derived from those same two numbers, so a stale .fam agrees with
      # itself and hands back a different variant's bytes as a confident
      # genotype. Only the .bed's own size settles it, and only the browser can
      # see it, so it is carried across and re-checked here.
      bed_check_dims(self$options$sourceDims, length(fam$iid))

      # The structural check that the payload matches the .bim it came with: the
      # byte count must be exactly one block of ceil(n/4) per .bim line.
      # A .bed cannot carry an unreadable call: two bits, and one of the four
      # values is missing. Only the allele *labels* in its .bim can be
      # unreadable, which geno_factor deals with where the labels are used.
      list(fam = fam, bim = bim,
           dose = bed_decode_slices(raw, length(fam$iid), length(bim$id)),
           raw_bytes = length(raw), skipped = 0L, unreadable = 0L,
           format = "PLINK binary (.bed)",
           dims_checked = nzchar(self$options$sourceDims))
    },

    # Fixed for a given set of options, so the report keeps its shape across
    # clicks. Rows whose value does not apply say so rather than disappearing.
    .provKeys = function() {
      k <- c("Files", "Samples", "SNPs requested", "SNPs received",
             "SNPs not found")
      if (self$options$filterSamples) k <- c(k, "Samples dropped")
      # .fillProvenance writes this row, but a key it did not create is skipped
      # rather than added, so without this the warning never appeared. Read from
      # the cache, which both callers have already filled.
      if (!is.null(private$.cache) && !isTRUE(private$.cache$fam$iid_unique))
        k <- c(k, "Sample IDs")
      if (!is.null(private$.cache) && isTRUE(private$.cache$unreadable > 0))
        k <- c(k, "Unreadable calls")
      # Only for a .bed loaded before the trio check existed. Saying so is the
      # honest position: the import is not wrong, it is unverified.
      if (identical(self$options$sourceFormat, "bed") &&
          !is.null(private$.cache) && !isTRUE(private$.cache$dims_checked))
        k <- c(k, "File consistency")
      if (self$options$applyFilters)  k <- c(k, "SNPs kept after filters")
      if (!is.null(self$options$covFile))
        k <- c(k, "Covariates", "Covariates matched")
      c(k, "Format", "Genotype data")
    },

    # Decode, compute and filter -- once per analysis object.
    #
    # Not named .data: jmvcore's Analysis base class keeps the dataset in a
    # private field of that name and assigns to it, which an R6 method binding
    # would block.
    #
    # .init() and .run() share one object (jamovi builds a fresh one per option
    # click), so caching here means a click pays for the decode once instead of
    # twice. NULL means the payload could not be read; .run() reports why.
    .prepared = function() {
      if (private$.cacheDone) return(private$.cache)
      private$.cacheDone <- TRUE

      d <- try({
        src <- private$.readSource()
        fam <- src$fam; bim <- src$bim; dose <- src$dose
        raw <- src$raw_bytes

        samp <- private$.sampleFilter(dose)
        if (!all(samp$keep)) {
          dose <- dose[samp$keep, , drop = FALSE]
          fam  <- private$.subsetFam(fam, samp$keep)
        }
        # Which samples the HWE test uses. Controls only by default, matching
        # PLINK's --hwe, but only when the .fam actually says who the controls
        # are -- otherwise the test would silently run on nobody.
        ctrl <- !is.na(fam$pheno) & fam$pheno == "control"
        use_ctrl <- self$options$hweGroup == "controls" && any(ctrl)
        st <- snp_stats_all(dose, if (use_ctrl) ctrl else NULL)
        st$hwe_in <- if (use_ctrl) "controls" else "all samples"
        st$hwe_fellback <- self$options$hweGroup == "controls" && !any(ctrl)

        # Covariates are matched onto the samples that survived filtering, so a
        # dropped sample cannot reappear through the covariate file.
        cov <- NULL
        covFile <- self$options$covFile
        if (!is.null(covFile)) {
          tab <- read_covariate_table(
            file_lines(covFile$path, covFile$filename, "covariate file"),
            self$options$covIdCol)
          cov <- merge_covariates(tab, fam)
        }

        # Every column this analysis writes shares one namespace, so the names
        # are settled here, once, and the table and the opened dataset cannot
        # disagree. A covariate named after a SNP is the way two columns end up
        # wanting the same name in practice.
        #
        # The four sample names are reserved whether or not they are emitted,
        # so toggling emitSamples cannot make a SNP column collide with one.
        samp_keys <- c("FID", "IID", "sex", "phenotype")
        cov_names <- if (is.null(cov)) character(0)
                     else cov_unique_names(names(cov$values), samp_keys)
        vids <- variant_ids(bim)
        ids  <- cov_unique_names(vids, c(samp_keys, cov_names))

        list(fam = fam, bim = bim, dose = dose, stats = st, samp = samp,
             cov = cov, cov_names = cov_names,
             ids = ids, n_received = length(bim$id),
             # against vids, not ids: a name that moved to avoid a covariate is
             # not a variant that had no ID
             n_renamed = sum(vids != sanitise_id(bim$id)),
             raw_bytes = raw, skipped = src$skipped,
             unreadable = src$unreadable, dims_checked = src$dims_checked,
             format = src$format, keep = private$.variantFilter(st))
      }, silent = TRUE)

      if (inherits(d, "try-error")) {
        private$.cacheErr <- conditionMessage(attr(d, "condition"))
        return(NULL)
      }
      private$.cache <- d
      d
    },

    .run = function() {

      # The notice is not cleared on every option change (see the clearWith note
      # in snpImport.r.yaml), so a stale warning would otherwise outlive the
      # condition that raised it.
      self$results$notice$setVisible(FALSE)

      if (!private$.haveData()) {
        private$.reportLoadProblem()
        return()
      }

      d <- private$.prepared()
      if (is.null(d)) {
        # .init has already created rows from the .bim; leaving them on screen
        # with empty statistics reads as a partial success.
        self$results$summary$setVisible(FALSE)
        self$results$provenance$setVisible(FALSE)
        private$.notice(paste0("<b>Import failed.</b> ",
                               .snpi_escape_html(private$.cacheErr)), "error")
        return()
      }

      # openNewState's row exists once .init has run; its cell says whether
      # perform() has already fired (see the r.yaml comment on openNewState).
      # Nothing is produced until the button is pressed, so say so rather than
      # letting a loaded import look like it did nothing.
      already_opened <- isTRUE(
        self$results$openNewState$getCell(col = "fired", rowNo = 1)$value == "yes")
      if (!already_opened)
        private$.notice(paste0(
          "<b>Genotypes are loaded.</b> Press <b>Open as new dataset</b> to ",
          "launch them in a new jamovi window \u2014 this panel does not write ",
          "into the sheet it runs in."),
          "warn")

      private$.fillSummary(d)
      private$.fillProvenance(d)

      # Every sample dropped is not a small import, it is an empty one: the
      # summary would list every SNP with N = 0, and nothing on screen would say
      # why. The threshold is the user's, so this reports rather than
      # overrides it -- but it reports loudly, and .performOpen refuses to open
      # a dataset with no rows rather than opening an empty one silently.
      if (length(d$fam$iid) == 0) {
        dropped <- d$samp$by_missing + d$samp$by_het
        private$.notice(sprintf(paste0(
          "<b>Every sample was dropped by the sample filter.</b> All %s of them, ",
          "so there is nothing left to import. Raise <b>Max missing per sample</b> ",
          "above %.0f%%%s, or untick <b>Drop low-quality samples</b>."),
          format(dropped, big.mark = ","), self$options$maxIndMissing,
          if (d$samp$by_het > 0)
            sprintf(" (%s of them went to the heterozygosity rule)",
                    format(d$samp$by_het, big.mark = ",")) else ""), "error")
        self$results$summary$setVisible(FALSE)
      }

      private$.performOpen(d)
    },

    # Per-individual QC. Returns the keep mask plus what each rule removed, so
    # the report can attribute the drops.
    .sampleFilter = function(dose) {
      n <- nrow(dose)
      out <- list(keep = rep(TRUE, n), by_missing = 0L, by_het = 0L,
                  het_unreliable = FALSE)
      if (!self$options$filterSamples || n == 0) return(out)

      ss <- sample_stats(dose)

      bad_miss <- ss$missing_rate > self$options$maxIndMissing / 100
      out$by_missing <- sum(bad_miss)
      out$keep <- !bad_miss

      sd_thresh <- self$options$hetSd
      if (sd_thresh > 0) {
        h <- ss$het_rate
        h[!out$keep] <- NA                    # judge spread on the survivors
        mu <- mean(h, na.rm = TRUE); sdv <- stats::sd(h, na.rm = TRUE)
        if (is.finite(sdv) && sdv > 0) {
          bad_het <- !is.na(h) & abs(h - mu) > sd_thresh * sdv
          out$by_het <- sum(bad_het)
          out$keep <- out$keep & !bad_het
        }
        # Heterozygosity over a few hundred SNPs is dominated by sampling noise;
        # say so rather than letting a tidy-looking count imply otherwise.
        out$het_unreliable <- ncol(dose) < 1000
      }
      out
    },

    .subsetFam = function(fam, keep) {
      fam$fid <- fam$fid[keep]; fam$iid <- fam$iid[keep]
      fam$sex <- fam$sex[keep]; fam$pheno <- fam$pheno[keep]
      fam$fid_informative <- length(unique(fam$fid)) > 1L
      fam$iid_unique      <- !anyDuplicated(fam$iid)
      fam
    },

    # Requested IDs, plus effect alleles when a weights file supplied them.
    .selection = function() {
      if (nzchar(self$options$snpListContent))
        return(parse_snp_selection(
          payload_lines(self$options$snpListContent, "selection file")))
      txt <- self$options$snpListText
      if (nzchar(trimws(txt))) {
        ids <- unlist(strsplit(txt, "[\r\n \t,;]+"))
        return(list(ids = ids[nzchar(ids)], effect_allele = NULL, weight = NULL))
      }
      list(ids = character(0), effect_allele = NULL, weight = NULL)
    },

    # Filters are evaluated on every SNP so the summary can show *why* each one
    # was dropped, rather than silently shortening the table.
    .variantFilter = function(st) {
      k <- length(st$n)
      if (!self$options$applyFilters) return(rep(TRUE, k))
      vapply(seq_len(k), function(j) !nzchar(private$.filterReason(st, j)),
             logical(1))
    },

    .filterReason = function(st, j) {
      if (!self$options$applyFilters) return("")
      tot <- st$n[j] + st$missing[j]
      if (tot > 0 && st$missing[j] / tot > self$options$maxMissing / 100)
        return("dropped: missingness")
      if (!is.na(st$maf[j]) && st$maf[j] < self$options$minMaf)
        return("dropped: MAF")
      if (!is.na(st$hwe_p[j]) && st$hwe_p[j] < self$options$hweP)
        return("dropped: HWE")
      ""
    },

    # -- output columns ----------------------------------------------------

    # Column names, titles, descriptions, measure types and values, in the
    # order the opened dataset gets them. Kept separate from .performOpen so
    # the "what to build" and "how to hand it to jamovi" halves stay apart.
    .buildColumns = function(d) {
      fam <- d$fam; bim <- d$bim; ids <- d$ids; dose <- d$dose
      n    <- length(fam$iid)
      dos  <- self$options$dosage
      # No samples means no columns, not columns of length zero: an empty key
      # set here is what makes .performOpen report an error instead of opening
      # a dataset with no data in it.
      idx  <- if (n == 0) integer(0) else which(d$keep)

      keys   <- character(0); titles <- character(0)
      descs  <- character(0); types  <- character(0)
      values <- list()

      add <- function(key, title, desc, type, value) {
        keys   <<- c(keys, key);     titles <<- c(titles, title)
        descs  <<- c(descs, desc);   types  <<- c(types, type)
        values[[key]] <<- value
      }

      # Sample columns first, and only when asked. Sex and phenotype are
      # emitted only when the .fam actually carries them.
      if (self$options$emitSamples && n > 0) {
        # A .fam row is identified by the (FID, IID) pair; IID alone is only
        # required to be unique within a family. FID is emitted when it carries
        # information -- all-identical FIDs are noise -- and always when IIDs
        # repeat, because then IID alone does not identify a sample.
        if (fam$fid_informative || !fam$iid_unique)
          add("FID", "FID", "Family ID from .fam", "nominal",
              factor(fam$fid, levels = unique(fam$fid)))
        add("IID", "IID", "Subject ID from .fam", "nominal",
            factor(fam$iid, levels = unique(fam$iid)))
        if (any(!is.na(fam$sex)))
          add("sex", "sex", "Sex from .fam (1 = male, 2 = female)",
              "nominal", fam$sex)
        if (any(!is.na(fam$pheno)))
          add("phenotype", "phenotype",
              "Phenotype from .fam (1 = control, 2 = case)",
              "nominal", fam$pheno)
      }

      # Covariates go between the phenotype and the SNPs.
      if (!is.null(d$cov) && n > 0) {
        nms <- d$cov_names
        for (i in seq_along(nms))
          add(nms[i], nms[i],
              sprintf("Covariate '%s' matched by %s",
                      names(d$cov$values)[i], d$cov$matched_on),
              if (is.factor(d$cov$values[[i]])) "nominal" else "continuous",
              d$cov$values[[i]])
      }

      for (j in idx) {
        key <- ids[j]
        if (dos)
          add(key, key,
              sprintf("Dosage of %s at %s", bim$a2[j], bim$id[j]),
              "continuous", as.numeric(dose[, j]))
        else
          add(key, key,
              sprintf("Genotype at %s (%s/%s)", bim$id[j], bim$a1[j], bim$a2[j]),
              "nominal", geno_factor(dose[, j], bim$a1[j], bim$a2[j]))
      }

      list(keys = keys, titles = titles, descs = descs, types = types,
           values = values)
    },

    # Hands .buildColumns' columns to jamovi as a new dataset: perform()
    # writes them straight to a temp .omv on disk and opens it through
    # jamovi's normal file-open pipeline, replacing the current window's data
    # if it is still blank and unedited, opening a new window otherwise
    # (client/main/backstage.ts's requestOpen -- not this module's doing).
    #
    # self$options$openNew stays TRUE forever once the button is clicked --
    # jamovi does not reset an Action option's value, it only disables the
    # button client-side -- so every later .run(), for any option change, sees
    # openNew == TRUE again. perform() must therefore fire at most once. The
    # guard is openNewState's cell (see the r.yaml comment on it for why that
    # table exists at all: the natural guard, the 'openNew' results item
    # itself, cannot be used -- declaring it crashes this jamovi/jmvcore).
    .performOpen = function(d) {
      if (!isTRUE(self$options$openNew)) return(invisible(FALSE))

      state <- self$results$openNewState
      if (isTRUE(state$getCell(col = "fired", rowNo = 1)$value == "yes"))
        return(invisible(FALSE))

      option <- self$options$option("openNew")
      if (is.null(option) || is.null(option$perform)) return(invisible(FALSE))

      option$perform(function(action) {
        cols <- private$.buildColumns(d)
        if (length(cols$keys) == 0)
          return(list(status = "error",
                      message = paste("Nothing to open: no samples or",
                                       "columns to write.")))
        df <- as.data.frame(cols$values[cols$keys], stringsAsFactors = FALSE,
                             check.names = FALSE)
        names(df) <- cols$keys
        list(data = df, title = "Imported genotypes")
      })
      state$setRow(rowNo = 1, values = list(fired = "yes"))
      invisible(TRUE)
    },

    # -- tables ------------------------------------------------------------

    .fillSummary = function(d) {
      tbl <- self$results$summary
      if (!self$options$showSummary) return()

      bim <- d$bim; ids <- d$ids; st <- d$stats
      sel <- private$.selection()
      idx <- which(d$keep)

      ea <- sel$effect_allele
      if (!is.null(ea)) {
        m <- match(bim$id, sel$ids)
        ea <- ea[m]
      }
      # A weights file names the allele each weight belongs to, so the useful
      # frequency is that allele's, not the rarer one's: a PGS weight for an
      # allele carried by 88% of samples is checked against 0.88, and "MAF 0.12"
      # invites reading it the wrong way round. The column changes what it holds
      # and says so, rather than quietly changing what MAF means.
      eaf <- !is.null(ea) && any(!is.na(ea) & nzchar(ea))
      freq_title <- if (eaf) "EAF" else "MAF"
      tbl$getColumn("maf")$setTitle(freq_title)
      if (!is.null(ea)) tbl$getColumn("allele")$setVisible(TRUE)

      tbl$setNote("note", if (eaf)
        paste("Alleles are (allele 1 / allele 2) as in the source file; ",
              "EAF is the frequency of the effect allele named in the weights file. ",             
              "Min MAF filter always uses the minor allele frequency.")
      else
        paste("Alleles are (allele 1 / allele 2) as recorded in the source file. ",
              "MAF is the frequency of whichever allele is rarer."))

      # This table and the spreadsheet order the two alleles by different rules,
      # both on purpose: here they are in the file's order so the counts can be
      # read against it, and there they are sorted so the same heterozygote
      # reads the same whichever format it came from. Nothing said so, and the
      # first column of "AA / AB / BB" is not always the first level of the
      # column next to it.
      #
      # None of which applies in dosage mode: the columns are numeric counts of
      # allele 2, they carry no labels to be in any order, and a note about
      # A/A and A/G describes something the spreadsheet does not contain.
      # Cleared rather than skipped, because `dosage` is not in this table's
      # clearWith, so a note set on a previous run would otherwise survive the
      # checkbox being ticked.
      tbl$setNote("labels", if (self$options$dosage) NULL else paste(
        "Genotype labels in the imported columns are written with the two",
        "alleles in alphabetical order, which is not always the order here:",
        "a variant recorded as G/A is counted as G/G, G/A, A/A in this table",
        "and labelled A/A, A/G, G/G in the spreadsheet." ))

      # Every other column describes all samples; only HWE may not, so the
      # table has to say which set it used and how big that set was.
      tbl$setNote("hwe", sprintf(
        "HWE is the exact test computed in %s (n = %s)%s. N, missing, %s and the genotype counts describe all %s imported samples.",
        st$hwe_in, format(st$hwe_n, big.mark = ","),
        if (isTRUE(st$hwe_fellback))
          " \u2014 'controls only' was asked for, but the .fam carries no phenotype" else "",
        freq_title, format(length(d$fam$iid), big.mark = ",")))
      # The rows already exist (see .init) and cover exactly the kept SNPs;
      # only the computed cells are written, so the table keeps its shape.
      if (tbl$rowCount != length(idx)) return()

      for (j in idx) {
        status <- ""

        if (!is.null(ea) && !is.na(ea[j])) {
          am <- allele_match(ea[j], bim$a1[j], bim$a2[j])
          if (am != "match")
            status <- paste(c(status, sprintf("allele %s", am))[nzchar(c(status, "x"))],
                            collapse = "; ")
        }

        # NA when the effect allele is not one of the two observed -- the
        # status column has already said which case that is.
        f <- if (eaf) effect_freq(ea[j], bim$a1[j], bim$a2[j], st$freq_a2[j])
             else NA_real_

        tbl$setRow(rowKey = j, values = list(
          n          = st$n[j],
          missing    = st$missing[j],
          maf        = if (is.na(f)) st$maf[j] else f,
          genoCounts = sprintf("%d / %d / %d", st$n11[j], st$n12[j], st$n22[j]),
          hwePval    = st$hwe_p[j],
          allele     = if (is.null(ea) || is.na(ea[j])) "" else ea[j],
          status     = if (nzchar(status)) status else "ok"
        ))
      }
    },

    .fillProvenance = function(d) {
      tbl <- self$results$provenance
      if (!self$options$showProvenance) return()

      # Rows come from .init; this only writes their values. A key that .init
      # did not create (because its option is off) is skipped rather than
      # silently changing the table's shape mid-run.
      keys <- private$.provKeys()
      add <- function(k, v) if (k %in% keys)
        tbl$setRow(rowKey = k, values = list(item = k, value = v))

      fam <- d$fam; bim <- d$bim; samp <- d$samp; keep <- d$keep
      sel <- private$.selection()
      n <- length(fam$iid)          # after sample filtering
      k <- d$n_received             # variants received from the browser
      requested <- length(sel$ids)
      matched   <- sum(sel$ids %in% bim$id)

      add("Files",     if (nzchar(self$options$genoFilename))
                         self$options$genoFilename else "(unnamed)")

      dropped <- samp$by_missing + samp$by_het
      add("Samples", if (dropped > 0)
            sprintf("%s kept of %s", format(n, big.mark = ","),
                    format(n + dropped, big.mark = ","))
          else format(n, big.mark = ","))

      if (self$options$filterSamples) {
        bits <- character(0)
        if (samp$by_missing > 0)
          bits <- c(bits, sprintf("%d over %.0f%% missing",
                                  samp$by_missing, self$options$maxIndMissing))
        if (samp$by_het > 0)
          bits <- c(bits, sprintf("%d beyond %.1f SD heterozygosity%s",
                                  samp$by_het, self$options$hetSd,
                                  if (samp$het_unreliable)
                                    " (unreliable under 1000 SNPs)" else ""))
        add("Samples dropped", if (length(bits)) paste(bits, collapse = "; ") else "none")
      }

      if (!fam$iid_unique)
        add("Sample IDs",
            "IIDs are not unique \u2014 FID is included, and the pair identifies a sample")

      add("SNPs requested",
          if (requested > 0) format(requested, big.mark = ",")
          else "(no list given \u2014 every variant in the .bim)")
      # A renamed column is one the user cannot find by the ID they asked for,
      # so it is stated rather than left to be discovered in the spreadsheet.
      add("SNPs received", paste0(
        format(k, big.mark = ","),
        if (isTRUE(d$n_renamed > 0))
          sprintf(" \u2014 %d had no usable ID and are named chr:bp", d$n_renamed)
        else ""))

      missing_ids <- if (requested > 0) setdiff(sel$ids, bim$id) else character(0)
      add("SNPs not found",
          if (length(missing_ids) == 0) "none"
          else sprintf("%d \u2014 %s%s", length(missing_ids),
                       paste(utils::head(missing_ids, 5), collapse = ", "),
                       if (length(missing_ids) > 5) ", \u2026" else ""))

      if (self$options$applyFilters)
        add("SNPs kept after filters", sprintf("%d of %d", sum(keep), k))

      if (!is.null(d$cov)) {
        cv <- d$cov
        add("Covariates", sprintf("%d: %s (matched on %s)",
                                  length(cv$values),
                                  paste(names(cv$values), collapse = ", "),
                                  cv$matched_on))
        # Both directions, because an ID convention mismatch shows up as
        # "nothing matched" and has to be diagnosable from this table alone.
        bits <- sprintf("%d of %d samples", cv$matched, cv$n_geno)
        if (cv$matched == 0)
          bits <- paste(bits, "\u2014 check the ID column matches the .fam IIDs")
        if (length(cv$unmatched_geno))
          bits <- paste0(bits, "; no covariates for ",
                         length(cv$unmatched_geno), " (",
                         paste(utils::head(cv$unmatched_geno, 3), collapse = ", "),
                         if (length(cv$unmatched_geno) > 3) ", \u2026" else "", ")")
        if (length(cv$unmatched_cov))
          bits <- paste0(bits, "; ", length(cv$unmatched_cov),
                         " covariate rows unused")
        if (cv$dup_cov)
          bits <- paste0(bits, "; duplicate IDs in the covariate file \u2014 the ",
                         "first match was used")
        add("Covariates matched", bits)
      }


      # Genotype calls the file wrote in a code this variant does not have.
      # Reported rather than left to be inferred from the missing rate: a
      # column of unexpected codes reads as a badly typed SNP, and the
      # difference matters when deciding whether to re-export the file.
      add("File consistency", paste(
        "not checked \u2014 this import was loaded by an older version of the",
        "module; press Load genotypes again to verify the .bed, .bim and .fam",
        "belong together"))

      if (isTRUE(d$unreadable > 0))
        add("Unreadable calls", sprintf(
          "%s set to missing \u2014 the code was neither allele nor a missing marker",
          format(d$unreadable, big.mark = ",")))

      add("Format", if (isTRUE(d$skipped > 0))
            sprintf("%s \u2014 %d line%s skipped (multi-allelic or indel)",
                    d$format, d$skipped, if (d$skipped == 1) "" else "s")
          else d$format)
      # "0.00 MB" is what every ordinary import used to report: a few hundred
      # SNPs is tens of kilobytes.
      add("Genotype data", if (d$raw_bytes < 1024^2)
            sprintf("%.0f KB", d$raw_bytes / 1024)
          else sprintf("%.2f MB", d$raw_bytes / 1024^2))
    },

    # A refusal in the browser leaves the carriers empty, so from R it looks
    # exactly like an import nobody has started yet: the instructions come back
    # and the reason is only in the status box, which is one line wide and shows
    # about four words of a message that runs to three sentences. Restate it
    # here, where there is room for it.
    .reportLoadProblem = function() {
      kind <- self$options$loadProblem
      msg  <- self$options$loadStatus
      if (!nzchar(kind) || !nzchar(msg)) return()
      private$.notice(
        paste0("<b>",
               if (identical(kind, "warn")) "Nothing to import."
               else "The genotypes were not loaded.",
               "</b> ", .snpi_escape_html(msg)),
        if (identical(kind, "warn")) "warn" else "error")
    },

    # -- presentation ------------------------------------------------------

    .notice = function(html, kind = "warn") {
      bg <- if (kind == "error") "#fdecea" else "#fff8e1"
      bd <- if (kind == "error") "#d93025" else "#f9a825"
      self$results$notice$setContent(sprintf(
        "<div style='background:%s;border-left:4px solid %s;padding:8px 12px;
         margin:4px 0;border-radius:3px'>%s</div>", bg, bd, html))
      self$results$notice$setVisible(TRUE)
    },

    .setInstructions = function() {
      if (nzchar(self$options$genoContent)) {
        self$results$instructions$setVisible(FALSE)
        return()
      }
      self$results$instructions$setContent(paste0(
        "<div style='padding:8px 12px'>",
        "<p><b>Load genotypes from PLINK and VCF files, and open them as a ",
        "new dataset.</b> This panel is a loader, not part of the analysis: ",
        "it never writes into the sheet it runs in. Once genotypes are ",
        "loaded, <b>Open as new dataset</b> launches a new jamovi window ",
        "with the samples, covariates and genotypes already in it, ready ",
        "for SNP Analysis and Polygenic Score \u2014 continue there.</p>",
        "<ol>",
        "<li><b>Genotype files</b> \u2014 <code>.bed</code>+<code>.bim</code>",
        "+<code>.fam</code>, <code>.ped</code>+<code>.map</code>, ",
        "<code>.tped</code>+<code>.tfam</code>, or a <code>.vcf</code> / ",
        "<code>.vcf.gz</code> on its own. Where a set is needed, select all of ",
        "it at once: the web page cannot open a file you did not pick.</li>",
        "<li><b>SNPs to import</b> (optional) \u2014 paste the SNP IDs, or load ",
        "an SNP ID list or a PGS-Catalog weights file. Leaving it empty reads ",
        "every SNP in the file, up to the transfer limit.</li>",
        "<li><b>Covariates</b> (optional) \u2014 a delimited file with one ",
        "row per sample, matched to the genotypes by sample ID. The two files ",
        "need not hold the same samples: a sample with no covariate row keeps ",
        "empty cells, and a covariate row matching no genotype is dropped.</li>",
        "<li><b>Load genotypes</b> \u2014 nothing is read until you press it.</li>",
        "<li><b>Open as new dataset</b> \u2014 builds the samples, covariates ",
        "and genotypes described in the import report and opens them ",
        "as a new dataset. </li>",
        "</ol>",
        "<p><b>Prefer <code>.bed/.bim/.fam</code>.</b> They load faster: only ",
        "the SNPs you ask for are read, so the file itself can be any size. ",
        "For a format that is not supported, convert it with plink first:<br>",
        "<code>plink2 --pfile mydata --make-bed --out mydata</code></p>",
        # Not a footnote. The payload is stored in the analysis options, so the
        # .omv is the genotypes -- and an .omv is a file people email each
        # other. Nothing else in the panel says so, and the consequence for
        # human genetic data is not one to leave to be discovered.
        "<p style='background:#fff8e1;border-left:4px solid #f9a825;",
        "padding:6px 10px;border-radius:3px'><b>Saving this panel's own ",
        "analysis saves the genotypes too.</b> The selected genotype bytes ",
        "and the sample IDs stay inside its <code>.omv</code> for as long as ",
        "the files are loaded here, whether or not you have opened them as a ",
        "dataset yet \u2014 so saving and sharing a copy of <b>this</b> ",
        "analysis (not the dataset it opens) shares individual-level genetic ",
        "data. Pick a different, empty selection of files before sharing one, ",
        "or delete this analysis.</p>",
        "</div>"))
      self$results$instructions$setVisible(TRUE)
    }
  )
)

# Is the private-API fast path in .addRows still equivalent to jmvcore's own
# addRow()?
#
# .addRows writes straight into a Table's private state to avoid jmvcore's
# quadratic row building. That is a dependency on internals that carry no
# compatibility promise, and testing that the *fields* exist proves nothing
# about what they now mean -- a release that kept `.rowKeys` and changed how it
# is indexed would fill every table with wrong rows and no error.
#
# So the two are compared on a throwaway table the first time a table is built,
# and the fast path is used only if they agree. The cost is one two-row table
# per session; the alternative is trusting a private interface indefinitely.
# The body of the fast path, out here so the self-check below exercises the
# same code the analysis runs rather than a copy of it that could drift.
.add_rows_fast <- function(tbl, keys, values = NULL) {
  pr <- tbl$.__enclos_env__$private
  n0 <- pr$.rowCount
  for (i in seq_along(keys)) {
    idx <- n0 + i
    pr$.rowKeys[idx] <- list(keys[[i]])
    pr$.rowCount <- idx
    for (column in pr$.columns) {
      nm <- column$name
      if (!is.null(values) && nm %in% names(values))
        column$addCell(values[[nm]][i], .key = keys[[i]], .index = idx)
      else
        column$addCell(.key = keys[[i]], .index = idx)
    }
  }
  # The whole point: one serialisation pass instead of one per row.
  tj <- getFromNamespace("toJSON", "jmvcore")
  pr$.rowNames <- vapply(pr$.rowKeys, tj, "", USE.NAMES = FALSE)
  invisible(TRUE)
}

.ADDROWS_OK <- new.env(parent = emptyenv())

.addRowsUsable <- function() {
  if (!is.null(.ADDROWS_OK$value)) return(.ADDROWS_OK$value)
  .ADDROWS_OK$value <- isTRUE(tryCatch({
    mk <- function() {
      t <- jmvcore::Table$new(options = jmvcore::Options$new(), name = "probe")
      t$addColumn(name = "v", type = "text")
      t$addColumn(name = "n", type = "integer")
      t
    }
    vals <- list(v = c("a", "b"), n = c(10L, 20L))

    ref <- mk()
    for (i in 1:2)
      ref$addRow(rowKey = i, values = lapply(vals, function(x) x[i]))

    fast <- mk()
    .add_rows_fast(fast, as.list(1:2), vals)

    # Same cells, same shape, same keys. rowCount is compared by value: the
    # fast path arrives at it by arithmetic and can hold a double where addRow
    # holds an integer, which is not a disagreement about anything.
    identical(as.data.frame(ref), as.data.frame(fast)) &&
      ref$rowCount == fast$rowCount &&
      identical(ref$.__enclos_env__$private$.rowNames,
                fast$.__enclos_env__$private$.rowNames)
  }, error = function(e) FALSE))
  .ADDROWS_OK$value
}

# Minimal escaping for text interpolated into the notice panel. Only used for
# error messages we generate ourselves, but those can contain a file name.
#
# Prefixed rather than called htmlEscape: everything in this package shares one
# namespace with SNPstats after the merge, and a generic utility name is the one
# most likely to already exist there. The domain names (snp_stats,
# MISSING_ALLELE, hwe_exact_p) are deliberately left alone -- they have to agree
# with SNPstats at merge time, and renaming them now would make that harder, not
# easier.
.snpi_escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
