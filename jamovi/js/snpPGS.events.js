"use strict";

// ─────────────────────────────────────────────────────────────────────────────
// snpPGS.events.js
//
// Handles two UI responsibilities:
//
//   1. onChange_snpCols  — keeps snpGrid in sync with the selected SNP columns
//      (same pattern as weightedmean.events.js):
//        • Remove rows for SNPs no longer selected.
//        • Add skeleton rows for newly added SNPs (empty alleles, weight = 1).
//        • Preserve any edits the user has already made to existing rows.
//
//   2. onChange_weightsPath / onChange_weightsSep  — reads the PGS Catalog
//      file, parses the metadata header block and the tab-delimited data
//      section, then merges catalog fields into the snpGrid rows that match
//      by rsID.  Rows with no catalog match keep their current values.
//
//   3. update  — required by jamovi for VariableSupplier; intentionally empty.
//
// PGS Catalog file format (v2.0) notes
// ──────────────────────────────────────
//   • Comment / metadata lines start with '#' (including '##' and '###').
//   • Key-value metadata lines look like:  #pgs_id=PGS001901
//     We extract: pgs_id, pgs_name, trait_reported, weight_type,
//                 genome_build, variants_number.
//   • The first non-comment line is the tab-delimited column header.
//   • Data columns of interest:
//       rsID            → rsid   (note capital "ID" in the real files)
//       chr_name        → chr
//       chr_position    → pos
//       effect_allele   → effect_allele
//       other_allele    → other_allele
//       effect_weight   → effect_weight  (may be negative beta — never coerce 0→1)
// ─────────────────────────────────────────────────────────────────────────────

// ── Helpers ──────────────────────────────────────────────────────────────────

/**
 * Build a blank grid row for a SNP column that has no catalog entry yet.
 */
function blankRow(colName) {
    return {
        rsid:          colName,
        effect_allele: "",
        other_allele:  "",
        effect_weight: null,   // null = "not yet set from catalog"
        chr:           "",
        pos:           "",
        matched:       false
    };
}

/**
 * Parse the PGS Catalog '#key=value' metadata header lines.
 * Returns a plain object with the extracted fields.
 *
 * Example input lines:
 *   #pgs_id=PGS001901
 *   #weight_type=beta
 *   #genome_build=GRCh37
 *   #variants_number=5
 */
function parseMetadata(lines) {
    const meta = {};
    const keys = ["pgs_id", "pgs_name", "trait_reported", "trait_mapped",
                  "weight_type", "genome_build", "variants_number"];
    for (const line of lines) {
        if (!line.startsWith("#")) break;           // stop at first data line
        const inner = line.replace(/^#+/, "");      // strip leading hashes
        const eq = inner.indexOf("=");
        if (eq === -1) continue;
        const k = inner.slice(0, eq).trim().toLowerCase();
        const v = inner.slice(eq + 1).trim();
        if (keys.includes(k)) meta[k] = v;
    }
    return meta;
}

/**
 * Detect delimiter from the first non-comment data line.
 * PGS Catalog files are always tab-delimited; this also handles CSV exports.
 */
function detectSep(lines) {
    const dataLine = lines.find(l => l.trim().length > 0 && !l.startsWith("#")) || "";
    const tabs   = (dataLine.match(/\t/g)  || []).length;
    const commas = (dataLine.match(/,/g)   || []).length;
    return tabs >= commas ? "\t" : ",";
}

/**
 * Resolve a column index from the header row using a list of candidate names.
 * Matching is case-insensitive.
 *
 * @param {string[]} header      - lowercased header tokens
 * @param {string[]} candidates  - preferred names in priority order
 * @returns {number}  column index, or -1 if not found
 */
function colIdx(header, candidates) {
    for (const c of candidates) {
        const idx = header.indexOf(c.toLowerCase());
        if (idx !== -1) return idx;
    }
    return -1;
}

/**
 * Parse a PGS Catalog-format file string.
 *
 * Returns { meta, catalogMap } where:
 *   meta       — object with PGS metadata fields (see parseMetadata)
 *   catalogMap — Map: rsid_lowercase → { effect_allele, other_allele,
 *                                         effect_weight, chr, pos }
 *
 * Handles the real PGS Catalog v2.0 format:
 *   - '###', '##', and '#' comment/metadata lines are all stripped
 *   - Column header:  rsID  chr_name  chr_position  effect_allele  other_allele  effect_weight
 *   - effect_weight may be a negative float (beta); we never replace a
 *     valid numeric (including 0) with a fallback of 1.
 */
function parseCatalogText(text, sep) {
    const lines = text.split(/\r?\n/);

    const meta = parseMetadata(lines);

    // Data lines = non-empty, non-comment lines
    const dataLines = lines.filter(l => l.trim().length > 0 && !l.startsWith("#"));
    if (dataLines.length === 0) return { meta, catalogMap: new Map() };

    const header = dataLines[0].split(sep).map(h => h.trim().toLowerCase());

    // Column resolution — list real PGS Catalog names first, then alternatives
    const iRsid   = colIdx(header, ["rsid",          "variant_id", "snp", "snp_id", "marker_name"]);
    const iEA     = colIdx(header, ["effect_allele"]);
    const iOA     = colIdx(header, ["other_allele",   "ref_allele", "non_effect_allele", "reference_allele"]);
    const iWeight = colIdx(header, ["effect_weight",  "beta",       "weight", "or"]);
    const iChr    = colIdx(header, ["chr_name",       "chromosome", "chr",    "chrom"]);
    const iPos    = colIdx(header, ["chr_position",   "position",   "pos",    "bp"]);

    if (iRsid === -1) return { meta, catalogMap: new Map() };

    const catalogMap = new Map();
    for (let i = 1; i < dataLines.length; i++) {
        const cells = dataLines[i].split(sep).map(c => c.trim());
        if (cells.length <= iRsid) continue;

        const rsid = cells[iRsid].toLowerCase();
        if (!rsid) continue;

        // Parse weight — use null (not 1) as "missing" so we can distinguish
        // "not in file" from a real zero/negative beta.
        let weight = null;
        if (iWeight >= 0 && cells[iWeight] !== "" && cells[iWeight] !== undefined) {
            const parsed = parseFloat(cells[iWeight]);
            if (!isNaN(parsed)) weight = parsed;
        }

        catalogMap.set(rsid, {
            effect_allele: iEA  >= 0 ? (cells[iEA]  || "") : "",
            other_allele:  iOA  >= 0 ? (cells[iOA]  || "") : "",
            effect_weight: weight,
            chr:           iChr >= 0 ? (cells[iChr] || "") : "",
            pos:           iPos >= 0 ? (cells[iPos] || "") : ""
        });
    }

    return { meta, catalogMap };
}

/**
 * Read a local file using jamovi's fs API.
 * Returns Promise<string>.
 */
function readFileAsync(path) {
    if (window.jamovi && window.jamovi.fs && typeof window.jamovi.fs.readFile === "function") {
        return window.jamovi.fs.readFile(path);
    }
    return Promise.reject(new Error("jamovi fs API not available"));
}

// ─────────────────────────────────────────────────────────────────────────────
// Module-level cache: last successfully parsed catalog, keyed by path+sep.
// Avoids re-reading the file when only snpCols changes.
// ─────────────────────────────────────────────────────────────────────────────
let _catalogCache = null;  // { cacheKey, meta, catalogMap }

// ─────────────────────────────────────────────────────────────────────────────
// Exported event handlers
// ─────────────────────────────────────────────────────────────────────────────

const events = {

    // ── 1. Required by jamovi for VariableSupplier; intentionally a no-op ──
    update: function(ui) {
        // Do NOT call onChange_snpCols here — fires before snpCols is
        // populated and would reset the grid on every load.
    },

    // ── 2. Sync grid rows when the SNP column list changes ──────────────────
    onChange_snpCols: function(ui) {
        const vars = ui.snpCols.value() || [];
        const grid = ui.snpGrid.value()  || [];

        // Preserve existing rows (keeps user edits), drop deselected SNPs
        let newGrid = grid.filter(row => vars.includes(row.rsid));

        // Add skeleton rows for newly added columns
        vars.forEach(v => {
            if (!newGrid.some(row => row.rsid === v)) {
                newGrid.push(blankRow(v));
            }
        });

        // Follow the variable-list order
        newGrid.sort((a, b) => vars.indexOf(a.rsid) - vars.indexOf(b.rsid));

        ui.snpGrid.setValue(newGrid);

        // Fill catalog data for any new rows using the cached catalog
        events._mergeFromCache(ui);
    },

    // ── 3. (Re)load catalog file when path or delimiter changes ─────────────
    onChange_weightsPath: function(ui) {
        _catalogCache = null;   // invalidate cache — path has changed
        events._applyCatalogToGrid(ui);
    },

    onChange_weightsSep: function(ui) {
        _catalogCache = null;   // delimiter change may produce different parse
        events._applyCatalogToGrid(ui);
    },

    // ── Internal: load file, parse, cache, merge into grid ──────────────────
    _applyCatalogToGrid: function(ui) {
        const path = (ui.weightsPath.value() || "").trim();
        if (!path) return;

        const sepOption = ui.weightsSep.value() || "auto";
        const cacheKey  = path + "|" + sepOption;

        // Use cache if valid
        if (_catalogCache && _catalogCache.cacheKey === cacheKey) {
            events._mergeFromCache(ui);
            return;
        }

        readFileAsync(path).then(function(text) {
            const lines = text.split(/\r?\n/);
            const sep   = (sepOption === "auto")  ? detectSep(lines.join("\n"))
                        : (sepOption === "comma") ? ","
                        : "\t";

            const { meta, catalogMap } = parseCatalogText(text, sep);
            if (catalogMap.size === 0) return;

            // Store in cache
            _catalogCache = { cacheKey, meta, catalogMap };

            events._mergeFromCache(ui);

        }).catch(function(err) {
            // File not readable — R will show a proper error in the results panel
            console.warn("snpPGS: could not read weights file:", err);
        });
    },

    // ── Internal: apply cached catalog data to the current grid ─────────────
    //
    // Merge rules (designed to preserve user edits):
    //   effect_allele / other_allele: fill only if the current value is empty
    //   effect_weight: fill only if the row has never been set from a catalog
    //                  (i.e. effect_weight === null, the blankRow sentinel)
    //                  A user who has manually typed a weight keeps it.
    //   chr / pos: always take from catalog (read-only display fields)
    //   matched: set true when a catalog entry is found
    // ────────────────────────────────────────────────────────────────────────
    _mergeFromCache: function(ui) {
        if (!_catalogCache) return;

        const { catalogMap } = _catalogCache;
        const grid = ui.snpGrid.value() || [];
        if (grid.length === 0) return;

        const updatedGrid = grid.map(row => {
            const entry = catalogMap.get(row.rsid) ||
                          catalogMap.get((row.rsid || "").toLowerCase());

            if (!entry) return row;   // no catalog match — keep as-is

            return {
                rsid:          row.rsid,
                effect_allele: row.effect_allele || entry.effect_allele,
                other_allele:  row.other_allele  || entry.other_allele,
                // null sentinel = weight was never loaded; replace with catalog value.
                // Any other value (including 0, negative) = user has set it; keep it.
                effect_weight: (row.effect_weight === null)
                                   ? entry.effect_weight
                                   : row.effect_weight,
                chr:           entry.chr || row.chr,
                pos:           entry.pos || row.pos,
                matched:       true
            };
        });

        ui.snpGrid.setValue(updatedGrid);
    }
};

module.exports = events;
