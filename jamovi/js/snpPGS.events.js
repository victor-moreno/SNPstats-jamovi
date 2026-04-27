"use strict";

// ─────────────────────────────────────────────────────────────────────────────
// snpPGS.events.js
//
// Handles UI synchronisation only.  File reading has been removed because
// jamovi's JS sandbox does NOT expose a reliable file-system API
// (window.jamovi.fs.readFile does not exist).  The PGS Catalog file is
// parsed exclusively in R (.buildWeightTable branch b), which writes the
// resulting rows back to the snpGrid option via self$options$.set().
// Subsequent R runs then enter branch (a) and honour any user edits.
//
// Responsibilities:
//   1. onChange_snpCols  — keeps snpGrid rows in sync with selected SNP
//      columns: removes deselected rows, adds blank skeletons for new ones,
//      and preserves any edits already made to existing rows.
//   2. onChange_weightsPath / onChange_weightsSep — no-ops; R re-runs
//      automatically when these options change and repopulates the grid.
//   3. update — required by jamovi for VariableSupplier; intentionally empty.
// ─────────────────────────────────────────────────────────────────────────────

// ── Helpers ──────────────────────────────────────────────────────────────────

/**
 * Build a blank grid row for a SNP column that has no catalog entry yet.
 * effect_weight is set to null (sentinel = "not yet from catalog").
 * Note: jamovi may serialise null→1 via the YAML default; R compensates
 * by treating (matched=false, weight=1) as unset.
 */
function blankRow(colName) {
    return {
        rsid:          String(colName || ""),
        effect_allele: "",
        other_allele:  "",
        effect_weight: null,
        chr:           "",
        pos:           "",
        matched:       false
    };
}

/**
 * Safely extract a plain string from an rsid cell value.
 * jamovi VariableLabel cells return objects like { name: "rs123", ... }.
 * Plain Label / String cells return the string directly.
 */
function rsidString(val) {
    if (val == null) return "";
    if (typeof val === "string") return val;
    if (typeof val === "object") {
        return String(val.name || val.label || val.value || "");
    }
    return String(val);
}

// ─────────────────────────────────────────────────────────────────────────────
// Exported event handlers
// ─────────────────────────────────────────────────────────────────────────────

const events = {

    // ── 1. Required by jamovi for VariableSupplier; intentionally a no-op ──
    update: function(ui) {
        // Do NOT sync snpCols here — fires before snpCols is populated
        // and would reset the grid on every load.
    },

    // ── 2. Sync grid rows when the SNP column list changes ──────────────────
    //
    // Rules:
    //   • Drop rows for SNPs that are no longer selected.
    //   • Add blank skeleton rows for newly added SNPs.
    //   • Preserve any edits the user has already made to existing rows.
    //   • Follow the variable-list order.
    //
    // We do NOT attempt to read the weights file here — that is R's job.
    // ────────────────────────────────────────────────────────────────────────
    onChange_snpCols: function(ui) {
        const vars = ui.snpCols.value() || [];
        const grid = ui.snpGrid.value()  || [];

        // Keep existing rows that are still selected (preserves user edits)
        let newGrid = grid.filter(row => vars.includes(rsidString(row.rsid)));

        // Add skeleton rows for newly added columns
        vars.forEach(v => {
            if (!newGrid.some(row => rsidString(row.rsid) === v)) {
                newGrid.push(blankRow(v));
            }
        });

        // Follow the variable-list order
        newGrid.sort((a, b) =>
            vars.indexOf(rsidString(a.rsid)) - vars.indexOf(rsidString(b.rsid)));

        ui.snpGrid.setValue(newGrid);
        // R will re-run automatically and merge catalog data into the grid
        // via self$options$.set(snpGrid = ...) in .buildWeightTable.
    },

    // ── 3. Weights-file options changed — R re-runs and repopulates grid ────
    //
    // No JS action needed.  Defining these handlers prevents jamovi from
    // raising "unknown event" warnings if the UI yaml wires them up.
    // ────────────────────────────────────────────────────────────────────────
    onChange_weightsPath: function(ui) { /* R handles file parsing */ },
    onChange_weightsSep:  function(ui) { /* R handles file parsing */ }

};

module.exports = events;
