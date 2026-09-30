'use strict';

// snpPGS.js — cosmetic-only view handler.
//
// weightsFile is jamovi's native File option (a FileSelector control); no
// data/option logic lives here at all -- doing that would defeat the point
// of switching to it (see snpPGS.b.R / snpPGS.a.yaml). This file only
// repositions the picked file beside the Browse button and gives it a boxed
// look, to match SNPstats' other file pickers (the still-custom ones in
// snpimport.js) instead of FileSelector's own default layout (button, then
// a plain list below it).
//
// Targets jamovi's internal FileSelector DOM
// (client/analysisui/fileselector.ts: 'body' and 'list' are public fields
// on the live control instance, exposed here as ui.weightsFile.body/.list) —
// this is not a documented/stable API, so a future jamovi client update
// could silently drop the styling. It cannot break the analysis itself:
// nothing here touches data or options, only cosmetics.

module.exports = {
    view_updated: function(ui) { _styleFileSelector(ui); },
    view_loaded:  function(ui) { _styleFileSelector(ui); }
};

var STYLE_ID = 'snppgs-inline-fileselector-css';
var ROW_CLASS = 'snppgs-inline-fs';

// Translucent neutral grey instead of hex greys: on jamovi's white panel these
// render as the old #bbb / #f0f0f0 exactly, and on a darker background (jamovi
// 28.4 has no dark theme, but a future one) they stay a subtle tint instead of
// a light box. OS high-contrast mode replaces them either way. Same values as
// snpimport.js.
var BORDER = 'rgba(128, 128, 128, 0.535)';
var FILL = 'rgba(128, 128, 128, 0.118)';

function _injectCss() {
    if (document.getElementById(STYLE_ID)) return;
    var style = document.createElement('style');
    style.id = STYLE_ID;
    style.textContent =
        '.' + ROW_CLASS + ' { display: flex; flex-direction: row; ' +
        'align-items: center; flex-wrap: wrap; gap: 6px; }' +
        '.' + ROW_CLASS + ' .jmv-file-selector-list { flex: 1 1 auto; min-width: 0; }' +
        '.' + ROW_CLASS + ' .jmv-file-selector-item { width: 100%; box-sizing: border-box; ' +
        'border: 1px solid ' + BORDER + '; border-radius: 3px; padding: 3px 8px; ' +
        'background: ' + FILL + '; min-height: 14px; }';
    document.head.appendChild(style);
}

// Re-run on every view_updated: FileSelector redraws its own list on every
// value change (client/analysisui/fileselector.ts:update()), but a class-
// based CSS rule (rather than per-element inline styles) keeps applying to
// whatever it redraws, so this only needs to make sure the class is present.
function _styleFileSelector(ui) {
    var ctrl = ui.weightsFile;
    if (!ctrl || !ctrl.body) return;
    _injectCss();
    ctrl.body.classList.add(ROW_CLASS);
}
