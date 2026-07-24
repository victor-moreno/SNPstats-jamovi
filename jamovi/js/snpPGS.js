'use strict';

module.exports = {

    view_updated: function(ui, event) {
        _hideCarrierControls(ui);
        _injectBrowseButton(ui);
    },

    view_loaded: function(ui, event) {
        _hideCarrierControls(ui);
        _injectBrowseButton(ui);
    }

};

// ── Hide the option "carrier" TextBoxes ───────────────────────────────────────
// weightsFilename / weightsContent are option carriers, not user inputs: the
// browse button writes the picked file's name and its base64 contents into them
// so the file travels to the R engine even in jamovi cloud (where the engine
// runs on a different machine than the browser and a local path is useless).
// They must exist as controls (so ui.<name>.setValue works) but never show.

function _hideCarrierControls(ui) {
    ['weightsFilename', 'weightsContent'].forEach(function(name) {
        var ctrl = ui[name];
        if (!ctrl || !ctrl.$input || ctrl.$input.length === 0) return;
        ctrl.$input.css({ display: 'none' });
        var $lbl = ctrl.$input.parent();
        if ($lbl && $lbl.length) $lbl.css({ display: 'none' });
    });
}

// ── File-browse button for the weights path TextBox ───────────────────────────

function _injectBrowseButton(ui) {

    var ctrl = ui.weightsPath;
    if (!ctrl) return;

    var $input = ctrl.$input;
    if (!$input || $input.length === 0) return;

    if ($input.next('.pgs-browse-btn').length > 0) return;

    var jq = $input.constructor;

    var $row = jq('<div class="pgs-browse-row"></div>').css({
        display: 'flex', alignItems: 'center', width: '100%', gap: '4px'
    });

    var $btn = jq('<button type="button" class="pgs-browse-btn" title="Browse…">📁</button>').css({
        flexShrink: '0', cursor: 'pointer', padding: '1px 7px',
        fontSize: '14px', lineHeight: '1.4', border: '1px solid #bbb',
        borderRadius: '3px', background: '#f0f0f0', whiteSpace: 'nowrap'
    });

    $input.css({ flex: '1 1 auto', minWidth: 0 });
    $input.wrap($row);
    $input.after($btn);

    $btn.on('click', function(e) {
        e.preventDefault();
        e.stopPropagation();

        // Always read the file's bytes in the browser (works on desktop and in
        // cloud) and embed them, rather than passing a path the R engine may not
        // be able to see.
        var fileInput = document.createElement('input');
        fileInput.type = 'file';
        fileInput.accept = '.csv,.tsv,.txt,.gz';
        fileInput.style.display = 'none';
        document.body.appendChild(fileInput);

        fileInput.addEventListener('change', function() {
            var file = fileInput.files && fileInput.files[0];
            if (file) _embedFile(ui, file);
            document.body.removeChild(fileInput);
        });

        fileInput.click();
    });
}

// Read `file` as base64 and stash it into the carrier options so the R backend
// can decode it (see R/snpPGS.b.R .weightsRawLines). weightsPath is set to the
// name for display; a manually-typed path (no content) still works on desktop.
function _embedFile(ui, file) {
    var reader = new FileReader();
    reader.onload = function() {
        var result = reader.result || '';
        var comma = result.indexOf(',');
        var b64 = comma >= 0 ? result.slice(comma + 1) : result;
        if (ui.weightsFilename) ui.weightsFilename.setValue(file.name);
        if (ui.weightsContent)  ui.weightsContent.setValue(b64);
        if (ui.weightsPath)     ui.weightsPath.setValue(file.name);
    };
    reader.onerror = function() {
        if (ui.weightsPath) ui.weightsPath.setValue(file.path || file.name);
    };
    reader.readAsDataURL(file);
}
