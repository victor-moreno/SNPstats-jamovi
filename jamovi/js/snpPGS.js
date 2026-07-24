'use strict';

module.exports = {

    view_updated: function(ui, event) {
        _injectBrowseButton(ui);
    },

    view_loaded: function(ui, event) {
        _injectBrowseButton(ui);
    }

};

// Set an option by name. weightsContent / weightsFilename are hidden options
// (hidden: true in snpPGS.a.yaml) so they render no control. Try, in order:
// a real control's setValue (visible options like weightsPath), the view-level
// setOptionValue(name, value), then resolving the option object via getOption().
// One of these is available depending on the jamovi version; hidden options have
// no control so they go through setOptionValue/getOption.
function _setOpt(ui, name, value) {
    if (ui[name] && typeof ui[name].setValue === 'function') {
        ui[name].setValue(value);
        return;
    }
    if (typeof ui.setOptionValue === 'function') {
        ui.setOptionValue(name, value);
        return;
    }
    if (typeof ui.getOption === 'function') {
        var opt = ui.getOption(name);
        if (opt && typeof opt.setValue === 'function') opt.setValue(value);
    }
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

// Read `file` as base64 and stash it into the hidden carrier options so the R
// backend can decode it (see R/snpPGS.b.R .weightsRawLines). weightsPath is set
// to the name for display; a manually-typed path (no content) still works on
// desktop.
function _embedFile(ui, file) {
    var reader = new FileReader();
    reader.onload = function() {
        var result = reader.result || '';
        var comma = result.indexOf(',');
        var b64 = comma >= 0 ? result.slice(comma + 1) : result;
        _setOpt(ui, 'weightsFilename', file.name);
        _setOpt(ui, 'weightsContent', b64);
        _setOpt(ui, 'weightsPath', file.name);
    };
    reader.onerror = function() {
        _setOpt(ui, 'weightsPath', file.path || file.name);
    };
    reader.readAsDataURL(file);
}
