'use strict';

// snpPGS.js — view event handlers for the PGS analysis.
//
// Naming matters here. The jamovi UI compiler emits
//     this.handlers = require('./snpPGS')
// into the generated snpPGS.src.js *because this file is named after the
// analysis*; view-level handlers (view_loaded / view_updated) and bare
// per-control handler names are then resolved against it. That is the same
// convention jmv's own anova.js uses. A `<name>.events.js` file is a different
// mechanism — it is only loaded when snpPGS.u.yaml declares an explicit
// `events:` block with `./snpPGS.events::<export>` references, which it does
// not. An empty snpPGS.events.js used to sit alongside this file and reads as
// if it were the live one; it was removed to stop that confusion.

module.exports = {

    view_updated: function(ui, event) {
        _injectBrowseButton(ui);
    },

    view_loaded: function(ui, event) {
        _injectBrowseButton(ui);
    }

};

// Set an option by name. weightsContent is a hidden option (hidden: true in
// snpPGS.a.yaml) so it renders no control. Try, in order: a real control's
// setValue (visible options like weightsFilename), the view-level
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

// ── File-browse button for the weights file ──────────────────────────────────
// The TextBox shows the chosen file's *name* only, and is disabled: there is no
// path option to type into any more. A typed path would be saved into the .omv
// and re-read on whoever opened it, so picking a file is the only way in and it
// embeds the file's bytes rather than its location.
//
// `disabled` rather than a hand-styled grey: a natively disabled input is
// greyed by the platform's own stylesheet, so it stays legible in jamovi's dark
// theme. It is re-asserted on every view_updated, so a jamovi refresh cannot
// quietly re-enable it. Disabling blocks user typing only — _embedFile still
// writes the file name into the option and into the field.
//
// Layout: [📁] [ filename ]. The button comes first because it is the control —
// the field is a passive readout of what the button loaded — and the field is
// hidden entirely until a file has actually been chosen, so the empty state is
// just the label and the button.

function _injectBrowseButton(ui) {

    var ctrl = ui.weightsFilename;
    if (!ctrl) return;

    var $input = ctrl.$input;
    if (!$input || $input.length === 0) return;

    $input.prop('readonly', true);
    $input.prop('disabled', true);
    $input.css('cursor', 'default');

    // The button sits BEFORE the input, so the already-injected test is prev(),
    // not next(). Getting this wrong re-injects a button on every view_updated.
    if ($input.prev('.pgs-browse-btn').length === 0) {

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
        $input.before($btn);

        $btn.on('click', function(e) {
            e.preventDefault();
            e.stopPropagation();

            // Always read the file's bytes in the browser (works on desktop and
            // in cloud) and embed them. No local path is ever written into an
            // option.
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

    // Outside the injection guard: runs on every view_updated, so the field
    // appears as soon as a file is chosen and after a saved analysis restores.
    _syncNameVisibility(ui);
}

// Show the file-name field only when there is a name to show.
function _syncNameVisibility(ui) {
    var ctrl = ui.weightsFilename;
    if (!ctrl || !ctrl.$input || ctrl.$input.length === 0) return;
    var name = _getName(ui);
    ctrl.$input.css('display', (name && String(name).length > 0) ? '' : 'none');
}

// Current value of weightsFilename. Read from the options model rather than the
// input's DOM value: on the first view_loaded of a restored analysis the field
// may not have been populated yet, and reading the DOM there would hide a field
// that does have a file.
function _getName(ui) {
    var ctrl = ui.weightsFilename;
    if (ctrl && typeof ctrl.value === 'function') {
        var v = ctrl.value();
        if (v !== undefined && v !== null) return v;
    }
    if (typeof ui.getOptionValue === 'function') {
        var v2 = ui.getOptionValue('weightsFilename');
        if (v2 !== undefined && v2 !== null) return v2;
    }
    if (typeof ui.getOption === 'function') {
        var opt = ui.getOption('weightsFilename');
        if (opt && typeof opt.value === 'function') return opt.value();
    }
    return (ctrl && ctrl.$input) ? ctrl.$input.val() : '';
}

// Read `file` as base64 and stash it into the carrier options so the R backend
// can decode it (see R/snpPGS.b.R .weightsRawLines).
function _embedFile(ui, file) {
    var reader = new FileReader();
    reader.onload = function() {
        var result = reader.result || '';
        var comma = result.indexOf(',');
        var b64 = comma >= 0 ? result.slice(comma + 1) : result;
        _setName(ui, file.name);
        _setOpt(ui, 'weightsContent', b64);
    };
    reader.onerror = function() {
        // Report the failure; never fall back to writing a filesystem path.
        _setOpt(ui, 'weightsContent', '');
        _setName(ui, '(could not read ' + file.name + ')');
    };
    reader.readAsDataURL(file);
}

// Set the displayed file name. The option is what the R backend reads; the
// direct field write is belt-and-braces, so the name still appears even if
// jamovi were to skip refreshing a disabled control's DOM. Reveal it here too
// rather than waiting for the next view_updated.
function _setName(ui, name) {
    _setOpt(ui, 'weightsFilename', name);
    var ctrl = ui.weightsFilename;
    if (ctrl && ctrl.$input && ctrl.$input.length > 0) ctrl.$input.val(name);
    _syncNameVisibility(ui);
}
