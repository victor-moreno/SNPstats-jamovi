'use strict';

module.exports = {

    view_updated: function(ui, event) {
        _injectBrowseButton(ui);
    },

    view_loaded: function(ui, event) {
        _injectBrowseButton(ui);
    }

};

function _injectBrowseButton(ui) {

    var ctrl = ui.weightsPath;
    if (!ctrl) return;

    // jamovi TextBox exposes $input directly; $el is not populated
    var $input = ctrl.$input;
    if (!$input || $input.length === 0) return;

    // Guard: only inject once
    if ($input.next('.pgs-browse-btn').length > 0) return;

    // jQuery is not available as $ in the jamovi sandbox — use the input's
    // own jQuery reference via its constructor, or fall back to vanilla DOM.
    var jq = $input.constructor;  // this IS jQuery — $input is already a jQuery object

    var $row = jq('<div class="pgs-browse-row"></div>').css({
        display: 'flex', alignItems: 'center', width: '100%', gap: '4px'
    });

    var $btn = jq('<button type="button" class="pgs-browse-btn" title="Browse\u2026">\uD83D\uDCC1</button>').css({
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

        // Try jamovi native file-picker first (desktop >= 1.6)
        var lm = (ui.view && (ui.view.layoutManager || (ui.view.model && ui.view.model.layoutManager)))
               || window.layoutManager;

        if (lm && typeof lm.openFileBrowser === 'function') {
            lm.openFileBrowser({ title: 'Select PGS weights file' }, function(result) {
                if (result && result.path) ctrl.setValue(result.path);
            });
            return;
        }

        // Fallback: vanilla DOM file input (Electron gives file.path)
        var fileInput = document.createElement('input');
        fileInput.type = 'file';
        fileInput.accept = '.csv,.tsv,.txt,.gz';
        fileInput.style.display = 'none';
        document.body.appendChild(fileInput);

        fileInput.addEventListener('change', function() {
            var file = fileInput.files && fileInput.files[0];
            if (file) ctrl.setValue(file.path || file.name);
            document.body.removeChild(fileInput);
        });

        fileInput.click();
    });
}
