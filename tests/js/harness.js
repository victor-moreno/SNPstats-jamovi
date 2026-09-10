'use strict';

// A jamovi-shaped stand-in for snpimport.js to run against.
//
// The regex tests in test-groupfiles.js check what the source says; this runs
// it. It exists because the Load button shipped dead twice — once from an
// exception while building it, once from a stray reference in its own body —
// and neither could be seen without actually injecting the panel and clicking.
//
// What is faked: just enough jQuery for the injection (wrap/before/after/on),
// a document whose file input can be handed files, and jamovi's option
// accessors. Files are real: Node's File over the fixture bytes, so
// File.slice/stream/arrayBuffer behave as the browser's do.

const fs = require('fs');

// ── a minimal jQuery ────────────────────────────────────────────────────────

function makeEl(tag, cls) {
    return { tag: tag, cls: cls || '', parent: null, children: [],
             handlers: {}, props: {}, style: {} };
}

function jq(html) {
    if (typeof html !== 'string') return wrap([html]);   // jQuery(element)
    const t = /^<(\w+)/.exec(html);
    const c = /class="([^"]*)"/.exec(html);
    return wrap([makeEl(t ? t[1] : 'div', c ? c[1] : '')]);
}

function insertAt(parent, node, index) {
    if (node.parent) {
        const at = node.parent.children.indexOf(node);
        if (at >= 0) node.parent.children.splice(at, 1);
    }
    node.parent = parent;
    parent.children.splice(index, 0, node);
}

function matches(el, sel) {
    return sel.charAt(0) === '.' ? el.cls.split(/\s+/).includes(sel.slice(1))
                                 : el.tag === sel;
}

function wrap(els) {
    const o = {
        length: els.length,
        els: els,
        prop: function(k, v) { els.forEach(e => { e.props[k] = v; }); return o; },
        css:  function(k, v) {
            els.forEach(e => {
                if (typeof k === 'object') Object.assign(e.style, k);
                else e.style[k] = v;
            });
            return o;
        },
        prev: function(sel) { return sibling(-1, sel); },
        next: function(sel) { return sibling(+1, sel); },
        parent: function() {
            const p = els[0] && els[0].parent;
            return p ? wrap([p]) : wrap([]);
        },
        text: function(t) { els.forEach(e => { e.text = t; }); return o; },
        attr: function(k, v) { els.forEach(e => { e.props[k] = v; }); return o; },
        wrap: function($w) {
            const w = $w.els[0], e = els[0];
            const p = e.parent, at = p ? p.children.indexOf(e) : 0;
            if (p) { p.children.splice(at, 1, w); w.parent = p; }
            e.parent = w; w.children.push(e);
            return o;
        },
        before: function($x) { return place($x, 0); },
        after:  function($x) { return place($x, 1); },
        on: function(ev, fn) {
            els.forEach(e => { (e.handlers[ev] = e.handlers[ev] || []).push(fn); });
            return o;
        }
    };
    function sibling(dir, sel) {
        const e = els[0];
        if (!e || !e.parent) return wrap([]);
        const at = e.parent.children.indexOf(e) + dir;
        const s = e.parent.children[at];
        return (s && (!sel || matches(s, sel))) ? wrap([s]) : wrap([]);
    }
    function place($x, offset) {
        const e = els[0];
        if (!e.parent) throw new Error('cannot insert next to a detached node');
        insertAt(e.parent, $x.els[0], e.parent.children.indexOf(e) + offset);
        return o;
    }
    o.constructor = jq;
    els.forEach((e, i) => { o[i] = e; });
    return o;
}

function find(root, cls, out) {
    out = out || [];
    if (root.cls && root.cls.split(/\s+/).includes(cls)) out.push(root);
    root.children.forEach(c => find(c, cls, out));
    return out;
}

// ── the panel ───────────────────────────────────────────────────────────────

const FIELDS = ['genoFilename', 'snpListFilename', 'covFilename', 'loadStatus'];

function makeUi(initial) {
    const values = Object.assign({
        genoContent: '', variantContent: '', sampleContent: '',
        snpListContent: '', snpListText: '', covContent: '',
        genoFilename: '', snpListFilename: '', covFilename: '',
        loadStatus: '', loadProblem: '', sourceFormat: 'bed',
        sourceDims: '',
        covIdCol: '', dosage: false, openNew: false
    }, initial || {});

    const root = makeEl('div', 'panel');
    const ui = { root: root, values: values };

    // jamovi's shape, not a convenient one: the handler is handed a resource
    // per option — a control where the .u.yaml declares one, a bare
    // `_hiddenOption` where it does not — and each carries value()/setValue().
    // The view object itself has no getOptionValue/setOptionValue, so offering
    // them here let the panel pass on a path that does not exist in jamovi.
    Object.keys(values).forEach(n => {
        ui[n] = { value: () => values[n],
                  setValue: (v) => { values[n] = v; } };
    });
    FIELDS.forEach(n => {
        const el = makeEl('input', n);
        insertAt(root, el, root.children.length);
        ui[n].$input = wrap([el]);
    });
    return ui;
}

// The file dialog: createElement('input') hands back something whose click()
// fires the change listener with whatever the test queued.
function installDocument(queue) {
    global.document = {
        createElement: function() {
            const fi = { type: '', accept: '', style: {}, files: null,
                         listeners: [] };
            fi.addEventListener = function(ev, fn) {
                if (ev === 'change') fi.listeners.push(fn);
            };
            fi.click = function() {
                fi.files = queue.shift() || [];
                fi.listeners.forEach(fn => fn());
            };
            return fi;
        },
        body: { appendChild: function() {}, removeChild: function() {} }
    };
}

function click(ui, cls) {
    const btns = find(ui.root, cls);
    if (btns.length === 0) throw new Error('no ' + cls + ' button in the panel');
    // A button nobody can see is as dead as one with no handler: the ✕ buttons
    // shipped display:none, waiting for a view_updated that never came.
    if (btns[0].style.display === 'none')
        throw new Error(cls + ' is display:none — the user cannot click it');
    const hs = btns[0].handlers.click || [];
    if (hs.length === 0)
        throw new Error(cls + ' has no click handler — the button is dead');
    hs.forEach(fn => fn({ preventDefault() {}, stopPropagation() {} }));
}

function file(path, name) {
    return new File([fs.readFileSync(path)], name || path.split('/').pop());
}

// The click handlers are fire-and-forget by design — a DOM handler cannot be
// awaited — so a test waits on the status line instead, which is the same
// signal the user gets. File reads are real I/O here, so a microtask drain is
// not enough; this polls a timer.
async function until(fn, what, ms) {
    const deadline = Date.now() + (ms || 120000);
    while (Date.now() < deadline) {
        if (fn()) return;
        await new Promise(r => setTimeout(r, 5));
    }
    throw new Error('timed out waiting for ' + (what || 'the condition'));
}

const settled = (ui) => until(
    () => ui.values.loadStatus !== '' && !/…$/.test(ui.values.loadStatus),
    'the load to finish (status: "' + ui.values.loadStatus + '")');

module.exports = { jq, wrap, makeEl, makeUi, installDocument, click, find,
                   file, until, settled };
