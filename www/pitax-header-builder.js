// PITAX Custom FASTA Header Builder: drag/reorder UI + browser profile persistence.
// Profiles only (never per-record metadata values) are stored in localStorage.
(function () {
  'use strict';

  var STORAGE_KEY = 'PITAX_FASTA_HEADER_PROFILES_V2';
  var state = {
    library: [],
    activeFields: [],
    draft: null,
    dirty: false,
    persistence: {},
    structureIssues: []
  };

  function shinyReady() {
    return window.Shiny && typeof Shiny.setInputValue === 'function';
  }

  function send(name, value) {
    if (!shinyReady()) return;
    Shiny.setInputValue(name, value, { priority: 'event' });
  }

  function reportSave(name, ok, error) {
    send('pitax_header_browser_save_result', {
      ok: Boolean(ok),
      name: name || '',
      error: error || '',
      nonce: Date.now()
    });
  }

  function restoreProfiles() {
    var saved = null;
    try {
      var raw = window.localStorage.getItem(STORAGE_KEY);
      if (raw) saved = JSON.parse(raw);
    } catch (e) {
      reportSave('', false, 'Browser storage unavailable: ' + String(e.message || e));
    }
    send('pitax_header_browser_profiles_restore', {
      data: saved,
      nonce: Date.now()
    });
  }

  var saveHandlerRegistered = false;
  var restoreSent = false;

  function ensureSaveHandler() {
    if (saveHandlerRegistered) return true;
    if (!window.Shiny || typeof Shiny.addCustomMessageHandler !== 'function') return false;
    Shiny.addCustomMessageHandler('pitaxStoreHeaderProfiles', function (payload) {
      var name = (payload && (payload.active_name || payload.active_profile_id)) || '';
      try {
        if (!payload || Number(payload.schema) < 2 || !Array.isArray(payload.profiles)) {
          throw new Error('Invalid profile payload from server.');
        }
        var raw = JSON.stringify(payload);
        window.localStorage.setItem(STORAGE_KEY, raw);
        if (window.localStorage.getItem(STORAGE_KEY) !== raw) {
          throw new Error('Browser did not retain saved profiles.');
        }
        // Read-back verification for honest save acknowledgement.
        var verified = JSON.parse(window.localStorage.getItem(STORAGE_KEY));
        if (!verified || !Array.isArray(verified.profiles) ||
            verified.profiles.length !== payload.profiles.length) {
          throw new Error('Read-back verification failed.');
        }
        reportSave(name, true, '');
      } catch (e) {
        reportSave(name, false, 'Could not save to this browser: ' + String(e.message || e));
      }
    });
    saveHandlerRegistered = true;
    return true;
  }

  function onConnected() {
    if (!ensureSaveHandler() || restoreSent) return;
    restoreSent = true;
    restoreProfiles();
  }

  function installConnectionHooks() {
    if (window.jQuery) {
      window.jQuery(document).on('shiny:connected', onConnected);
      window.jQuery(document).on('shiny:disconnected', function () {
        restoreSent = false;
      });
    } else {
      document.addEventListener('shiny:connected', onConnected);
      document.addEventListener('shiny:disconnected', function () {
        restoreSent = false;
      });
    }
  }

  installConnectionHooks();
  var attempts = 0;
  var bootTimer = window.setInterval(function () {
    attempts += 1;
    if (ensureSaveHandler()) {
      var app = window.Shiny && window.Shiny.shinyapp;
      var connected = app && app.$socket && app.$socket.readyState === 1;
      if (connected) onConnected();
    }
    if ((saveHandlerRegistered && restoreSent) || attempts > 300) {
      window.clearInterval(bootTimer);
    }
  }, 100);

  function byId(id) { return document.getElementById(id); }

  function fieldMeta(fieldId) {
    for (var i = 0; i < state.library.length; i++) {
      if (state.library[i].field_id === fieldId) return state.library[i];
    }
    return {
      field_id: fieldId,
      display: fieldId,
      canonical_tag: fieldId,
      destination: 'Custom / Extras',
      destination_class: 'extras',
      value_source: 'user',
      description: '',
      group: 'C'
    };
  }

  function emitOrder() {
    send('pitax_header_builder_order', {
      fields: state.activeFields.slice(),
      nonce: Date.now()
    });
  }

  function renderLibrary() {
    var root = byId('fhb_library');
    if (!root) return;
    var groups = { A: [], B: [], C: [] };
    state.library.forEach(function (f) {
      if (f.field_id === 'final_name') return;
      var g = f.group || 'C';
      if (!groups[g]) groups[g] = [];
      groups[g].push(f);
    });
    var html = '';
    [['A', 'Primary biological / record'],
     ['B', 'Experimental metadata (Notes downstream)'],
     ['C', 'Custom / Extras']].forEach(function (pair) {
      var key = pair[0], title = pair[1];
      if (!groups[key] || !groups[key].length) return;
      html += '<div class="fhb-group-label">' + title + '</div>';
      groups[key].forEach(function (f) {
        var used = state.activeFields.indexOf(f.field_id) >= 0;
        var cls = 'fhb-field-card';
        if (used) cls += ' is-used';
        if (f.destination_class === 'notes_folded') cls += ' is-notes';
        else if (f.destination_class === 'extras') cls += ' is-extras';
        else cls += ' is-first';
        html += '<div class="' + cls + '" draggable="' + (!used) + '" data-field-id="' +
          f.field_id + '" tabindex="0" role="button" aria-label="Add ' + f.display + '">' +
          '<div class="fhb-field-title">' + f.display + '</div>' +
          '<div class="fhb-field-tag">' + (f.canonical_tag || '') + '</div>' +
          '<div class="fhb-field-dest">' + (f.destination || '') + '</div>' +
          '</div>';
      });
    });
    root.innerHTML = html || '<div class="fhb-panel-hint">No fields available.</div>';
    root.querySelectorAll('.fhb-field-card').forEach(function (el) {
      el.addEventListener('dragstart', function (ev) {
        if (el.classList.contains('is-used')) {
          ev.preventDefault();
          return;
        }
        ev.dataTransfer.setData('text/plain', el.getAttribute('data-field-id'));
        ev.dataTransfer.effectAllowed = 'copy';
      });
      el.addEventListener('dblclick', function () {
        addField(el.getAttribute('data-field-id'));
      });
      el.addEventListener('keydown', function (ev) {
        if (ev.key === 'Enter' || ev.key === ' ') {
          ev.preventDefault();
          addField(el.getAttribute('data-field-id'));
        }
      });
    });
  }

  function addField(fieldId) {
    if (!fieldId || state.activeFields.indexOf(fieldId) >= 0) return;
    send('pitax_header_builder_add', { field_id: fieldId, nonce: Date.now() });
  }

  function removeField(fieldId) {
    send('pitax_header_builder_remove', { field_id: fieldId, nonce: Date.now() });
  }

  function moveField(fieldId, delta) {
    var idx = state.activeFields.indexOf(fieldId);
    if (idx < 0) return;
    var next = idx + delta;
    if (next < 0 || next >= state.activeFields.length) return;
    var copy = state.activeFields.slice();
    var tmp = copy[idx];
    copy[idx] = copy[next];
    copy[next] = tmp;
    state.activeFields = copy;
    renderConstruction();
    emitOrder();
  }

  function renderConstruction() {
    var root = byId('fhb_construction');
    if (!root) return;
    var html = '<div class="fhb-primary-chip"><span>Primary identifier</span><span>{final_name}</span></div>';
    state.activeFields.forEach(function (fid, i) {
      var meta = fieldMeta(fid);
      var sourceCls = meta.value_source === 'system' ? 'is-system' : 'is-user';
      var sourceLabel = meta.value_source === 'system' ? 'auto' : 'manual';
      html += '<div class="fhb-chip" draggable="true" data-field-id="' + fid + '" data-index="' + i + '">' +
        '<div class="fhb-chip-meta">' +
        '<div class="fhb-field-title">' + meta.display +
        '<span class="fhb-source ' + sourceCls + '">' + sourceLabel + '</span></div>' +
        '<div class="fhb-field-tag">' + (meta.canonical_tag || fid) + '</div>' +
        '</div>' +
        '<div class="fhb-chip-actions">' +
        '<button type="button" data-act="up" title="Move up" aria-label="Move up">↑</button>' +
        '<button type="button" data-act="down" title="Move down" aria-label="Move down">↓</button>' +
        '<button type="button" data-act="remove" title="Remove" aria-label="Remove">×</button>' +
        '</div></div>';
    });
    if (!state.activeFields.length) {
      html += '<div class="fhb-panel-hint">Drag fields here, double-click a library card, or press Enter on a focused card.</div>';
    }
    root.innerHTML = html;

    root.querySelectorAll('.fhb-chip').forEach(function (chip) {
      chip.addEventListener('dragstart', function (ev) {
        ev.dataTransfer.setData('text/plain', 'reorder:' + chip.getAttribute('data-field-id'));
        ev.dataTransfer.effectAllowed = 'move';
      });
      chip.addEventListener('dragover', function (ev) {
        ev.preventDefault();
      });
      chip.addEventListener('drop', function (ev) {
        ev.preventDefault();
        var payload = ev.dataTransfer.getData('text/plain') || '';
        var target = chip.getAttribute('data-field-id');
        if (payload.indexOf('reorder:') === 0) {
          var src = payload.slice(8);
          var from = state.activeFields.indexOf(src);
          var to = state.activeFields.indexOf(target);
          if (from < 0 || to < 0 || from === to) return;
          var copy = state.activeFields.slice();
          copy.splice(from, 1);
          copy.splice(to, 0, src);
          state.activeFields = copy;
          renderConstruction();
          emitOrder();
        } else if (payload) {
          addField(payload);
        }
      });
      chip.querySelectorAll('button').forEach(function (btn) {
        btn.addEventListener('click', function () {
          var act = btn.getAttribute('data-act');
          var fid = chip.getAttribute('data-field-id');
          if (act === 'up') moveField(fid, -1);
          if (act === 'down') moveField(fid, 1);
          if (act === 'remove') removeField(fid);
        });
      });
    });
  }

  function renderPreview(template) {
    var box = byId('fhb_preview_template');
    if (box) box.textContent = '>' + (template || '');
  }

  function renderValidation() {
    var box = byId('fhb_structure_validation');
    if (!box) return;
    var issues = state.structureIssues || [];
    if (!issues.length) {
      box.innerHTML = '<div class="ok">Profile structure is valid.</div>';
      return;
    }
    box.innerHTML = issues.map(function (msg) {
      var cls = String(msg).indexOf('Warning:') === 0 ? 'warn' : 'err';
      return '<div class="' + cls + '">' + msg + '</div>';
    }).join('');
  }

  function renderStatus() {
    var box = byId('fhb_persistence_status');
    if (!box) return;
    var p = state.persistence || {};
    var pills = [];
    if (state.dirty) pills.push('<span class="fhb-status-pill is-dirty">Unsaved draft</span>');
    else pills.push('<span class="fhb-status-pill">Draft matches saved profile</span>');
    if (p.browser_saved) pills.push('<span class="fhb-status-pill is-browser">Saved to this browser</span>');
    else pills.push('<span class="fhb-status-pill">Not saved to this browser</span>');
    if (p.project_included) pills.push('<span class="fhb-status-pill is-project">Included in project state</span>');
    if (p.portable_exported) pills.push('<span class="fhb-status-pill">Portable export available</span>');
    if (p.last_status) pills.push('<span class="fhb-status-pill">' + p.last_status + '</span>');
    box.innerHTML = pills.join('');
  }

  function applyState(payload) {
    if (!payload) return;
    state.library = payload.library || [];
    // library may arrive as object-of-columns from R; normalize to row objects.
    if (!Array.isArray(state.library) && state.library.field_id) {
      var rows = [];
      var n = state.library.field_id.length;
      for (var i = 0; i < n; i++) {
        rows.push({
          field_id: state.library.field_id[i],
          display: state.library.display[i],
          canonical_tag: state.library.canonical_tag[i],
          group: state.library.group[i],
          destination: state.library.destination[i],
          destination_class: state.library.destination_class[i],
          value_source: state.library.value_source[i],
          description: state.library.description[i]
        });
      }
      state.library = rows;
    }
    state.activeFields = (payload.active_fields || []).slice();
    state.draft = payload.draft || null;
    state.dirty = !!payload.dirty;
    state.persistence = payload.persistence || {};
    state.structureIssues = payload.structure_issues || [];
    renderLibrary();
    renderConstruction();
    renderPreview(payload.preview_template || '');
    renderValidation();
    renderStatus();
  }

  function bindConstructionDropZone() {
    var root = byId('fhb_construction');
    if (!root || root._fhbBound) return;
    root._fhbBound = true;
    root.addEventListener('dragover', function (ev) {
      ev.preventDefault();
      root.classList.add('is-dragover');
    });
    root.addEventListener('dragleave', function () {
      root.classList.remove('is-dragover');
    });
    root.addEventListener('drop', function (ev) {
      ev.preventDefault();
      root.classList.remove('is-dragover');
      var payload = ev.dataTransfer.getData('text/plain') || '';
      if (!payload || payload.indexOf('reorder:') === 0) return;
      addField(payload);
    });
  }

  if (window.jQuery) {
    window.jQuery(document).on('shiny:sessioninitialized', function () {
      bindConstructionDropZone();
    });
  }
  document.addEventListener('DOMContentLoaded', bindConstructionDropZone);

  if (window.Shiny && typeof Shiny.addCustomMessageHandler === 'function') {
    Shiny.addCustomMessageHandler('pitaxHeaderBuilderState', applyState);
  } else {
    var wait = window.setInterval(function () {
      if (window.Shiny && typeof Shiny.addCustomMessageHandler === 'function') {
        Shiny.addCustomMessageHandler('pitaxHeaderBuilderState', applyState);
        window.clearInterval(wait);
      }
    }, 50);
  }

  // Expose for Playwright tests.
  window.PITAXHeaderBuilder = {
    getState: function () { return state; },
    addField: addField,
    removeField: removeField
  };
})();
