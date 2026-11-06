// FASTA header profiles are portable layout definitions only.
// Browser localStorage preserves them across Shiny sessions in this browser.
// Per-read metadata values remain inside .sangerproject; never store those here.
(function () {
  'use strict';
  var key = 'PITAX_FASTA_HEADER_PROFILES_V1';

  function report(name, ok, error) {
    if (!window.Shiny || typeof Shiny.setInputValue !== 'function') return;
    Shiny.setInputValue('pitax_header_browser_save_result', {
      ok: Boolean(ok), name: name || '', error: error || '', nonce: Date.now()
    }, { priority: 'event' });
  }

  function restore() {
    if (!window.Shiny || typeof Shiny.setInputValue !== 'function') return;
    var saved = null;
    try {
      var raw = window.localStorage.getItem(key);
      if (raw) saved = JSON.parse(raw);
    } catch (e) {
      report('', false, 'Browser storage unavailable: ' + String(e.message || e));
    }
    Shiny.setInputValue('pitax_header_browser_profiles_restore', {
      data: saved, nonce: Date.now()
    }, { priority: 'event' });
  }

  var saveHandlerRegistered = false;
  var restoreSent = false;

  function ensureSaveHandler() {
    if (saveHandlerRegistered) return true;
    if (!window.Shiny || typeof Shiny.addCustomMessageHandler !== 'function') return false;
    Shiny.addCustomMessageHandler('pitaxStoreHeaderProfiles', function (payload) {
      var name = payload && payload.active_template || '';
      try {
        if (!payload || payload.schema !== 1 ||
            !Array.isArray(payload.fields) || !Array.isArray(payload.templates)) {
          throw new Error('Invalid profile data received from the server.');
        }
        var raw = JSON.stringify(payload);
        window.localStorage.setItem(key, raw);
        if (window.localStorage.getItem(key) !== raw) {
          throw new Error('Browser did not retain saved profiles.');
        }
        report(name, true, '');
      } catch (e) {
        report(name, false, 'Could not save to this browser: ' + String(e.message || e));
      }
    });
    saveHandlerRegistered = true;
    return true;
  }

  function onConnected() {
    if (!ensureSaveHandler() || restoreSent) return;
    restoreSent = true;
    restore();
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
  // Custom scripts can load before Shiny finishes initializing. Registering
  // only at script-evaluation time silently loses the save message handler.
  // Wait for the Shiny API, then restore once its websocket is connected.
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
})();
