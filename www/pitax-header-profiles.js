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

  if (window.Shiny && typeof Shiny.addCustomMessageHandler === 'function') {
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
  }

  // Shiny fires this jQuery event when each new browser session connects.
  // It runs again on page reload; no server-side persistence is assumed.
  if (window.jQuery) {
    window.jQuery(document).on('shiny:connected', restore);
  } else {
    document.addEventListener('shiny:connected', restore);
  }
})();
