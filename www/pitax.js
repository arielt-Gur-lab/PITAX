      Shiny.addCustomMessageHandler('openUrl', function(message) {
        window.open(message.url, '_blank');
      });
      Shiny.addCustomMessageHandler('copyText', function(message) {
        navigator.clipboard.writeText(message.text);
      });
      (function() {
        var owner = null;
        var pending = null;
        var fields = ['assay_name', 'assay_name_source', 'target', 'forward_primer',
          'reverse_primer', 'forward_primer_seq', 'reverse_primer_seq',
          'expected_amplicon_len', 'absolute_max_base_index', 'sequencing_primer',
          'window', 'min_peak_ratio', 'min_relative_signal', 'min_len_before_collapse',
          'bad_run_windows', 'min_usable_len', 'enable_primer_mapping'];
        function publishEditor() {
          if (!owner || document.getElementById('assay_editor_select').value !== owner.assay_id) return;
          var values = {};
          fields.forEach(function(id) {
            var radio = document.querySelector('input[name="' + id + '"]:checked');
            var el = document.getElementById(id);
            if (radio) values[id] = radio.value;
            else if (el && el.type === 'checkbox') values[id] = el.checked;
            else if (el && el.type === 'number') values[id] = el.value === '' || !isFinite(Number(el.value)) ? null : Number(el.value);
            else if (el) values[id] = el.value;
          });
          Shiny.setInputValue('assay_editor_snapshot', {
            assay_id: owner.assay_id, revision: owner.revision, values: values
          }, {priority: 'event'});
        }
        function applyPendingEditor() {
          var message = pending;
          if (!message) return;
          // Apply all profile fields and their owner together; no input from
          // the previous profile can be published under the new profile ID.
          owner = null;
          setTimeout(function() {
            if (pending !== message) return;
            var active = document.getElementById('assay_editor_select');
            if (!active || active.value !== message.assay_id) return;
            var locus = document.getElementById('target');
            if (locus && locus.selectize && message.locus_choices) {
              locus.selectize.clear(true);
              locus.selectize.clearOptions();
              message.locus_choices.forEach(function(choice) {
                locus.selectize.addOption(choice);
              });
            }
            Object.keys(message.values || {}).forEach(function(id) {
              var value = message.values[id];
              var radios = document.querySelectorAll('input[name="' + id + '"][type="radio"]');
              if (radios.length) {
                Array.prototype.forEach.call(radios, function(el) { el.checked = el.value === value; });
                $(radios).filter(':checked').trigger('change');
                return;
              }
              var el = document.getElementById(id);
              if (!el) return;
              if (el.selectize) el.selectize.setValue(value, true);
              else el.value = value == null ? '' : value;
              $(el).trigger('change');
            });
            owner = message;
            pending = null;
            publishEditor();
          }, 0);
        }
        Shiny.addCustomMessageHandler('assayEditorLoaded', function(message) {
          // Shiny can deliver this before its selector update; retain it
          // until that selector has reached the same profile.
          pending = message;
          applyPendingEditor();
        });
        $(document).on('change.pitaxAssaySelection', '#assay_editor_select', applyPendingEditor);
        Shiny.addCustomMessageHandler('assayEditorChoices', function(message) {
          var el = document.getElementById('assay_editor_select');
          if (!el || !el.selectize) return;
          (message.choices || []).forEach(function(choice) {
            if (el.selectize.options[choice.value]) el.selectize.updateOption(choice.value, choice);
          });
        });
        $(document).on('input.pitaxAssay change.pitaxAssay', function(event) {
          var el = event.target;
          if (fields.indexOf(el.id) >= 0 || fields.indexOf(el.name) >= 0) publishEditor();
        });
      })();

      // A Shiny modal with easyClose=true can disappear via ESC, backdrop
      // or its X control without firing the server's Close actionButton.
      // Capture the currently open guided-review shell on hide and notify
      // the server once the Bootstrap transition finishes. The generation
      // token rejects late dismissals when a user quickly reopens a row.
      (function() {
        $(document).on('hide.bs.modal.pitaxPeakReview', '#shiny-modal', function() {
          var marker = $(this).find('.pitax-peak-review-shell-marker').first();
          var generation = marker.length ? Number(marker.attr('data-generation')) : NaN;
          $(this).data('pitaxPeakReviewHideGeneration', Number.isInteger(generation) ? generation : null);
        });
        $(document).on('hidden.bs.modal.pitaxPeakReview', '#shiny-modal', function() {
          var generation = $(this).data('pitaxPeakReviewHideGeneration');
          $(this).removeData('pitaxPeakReviewHideGeneration');
          if (!Number.isInteger(generation) || !window.Shiny || !Shiny.setInputValue) return;
          Shiny.setInputValue('peak_review_modal_dismissed',
            {generation: generation, nonce: Date.now() + Math.random()},
            {priority: 'event'});
        });
      })();

      // Developer-only browser latency probe. Dormant until Enable diagnostics
      // is checked, or until window.__PITAX_PERF_CONSOLE__ / localStorage
      // pitaxPerfConsole=1 enables a console-only stopwatch (for OFF vs ON
      // comparisons without server table updates). Measures click -> first
      // server message -> tracked output values -> shiny:idle. Never transmits
      // DNA, filenames, samples or URLs.
      (function() {
        var pending = null;
        var TRACKED = {
          'peak_review_local_plot': true,
          'peak_review_side_ui': true,
          'peak_review_modal_title': true,
          'peak_review_advanced_ui': true,
          'ambiguous_peak_table': true,
          'ambiguous_peak_summary': true,
          'pitax_perf_summary': true,
          'pitax_perf_recent': true
        };
        function recordingEnabled() {
          var checkbox = document.getElementById('pitax_perf_enabled');
          return !!(checkbox && checkbox.checked);
        }
        function consoleEnabled() {
          if (window.__PITAX_PERF_CONSOLE__) return true;
          try { return localStorage.getItem('pitaxPerfConsole') === '1'; }
          catch (e) { return false; }
        }
        function active() { return recordingEnabled() || consoleEnabled(); }
        function markAction(action) {
          pending = {
            action: action,
            started: performance.now(),
            firstMessageAt: null,
            lastValueAt: null,
            values: {},
            finalizing: false,
            idleAt: null
          };
        }

        document.addEventListener('click', function(event) {
          if (!active()) return;
          var target = event.target;
          if (!target || !target.closest) return;
          var action = null;
          if (target.closest('#ambiguous_peak_table tbody tr')) action = 'guided.row';
          else if (target.closest('#peak_review_start')) action = 'guided.start';
          else if (target.closest('#peak_review_continue')) action = 'guided.continue';
          else if (target.closest('#peak_review_flag')) action = 'guided.flag';
          else if (target.closest('#confirm_curation_action')) action = 'guided.confirm';
          else if (target.closest('#cancel_curation_action')) action = 'guided.cancel';
          if (action) markAction(action);
        }, true);

        $(document).on('shiny:message.pitaxDevPerf', function() {
          if (!pending || pending.firstMessageAt != null) return;
          pending.firstMessageAt = performance.now();
        });
        $(document).on('shiny:value.pitaxDevPerf', function(event) {
          if (!pending) return;
          var name = event && event.name;
          if (!name || !TRACKED[name]) return;
          var now = performance.now();
          pending.values[name] = Math.round(now - pending.started);
          pending.lastValueAt = now;
        });
        $(document).on('shiny:idle.pitaxDevPerf', function() {
          if (!pending || !active() || pending.finalizing) return;
          // shiny:idle can fire in the same turn as shiny:value; defer so
          // value handlers run, then briefly wait for late plot/DT bindings.
          var current = pending;
          current.finalizing = true;
          current.idleAt = performance.now();
          setTimeout(function() {
            if (pending !== current) return;
            setTimeout(function() {
              if (pending !== current) return;
              pending = null;
              var idleElapsed = current.idleAt - current.started;
              var settledElapsed = performance.now() - current.started;
              if (!isFinite(idleElapsed) || idleElapsed < 0) return;
              var toMsg = current.firstMessageAt == null ? null :
                Math.round(current.firstMessageAt - current.started);
              var toVal = current.lastValueAt == null ? null :
                Math.round(current.lastValueAt - current.started);
              var valueKeys = Object.keys(current.values);
              var details = 'msg=' + (toMsg == null ? 'na' : toMsg) +
                ';ui=' + (toVal == null ? 'na' : toVal) +
                ';outs=' + valueKeys.length +
                ';settled=' + Math.round(settledElapsed);
              var record = {
                action: current.action,
                duration_ms: Math.round(idleElapsed * 100) / 100,
                settled_ms: Math.round(settledElapsed * 100) / 100,
                to_first_message_ms: toMsg,
                to_last_value_ms: toVal,
                values: current.values,
                details: details,
                recording: recordingEnabled(),
                t: new Date().toISOString()
              };
              if (!window.__pitaxPerfLog) window.__pitaxPerfLog = [];
              window.__pitaxPerfLog.push(record);
              if (window.__pitaxPerfLog.length > 100) window.__pitaxPerfLog.shift();
              if (consoleEnabled() && window.console && console.info) {
                console.info('[PITAX perf]', current.action, record);
              }
              if (!recordingEnabled() || !window.Shiny || !Shiny.setInputValue) return;
              Shiny.setInputValue('pitax_perf_browser_event', {
                action: current.action,
                duration_ms: idleElapsed,
                to_first_message_ms: toMsg,
                to_last_value_ms: toVal,
                details: details,
                nonce: Date.now() + Math.random()
              }, {priority: 'event'});
            }, 50);
          }, 0);
        });
        $(document).on('change.pitaxDevPerf', '#pitax_perf_enabled', function() {
          if (!this.checked && !consoleEnabled()) pending = null;
        });
      })();

      (function() {
        var loaderShownAt = 0;
        var loaderFallback = null;
        var hideTimer = null;
        var tableAdjustTimer = null;

        function adjustVisibleDataTables() {
          if (!$.fn.dataTable) return;
          var api = $.fn.dataTable.tables({ visible: true, api: true });
          if (api && api.columns) api.columns.adjust();
          $('.dataTables_scrollBody').each(function() {
            var body = this;
            var head = $(body).closest('.dataTables_scroll').find('.dataTables_scrollHead').get(0);
            if (head) head.scrollLeft = body.scrollLeft;
            $(body).off('scroll.pitaxAlignment').on('scroll.pitaxAlignment', function() {
              if (head) head.scrollLeft = body.scrollLeft;
            });
          });
        }
        function scheduleTableAdjust(delay) {
          clearTimeout(tableAdjustTimer);
          tableAdjustTimer = setTimeout(adjustVisibleDataTables, delay || 80);
        }
        function clearLoaderTimers() {
          clearTimeout(loaderFallback);
          clearTimeout(hideTimer);
          loaderFallback = null;
          hideTimer = null;
        }
        function showStepLoader(text) {
          clearLoaderTimers();
          loaderShownAt = Date.now();
          $('#app_loading_text').text(text || 'Loading workspace...');
          $('#app_loading_overlay').addClass('visible').attr('aria-hidden', 'false');
          loaderFallback = setTimeout(function() { hideStepLoader(true); }, 12000);
        }
        function hideStepLoader(force) {
          clearTimeout(loaderFallback);
          loaderFallback = null;
          var doHide = function() {
            $('#app_loading_overlay').removeClass('visible').attr('aria-hidden', 'true');
            hideTimer = null;
          };
          clearTimeout(hideTimer);
          if (force) {
            doHide();
            return;
          }
          var elapsed = Date.now() - loaderShownAt;
          var wait = Math.max(0, 280 - elapsed);
          hideTimer = setTimeout(doHide, wait);
        }
        Shiny.addCustomMessageHandler('showLoader', function(message) {
          showStepLoader(message && message.text ? message.text : 'Loading workspace...');
        });
        Shiny.addCustomMessageHandler('hideLoader', function(message) { hideStepLoader(true); });
        $(document).on('shiny:busy', function() { $('#shiny_activity_bar').addClass('visible'); });
        $(document).on('shiny:idle', function() {
          $('#shiny_activity_bar').removeClass('visible');
          if ($('#app_loading_overlay').hasClass('visible')) hideStepLoader(true);
          scheduleTableAdjust(40);
        });
        $(document).on('shown.bs.tab', 'a[data-toggle="tab"]', function() { scheduleTableAdjust(60); });
        $(document).on('init.dt draw.dt', function() { scheduleTableAdjust(40); });
        $(document).on('shiny:value', function() { scheduleTableAdjust(90); });
        $(window).on('resize', function() { scheduleTableAdjust(100); });
        Shiny.addCustomMessageHandler('adjustDataTables', function(message) { scheduleTableAdjust(30); });

        Shiny.addCustomMessageHandler('setHelpActive', function(message) {
          var btn = $('#workflow_open_help');
          if (!btn.length) return;
          btn.toggleClass('current', !!(message && message.active));
          btn.toggleClass('available', !(message && message.active));
        });

        // Step chips are Shiny actionButtons (workflow_goto_*). Do not put a
        // blocking overlay on chip clicks - that freezes the entire UI.
        $(document).on('click', '#workflow_open_help', function() {
          hideStepLoader(true);
        });
      })();
