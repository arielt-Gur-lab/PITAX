      Shiny.addCustomMessageHandler('openUrl', function(message) {
        window.open(message.url, '_blank');
      });
      Shiny.addCustomMessageHandler('copyText', function(message) {
        navigator.clipboard.writeText(message.text);
      });

      (function() {
        var loaderShownAt = 0;
        var loaderFallback = null;
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
        function showStepLoader(text) {
          loaderShownAt = Date.now();
          $('#app_loading_text').text(text || 'Loading workspace...');
          $('#app_loading_overlay').addClass('visible').attr('aria-hidden', 'false');
          clearTimeout(loaderFallback);
          loaderFallback = setTimeout(hideStepLoader, 20000);
        }
        function hideStepLoader() {
          var elapsed = Date.now() - loaderShownAt;
          var wait = Math.max(0, 280 - elapsed);
          setTimeout(function() {
            $('#app_loading_overlay').removeClass('visible').attr('aria-hidden', 'true');
          }, wait);
          clearTimeout(loaderFallback);
        }
        Shiny.addCustomMessageHandler('showLoader', function(message) {
          showStepLoader(message && message.text ? message.text : 'Loading workspace...');
        });
        Shiny.addCustomMessageHandler('hideLoader', function() { hideStepLoader(); });
        $(document).on('shiny:busy', function() { $('#shiny_activity_bar').addClass('visible'); });
        $(document).on('shiny:idle', function() {
          $('#shiny_activity_bar').removeClass('visible');
          if ($('#app_loading_overlay').hasClass('visible')) hideStepLoader();
          scheduleTableAdjust(40);
        });
        $(document).on('shown.bs.tab', 'a[data-toggle="tab"]', function() { scheduleTableAdjust(60); });
        $(document).on('init.dt draw.dt', function() { scheduleTableAdjust(40); });
        $(document).on('shiny:value', function() { scheduleTableAdjust(90); });
        $(window).on('resize', function() { scheduleTableAdjust(100); });
        Shiny.addCustomMessageHandler('adjustDataTables', function() { scheduleTableAdjust(30); });

        Shiny.addCustomMessageHandler('setHelpActive', function(message) {
          var btn = $('#workflow_open_help');
          if (!btn.length) return;
          btn.toggleClass('current', !!(message && message.active));
          btn.toggleClass('available', !(message && message.active));
        });

        // Full workflow stepper chips (Bootstrap tab strip is hidden).
        $(document).on('click', '.workflow-chip[data-step]', function(evt) {
          var btn = $(this);
          var step = btn.attr('data-step');
          if (!step) return;
          if (btn.is(':disabled') || btn.hasClass('locked') || btn.attr('aria-disabled') === 'true') {
            evt.preventDefault();
            evt.stopPropagation();
            return;
          }
          showStepLoader('Loading step...');
          Shiny.setInputValue('workflow_nav_step', step, {priority: 'event'});
        });
        $(document).on('click', '#workflow_open_help', function() {
          showStepLoader('Opening Help...');
        });
      })();
