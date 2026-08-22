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
          $('#app_loading_text').text(text || 'Loading workspace\u2026');
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
          showStepLoader(message && message.text ? message.text : 'Loading workspace\u2026');
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
        $(document).on('click', '#pipeline_step > li > a', function() {
          showStepLoader('Loading step\u2026');
        });
      })();

      // The source keeps the large QC panel before Rename for maintainability,
      // while the user-facing workflow is Rename -> QC.
      function pitaxOrderWorkflowTabs() {
        var nav = $('#pipeline_step');
        var renameTab = nav.find('a[data-value="rename"]').parent();
        var qcTab = nav.find('a[data-value="qc"]').parent();
        if (renameTab.length && qcTab.length) renameTab.insertBefore(qcTab);
      }
      $(pitaxOrderWorkflowTabs);
      $(document).on('shiny:connected', pitaxOrderWorkflowTabs);
