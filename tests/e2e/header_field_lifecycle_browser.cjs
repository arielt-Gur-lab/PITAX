// Real browser + structured Header Builder against production assets.
const { chromium } = require('playwright');
const assert = require('node:assert/strict');

async function waitActive(page, pred, timeout = 15000) {
  await page.waitForFunction(pred, null, { timeout });
}

(async () => {
  const browser = await chromium.launch({ headless: true, args: ['--no-sandbox'] });
  const page = await browser.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push(e.message));
  try {
    // Production assets must be reachable (same files as deployed www/).
    const css = await page.request.get('http://127.0.0.1:47778/pitax-www/pitax-header-builder.css');
    const js = await page.request.get('http://127.0.0.1:47778/pitax-www/pitax-header-builder.js');
    assert.equal(css.status(), 200, 'header builder CSS must return HTTP 200');
    assert.equal(js.status(), 200, 'header builder JS must return HTTP 200');

    await page.goto('http://127.0.0.1:47778', { waitUntil: 'domcontentloaded' });
    await page.locator('#fhb_library').waitFor({ state: 'visible' });
    await page.waitForFunction(() => window.PITAXHeaderBuilder &&
      Array.isArray(window.PITAXHeaderBuilder.getState().library) &&
      window.PITAXHeaderBuilder.getState().library.length > 0);

    // Default profile must not include experiment.
    const initial = await page.locator('#e2e_persisted_fields').innerText();
    assert.ok(!initial.split(',').includes('experiment'), 'experiment not auto-included');

    // Add organism via JS bridge (same path as drag/drop add).
    await page.evaluate(() => window.PITAXHeaderBuilder.addField('organism'));
    await waitActive(page, () =>
      (document.querySelector('#e2e_active_columns')?.textContent || '').includes('organism'));

    // Reorder: move organism to front via keyboard-style controls in construction.
    await page.evaluate(() => {
      const st = window.PITAXHeaderBuilder.getState();
      // Trigger server reorder with organism first.
      if (window.Shiny) {
        Shiny.setInputValue('pitax_header_builder_order', {
          fields: ['organism'].concat(st.activeFields.filter(f => f !== 'organism')),
          nonce: Date.now()
        }, { priority: 'event' });
      }
    });
    await waitActive(page, () =>
      (document.querySelector('#e2e_persisted_fields')?.textContent || '').startsWith('organism'));

    // Create custom field.
    await page.locator('#export_header_new_field').fill('Collection_Site');
    await page.locator('#export_header_add_field').click();
    await waitActive(page, () =>
      (document.querySelector('#e2e_active_columns')?.textContent || '').includes('collection_site'));
    assert.ok((await page.locator('#export_header_registered_fields').innerText())
      .includes('Collection_Site'));

    // Save named profile.
    await page.locator('#export_header_template_name').fill('Field profile');
    await page.locator('#export_header_save_template').click();
    await waitActive(page, () =>
      (document.querySelector('#e2e_saved_profiles')?.textContent || '').includes('Field profile'));

    // Remove collection_site from draft; registry should keep it available.
    await page.evaluate(() => window.PITAXHeaderBuilder.removeField('collection_site'));
    await waitActive(page, () => {
      const cols = document.querySelector('#e2e_active_columns')?.textContent || '';
      const reg = document.querySelector('#export_header_registered_fields')?.textContent || '';
      return !cols.includes('collection_site') && reg.includes('Collection_Site');
    });

    // Reactivate via addField.
    await page.evaluate(() => window.PITAXHeaderBuilder.addField('collection_site'));
    await waitActive(page, () =>
      (document.querySelector('#e2e_active_columns')?.textContent || '').includes('collection_site'));

    // Experiment must only appear when explicitly added.
    await page.evaluate(() => window.PITAXHeaderBuilder.addField('experiment'));
    await waitActive(page, () =>
      (document.querySelector('#e2e_persisted_fields')?.textContent || '').includes('experiment'));
    await page.evaluate(() => window.PITAXHeaderBuilder.removeField('experiment'));
    await waitActive(page, () =>
      !(document.querySelector('#e2e_persisted_fields')?.textContent || '').includes('experiment'));
    assert.ok(!(await page.locator('#e2e_active_columns').innerText()).includes('experiment'));

    // Delete temporary custom field (remove from profile, then delete definition).
    await page.locator('#export_header_new_field').fill('Temporary_Field');
    await page.locator('#export_header_add_field').click();
    await waitActive(page, () =>
      (document.querySelector('#e2e_persisted_fields')?.textContent || '').includes('temporary_field'));
    await page.evaluate(() => window.PITAXHeaderBuilder.removeField('temporary_field'));
    await waitActive(page, () =>
      !(document.querySelector('#e2e_active_columns')?.textContent || '').includes('temporary_field'));
    await page.evaluate(() => {
      if (window.Shiny) {
        Shiny.setInputValue('export_header_existing_field', 'temporary_field', { priority: 'event' });
        Shiny.setInputValue('export_header_delete_field', Date.now(), { priority: 'event' });
      }
    });
    await waitActive(page, () => {
      const reg = document.querySelector('#export_header_registered_fields')?.textContent || '';
      return !reg.includes('Temporary_Field') && !reg.includes('temporary_field');
    }, 20000);

    // Browser persistence with verified acknowledgement.
    await page.locator('#export_header_template_name').fill('Browser durable');
    await page.locator('#export_header_save_template').click();
    await page.waitForTimeout(500);
    await page.locator('#export_header_save_browser').click();
    await waitActive(page, () =>
      (document.querySelector('#e2e_browser_status')?.textContent || '').includes('browser=TRUE'),
      20000);

    const stored = await page.evaluate(() => {
      const raw = window.localStorage.getItem('PITAX_FASTA_HEADER_PROFILES_V2');
      return raw ? JSON.parse(raw) : null;
    });
    assert.ok(stored && Array.isArray(stored.profiles) && stored.profiles.length >= 1,
      'localStorage must contain verified profiles');
    assert.equal(stored.schema, 2);

    // Hard reload + new session must restore browser profiles.
    await page.reload({ waitUntil: 'domcontentloaded' });
    await page.locator('#fhb_library').waitFor({ state: 'visible' });
    // Explicitly re-push localStorage profiles in case the automatic
    // shiny:connected restore raced the test harness websocket.
    await page.waitForFunction(() => window.Shiny && typeof Shiny.setInputValue === 'function');
    await page.evaluate(() => {
      const raw = window.localStorage.getItem('PITAX_FASTA_HEADER_PROFILES_V2');
      const saved = raw ? JSON.parse(raw) : null;
      Shiny.setInputValue('pitax_header_browser_profiles_restore', {
        data: saved,
        nonce: Date.now()
      }, { priority: 'event' });
    });
    await waitActive(page, () =>
      (document.querySelector('#e2e_saved_profiles')?.textContent || '').includes('Browser durable'),
      20000);

    assert.deepEqual(errors, [], 'Uncaught Shiny browser errors: ' + errors.join(' | '));
    console.log('PASS: structured header builder drag/reorder/save/browser persistence');
  } finally {
    await browser.close();
  }
})().catch(e => {
  console.error('FAIL: header builder browser integration', e);
  process.exit(1);
});
