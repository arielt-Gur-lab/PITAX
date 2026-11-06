// Real browser + real PITAX header-builder observers, not a simulated model UI.
const { chromium } = require('playwright');
const assert = require('node:assert/strict');

async function chooseShinySelect(page, id, value) {
  await page.locator('#' + id).evaluate((el, val) => {
    // Shiny selectInput defaults to the Selectize widget; dynamically added
    // choices live in selectize.options, not necessarily DOM <option> nodes.
    if (el.selectize) el.selectize.setValue(val);
    else {
      el.value = val;
      el.dispatchEvent(new Event('change', { bubbles: true }));
    }
  }, value);
}

(async () => {
  const browser = await chromium.launch({ headless: true, args: ['--no-sandbox'] });
  const page = await browser.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push(e.message));
  const active = page.locator('#e2e_active_columns');
  const registry = page.locator('#export_header_registered_fields');
  const template = page.locator('#export_header_template');
  try {
    await page.goto('http://127.0.0.1:47778', { waitUntil: 'domcontentloaded' });
    await registry.waitFor({ state: 'visible' });

    // Case 1: old Experiment removed -> no column but still registered.
    await template.fill('{Final_Name} organism={Organism}');
    await template.press('Tab');
    await page.waitForFunction(() =>
      document.querySelector('#export_header_active_template')?.textContent?.includes(
        'Active export template: >{Final_Name} organism={Organism}'));
    await page.waitForFunction(() => {
      const cols = document.querySelector('#e2e_active_columns')?.textContent || '';
      const reg = document.querySelector('#export_header_registered_fields')?.textContent || '';
      return !cols.includes('Experiment') && reg.includes('Experiment [available]');
    });
    assert.ok(!(await active.innerText()).includes('Experiment'));

    // Case 2: create new field from UI. Must IMMEDIATELY appear in
    // registry, template, active editor and bulk selector.
    await page.locator('#export_header_new_field').fill('Collection_Site');
    await page.locator('#export_header_add_field').click();
    await page.waitForFunction(() =>
      document.querySelector('#export_header_template')?.value.includes(
        'collection_site={Collection_Site}'));
    try {
      await page.waitForFunction(() => {
        const cols = document.querySelector('#e2e_active_columns')?.textContent || '';
        const reg = document.querySelector('#export_header_registered_fields')?.textContent || '';
        const select = document.getElementById('export_header_bulk_field');
        const isOption = Boolean(select?.selectize?.options?.Collection_Site) ||
          Boolean(select?.querySelector('option[value="Collection_Site"]'));
        return cols.includes('Collection_Site') &&
          reg.includes('Collection_Site [active]') && isOption;
      }, null, { timeout: 12000 });
    } catch (e) {
      console.error('FIELD DIAGNOSTIC', {
        template: await template.inputValue(),
        cols: await active.innerText(),
        registry: await registry.innerText(),
        bulkOptions: await page.locator('#export_header_bulk_field option').allTextContents(),
        existingOptions: await page.locator('#export_header_existing_field option').allTextContents(),
        persistedFields: await page.locator('#e2e_persisted_fields').innerText(),
        pageErrors: errors
      });
      throw e;
    }
    assert.ok(!(await active.innerText()).includes('Experiment'),
      'Adding a new field must not reactivate an unwanted Experiment');

    // A real metadata edit through the bulk action must update the table.
    await chooseShinySelect(page, 'export_header_bulk_field', 'Collection_Site');
    await page.locator('#export_header_bulk_value').fill('Golan');
    await page.locator('#export_header_apply_bulk').click();
    await page.waitForFunction(() =>
      document.querySelector('#export_header_values_table')?.textContent?.includes('Golan'));

    // Case 3: Save, remove from draft (not registry), reactivate, then
    // try adding with different capitalization. No duplicate field or token.
    await page.locator('#export_header_template_name').fill('Field profile');
    await page.locator('#export_header_save_template').click();
    await page.waitForFunction(() =>
      document.querySelector('#e2e_saved_templates')?.textContent?.includes('Field profile'));

    await template.fill('{Final_Name} organism={Organism}');
    await template.press('Tab');
    await page.waitForFunction(() => {
      const cols = document.querySelector('#e2e_active_columns')?.textContent || '';
      const reg = document.querySelector('#export_header_registered_fields')?.textContent || '';
      return !cols.includes('Collection_Site') && reg.includes('Collection_Site [available]');
    });
    await chooseShinySelect(page, 'export_header_existing_field', 'Collection_Site');
    await page.locator('#export_header_activate_field').click();
    await page.waitForFunction(() => {
      const reg = document.querySelector('#export_header_registered_fields')?.textContent || '';
      return reg.includes('Collection_Site [active]');
    });
    await page.locator('#export_header_new_field').fill('collection_site');
    await page.locator('#export_header_add_field').click();
    await page.waitForFunction(() =>
      !document.querySelector('#export_header_new_field')?.value);
    const fields = (await page.locator('#e2e_persisted_fields').innerText()).split(',');
    assert.equal(fields.filter(f => f.toLowerCase() === 'collection_site').length, 1);
    const header = await template.inputValue();
    assert.equal((header.match(/\{Collection_Site\}/g) || []).length, 1);

    // Reload the saved profile; both the field and template must survive.
    await chooseShinySelect(page, 'export_header_saved_template', 'Field profile');
    await page.locator('#export_header_load_template').click();
    await page.waitForFunction(() =>
      document.querySelector('#export_header_template')?.value.includes(
        'collection_site={Collection_Site}'));
    assert.ok((await registry.innerText()).includes('Collection_Site [active]'));
    assert.ok((await active.innerText()).includes('Collection_Site'));

    // Explicitly remove Experiment from the active header, preserving its
    // registered definition and avoiding accidental deletion of older templates.
    await chooseShinySelect(page, 'export_header_existing_field', 'Experiment');
    await page.locator('#export_header_activate_field').click();
    await page.waitForFunction(() =>
      document.querySelector('#export_header_template')?.value.includes('{Experiment}'));
    await page.locator('#export_header_remove_from_header').click();
    await page.waitForFunction(() =>
      !document.querySelector('#export_header_template')?.value.includes('{Experiment}'));
    assert.ok((await registry.innerText()).includes('Experiment [available]'));

    // Delete a new unsaved field: draft token, registry entry and values
    // must be removed, without touching the other registered fields.
    await page.locator('#export_header_new_field').fill('Temporary_Field');
    await page.locator('#export_header_add_field').click();
    await page.waitForFunction(() =>
      document.querySelector('#export_header_template')?.value.includes('{Temporary_Field}'));
    await chooseShinySelect(page, 'export_header_existing_field', 'Temporary_Field');
    await page.locator('#export_header_delete_field').click();
    await page.waitForFunction(() => {
      const reg = document.querySelector('#export_header_registered_fields')?.textContent || '';
      const tmpl = document.querySelector('#export_header_template')?.value || '';
      return !reg.includes('Temporary_Field') && !tmpl.includes('{Temporary_Field}');
    });
    assert.ok((await registry.innerText()).includes('Experiment [available]'));
    assert.deepEqual(errors, [], 'Uncaught Shiny browser errors');
    console.log('PASS: create/activate/edit/reuse/save/load metadata field through real Shiny UI');
  } finally {
    await browser.close();
  }
})().catch(e => { console.error('FAIL: metadata field UI integration', e); process.exit(1); });
