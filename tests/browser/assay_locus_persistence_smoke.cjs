// With PITAX running locally and Playwright available:
// node tests/browser/assay_locus_persistence_smoke.cjs [http://127.0.0.1:3838]
// Optional: PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH for an existing Chromium binary.
const { chromium } = require('playwright');
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');

(async () => {
  const browser = await chromium.launch({
    headless: true,
    ...(process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH
      ? { executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH } : {})
  });
  try {
    const url = process.argv[2] || 'http://127.0.0.1:3838';
    const page = await browser.newPage({ acceptDownloads: true });
    const errors = [];
    const trackErrors = p => p.on('pageerror', error => errors.push(error.message));
    trackErrors(page);
    const choose = (p, id, value) => p.evaluate(({ id, value }) => {
      document.getElementById(id).selectize.setValue(value);
    }, { id, value });
    const expectProfile = (p, name, locus) => p.waitForFunction(({ name, locus }) =>
      document.getElementById('assay_name').value === name &&
      document.getElementById('target').value === locus,
    { name, locus }, { timeout: 10000 }).catch(async error => {
      console.error('Unexpected assay state:', await p.evaluate(() => ({
        selected: document.getElementById('assay_editor_select').value,
        name: document.getElementById('assay_name').value,
        locus: document.getElementById('target').value,
        snapshot: Shiny.shinyapp.$inputValues.assay_editor_snapshot,
        choices: document.getElementById('assay_editor_select').selectize.options
      })));
      throw error;
    });
    const addLocus = async (id, label) => {
      await page.click('#locus_add');
      await page.fill('#locus_new_id', id);
      await page.fill('#locus_display_name', label);
      await page.click('#locus_save');
      await page.waitForFunction(id => document.getElementById('target').value === id, id);
    };
    await page.goto(url);
    await page.waitForSelector('#workflow_goto_settings');
    const raw = fs.readFileSync(path.join(__dirname, '../fixtures/Stage3_synthetic_pair/TESTPAIR001_forward.ab1'));
    await page.locator('#ab1_files').setInputFiles([
      { name: 'BAR.1.ab1', mimeType: 'application/octet-stream', buffer: raw }
    ]);
    await page.click('#workflow_goto_settings');
    const first = await page.locator('#assay_editor_select').inputValue();
    await addLocus('VTC1', 'Custom locus one');
    await page.fill('#assay_name', 'First custom assay');
    await page.click('#assay_add_profile');
    await page.waitForFunction(first => document.getElementById('assay_editor_select').value !== first, first);
    await expectProfile(page, 'New assay', 'ITS');
    const second = await page.locator('#assay_editor_select').inputValue();
    await addLocus('VTC2', 'Custom locus two');
    await page.fill('#assay_name', 'Second custom assay');
    await choose(page, 'assay_editor_select', first);
    await expectProfile(page, 'First custom assay', 'VTC1');
    await choose(page, 'assay_editor_select', second);
    await expectProfile(page, 'Second custom assay', 'VTC2');

    const downloadPromise = page.waitForEvent('download');
    await page.click('#save_project');
    const project = await (await downloadPromise).path();
    const loaded = await browser.newPage();
    trackErrors(loaded);
    await loaded.goto(url);
    await loaded.waitForSelector('#workflow_goto_settings');
    await loaded.locator('#load_project').setInputFiles(project);
    await loaded.click('#workflow_goto_settings');
    await expectProfile(loaded, 'Second custom assay', 'VTC2');
    await choose(loaded, 'assay_editor_select', first);
    await expectProfile(loaded, 'First custom assay', 'VTC1');
    const catalog = await loaded.locator('#target').evaluate(el => el.selectize.options);
    assert.equal(catalog.VTC1.label || catalog.VTC1.text, 'Custom locus one');
    assert.equal(catalog.VTC2.label || catalog.VTC2.text, 'Custom locus two');
    assert.deepEqual(errors, []);
    console.log('PASS custom loci and profiles survive edits, switching, save and fresh-session load');
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
