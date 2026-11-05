// Real browser integration: click visible Shiny download link after preparing.
const browsers = require('playwright');
const engine = process.env.PITAX_BROWSER || 'chromium';
const assert = require('node:assert/strict');
const fs = require('node:fs');

(async () => {
  assert.ok(['chromium', 'firefox'].includes(engine), 'Unsupported test browser');
  const browser = await browsers[engine].launch({ headless: true,
    args: engine === 'chromium' ? ['--no-sandbox'] : [] });
  const context = await browser.newContext({ acceptDownloads: true });
  const page = await context.newPage();
  try {
    await page.goto('http://127.0.0.1:47777', { waitUntil: 'domcontentloaded' });
    await page.getByRole('button', { name: '1. Prepare ZIP' }).waitFor();
    assert.equal(await page.locator('#download_export_package').count(), 0,
      'Download link should not appear before an archive exists');
    await page.getByRole('button', { name: '1. Prepare ZIP' }).click();
    const link = page.locator('#download_export_package');
    await link.waitFor({ state: 'visible', timeout: 30000 });
    assert.match(await page.locator('#status').innerText(), /ZIP ready/);

    const waitForDownload = page.waitForEvent('download', { timeout: 30000 });
    await link.click(); // native, real user click; no synthetic Shiny JS message
    const download = await waitForDownload;
    assert.equal(download.suggestedFilename(), 'PITAX_browser_test.zip');
    const downloadedPath = await download.path();
    const content = fs.readFileSync(downloadedPath);
    assert.ok(content.length > 100, 'ZIP should contain real bytes');
    assert.equal(content.toString('ascii', 0, 2), 'PK');
    console.log('PASS:', engine, 'visible download button triggered a real browser ZIP download', content.length, 'bytes');
  } finally {
    await browser.close();
  }
})().catch(error => {
  console.error('FAIL: real browser ZIP download:', error);
  process.exit(1);
});
