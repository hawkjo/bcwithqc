# Config builder: notes for maintainers

`barcode_configurator.html` is the browser config builder linked from the main
README. It is one self-contained file with no build step: open it from disk or
from GitHub Pages and it works the same way, with no network access.

## Layout of the file

In order:

1. `<meta charset="utf-8">`, styles and the page skeleton.
2. `<script id="barcode_core">`: pure functions with no DOM access, exposed
   as `window.barcode_core`. All config logic lives here: reading and writing
   configs, migration of older formats, import, validation, decoder labels and
   the descriptions of bcwithqc's outputs.
3. `<script id="barcode_ui">`: state, rendering and event handlers. It calls
   the core and draws the result.
4. `<script type="application/json" id="template_data">`: the templates.
5. `<script id="vendor_sheetjs">`: SheetJS Community Edition (Apache-2.0),
   unmodified, for spreadsheet import. It must stay the last script.

The written config is exactly Python's `json.dumps(config, indent=4)` plus a
newline, in ASCII and in bcwithqc's key order.

## Running the tests

```
pip install -e ".[test]"            # bcwithqc itself, needed only for tier F
pip install --group gui-test        # pytest-playwright and spreadsheet writers (pip 25.1+)
playwright install chromium firefox webkit
pytest tests/gui                    # every engine; add --browser chromium for one
```

Without Playwright, `tests/gui` is skipped, so bcwithqc's own test run is
unaffected. `tests/gui/fixtures/SOURCES.md` records where each fixture came
from.

## Common changes

**Templates.** The templates are copies of `examples/cDNA_paired_end_config.json`,
`examples/gDNA_paired_end_config.json` and the example in
`doc/config_format.md`, without per-sequence names. After changing a source, or
the list in `update_template_data.py`, run `python gui/update_template_data.py`.
A test fails while the page's copy is out of date.

**SheetJS.** To update, replace the contents of the `vendor_sheetjs` script
with the new `xlsx.full.min.js`, unmodified, and update the size, sha256 and
version that `tests/gui/test_gui_spreadsheets.py` pins.

**Temporary warnings.** Some warnings describe what current bcwithqc can't
handle. They are listed in `TEMPORARY_WARNINGS` in the core script, each with
the bcwithqc function involved, its item in the project's future-work notes,
and the tier F test that reproduces the limitation. That test is marked
`xfail(strict=True)`, so once bcwithqc is fixed it passes unexpectedly and
fails the run. Then delete the warning, its registry entry and its test
together. Messages shown to users name neither notes items nor line numbers.
