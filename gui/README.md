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
newline, in ASCII. Its first key, `generated_with`, records the name and
version of the tool the page writes for, as `pyproject.toml` gives them;
bcwithqc's own keys follow in bcwithqc's order.

## Wording

Everything the page shows calls bcwithqc's blocks *regions*. The config keeps
bcwithqc's keys (`blocks`, `blockname`, `keepblock`), and the code may say
block. The keep/drop switch says what changes in the reads bcwithqc writes,
"In output reads: Keep | Leave out", and doesn't name a next step: the output
reads may be aligned, pseudoaligned or counted. A region that is left out is
still written as its tag (CB, UR, …). The Outputs panel lists the output reads
under "Output reads".
The switch for what follows the regions always says "Rest of Read N", also
while a read has no regions. A read without regions is always kept whole,
so its "Leave out" is disabled.
Text doesn't assume the data exists yet: users may be planning an experiment.
Steps not done yet (no read layout, no regions, an empty list) are *to do*
items. They block the download like errors, but are not shown as mistakes.

In a question, Return presses the filled button. When the focus starts on a
button, that button is the filled one, and its focus ring always shows. Where
a question can replace or delete work, the filled button is Cancel. In a
question that starts on a choice (single- or paired-end, which column), the
focus is on the choice, and Return presses the filled button once a choice is
made.

## Browsers

The page needs Chrome or Edge 86+, Firefox 113+ or Safari 16.4+ (gzip
streams, `<dialog>` and `replaceChildren`); `missing_features()` names what an
older browser lacks in a banner. Colours mixed with `color-mix()` have plain
fallbacks in an `@supports not` block for Chrome and Edge before 111.

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
the list in `update_page.py`, run `python gui/update_page.py`. A test fails
while the page's copy is out of date.

**Version bumps.** The `generated_with` stamp comes from `GENERATED_WITH` in
the core script, which `python gui/update_page.py` writes from
`pyproject.toml`: the name from `[project] name`, the version from
`[project] version` or the file `[tool.hatch.version]` points to. After a
bump, run the script and commit the page with the bump; a test fails, in CI
too, until you do. Bump the version whenever the config format changes, so
each config says which format it follows.

**Renaming the project.** The rename itself is a find-and-replace across the
repo. For the config builder:

1. Change `[project] name` in `pyproject.toml` (and `[tool.hatch.version]`
   `path`, if the package folder moves), then run `python gui/update_page.py`.
   New configs then carry the new name. Configs written before the rename
   keep the old one, and the page reports it when it opens them.
2. Replace the old name in the messages of the core and UI scripts, and in
   this file.
3. Tier F runs the package by name (`python -m bcwithqc`, imports from
   `bcwithqc.*`): update `tests/gui/test_gui_pipeline_contract.py`.
4. GitHub redirects a renamed repository, but not its GitHub Pages site: the
   page moves to `https://<owner>.github.io/<new name>/gui/barcode_configurator.html`
   and the old address stops working. Update the link in the main README and
   `PAGES_ROOT` in `tests/gui/test_gui_docs.py`. A custom domain for Pages,
   set up before the rename, keeps the address stable.

**SheetJS.** To update, replace the contents of the `vendor_sheetjs` script
with the new `xlsx.full.min.js`, unmodified, and update the size, sha256 and
version that `tests/gui/test_gui_spreadsheets.py` pins.

**Temporary warnings.** Some warnings describe what current bcwithqc can't
handle, and some notes in the Outputs panel say where current bcwithqc writes
something else. They are listed in `TEMPORARY_WARNINGS` and
`TEMPORARY_OUTPUT_NOTES` in the core script, each with the bcwithqc function
involved, its item in the project's future-work notes, and the tier F test that
reproduces the limitation. That test is marked
`xfail(strict=True)`, so once bcwithqc is fixed it passes unexpectedly and
fails the run. Then delete the warning, its registry entry and its test
together. Messages shown to users name neither notes items nor line numbers.
