"""Helpers shared by the config builder's browser tests.

Expected values are computed here in Python, independently of the page's
JavaScript: a config's canonical text is Python's own json.dumps output, and
the name and version it records are read from pyproject.toml.
"""

import copy
import gzip
import json
import re
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = REPO_ROOT / "examples"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
GUI_PATH = REPO_ROOT / "gui" / "barcode_configurator.html"
DESIGN_DOC_PATH = REPO_ROOT / "doc" / "config_format.md"
READ_KEYS = ("barcode_struct_r1", "barcode_struct_r2")


def example_config_names() -> list[str]:
    """Return the file names of bcwithqc's example configs, sorted."""
    return sorted(path.name for path in EXAMPLES_DIR.glob("*.json"))


def read_bytes_text(path: Path) -> str:
    """Return a config file's text, decompressing .gz files."""
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode("utf-8")


def read_config(path: Path) -> dict:
    """Return a config file parsed as JSON, decompressing .gz files."""
    return json.loads(read_bytes_text(path))


def generated_with() -> dict:
    """Return the name and version pyproject.toml gives, which the page must stamp.

    The version is [project] version or, when it is dynamic, the __version__
    in the file [tool.hatch.version] points to. Read when called, not at
    import, so a problem here fails the tests that need it instead of
    stopping pytest's collection.
    """
    pyproject = tomllib.loads(
        (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )
    project = pyproject["project"]
    version = project.get("version")
    if version is None:
        version_file = REPO_ROOT / pyproject["tool"]["hatch"]["version"]["path"]
        match = re.search(
            r"^__version__\s*=\s*['\"]([^'\"]+)['\"]",
            version_file.read_text(encoding="utf-8"),
            flags=re.MULTILINE,
        )
        assert match is not None, f"{version_file} sets no __version__"
        version = match.group(1)
    return {"name": project["name"], "version": version}


def canonical_text(config: dict) -> str:
    """Return the exact text the GUI must write for a config.

    That is json.dumps(indent=4) plus a newline, with "generated_with" as the
    first key, set to pyproject.toml's name and version whatever the config
    says.
    """
    written = {"generated_with": generated_with()}
    written.update(
        (key, value) for key, value in config.items() if key != "generated_with"
    )
    return json.dumps(written, indent=4) + "\n"


def without_names(config: dict) -> dict:
    """Return a copy of a config with every block's per-sequence "name" list removed."""
    result = copy.deepcopy(config)
    for key in READ_KEYS:
        for block in result.get(key, {}).get("blocks", []):
            block.pop("name", None)
    return result


def with_sorted_lists(config: dict) -> dict:
    """Return a copy of a config with each barcode list sorted A→Z; constant regions unchanged."""
    result = copy.deepcopy(config)
    for key in READ_KEYS:
        for block in result.get(key, {}).get("blocks", []):
            if block["blocktype"] == "barcodeList":
                block["sequence"] = sorted(block["sequence"])
    return result


def design_doc_example() -> dict:
    """Return the example config under "Example Config.json" in doc/config_format.md."""
    text = DESIGN_DOC_PATH.read_text(encoding="utf-8")
    match = re.search(
        r"^## Example Config\.json\s*```json\n(.*?)\n```",
        text,
        flags=re.DOTALL | re.MULTILINE,
    )
    assert match is not None, "doc/config_format.md has no example config"
    return json.loads(match.group(1))


OLD_GUI_BARCODES = ["CCTTAACAT", "AGCGTAGAA", "TATAGGTGT", "GGCGACTCA"]
OLD_GUI_SAMPLE_BARCODES = ["TCGCCTTA", "CTAGTACG", "TTCTGCCT", "GCTCAGGA"]


def old_gui_config() -> dict:
    """Return a config laid out exactly as barcodinggui 2b473d1 saves one.

    barcodinggui writes maxerrors and blockfunction, then the sequence table
    (a "name" list, left out for barcode lists whose names are all empty, then
    "sequence") or the length, then blocktype, blockname and, for linked
    blocks, linkgroup. It writes no keepblock and no single_end_reads, and
    names new blocks "sequence 0", "sequence 1", ...
    """
    return {
        "unknown_read_orientation": True,
        "barcode_struct_r1": {
            "keep_nonbarcode": False,
            "blocks": [
                {
                    "maxerrors": 1,
                    "blockfunction": "barcode",
                    "sequence": OLD_GUI_BARCODES,
                    "blocktype": "barcodeList",
                    "blockname": "sequence 0",
                    "linkgroup": 0,
                },
                {
                    "maxerrors": 1,
                    "blockfunction": "barcode",
                    "name": ["linker"],
                    "sequence": ["GGGTCAGTACGTACGAGTCCC"],
                    "blocktype": "constantRegion",
                    "blockname": "sequence 1",
                },
                {
                    "maxerrors": 1,
                    "blockfunction": "discard",
                    "name": ["s1", "s2", "s3", "s4"],
                    "sequence": OLD_GUI_SAMPLE_BARCODES,
                    "blocktype": "barcodeList",
                    "blockname": "sequence 2",
                    "linkgroup": 0,
                },
                {
                    "maxerrors": 0,
                    "blockfunction": "UMI",
                    "length": 8,
                    "blocktype": "randomBarcode",
                    "blockname": "sequence 3",
                },
            ],
        },
        "barcode_struct_r2": {"keep_nonbarcode": True, "blocks": []},
    }


def write_old_gui_config(path: Path) -> None:
    """Write old_gui_config() the way barcodinggui saves: gzip level 9, json.dump(indent=4), no final newline."""
    with gzip.open(path, "wt", compresslevel=9, encoding="utf8") as handle:
        json.dump(old_gui_config(), handle, indent=4)


# Small hand-built configs. Layout letters: L barcode list, C constant region, R random block.
SAMPLE_BARCODES = ["CCTTAACAT", "AGCGTAGAA", "TATAGGTGT", "GGCGACTCA"]
SAMPLE_CONSTANT = "GGGTCAGTACGTACGAGTCCC"


def barcode_list(name, sequences=None, maxerrors=1, function="barcode", keep=False):
    """Return a barcodeList block in bcwithqc's key order."""
    return {
        "maxerrors": maxerrors,
        "blockfunction": function,
        "sequence": list(SAMPLE_BARCODES if sequences is None else sequences),
        "blocktype": "barcodeList",
        "blockname": name,
        "keepblock": keep,
    }


def constant_region(name, sequences=None, maxerrors=1, function="barcode", keep=False):
    """Return a constantRegion block in bcwithqc's key order."""
    return {
        "maxerrors": maxerrors,
        "blockfunction": function,
        "sequence": list([SAMPLE_CONSTANT] if sequences is None else sequences),
        "blocktype": "constantRegion",
        "blockname": name,
        "keepblock": keep,
    }


def random_block(name, length=8, maxerrors=0, function="UMI", keep=False):
    """Return a randomBarcode block in bcwithqc's key order."""
    return {
        "maxerrors": maxerrors,
        "blockfunction": function,
        "length": length,
        "blocktype": "randomBarcode",
        "blockname": name,
        "keepblock": keep,
    }


def make_config(
    read_1,
    read_2=None,
    single_end=None,
    orientation=False,
    keep_r1=False,
    keep_r2=True,
):
    """Return a config: single-end unless read_2 is given (a list, possibly empty)."""
    if single_end is None:
        single_end = read_2 is None
    config = {
        "unknown_read_orientation": orientation,
        "single_end_reads": single_end,
        "barcode_struct_r1": {"keep_nonbarcode": keep_r1, "blocks": list(read_1)},
    }
    if read_2 is not None:
        config["barcode_struct_r2"] = {
            "keep_nonbarcode": keep_r2,
            "blocks": list(read_2),
        }
    return config


def layout_blocks(layout):
    """Return blocks for a layout string such as "LRCL", named L1, R1, C1, L2, ..."""
    counts = {"L": 0, "C": 0, "R": 0}
    blocks = []
    for letter in layout:
        counts[letter] += 1
        name = f"{letter}{counts[letter]}"
        if letter == "L":
            blocks.append(barcode_list(name))
        elif letter == "C":
            blocks.append(constant_region(name))
        else:
            blocks.append(random_block(name))
    return blocks


# ---------------------------------------------------------------------------
# Acting like a user, for the UI tests. The page calls blocks "regions".

READY = "document.getElementById('app').dataset.busy === 'false'"
SIMPLE_EXAMPLE = "Simple example: barcode, linker, UMI"
ADD_BUTTONS = {
    "Barcode list": "#add_list_button",
    "Constant region": "#add_constant_button",
    "Random (UMI)": "#add_random_button",
}


def wait_ready(app):
    """Wait until the page reports it is idle (data-busy="false")."""
    app.wait_for_function(READY)


def click_button(app, name):
    """Click the button with exactly this accessible name."""
    app.get_by_role("button", name=name, exact=True).click()


def add_region(app, kind):
    """Add a region ("Barcode list", "Constant region" or "Random (UMI)") after the selected one."""
    app.click(ADD_BUTTONS[kind])


def select_region(app, name):
    """Select a region in the read diagram by its name."""
    app.get_by_role("button", name=f"Region {name}", exact=True).click()


def fill(app, label, value):
    """Fill the field with exactly this label."""
    app.get_by_label(label, exact=True).fill(str(value))


def type_sequence(app, sequence):
    """Type one sequence into the selected list's add box and press Enter."""
    box = app.get_by_label("Add a sequence", exact=True)
    box.fill(sequence)
    box.press("Enter")


def paste_sequences(app, text):
    """Paste text into the selected list's add box, as the browser delivers a paste.

    A synthetic paste event carries the text, so this works in every engine
    without clipboard permissions. Synthetic events never insert text
    themselves, so a paste the page lets through leaves the box unchanged.

    Returns:
        True when the page took the paste as a list (it prevented the
        browser's own paste), False when it let a single value through.
    """
    prevented = app.get_by_label("Add a sequence", exact=True).evaluate(
        """(box, text) => {
            const event = new Event('paste', { bubbles: true, cancelable: true });
            Object.defineProperty(event, 'clipboardData', {
                value: { getData: (type) => (type === 'text/plain' || type === 'text' ? text : '') },
            });
            box.dispatchEvent(event);
            return event.defaultPrevented;
        }""",
        text,
    )
    wait_ready(app)
    return prevented


def set_order(app, order):
    """Set the selected list's order: "A→Z" or "As entered"."""
    app.get_by_role("group", name="Order").get_by_role(
        "button", name=order, exact=True
    ).click()


def show_read(app, number):
    """Show Read 1 or Read 2 (paired-end only)."""
    app.get_by_role(
        "tab", name=f"Read {number} · {'first' if number == 1 else 'second'} file"
    ).click()


def set_rest_in_output_reads(app, choice):
    """Set the rest of the shown read to "Keep" or "Leave out" in the output reads (a read without regions can't be left out)."""
    app.locator("#rest_switch").get_by_role("button", name=choice, exact=True).click()


def set_region_in_output_reads(app, choice):
    """Set the selected region to "Keep" or "Leave out" in the output reads."""
    app.locator("#region_switch").get_by_role("button", name=choice, exact=True).click()


def answer_replace(app, asks):
    """Answer "Replace" to the replace question if it is expected (the page wasn't the welcome panel)."""
    if asks:
        app.get_by_role("dialog").get_by_role(
            "button", name="Replace", exact=True
        ).click()


def on_welcome(app):
    """True while the page shows the welcome panel."""
    return app.locator("#welcome").count() == 1


def open_template(app, title):
    """Open a template from the Templates menu, replacing what the page holds."""
    asks = not on_welcome(app)
    app.click("#templates_button")
    app.get_by_role("menuitem", name=title, exact=True).click()
    answer_replace(app, asks)
    wait_ready(app)


def open_config(app, path):
    """Choose a config file with Open…, answering the replace question if it comes.

    Questions that follow (older format, unknown keys) are left to the caller.
    """
    asks = not on_welcome(app)
    app.set_input_files("#open_input", str(path))
    answer_replace(app, asks)


def download_text(app, gzipped=False):
    """Download the config (.json, or .json.gz unzipped) and return its text."""
    with app.expect_download() as info:
        if gzipped:
            app.click("#download_menu_button")
            app.click("#download_gz_button")
        else:
            app.click("#download_button")
    data = Path(info.value.path()).read_bytes()
    return gzip.decompress(data).decode("ascii") if gzipped else data.decode("ascii")


def rule_ids(app):
    """The rule ids listed under Checks, in order."""
    return app.locator("#checks .issue").evaluate_all(
        "(items) => items.map((item) => item.dataset.rule)"
    )


def table_sequences(app):
    """The sequences shown in the selected list's table, top to bottom."""
    return app.locator("table.seqs td.seq").all_inner_texts()
