"""Tier B: the user's notes on the first published GUI.

notes_on_gui_v1.txt holds the user's notes on the first published config
builder, ledgered U01–U19 in the plan. Each test names the note it covers in
its docstring. Expected texts and values follow the behaviour agreed for each
note and are written out here, not read from the page. U11 and U19 ask for
example files and a verification checklist rather than page changes, so no
test covers them.

The user's notes on the second version (notes_on_gui_v2.txt, ledgered
V01–V15, with the review findings R01–R25) superseded U04: a read without
regions is always kept whole (V03), so its "Leave out" is disabled instead of
writing keep_nonbarcode false. V04 replaced U03's "For alignment" with "In
output reads: Keep | Leave out". Where the second version changed the texts
or controls of another note, such as the switch label (V02, V04) or the
questions' filled Cancel button (V11), these tests expect the second version.
test_gui_feedback_v2.py covers the second version's notes themselves.

Drags from outside the page are plain events carrying a dataTransfer object
(types, files, getData and a writable dropEffect), which every engine
dispatches the same way. Texts that name the tool take its name from
pyproject.toml.
"""

import base64
import json
import re
from pathlib import Path

import pytest
from gui_support import (
    ADD_BUTTONS,
    EXAMPLES_DIR,
    FIXTURES_DIR,
    SAMPLE_BARCODES,
    SAMPLE_CONSTANT,
    SIMPLE_EXAMPLE,
    add_region,
    barcode_list,
    canonical_text,
    constant_region,
    design_doc_example,
    download_text,
    fill,
    generated_with,
    make_config,
    on_welcome,
    open_template,
    paste_sequences,
    rule_ids,
    select_region,
    set_order,
    set_region_in_output_reads,
    set_rest_in_output_reads,
    show_read,
    table_sequences,
    type_sequence,
    wait_ready,
)

PAIRED = "Paired-end"
SINGLE = "Single-end"
NAME_HELP = (
    'Use letters, digits, "_", "." and "-", up to 64 characters. Two regions in '
    "the same read can't share a name, even in different case: CB1 and cb1 clash."
)
CONTIGUITY_NOTE = (
    "Regions start at base 1 and follow each other without gaps. "
    "Anything after the last region is the rest of the read."
)
# The data may not exist yet: users may be planning an experiment.
PAST_TENSE_PHRASES = (
    "were sequenced",
    "was sequenced",
    "your reads were",
    "have been sequenced",
    "has been sequenced",
    "you sequenced",
)
# U03 replaced "send ... to alignment"; V04 replaced "for alignment".
OLD_KEEP_WORDING = ("send the rest", "to alignment", "for alignment")
ALWAYS_KEPT = "A read without regions is always kept."
BARCODES_TXT = ("\n".join(SAMPLE_BARCODES) + "\n").encode("ascii")
# Two 9-nt barcodes and one 7-nt barcode.
MIXED_LENGTHS = ["CCTTAACAT", "AGCGTAGAA", "ACCTGGA"]


# ---------------------------------------------------------------------------
# Helpers


def tool_name() -> str:
    """Return the tool's name as pyproject.toml gives it, for texts that name it."""
    return generated_with()["name"]


def start(app, layout: str) -> None:
    """Choose the read layout ("Paired-end" or "Single-end") on the welcome panel."""
    app.locator("#welcome").get_by_role("button", name=layout, exact=True).click()


def visible_text(app) -> str:
    """Return the text the page shows (document.body.innerText)."""
    return app.evaluate("document.body.innerText")


def checks(app) -> list[tuple[str, str, str]]:
    """Return (rule, level, message) for each item under Checks, top to bottom."""
    items = app.locator("#checks .issue").evaluate_all(
        "(items) => items.map((item) => "
        "[item.dataset.rule, item.dataset.level, item.lastElementChild.textContent])"
    )
    return [tuple(item) for item in items]


def region_tile(app, name: str):
    """Return the read diagram's button for the region with this name."""
    return app.get_by_role("button", name=f"Region {name}", exact=True)


def diagram_regions(app) -> list[str]:
    """Return the region names in the shown read's diagram, left to right."""
    return app.locator(".track button.blk > b").all_text_contents()


def selected_region_card(app):
    """Return the selected region's card (a section labelled "Selected region")."""
    return app.get_by_role("region", name="Selected region")


def classes(locator) -> list[str]:
    """Return the class names of the element a locator finds."""
    return (locator.get_attribute("class") or "").split()


def pressed_state(group) -> dict[str, str]:
    """Return {button text: aria-pressed} for the buttons in a group, in order."""
    pairs = group.locator("button").evaluate_all(
        "(buttons) => buttons.map((button) => "
        "[button.textContent, button.getAttribute('aria-pressed')])"
    )
    return dict(pairs)


def json_preview(app) -> dict:
    """Return the config the JSON panel shows, opening the panel and closing it."""
    app.get_by_role("button", name="▸ show", exact=True).click()
    config = json.loads(app.locator("#json_preview").inner_text())
    app.get_by_role("button", name="▾ hide", exact=True).click()
    return config


def json_region_names(app, read_key: str = "barcode_struct_r1") -> list[str]:
    """Return the region names of one read in the JSON panel, in file order."""
    return [block["blockname"] for block in json_preview(app)[read_key]["blocks"]]


def output_files(app) -> list[str]:
    """Return the Outputs panel's lines about the output reads."""
    return app.locator("#out_files > span:not(.hint)").all_inner_texts()


def add_list_and_constant(app) -> None:
    """Add bc1 (SAMPLE_BARCODES) and con1 (SAMPLE_CONSTANT) to the shown read."""
    add_region(app, "Barcode list")
    assert paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    add_region(app, "Constant region")
    type_sequence(app, SAMPLE_CONSTANT)


def choose_template(app, title: str) -> None:
    """Pick a template from the Templates menu, leaving any question to the caller."""
    app.click("#templates_button")
    app.get_by_role("menuitem", name=title, exact=True).click()


def replace_question(app) -> tuple:
    """Return the open replace question and its text, checking heading and buttons.

    Cancel, the default answer, is the filled button (V11).
    """
    dialog = app.get_by_role("dialog")
    assert dialog.get_by_role("heading").inner_text() == "Replace the current config?"
    assert dialog.get_by_role("button").all_inner_texts() == ["Cancel", "Replace"]
    assert dialog.locator("button.primary").all_inner_texts() == ["Cancel"]
    return dialog, dialog.locator("p").inner_text()


def validate(core, config: dict) -> list[dict]:
    """Return validate_config's issues for a config opened as the page opens a file."""
    return core.run(
        "(core, config) => core.validate_config(core.config_from_json(config).model)",
        config,
    )


COMPUTED_STYLE = """(element, names) => {
    const style = getComputedStyle(element);
    return Object.fromEntries(names.map((name) => [name, style.getPropertyValue(name)]));
}"""

COMPUTED_RGBA = """(element, property) => {
    const value = getComputedStyle(element).getPropertyValue(property);
    const context = document.createElement('canvas').getContext('2d');
    context.fillStyle = '#010203';
    context.fillStyle = value;
    if (context.fillStyle === '#010203') throw new Error('Unreadable colour: ' + value);
    context.fillRect(0, 0, 1, 1);
    return Array.from(context.getImageData(0, 0, 1, 1).data);
}"""


def computed_style(locator, *properties: str) -> dict[str, str]:
    """Return the computed values of CSS properties on the element a locator finds."""
    return locator.evaluate(COMPUTED_STYLE, list(properties))


def computed_rgba(locator, css_property: str) -> list[int]:
    """Return a computed colour of an element as 8-bit [red, green, blue, alpha].

    Engines write computed colours differently (rgb(), rgba(), or
    color(srgb ...) for color-mix() results), so the browser resolves the
    value by painting it onto a 1 x 1 canvas.
    """
    return locator.evaluate(COMPUTED_RGBA, css_property)


DRAG_FROM_OUTSIDE = """(target, drag) => {
    const files = drag.files.map(([name, base64]) => new File(
        [Uint8Array.from(atob(base64), (character) => character.charCodeAt(0))], name));
    const types = (files.length ? ['Files'] : []).concat(drag.text === null ? [] : ['text/plain']);
    const data_transfer = {
        types,
        files,
        dropEffect: 'none',
        effectAllowed: 'all',
        getData: (type) => (drag.text !== null && (type === 'text/plain' || type === 'text') ? drag.text : ''),
    };
    return drag.events.map((name) => {
        const event = new Event(name, { bubbles: true, cancelable: true });
        Object.defineProperty(event, 'dataTransfer', { value: data_transfer });
        target.dispatchEvent(event);
        return { event: name, prevented: event.defaultPrevented, drop_effect: data_transfer.dropEffect };
    });
}"""


def drag_from_outside(
    locator, files=(), text=None, events=("dragenter", "dragover", "drop")
) -> dict:
    """Drag files or text from outside the page onto an element, as a user does.

    Args:
        locator: The element the drag ends on; the events bubble up from it.
        files: (file name, contents as bytes) pairs the drag carries.
        text: Plain text the drag carries, or None.
        events: Drag events to dispatch in order, sharing one dataTransfer.

    Returns:
        {event name: {"prevented": bool, "drop_effect": str}}, read right after
        each event: whether a handler called preventDefault, and the
        dataTransfer's dropEffect by then.
    """
    drag = {
        "files": [
            [name, base64.b64encode(data).decode("ascii")] for name, data in files
        ],
        "text": text,
        "events": list(events),
    }
    return {
        result["event"]: result for result in locator.evaluate(DRAG_FROM_OUTSIDE, drag)
    }


def is_busy(app) -> bool:
    """True while the page reports work in progress (data-busy="true")."""
    return app.locator("#app").get_attribute("data-busy") == "true"


# ---------------------------------------------------------------------------
# U01: a welcoming first page; steps not done yet are to-do items, not errors


def welcome_steps() -> list[str]:
    """Return the welcome panel's three steps."""
    return [
        "Choose the read layout.",
        (
            "Add the regions of each read: barcode lists, constant regions and "
            "random sequences such as UMIs."
        ),
        f"Download the .json and give it to {tool_name()} with --config.",
    ]


def test_u01_fresh_page_shows_the_welcome_panel(app):
    """U01: a fresh page shows a welcome panel, with no error anywhere.

    The panel says what the tool does, lists three steps and offers the first
    choice; the side rail stays hidden and Download waits for a read layout.
    """
    welcome = app.locator("#main section#welcome")
    assert welcome.is_visible()
    assert (
        welcome.get_by_role("heading", level=2).inner_text()
        == "Welcome to the barcode configurator"
    )
    text = welcome.inner_text()
    lead = (
        f"This page builds the config file that tells {tool_name()} where the "
        "barcodes, UMIs and constant sequences sit in each read."
    )
    assert lead in text
    assert welcome.locator("ol > li").all_inner_texts() == welcome_steps()
    assert "What you enter stays in this browser tab." in text
    assert "Choose the read layout to start:" in text
    assert "Or start from a template, or open a config you already have." in text
    assert "block" not in text.lower()

    assert welcome.locator("button.choice").count() == 2
    for layout in (PAIRED, SINGLE):
        choice = welcome.get_by_role("button", name=layout, exact=True)
        assert "choice" in classes(choice)

    assert app.locator("#rail").is_hidden()
    assert "welcome-mode" in classes(app.locator("#layout"))
    download = app.locator("#download_button")
    assert download.is_disabled()
    assert download.get_attribute("title") == "Choose the read layout first"
    # No red outline and no error: nothing has the old .need class.
    assert app.locator(".need, .chip.err, .has-error").count() == 0


def test_u01_welcome_text_does_not_assume_the_reads_exist(app):
    """U01: the first page doesn't assume the reads exist yet."""
    text = visible_text(app).lower()
    assert [phrase for phrase in PAST_TENSE_PHRASES if phrase in text] == []


@pytest.mark.parametrize(
    ("layout", "tabs", "to_do"),
    [
        (
            PAIRED,
            ["Read 1 · first file", "Read 2 · second file"],
            "Add at least one region to Read 1 or Read 2.",
        ),
        (SINGLE, [], "Add at least one region to Read 1."),
    ],
)
def test_u01_choosing_a_layout_leaves_only_a_to_do(app, layout, tabs, to_do):
    """U01: after the first choice, the next step is a to-do, not an error."""
    start(app, layout)
    assert not on_welcome(app)
    assert app.locator("#rail").is_visible()
    assert app.get_by_role("tab").all_inner_texts() == tabs
    layout_choice = pressed_state(app.get_by_role("group", name="Read layout"))
    assert [name for name, pressed in layout_choice.items() if pressed == "true"] == [
        layout
    ]
    assert checks(app) == [("E14", "todo", to_do)]
    chip = app.locator("#checks .issue .chip")
    assert chip.inner_text() == "to do"
    assert "todo" in classes(chip)
    assert app.locator("#checks_summary").inner_text() == "1 to do"
    assert app.locator(".need, .chip.err, .has-error").count() == 0
    download = app.locator("#download_button")
    assert download.is_disabled()
    assert download.get_attribute("title") == "Finish the to-do items in Checks first"


def test_u01_welcome_links_open_the_templates_menu_and_a_file(app):
    """U01: the welcome panel's links open the Templates menu and a file chooser."""
    app.click("#welcome_templates")
    assert app.locator("#templates_menu").is_visible()
    assert app.locator("#templates_button").get_attribute("aria-expanded") == "true"
    app.keyboard.press("Escape")
    assert app.locator("#templates_menu").is_hidden()

    with app.expect_file_chooser() as chooser_info:
        app.click("#welcome_open")
    chooser_info.value.set_files(str(EXAMPLES_DIR / "se_mini_config.json"))
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    region_tile(app, "cb1").wait_for()
    assert app.get_by_role("tab").count() == 0


def test_u01_u17_new_empty_config_asks_then_shows_the_welcome_again(app):
    """U01, U17: "New empty config" asks first, then brings back the welcome panel."""
    open_template(app, "SDR-seq RNA")
    choose_template(app, "New empty config")
    dialog, text = replace_question(app)
    assert text == "Starting a new empty config replaces the config on this page."
    dialog.get_by_role("button", name="Replace", exact=True).click()
    app.locator("#main section#welcome").wait_for()
    assert app.locator("#rail").is_hidden()
    assert (
        app.locator("#download_button").get_attribute("title")
        == "Choose the read layout first"
    )


def test_u01_a_file_without_a_read_layout_shows_it_as_a_to_do(app, tmp_path):
    """U01: a file that leaves the read layout open shows it as a to-do, not in red."""
    config = make_config([barcode_list("bc1"), constant_region("con1")], [])
    config["single_end_reads"] = None
    path = tmp_path / "no_layout.json"
    path.write_text(json.dumps(config, indent=4))
    app.set_input_files("#open_input", str(path))
    wait_ready(app)
    assert not on_welcome(app)
    assert checks(app) == [("E01", "todo", "Choose paired-end or single-end reads.")]
    assert pressed_state(app.get_by_role("group", name="Read layout")) == {
        PAIRED: "false",
        SINGLE: "false",
    }
    assert app.locator(".need, .chip.err, .has-error").count() == 0
    assert (
        app.locator("#download_button").get_attribute("title")
        == "Finish the to-do items in Checks first"
    )


def test_u01_an_empty_list_is_a_to_do_and_a_duplicate_is_an_error(app):
    """U01: an empty list is a to-do and a duplicate barcode an error.

    Both block the download; the button's title says which kind stands in the way.
    """
    start(app, SINGLE)
    add_region(app, "Barcode list")
    assert [(rule, level) for rule, level, _ in checks(app)] == [
        ("E03", "todo"),
        ("W11", "warning"),
    ]
    issue = app.locator('#checks .issue[data-rule="E03"]')
    assert issue.get_attribute("data-level") == "todo"
    assert issue.locator(".chip").inner_text() == "to do"
    assert "todo" in classes(issue.locator(".chip"))
    assert app.locator("#checks_summary").inner_text() == "1 to do · 1 warning"
    tile_classes = classes(region_tile(app, "bc1"))
    assert "has-todo" in tile_classes
    assert "has-error" not in tile_classes
    download = app.locator("#download_button")
    assert download.is_disabled()
    assert download.get_attribute("title") == "Finish the to-do items in Checks first"

    type_sequence(app, "CCTTAACAT")
    type_sequence(app, "CCTTAACAT")
    assert [(rule, level) for rule, level, _ in checks(app)] == [
        ("E04", "error"),
        ("W11", "warning"),
    ]
    assert app.locator("#checks_summary").inner_text() == "1 error · 1 warning"
    tile_classes = classes(region_tile(app, "bc1"))
    assert "has-error" in tile_classes
    assert "has-todo" not in tile_classes
    assert download.is_disabled()
    assert download.get_attribute("title") == "Fix the errors in Checks first"


def test_u01_checks_list_errors_then_to_dos_then_warnings(app):
    """U01: Checks lists errors, then to-dos, then warnings, and counts each kind."""
    start(app, SINGLE)
    add_region(app, "Random (UMI)")
    app.get_by_role("combobox", name="Function").select_option("barcode")
    add_region(app, "Constant region")
    add_region(app, "Barcode list")
    type_sequence(app, "CCTTAACAT")
    type_sequence(app, "CCTTAACAT")
    assert diagram_regions(app) == ["umi1", "con1", "bc1"]
    assert [(rule, level) for rule, level, _ in checks(app)] == [
        ("E04", "error"),
        ("E03", "todo"),
        ("W06", "warning"),
    ]
    assert (
        app.locator("#checks_summary").inner_text() == "1 error · 1 to do · 1 warning"
    )
    assert "has-warning" in classes(region_tile(app, "umi1"))
    assert "has-todo" in classes(region_tile(app, "con1"))
    assert "has-error" in classes(region_tile(app, "bc1"))


# ---------------------------------------------------------------------------
# U02: both reads start with the same default


def test_u02_both_reads_default_to_keep(app):
    """U02: right after choosing paired-end, the rest of each read is set to "Keep"."""
    start(app, PAIRED)
    for number in (1, 2):
        show_read(app, number)
        state = pressed_state(app.locator("#rest_switch"))
        assert list(state) == ["Keep", "Leave out"]
        assert state == {"Keep": "true", "Leave out": "false"}


def test_u02_a_config_built_from_a_fresh_page_outputs_both_reads(app):
    """U02: unless changed, a config from a fresh page outputs the rest of both reads.

    The switch changes only the read it belongs to.
    """
    start(app, PAIRED)
    add_list_and_constant(app)
    config = json.loads(download_text(app))
    assert config["barcode_struct_r1"]["keep_nonbarcode"] is True
    assert config["barcode_struct_r2"]["keep_nonbarcode"] is True

    set_rest_in_output_reads(app, "Leave out")
    config = json.loads(download_text(app))
    assert config["barcode_struct_r1"]["keep_nonbarcode"] is False
    assert config["barcode_struct_r2"]["keep_nonbarcode"] is True


# ---------------------------------------------------------------------------
# U03, V04: "In output reads", not "send ... to alignment" or "for alignment"


def test_u03_v04_keep_controls_say_in_output_reads(app):
    """U03, V04: the page says "In output reads:", never "for alignment" or "send ... to alignment"."""
    start(app, PAIRED)
    text = visible_text(app)
    assert text.count("In output reads:") == 1
    assert [wording for wording in OLD_KEEP_WORDING if wording in text.lower()] == []

    add_region(app, "Barcode list")
    text = visible_text(app)
    assert text.count("In output reads:") == 2
    assert [wording for wording in OLD_KEEP_WORDING if wording in text.lower()] == []

    open_template(app, SIMPLE_EXAMPLE)
    text = visible_text(app)
    assert [wording for wording in OLD_KEEP_WORDING if wording in text.lower()] == []


def test_u03_e11_names_the_flags_by_what_they_output(core):
    """U03, V04: E11 names the two flags by their label, "In output reads"."""
    config = make_config([barcode_list("bc1", keep="no"), constant_region("con1")])
    config["barcode_struct_r1"]["keep_nonbarcode"] = "yes"
    messages = [
        issue["message"] for issue in validate(core, config) if issue["rule"] == "E11"
    ]
    assert len(messages) == 2
    rest_message, region_message = messages
    assert rest_message.startswith("Read 1: ")
    assert (
        '"In output reads" for the rest of the read must be Keep or Leave out'
        in rest_message
    )
    assert region_message.startswith("bc1: ")
    assert '"In output reads" must be Keep or Leave out' in region_message
    assert [
        message
        for message in messages
        if any(wording in message.lower() for wording in OLD_KEEP_WORDING)
    ] == []


# ---------------------------------------------------------------------------
# U04, U05: a read without regions, the same on either read


@pytest.mark.parametrize("number", [1, 2])
def test_u04_u05_a_read_without_regions_reads_the_same_on_either_read(app, number):
    """U04, U05: a read without regions reads the same on Read 1 and Read 2.

    Both say all of the read is kept. The switch reads "Rest of Read N", as
    on a read with regions (V02), and its "Leave out" is disabled and says
    why: a read without regions is always kept (V03).
    """
    start(app, PAIRED)
    show_read(app, number)
    empty = app.locator("#empty_read")
    assert empty.locator("b").inner_text() == "No regions yet"
    assert (
        empty.locator("span").inner_text()
        == f"All of Read {number} is kept in the output reads."
    )
    align_row = app.locator(".align-row", has=app.locator("#rest_switch"))
    assert align_row.locator(".label").inner_text() == f"Rest of Read {number}"
    switch = app.get_by_role(
        "group", name=f"Rest of Read {number} in output reads", exact=True
    )
    assert switch.is_visible()
    assert pressed_state(switch) == {"Keep": "true", "Leave out": "false"}
    leave_out = switch.get_by_role("button", name="Leave out", exact=True)
    assert leave_out.is_disabled()
    assert leave_out.get_attribute("title") == ALWAYS_KEPT
    assert align_row.locator("#rest_always_kept").inner_text() == ALWAYS_KEPT
    assert output_files(app) == ["Read 1: the whole read", "Read 2: the whole read"]


def test_u04_v03_a_read_without_regions_is_always_kept(app):
    """U04, V03: a read without regions can't be left out, and the file keeps all of it.

    U04 let "Leave out" write keep_nonbarcode false for such a read; V03
    replaced that. Read 2's rest is left out while it holds a region; once
    that region is deleted, Read 2 is kept whole again, its "Leave out" is
    disabled, and the file says keep_nonbarcode true. A read with regions that
    keeps nothing still gets the note that nothing is kept.
    """
    start(app, PAIRED)
    add_list_and_constant(app)
    show_read(app, 2)
    add_region(app, "Random (UMI)")
    set_rest_in_output_reads(app, "Leave out")
    assert json_preview(app)["barcode_struct_r2"]["keep_nonbarcode"] is False

    selected_region_card(app).get_by_role(
        "button", name="Delete region", exact=True
    ).click()
    app.get_by_role("dialog").get_by_role(
        "button", name="Delete region", exact=True
    ).click()
    assert diagram_regions(app) == []
    switch = app.locator("#rest_switch")
    assert pressed_state(switch) == {"Keep": "true", "Leave out": "false"}
    assert switch.get_by_role("button", name="Leave out", exact=True).is_disabled()
    assert output_files(app) == ["Read 1: rest of read", "Read 2: the whole read"]
    config = json.loads(download_text(app))
    assert config["barcode_struct_r1"]["keep_nonbarcode"] is True
    assert config["barcode_struct_r2"] == {"keep_nonbarcode": True, "blocks": []}

    # A read with regions that keeps nothing at all gets the other note.
    show_read(app, 1)
    set_rest_in_output_reads(app, "Leave out")
    assert output_files(app)[0] == (
        f"Read 1: no file (nothing is kept); current {tool_name()} writes a file of empty reads"
    )


# ---------------------------------------------------------------------------
# U06: adding regions is obvious


def test_u06_add_buttons_show_without_a_menu_in_their_region_colours(app):
    """U06: the three add buttons are always shown, each in its region type's colour."""
    start(app, SINGLE)
    assert app.locator(".menu:not([hidden])").count() == 0
    assert app.locator("#add_block_button").count() == 0
    expected = {
        "Barcode list": ("+ Barcode list", "list", "bc1"),
        "Constant region": ("+ Constant region", "const", "con1"),
        "Random (UMI)": ("+ Random (UMI)", "rand", "umi1"),
    }
    for kind, (label, type_class, _) in expected.items():
        button = app.locator(ADD_BUTTONS[kind])
        assert button.is_visible()
        assert button.inner_text() == label
        assert {"add", type_class} <= set(classes(button))

    for kind in expected:
        add_region(app, kind)
    backgrounds = []
    for kind, (_, _, name) in expected.items():
        button_colour = computed_rgba(
            app.locator(ADD_BUTTONS[kind]), "background-color"
        )
        assert button_colour == computed_rgba(
            region_tile(app, name), "background-color"
        )
        backgrounds.append(tuple(button_colour))
    assert len(set(backgrounds)) == 3


def test_u06_add_label_says_where_the_new_region_goes(app):
    """U06: the label says where a new region goes, and it lands there, selected."""
    start(app, PAIRED)
    where = app.locator("#add_region_where")
    assert where.inner_text() == "Add a region to Read 1:"
    group = app.get_by_role("group", name="Add a region to Read 1:", exact=True)
    assert group.get_by_role("button").all_inner_texts() == [
        "+ Barcode list",
        "+ Constant region",
        "+ Random (UMI)",
    ]

    add_region(app, "Barcode list")
    assert where.inner_text() == "Add a region after bc1:"
    add_region(app, "Constant region")
    select_region(app, "bc1")
    assert where.inner_text() == "Add a region after bc1:"
    add_region(app, "Random (UMI)")
    assert diagram_regions(app) == ["bc1", "umi1", "con1"]
    assert region_tile(app, "umi1").get_attribute("aria-pressed") == "true"
    assert app.get_by_label("Name", exact=True).input_value() == "umi1"

    show_read(app, 2)
    assert where.inner_text() == "Add a region to Read 2:"
    show_read(app, 1)
    assert where.inner_text() == "Add a region at the end of Read 1:"
    add_region(app, "Barcode list")
    assert diagram_regions(app) == ["bc1", "umi1", "con1", "bc2"]

    # A region without a name is called by its position.
    fill(app, "Name", "")
    assert where.inner_text() == "Add a region after region 4:"


# ---------------------------------------------------------------------------
# U07: the simple example template


def test_u07_templates_menu_names_the_simple_example(app):
    """U07: the simple example has a plain name; "Design-doc example" is gone."""
    app.click("#templates_button")
    items = app.locator("#templates_menu").get_by_role("menuitem").all_inner_texts()
    assert items == [
        "New empty config",
        "SDR-seq RNA",
        "SDR-seq gDNA",
        "Simple example: barcode, linker, UMI",
    ]
    assert re.search(r"design[- ]doc", visible_text(app), flags=re.IGNORECASE) is None
    templates = app.evaluate(
        "window.barcode_templates.map((template) => [template.id, template.title])"
    )
    assert templates == [
        ["sdr_seq_rna", "SDR-seq RNA"],
        ["sdr_seq_gdna", "SDR-seq gDNA"],
        ["simple_example", "Simple example: barcode, linker, UMI"],
    ]


def test_u07_simple_example_downloads_the_documented_example(app):
    """U07: the simple example downloads as simple_example.json.

    The file holds the example config in doc/config_format.md.
    """
    open_template(app, SIMPLE_EXAMPLE)
    with app.expect_download() as info:
        app.click("#download_button")
    download = info.value
    assert download.suggested_filename == "simple_example.json"
    assert Path(download.path()).read_text(encoding="ascii") == canonical_text(
        design_doc_example()
    )


# ---------------------------------------------------------------------------
# U08: menus stay inside the window

MENUS = [
    ("#templates_button", "#templates_menu"),
    ("#download_menu_button", "#download_menu"),
    ("#list_more_button", "#list_more_menu"),
]


@pytest.mark.parametrize("width", [320, 420, 1280])
def test_u08_menus_open_inside_the_window(app, width):
    """U08: every menu opens fully inside the window, on narrow windows too."""
    app.set_viewport_size({"width": width, "height": 900})
    start(app, SINGLE)
    add_region(app, "Barcode list")
    window_width = app.evaluate("document.documentElement.clientWidth")
    outside = []
    for button, menu in MENUS:
        app.click(button)
        box = app.locator(menu).bounding_box()
        assert box is not None, f"{menu} did not open"
        left, right = box["x"], box["x"] + box["width"]
        if left < 0 or right > window_width:
            outside.append(f"{menu} spans {left:.0f}..{right:.0f} px of {window_width}")
        app.keyboard.press("Escape")
        app.locator(menu).wait_for(state="hidden")
    assert outside == []


# ---------------------------------------------------------------------------
# U09: drag and drop


def test_u09_an_empty_list_offers_a_drop_area(app, tmp_path):
    """U09: an empty list shows where to drop a file, with a link that chooses one."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    dropzone = app.locator("#list_dropzone")
    lines = dropzone.inner_text().splitlines()
    assert (
        lines[0] == "Drop a file here (CSV, TSV, TXT, Excel or ODS), or choose a file."
    )
    assert "paste" in lines[1]
    path = tmp_path / "barcodes.txt"
    path.write_bytes(BARCODES_TXT)
    with app.expect_file_chooser() as chooser_info:
        dropzone.get_by_role("button", name="choose a file", exact=True).click()
    chooser_info.value.set_files(str(path))
    wait_ready(app)
    assert table_sequences(app) == sorted(SAMPLE_BARCODES)
    assert dropzone.count() == 0


def test_u09_a_dropped_text_file_fills_an_empty_list_without_a_question(app):
    """U09: a .txt dropped on an empty list is imported without a question."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    results = drag_from_outside(
        app.locator("#list_dropzone"), files=[("barcodes.txt", BARCODES_TXT)]
    )
    assert results["drop"]["prevented"]
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    assert table_sequences(app) == sorted(SAMPLE_BARCODES)


@pytest.mark.parametrize(
    ("answer", "expected"),
    [
        ("Add to the existing 2", SAMPLE_BARCODES),
        ("Replace them", SAMPLE_BARCODES[2:]),
    ],
)
def test_u09_a_csv_dropped_on_a_list_with_entries_asks_add_or_replace(
    app, answer, expected
):
    """U09: a file dropped on a list with entries asks: add to them or replace them."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    for sequence in SAMPLE_BARCODES[:2]:
        type_sequence(app, sequence)
    csv = "name,barcode\n" + "".join(
        f"x{index},{sequence}\n" for index, sequence in enumerate(SAMPLE_BARCODES[2:])
    )
    drag_from_outside(
        selected_region_card(app), files=[("more.csv", csv.encode("ascii"))]
    )
    dialog = app.get_by_role("dialog")
    assert (
        dialog.get_by_role("heading").inner_text() == "Add to the list or replace it?"
    )
    dialog.get_by_role("button", name=answer, exact=True).click()
    wait_ready(app)
    assert table_sequences(app) == sorted(expected)


def test_u09_a_dropped_spreadsheet_fills_an_empty_list(app, tmp_path):
    """U09: an .xlsx dropped on an empty list is imported like Import file…."""
    from openpyxl import Workbook

    workbook = Workbook()
    sheet = workbook.active
    sheet.append(["barcode"])
    for sequence in SAMPLE_BARCODES:
        sheet.append([sequence])
    path = tmp_path / "plate.xlsx"
    workbook.save(path)

    start(app, SINGLE)
    add_region(app, "Barcode list")
    drag_from_outside(
        app.locator("#list_dropzone"), files=[("plate.xlsx", path.read_bytes())]
    )
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    assert table_sequences(app) == sorted(SAMPLE_BARCODES)


def test_u09_a_constant_region_takes_dropped_variants(app):
    """U09: a constant region's card takes a dropped file; variants keep their order."""
    variants = ["GGGTCAGTACGTACGAGTCCC", "GGGTCAGTACGTACGAGTCCA"]
    start(app, SINGLE)
    add_region(app, "Constant region")
    drag_from_outside(
        app.locator("#list_dropzone"),
        files=[("variants.txt", "\n".join(variants).encode())],
    )
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    assert table_sequences(app) == variants


def test_u09_two_dropped_files_import_nothing(app):
    """U09: dropping two files at once asks for one at a time and imports nothing."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    drag_from_outside(
        app.locator("#list_dropzone"),
        files=[("a.txt", BARCODES_TXT), ("b.txt", b"GGGGAAAA\nTTTTCCCC\n")],
    )
    dialog = app.get_by_role("dialog")
    assert dialog.get_by_role("heading").inner_text() == "Nothing imported"
    assert dialog.locator("p").inner_text() == "Drop one file at a time."
    dialog.get_by_role("button", name="OK", exact=True).click()
    wait_ready(app)
    assert app.locator("#sequence_count").inner_text() == "0 barcodes"
    assert app.locator("#list_dropzone").is_visible()


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("AAAACCCC\nCCCCAAAA\nGGGGTTTT\n", ["AAAACCCC", "CCCCAAAA", "GGGGTTTT"]),
        ("AAAACCCC, CCCCAAAA GGGGTTTT", ["AAAACCCC", "CCCCAAAA", "GGGGTTTT"]),
        ("AAAACCCC", ["AAAACCCC"]),
    ],
)
def test_u09_dropped_text_fills_an_empty_list(app, text, expected):
    """U09: text dropped on an empty list is imported without a question."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    results = drag_from_outside(selected_region_card(app), text=text)
    assert results["drop"]["prevented"]
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    assert table_sequences(app) == expected


def test_u09_dropped_text_asks_add_or_replace_when_the_list_has_entries(app):
    """U09: text dropped on a list with entries asks, unlike a paste into the box."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    type_sequence(app, SAMPLE_BARCODES[0])
    drag_from_outside(selected_region_card(app), text="\n".join(SAMPLE_BARCODES[1:]))
    dialog = app.get_by_role("dialog")
    dialog.get_by_role("button", name="Add to the existing 1", exact=True).click()
    wait_ready(app)
    assert table_sequences(app) == sorted(SAMPLE_BARCODES)


def test_u09_the_card_shows_it_takes_a_drop_while_dragged_over(app):
    """U09: a list's card is marked as the drop target while a file is over it."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    card = app.locator("section.block-card")
    assert "drop-target" not in classes(card)
    results = drag_from_outside(
        app.locator("#list_dropzone"),
        files=[("barcodes.txt", BARCODES_TXT)],
        events=("dragenter", "dragover"),
    )
    assert results["dragover"]["prevented"]
    assert results["dragover"]["drop_effect"] != "none"
    assert "drop-target" in classes(card)
    drag_from_outside(
        app.locator("#list_dropzone"),
        files=[("barcodes.txt", BARCODES_TXT)],
        events=("dragleave",),
    )
    assert "drop-target" not in classes(card)


@pytest.mark.parametrize("where", ["welcome panel", "read diagram", "random region"])
def test_u09_a_file_dropped_outside_a_list_is_refused(app, where):
    """U09: a file dropped anywhere but a list's card is refused and imports nothing.

    Refusing the drop keeps the browser from leaving the page to show the file.
    """
    if where == "welcome panel":
        target = app.locator("#welcome")
    else:
        start(app, SINGLE)
        add_region(app, "Barcode list")
        if where == "read diagram":
            target = region_tile(app, "bc1")
        else:
            add_region(app, "Random (UMI)")
            target = selected_region_card(app)
    results = drag_from_outside(target, files=[("barcodes.txt", BARCODES_TXT)])
    assert results["dragover"]["prevented"]
    assert results["dragover"]["drop_effect"] == "none"
    assert results["drop"]["prevented"]
    assert not is_busy(app)
    assert app.get_by_role("dialog").count() == 0
    assert app.locator(".drop-target").count() == 0
    if where != "welcome panel":
        assert json_preview(app)["barcode_struct_r1"]["blocks"][0]["sequence"] == []


def test_u09_dragging_a_region_still_reorders_and_imports_nothing(app):
    """U09: a dragged region still reorders, and a list's card doesn't import it."""
    start(app, SINGLE)
    add_list_and_constant(app)
    add_region(app, "Random (UMI)")
    region_tile(app, "umi1").drag_to(region_tile(app, "bc1"))
    assert json_region_names(app) == ["umi1", "bc1", "con1"]

    select_region(app, "bc1")
    region_tile(app, "con1").drag_to(app.locator("table.seqs"))
    assert app.get_by_role("dialog").count() == 0
    assert "drop-target" not in classes(selected_region_card(app))
    assert table_sequences(app) == sorted(SAMPLE_BARCODES)
    assert json_region_names(app) == ["umi1", "bc1", "con1"]


# ---------------------------------------------------------------------------
# U10: E05 says how to fix it


@pytest.mark.parametrize("followed_by", [None, "Random (UMI)"])
def test_u10_e05_says_to_add_a_constant_region_after_the_list(app, followed_by):
    """U10: E05 says to add a constant region after the list; doing so clears it."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    for sequence in MIXED_LENGTHS:
        type_sequence(app, sequence)
    if followed_by:
        add_region(app, followed_by)
        select_region(app, "bc1")
    e05 = [(level, message) for rule, level, message in checks(app) if rule == "E05"]
    assert e05 == [
        (
            "error",
            (
                "bc1 has barcodes of different lengths (9 nt ×2 · 7 nt ×1). A list "
                "like that needs a constant region right after it, to its right. Add "
                "a constant region after bc1, or give all its barcodes the same length."
            ),
        )
    ]

    add_region(app, "Constant region")
    assert diagram_regions(app)[:2] == ["bc1", "con1"]
    assert "E05" not in rule_ids(app)
    assert [rule for rule, level, _ in checks(app) if level == "error"] == []


# ---------------------------------------------------------------------------
# U12: pasting a list into the add box


def test_u12_lists_have_import_file_and_no_add_menu(app, tmp_path):
    """U12: a list has "Import file…" and a paste-aware add box, not an add menu."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    assert app.locator("#add_sequences_button").count() == 0
    assert app.locator("#import_button").is_visible()
    assert app.locator("#import_button").inner_text() == "Import file…"
    box = app.get_by_label("Add a sequence", exact=True)
    assert box.get_attribute("placeholder") == "Type a barcode, or paste a list"
    path = tmp_path / "barcodes.txt"
    path.write_bytes(BARCODES_TXT)
    with app.expect_file_chooser() as chooser_info:
        app.click("#import_button")
    chooser_info.value.set_files(str(path))
    wait_ready(app)
    assert table_sequences(app) == sorted(SAMPLE_BARCODES)


def test_u12_pasting_lines_into_an_empty_list_imports_them_without_a_question(app):
    """U12: pasting several lines into the add box imports each line as a barcode."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    assert paste_sequences(app, "\n".join(SAMPLE_BARCODES)) is True
    assert app.get_by_role("dialog").count() == 0
    assert table_sequences(app) == sorted(SAMPLE_BARCODES)


def test_u12_pasting_into_a_list_with_entries_appends_without_asking(app):
    """U12: a paste adds after the existing entries without asking add or replace."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    type_sequence(app, SAMPLE_BARCODES[0])
    assert paste_sequences(app, "\n".join(SAMPLE_BARCODES[1:])) is True
    assert app.get_by_role("dialog").count() == 0
    set_order(app, "As entered")
    assert table_sequences(app) == SAMPLE_BARCODES


@pytest.mark.parametrize(
    "text",
    [
        "AAAACCCC, CCCCAAAA GGGGTTTT",
        "AAAACCCC\tCCCCAAAA\tGGGGTTTT",
        "AAAACCCC;CCCCAAAA;GGGGTTTT",
    ],
)
def test_u12_one_pasted_line_of_several_barcodes_imports_each(app, text):
    """U12: one pasted line of several barcodes imports each, whatever separates them."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    assert paste_sequences(app, text) is True
    assert app.get_by_role("dialog").count() == 0
    assert table_sequences(app) == ["AAAACCCC", "CCCCAAAA", "GGGGTTTT"]


@pytest.mark.parametrize("text", ["CCTTAACAT", "CCTTAACAT\n", "AAAACCCC hello"])
def test_u12_a_paste_that_is_not_a_list_is_left_to_the_browser(app, text):
    """U12: a single value, or a line that isn't all barcodes, pastes as typed."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    assert paste_sequences(app, text) is False
    assert app.get_by_role("dialog").count() == 0
    assert app.locator("#sequence_count").inner_text() == "0 barcodes"


@pytest.mark.parametrize(
    ("answer", "expected"),
    [("A header", SAMPLE_BARCODES), ("A barcode", ["TAG", *SAMPLE_BARCODES])],
)
def test_u12_a_pasted_list_starting_with_tag_asks_about_a_header(app, answer, expected):
    """U12: a pasted list whose first entry may be a header still asks about it."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    assert paste_sequences(app, "TAG\n" + "\n".join(SAMPLE_BARCODES)) is True
    dialog = app.get_by_role("dialog")
    assert (
        dialog.get_by_role("heading").inner_text() == 'Is "TAG" a header or a barcode?'
    )
    dialog.get_by_role("button", name=answer, exact=True).click()
    set_order(app, "As entered")
    assert table_sequences(app) == expected


def test_u12_a_pasted_table_with_two_dna_columns_asks_for_the_column(app):
    """U12: a pasted table with two columns of DNA still asks which column to import."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    assert paste_sequences(app, "CCTTAACAT,ATGTTAAGG\nAGCGTAGAA,TTCTACGCT") is True
    dialog = app.get_by_role("dialog")
    import_button = dialog.get_by_role("button", name="Import", exact=True)
    assert import_button.is_disabled()
    dialog.get_by_label("Column 2").check()
    import_button.click()
    app.locator("table.seqs").wait_for()
    assert table_sequences(app) == ["ATGTTAAGG", "TTCTACGCT"]


def test_u12_pasting_several_variants_into_a_constant_region(app):
    """U12: several pasted lines become a constant region's variants, in order."""
    variants = ["GGGTCAGTACGTACGAGTCCC", "GGGTCAGTACGTACGAGTCCA"]
    start(app, SINGLE)
    add_region(app, "Constant region")
    assert paste_sequences(app, "\n".join(variants)) is True
    assert table_sequences(app) == variants


# ---------------------------------------------------------------------------
# U13: the keep/drop choice is a visible two-way switch


def test_u13_new_regions_start_left_out_on_a_named_switch(app):
    """U13: every new region shows "Keep" and "Leave out", with "Leave out" chosen."""
    start(app, SINGLE)
    for kind, name in (
        ("Barcode list", "bc1"),
        ("Constant region", "con1"),
        ("Random (UMI)", "umi1"),
    ):
        add_region(app, kind)
        assert app.get_by_role(
            "group", name=f"{name} in output reads", exact=True
        ).is_visible()
        state = pressed_state(app.locator("#region_switch"))
        assert list(state) == ["Keep", "Leave out"]
        assert state == {"Keep": "false", "Leave out": "true"}
    blocks = json_preview(app)["barcode_struct_r1"]["blocks"]
    assert [block["keepblock"] for block in blocks] == [False, False, False]


def test_u13_keep_sets_keepblock_and_marks_the_tile(app):
    """U13: "Keep" sets keepblock and labels the tile; the chosen side stands out."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    add_region(app, "Constant region")
    tile = region_tile(app, "con1")
    assert "in output reads" not in tile.inner_text()

    set_region_in_output_reads(app, "Keep")
    switch = app.locator("#region_switch")
    assert pressed_state(switch) == {"Keep": "true", "Leave out": "false"}
    assert "in output reads" in tile.inner_text().splitlines()
    blocks = json_preview(app)["barcode_struct_r1"]["blocks"]
    assert [block["keepblock"] for block in blocks] == [False, True]
    chosen = switch.get_by_role("button", name="Keep", exact=True)
    other = switch.get_by_role("button", name="Leave out", exact=True)
    assert computed_rgba(chosen, "background-color") != computed_rgba(
        other, "background-color"
    )

    set_region_in_output_reads(app, "Leave out")
    assert "in output reads" not in tile.inner_text()
    blocks = json_preview(app)["barcode_struct_r1"]["blocks"]
    assert [block["keepblock"] for block in blocks] == [False, False]


def test_u13_rest_switch_updates_the_rest_tile_and_the_file(app):
    """U13: the rest-of-read switch shows its choice on the rest tile and in the file."""
    start(app, SINGLE)
    add_list_and_constant(app)
    assert app.get_by_role(
        "group", name="Rest of Read 1 in output reads", exact=True
    ).is_visible()
    rest_tile = app.locator(".track .blk.rest")
    assert rest_tile.locator("b").inner_text() == "rest of read"
    assert rest_tile.locator("span").inner_text() == "in output reads"
    switch = app.locator("#rest_switch")
    assert computed_rgba(
        switch.get_by_role("button", name="Keep", exact=True), "background-color"
    ) != computed_rgba(
        switch.get_by_role("button", name="Leave out", exact=True), "background-color"
    )

    set_rest_in_output_reads(app, "Leave out")
    assert pressed_state(switch) == {"Keep": "false", "Leave out": "true"}
    assert rest_tile.locator("span").inner_text() == "left out"
    assert json_preview(app)["barcode_struct_r1"]["keep_nonbarcode"] is False


def test_u13_a_value_that_is_neither_presses_neither_side(app, tmp_path):
    """U13: a keepblock that is neither true nor false presses neither side."""
    config = make_config([barcode_list("bc1", keep="yes"), constant_region("con1")])
    path = tmp_path / "odd_keepblock.json"
    path.write_text(json.dumps(config, indent=4))
    app.set_input_files("#open_input", str(path))
    wait_ready(app)
    assert pressed_state(app.locator("#region_switch")) == {
        "Keep": "false",
        "Leave out": "false",
    }
    assert ("E11", "error") in [(rule, level) for rule, level, _ in checks(app)]
    set_region_in_output_reads(app, "Keep")
    assert "E11" not in rule_ids(app)
    assert json_preview(app)["barcode_struct_r1"]["blocks"][0]["keepblock"] is True


# ---------------------------------------------------------------------------
# U14: "Delete region" looks like a button


def test_u14_delete_region_is_an_outlined_red_button(app):
    """U14: "Delete region" is an outlined red button with a light red hover.

    It asks before deleting (V13); answering "Delete region" removes the region.
    """
    start(app, SINGLE)
    add_region(app, "Barcode list")
    delete = selected_region_card(app).get_by_role(
        "button", name="Delete region", exact=True
    )
    style = computed_style(delete, "border-top-width", "border-top-style")
    assert float(style["border-top-width"].removesuffix("px")) >= 1
    assert style["border-top-style"] not in ("none", "hidden")
    red, green, blue, alpha = computed_rgba(delete, "border-top-color")
    assert alpha > 0
    assert red > green and red > blue

    resting = computed_rgba(delete, "background-color")
    delete.hover()
    hovered = computed_rgba(delete, "background-color")
    assert hovered != resting
    red, green, blue, _ = hovered
    assert red > green and red > blue

    delete.click()
    app.get_by_role("dialog").get_by_role(
        "button", name="Delete region", exact=True
    ).click()
    assert diagram_regions(app) == []


# ---------------------------------------------------------------------------
# U15: "region", never "block", and name rules that say which names


def test_u15_the_page_says_region_and_the_name_help_is_short(app):
    """U15: with a region and its name help open, nothing on the page says "block"."""
    open_template(app, "SDR-seq RNA")
    assert selected_region_card(app).is_visible()
    app.get_by_role("button", name="Rules for region names", exact=True).click()
    help_text = app.locator("#help_name")
    assert help_text.is_visible()
    assert help_text.inner_text() == NAME_HELP
    assert re.search(r"\bQC\b", help_text.inner_text(), flags=re.IGNORECASE) is None
    assert (
        selected_region_card(app)
        .get_by_role("button", name="Delete region", exact=True)
        .is_visible()
    )
    assert app.locator("#structure_caption").inner_text().endswith(CONTIGUITY_NOTE)
    assert app.locator("#json_preview").count() == 0
    assert "block" not in visible_text(app).lower()
    labels = app.evaluate(
        """() => Array.from(document.querySelectorAll('[aria-label], [title], [placeholder]'))
            .flatMap((element) => ['aria-label', 'title', 'placeholder']
                .map((name) => element.getAttribute(name)).filter((value) => value))"""
    )
    assert [label for label in labels if "block" in label.lower()] == []


def test_u15_outputs_without_a_umi_region_say_region(app):
    """U15: without a UMI region, the Outputs panel's UB and UR notes don't say "block"."""
    open_template(app, "SDR-seq gDNA")
    for selector in ("#out_ub", "#out_ur"):
        text = app.locator(selector).inner_text()
        assert text.startswith("not written")
        assert "block" not in text.lower()
    assert "block" not in visible_text(app).lower()


def test_u15_duplicate_names_name_the_regions(app):
    """U15: the duplicate-name message says it is about region names in a read."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    add_region(app, "Barcode list")
    fill(app, "Name", "BC1")
    assert [
        (level, message) for rule, level, message in checks(app) if rule == "E08"
    ] == [
        (
            "error",
            "bc1 and BC1: two regions in Read 1 can't share a name, even in different case.",
        )
    ]


# ---------------------------------------------------------------------------
# U16: CB says it is corrected


def test_u16_cb_says_it_is_corrected(app):
    """U16: the Outputs panel says CB is corrected, as it says for UB.

    CR holds the raw bases of the CB regions (V09).
    """
    open_template(app, SIMPLE_EXAMPLE)
    assert app.locator("#out_cb").inner_text() == "cb1 corrected"
    assert app.locator("#out_cr").inner_text() == "raw bases of the CB regions"
    assert app.locator("#out_ub").inner_text() == "UR, corrected"


def test_u16_cb_names_random_regions_that_go_in_uncorrected(app):
    """U16: a random region set to barcode goes into CB as read, and CB says so."""
    start(app, SINGLE)
    add_list_and_constant(app)
    add_region(app, "Random (UMI)")
    app.get_by_role("combobox", name="Function").select_option("barcode")
    assert (
        app.locator("#out_cb").inner_text()
        == "bc1.umi1 barcode lists corrected; umi1 as read"
    )
    assert [
        (level, message) for rule, level, message in checks(app) if rule == "W06"
    ] == [
        (
            "warning",
            "umi1 is a random region set to barcode, so its bases go into CB without error correction.",
        )
    ]


# ---------------------------------------------------------------------------
# U17: ask before replacing the page


@pytest.mark.parametrize("action", ["template", "file", "new empty config"])
def test_u17_the_welcome_panel_is_replaced_without_a_question(app, action):
    """U17: the untouched welcome panel holds nothing to lose, so nothing is asked."""
    if action == "file":
        app.set_input_files("#open_input", str(EXAMPLES_DIR / "se_mini_config.json"))
    elif action == "template":
        choose_template(app, "SDR-seq RNA")
    else:
        choose_template(app, "New empty config")
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    if action == "new empty config":
        assert on_welcome(app)
    else:
        region_tile(app, "cb1").wait_for()


@pytest.mark.parametrize("how", ["Cancel", "Escape"])
def test_u17_a_template_over_another_asks_first(app, how):
    """U17: opening a template over another asks first.

    Cancel is the filled button (V11) and has the focus; Cancel or Escape keeps
    the page, Replace loads the template.
    """
    open_template(app, "SDR-seq RNA")
    choose_template(app, "SDR-seq gDNA")
    dialog, text = replace_question(app)
    assert (
        text == 'Opening the template "SDR-seq gDNA" replaces the config on this page.'
    )
    assert app.evaluate("document.activeElement.textContent") == "Cancel"
    if how == "Escape":
        app.keyboard.press("Escape")
    else:
        dialog.get_by_role("button", name="Cancel", exact=True).click()
    assert app.get_by_role("dialog").count() == 0
    assert region_tile(app, "cso").is_visible()
    assert region_tile(app, "commonseq1").count() == 0

    choose_template(app, "SDR-seq gDNA")
    replace_question(app)[0].get_by_role("button", name="Replace", exact=True).click()
    region_tile(app, "commonseq1").wait_for()
    assert region_tile(app, "cso").count() == 0


def test_u17_the_question_says_when_changes_are_not_downloaded(app):
    """U17: the question mentions changes not yet downloaded, only while there are some."""
    open_template(app, "SDR-seq RNA")
    fill(app, "Name", "cb_renamed")
    choose_template(app, "SDR-seq gDNA")
    dialog, text = replace_question(app)
    assert text == (
        'Opening the template "SDR-seq gDNA" replaces the config on this page. '
        "It has changes you have not downloaded."
    )
    dialog.get_by_role("button", name="Cancel", exact=True).click()

    download_text(app)
    choose_template(app, "SDR-seq gDNA")
    dialog, text = replace_question(app)
    assert (
        text == 'Opening the template "SDR-seq gDNA" replaces the config on this page.'
    )
    dialog.get_by_role("button", name="Cancel", exact=True).click()
    assert region_tile(app, "cb_renamed").is_visible()


def test_u17_open_asks_first_naming_the_file_then_the_older_format_question(app):
    """U17: Open… asks first, naming the file; the older-format question comes after."""
    open_template(app, "SDR-seq RNA")
    path = FIXTURES_DIR / "legacy_master_gDNA.json.gz"
    for answer in ("Cancel", "Replace"):
        with app.expect_file_chooser() as chooser_info:
            app.click("#open_button")
        chooser_info.value.set_files(str(path))
        dialog, text = replace_question(app)
        assert (
            text
            == "Opening legacy_master_gDNA.json.gz replaces the config on this page."
        )
        dialog.get_by_role("button", name=answer, exact=True).click()
        if answer == "Cancel":
            assert app.get_by_role("dialog").count() == 0
            assert region_tile(app, "cso").is_visible()

    older_format_question = app.get_by_role("dialog")
    assert (
        older_format_question.get_by_role("heading").inner_text()
        == "Single- or paired-end?"
    )
    older_format_question.get_by_role("button", name="Open", exact=True).click()
    wait_ready(app)
    region_tile(app, "commonseq1").wait_for()


def test_u17_control_o_asks_first_too(app):
    """U17: Ctrl/Cmd+O asks the same question before replacing the page."""
    open_template(app, SIMPLE_EXAMPLE)
    with app.expect_file_chooser() as chooser_info:
        app.keyboard.press("ControlOrMeta+o")
    chooser_info.value.set_files(str(EXAMPLES_DIR / "se_mini_config.json"))
    dialog, text = replace_question(app)
    assert text == "Opening se_mini_config.json replaces the config on this page."
    dialog.get_by_role("button", name="Replace", exact=True).click()
    wait_ready(app)
    # The simple example is paired-end; se_mini is single-end.
    assert app.get_by_role("tab").count() == 0


# ---------------------------------------------------------------------------
# U18: "As entered"


def test_u18_order_buttons_are_a_to_z_and_as_entered(app):
    """U18: the order choice reads "A→Z" or "As entered", which keeps the typed order."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    order = app.get_by_role("group", name="Order")
    assert pressed_state(order) == {"A→Z": "true", "As entered": "false"}
    assert list(pressed_state(order)) == ["A→Z", "As entered"]
    assert app.get_by_role("button", name="Entered", exact=True).count() == 0

    entered = ["TTTTCCCC", "AAAACCCC", "GGGGCCCC"]
    for sequence in entered:
        type_sequence(app, sequence)
    set_order(app, "As entered")
    assert pressed_state(order) == {"A→Z": "false", "As entered": "true"}
    assert table_sequences(app) == entered
    config = json.loads(download_text(app))
    assert config["barcode_struct_r1"]["blocks"][0]["sequence"] == entered
