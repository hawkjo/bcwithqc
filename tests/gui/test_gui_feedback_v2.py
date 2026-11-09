"""Tier B: changes from the second round of hands-on testing, and review fixes.

A second round of hands-on testing, of the page as changed after the first
round (test_gui_feedback_v1.py), produced notes numbered V01–V15. Code
reviews found R01–R25 in the first round's changes, and R26–R46 in the
second round's, together with a trial run of a manual checklist. The notes
and findings are not part of the repo: each test is named after the item it
covers (test_v04_…, test_r26_…), and its docstring states the behaviour it
expects. Expected texts and values are written out here, not read from the
page; read lengths are worked out by hand from the regions' lengths and
allowed errors. V01, V05 and V08 asked for no change to the page (a remark,
a question, and an item for bcwithqc's own to-do list), so no test covers
them. V04's wording, "In output reads: Keep | Leave out", is asserted
wherever a test reads the switches. V15 is fixed with V11, and R02 with V12.

Drags from outside the page are plain events carrying a dataTransfer object,
as in test_gui_feedback_v1.py. Keyboard use is simulated by focusing an
element and pressing keys, since WebKit neither focuses a clicked button nor
moves to one with Tab. Texts that name the tool take its name from
pyproject.toml.
"""

import base64
import json
from pathlib import Path

import pytest
from gui_support import (
    ADD_BUTTONS,
    EXAMPLES_DIR,
    FIXTURES_DIR,
    GUI_PATH,
    OLD_GUI_BARCODES,
    READ_KEYS,
    SAMPLE_BARCODES,
    SAMPLE_CONSTANT,
    SIMPLE_EXAMPLE,
    add_region,
    barcode_list,
    canonical_text,
    constant_region,
    download_text,
    fill,
    generated_with,
    make_config,
    old_gui_config,
    on_welcome,
    open_template,
    paste_sequences,
    random_block,
    read_config,
    rule_ids,
    select_region,
    set_region_in_output_reads,
    set_rest_in_output_reads,
    show_read,
    table_sequences,
    type_sequence,
    wait_ready,
)

PAIRED = "Paired-end"
SINGLE = "Single-end"
ALWAYS_KEPT = "A read without regions is always kept."
CONTIGUITY_NOTE = (
    "Regions start at base 1 and follow each other without gaps. "
    "Anything after the last region is the rest of the read."
)
BARCODES_TXT = ("\n".join(SAMPLE_BARCODES) + "\n").encode("ascii")
# Eight distinct 9-nt barcodes.
EIGHT_BARCODES = [*SAMPLE_BARCODES, "CATAGGTGT", "CACGACTCA", "AGCATTCCA", "TCCACCGGT"]
# SDR-seq's linker variants, 17, 16, 15 and 14 nt long.
CSO_VARIANTS = [
    "GTCAGTACGTACGAGTC",
    "TCAGTACGTACGAGTC",
    "CAGTACGTACGAGTC",
    "AGTACGTACGAGTC",
]
# Two more 21-nt variants of SAMPLE_CONSTANT.
OTHER_VARIANTS = ["GGGTCAGTACGTACGAGTCCA", "GGGTCAGTACGTACGAGTCCG"]
QUESTION_FOCUS_RING = {"outline-style": "solid", "outline-width": "2px"}
# The caption joins each read length to its "≥" and "nt" with a no-break space.
NBSP = chr(0xA0)
SINGLE_END_KEEPS_NOTHING = (
    "This config keeps nothing from Read 1 in the output reads. With single-end "
    'reads they would be empty: set the rest of Read 1, or a region, to "Keep" '
    "if you need them."
)
PAIRED_END_KEEPS_NOTHING = (
    "This config keeps nothing from Read 1, the barcode read, in the output reads. "
    "As with the old pipeline, use only the other read's file."
)


# ---------------------------------------------------------------------------
# Helpers


def tool_name() -> str:
    """Return the tool's name as pyproject.toml gives it, for texts that name it."""
    return generated_with()["name"]


def start(app, layout: str) -> None:
    """Choose the read layout ("Paired-end" or "Single-end") on the welcome panel."""
    app.locator("#welcome").get_by_role("button", name=layout, exact=True).click()


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


def is_focused(locator) -> bool:
    """True when the element a locator finds has the keyboard focus."""
    return locator.evaluate("(element) => element === document.activeElement")


def json_preview(app) -> dict:
    """Return the config the JSON panel shows, opening the panel and closing it."""
    app.get_by_role("button", name="▸ show", exact=True).click()
    config = json.loads(app.locator("#json_preview").inner_text())
    app.get_by_role("button", name="▾ hide", exact=True).click()
    return config


def output_files(app) -> list[str]:
    """Return the Outputs panel's lines about the output reads."""
    return app.locator("#out_files > span:not(.hint)").all_inner_texts()


def rest_row_parts(app) -> list[str]:
    """Return the texts of the rest-of-read row's parts, left to right."""
    return app.locator(".align-row", has=app.locator("#rest_switch")).evaluate(
        "(row) => Array.from(row.children, (child) => child.textContent)"
    )


def add_list_and_constant(app) -> None:
    """Add bc1 (SAMPLE_BARCODES) and con1 (SAMPLE_CONSTANT) to the shown read."""
    add_region(app, "Barcode list")
    assert paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    add_region(app, "Constant region")
    type_sequence(app, SAMPLE_CONSTANT)


def choose_template(app, title: str) -> None:
    """Pick a template from the Templates menu with the mouse, leaving any question to the caller."""
    app.click("#templates_button")
    app.get_by_role("menuitem", name=title, exact=True).click()


def write_json(path: Path, config: dict) -> Path:
    """Write a config as indented JSON and return its path."""
    path.write_text(json.dumps(config, indent=4), encoding="utf-8")
    return path


def open_file(app, path: Path) -> None:
    """Choose a config file with the Open… button."""
    with app.expect_file_chooser() as chooser_info:
        app.click("#open_button")
    chooser_info.value.set_files(str(path))


def import_file(app, path: Path) -> None:
    """Choose a file with the selected region's "Import file…" button."""
    with app.expect_file_chooser() as chooser_info:
        app.click("#import_button")
    chooser_info.value.set_files(str(path))


def question(app) -> tuple:
    """Return the open question, its heading and the text of its first paragraph."""
    dialog = app.get_by_role("dialog")
    heading = dialog.get_by_role("heading").inner_text()
    return dialog, heading, dialog.locator("p").first.inner_text()


def wait_closed(app) -> None:
    """Wait until no question is open."""
    app.locator("#dialog").wait_for(state="hidden")


def delete_selected_region(app, answer: str = "Delete region") -> None:
    """Press the selected region's "Delete region" and answer the question that follows.

    Args:
        app: The page.
        answer: "Delete region" or "Cancel" to click that button, or "Escape"
            to press the key.
    """
    app.click("#delete_region_button")
    dialog = app.get_by_role("dialog")
    if answer == "Escape":
        dialog.wait_for()
        app.keyboard.press("Escape")
    else:
        dialog.get_by_role("button", name=answer, exact=True).click()
    wait_closed(app)


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

    Engines write computed colours differently, so the browser resolves the
    value by painting it onto a 1 x 1 canvas.
    """
    return locator.evaluate(COMPUTED_RGBA, css_property)


def assert_cancel_is_the_filled_default(app, others: list[str]) -> None:
    """Check the open question's buttons: Cancel, then the others.

    Cancel is the filled button ("btn primary"), has the focus and shows a
    solid 2-px focus ring, so Return cancels; no other button is filled.
    """
    dialog = app.get_by_role("dialog")
    assert dialog.get_by_role("button").all_inner_texts() == ["Cancel", *others]
    cancel = dialog.get_by_role("button", name="Cancel", exact=True)
    assert classes(cancel) == ["btn", "primary"]
    assert is_focused(cancel)
    assert (
        computed_style(cancel, "outline-style", "outline-width") == QUESTION_FOCUS_RING
    )
    filled = [
        name
        for name in others
        if "primary" in classes(dialog.get_by_role("button", name=name, exact=True))
    ]
    assert filled == []


NOTICE_ITEMS = """(items) => items.map((item) => ({
    class: item.className,
    display: getComputedStyle(item).display,
    width: item.getBoundingClientRect().width,
    overflows: item.scrollWidth > item.clientWidth + 1,
    text: item.textContent,
}))"""


def notice_items(app) -> list[dict]:
    """Return each notice's class, computed display, width, overflow and text."""
    return app.locator("#notices li").evaluate_all(NOTICE_ITEMS)


def unreadable_notices(app, items: list[dict]) -> list[str]:
    """Return what keeps a notice from reading as a full-width list item, one line each.

    A notice is a list item of class notice-info or notice-warning, at least
    90% as wide as the notices box, whose text doesn't overflow it.
    """
    box_width = app.locator("#notices").bounding_box()["width"]
    problems = []
    for item in items:
        if item["class"] not in ("notice-info", "notice-warning"):
            problems.append(f"class {item['class']!r}: {item['text']}")
        if item["display"] != "list-item":
            problems.append(f"display {item['display']}: {item['text']}")
        if item["width"] < 0.9 * box_width:
            problems.append(
                f"{item['width']:.0f} px wide in a {box_width:.0f}-px box: {item['text']}"
            )
        if item["overflows"]:
            problems.append(f"text overflows: {item['text']}")
    return problems


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


def focus_from_keyboard(app, locator) -> None:
    """Focus an element the way keyboard navigation does.

    One Tab press puts the browser in keyboard mode, so :focus-visible
    applies to the focus that follows; the element is then focused directly,
    since WebKit's Tab skips buttons.
    """
    app.keyboard.press("Tab")
    locator.focus()


def next_frames(app) -> None:
    """Wait two animation frames, so resize handlers and layout have run."""
    app.evaluate(
        "() => new Promise((resolve) => "
        "requestAnimationFrame(() => requestAnimationFrame(resolve)))"
    )


def horizontal_extent(app, selector: str) -> tuple[float, float, float]:
    """Return an element's left and right edges and the window's width, in CSS px."""
    box = app.locator(selector).bounding_box()
    assert box is not None, f"{selector} is not shown"
    width = app.evaluate("document.documentElement.clientWidth")
    return box["x"], box["x"] + box["width"], width


def relative_luminance(rgba: list[int]) -> float:
    """Return the relative luminance of an 8-bit sRGB colour (WCAG 2)."""

    def linear(channel: int) -> float:
        value = channel / 255
        return value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4

    red, green, blue = (linear(channel) for channel in rgba[:3])
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


def contrast_ratio(first: list[int], second: list[int]) -> float:
    """Return the WCAG contrast ratio of two opaque colours, from 1 to 21."""
    lighter, darker = sorted(
        (relative_luminance(first), relative_luminance(second)), reverse=True
    )
    return (lighter + 0.05) / (darker + 0.05)


# ---------------------------------------------------------------------------
# V02, V03: "Rest of Read N" from the start; a read without regions is always kept


@pytest.mark.parametrize("number", [1, 2])
def test_v02_v03_a_fresh_read_says_rest_of_read_and_is_always_kept(app, number):
    """V02, V03: on a fresh paired-end page, either read's row says "Rest of Read N".

    The read has no regions, so it is always kept: "Keep" is pressed and
    "Leave out" is disabled, with the reason as its title and beside the
    switch. The diagram says the same for either read.
    """
    start(app, PAIRED)
    show_read(app, number)
    assert rest_row_parts(app) == [
        f"Rest of Read {number}",
        "In output reads:",
        "KeepLeave out",
        ALWAYS_KEPT,
    ]
    assert app.get_by_role(
        "group", name=f"Rest of Read {number} in output reads", exact=True
    ).is_visible()
    switch = app.locator("#rest_switch")
    assert pressed_state(switch) == {"Keep": "true", "Leave out": "false"}
    leave_out = switch.get_by_role("button", name="Leave out", exact=True)
    assert leave_out.is_disabled()
    assert leave_out.get_attribute("title") == ALWAYS_KEPT
    assert app.locator("#rest_always_kept").inner_text() == ALWAYS_KEPT
    empty = app.locator("#empty_read")
    assert classes(empty) == ["blk", "rest", "kept"]
    assert empty.locator("b").inner_text() == "No regions yet"
    assert (
        empty.locator("span").inner_text()
        == f"All of Read {number} is kept in the output reads."
    )


def test_v03_the_download_outputs_a_read_without_regions_whole(app):
    """V03: the file says keep_nonbarcode true for a read without regions.

    Leaving out the rest of the other read changes only that read.
    """
    start(app, PAIRED)
    add_list_and_constant(app)
    set_rest_in_output_reads(app, "Leave out")
    config = json.loads(download_text(app))
    assert config["barcode_struct_r1"]["keep_nonbarcode"] is False
    assert config["barcode_struct_r2"] == {"keep_nonbarcode": True, "blocks": []}


@pytest.mark.parametrize("number", [1, 2])
def test_v03_deleting_the_last_region_outputs_the_whole_read_again(app, number):
    """V03: a region enables "Leave out"; deleting the read's last one returns to Keep.

    Deleting a region that isn't the last keeps the choice.
    """
    start(app, PAIRED)
    show_read(app, number)
    switch = app.locator("#rest_switch")
    leave_out = switch.get_by_role("button", name="Leave out", exact=True)
    add_region(app, "Barcode list")
    assert leave_out.is_enabled()
    assert leave_out.get_attribute("title") is None
    assert app.locator("#rest_always_kept").count() == 0
    add_region(app, "Random (UMI)")
    set_rest_in_output_reads(app, "Leave out")
    assert pressed_state(switch) == {"Keep": "false", "Leave out": "true"}

    delete_selected_region(app)
    assert diagram_regions(app) == ["bc1"]
    assert pressed_state(switch) == {"Keep": "false", "Leave out": "true"}
    assert leave_out.is_enabled()

    select_region(app, "bc1")
    delete_selected_region(app)
    assert diagram_regions(app) == []
    assert pressed_state(switch) == {"Keep": "true", "Leave out": "false"}
    assert leave_out.is_disabled()
    assert app.locator("#rest_always_kept").inner_text() == ALWAYS_KEPT
    assert json_preview(app)[READ_KEYS[number - 1]]["keep_nonbarcode"] is True


def test_v03_a_file_that_leaves_out_a_read_without_regions_opens_as_output(
    app, tmp_path
):
    """V03: a region-less Read 2 that says false opens as Keep, with a notice."""
    config = make_config(
        [barcode_list("bc1"), constant_region("con1")], [], keep_r1=True, keep_r2=False
    )
    app.set_input_files(
        "#open_input", str(write_json(tmp_path / "r2_off.json", config))
    )
    wait_ready(app)
    assert [item["text"] for item in notice_items(app)] == [
        (
            "Read 2 has no regions, so all of it is kept in the output reads. The file "
            'had "keep_nonbarcode": false; the page writes true, which '
            f"{tool_name()} does for a read without regions anyway."
        )
    ]
    assert output_files(app) == ["Read 1: rest of read", "Read 2: the whole read"]
    show_read(app, 2)
    switch = app.locator("#rest_switch")
    assert pressed_state(switch) == {"Keep": "true", "Leave out": "false"}
    assert switch.get_by_role("button", name="Leave out", exact=True).is_disabled()
    config = json.loads(download_text(app))
    assert config["barcode_struct_r2"] == {"keep_nonbarcode": True, "blocks": []}


@pytest.mark.parametrize(
    ("value", "file_had"),
    [
        (False, '"keep_nonbarcode": false'),
        ("no", '"keep_nonbarcode": "no"'),
        (None, '"keep_nonbarcode": null'),
    ],
    ids=["false", "text", "null"],
)
def test_v03_opening_sets_a_read_without_regions_to_output(core, value, file_had):
    """V03: config_from_json sets keep_nonbarcode true on a read without regions and says so."""
    config = make_config(
        [barcode_list("bc1"), constant_region("con1")], [], keep_r2=value
    )
    opened = core.call("config_from_json", config)
    assert [read["keep_nonbarcode"] for read in opened["model"]["reads"]] == [
        False,
        True,
    ]
    assert opened["notices"] == [
        {
            "level": "info",
            "text": (
                "Read 2 has no regions, so all of it is kept in the output reads. "
                f"The file had {file_had}; the page writes true, which "
                f"{tool_name()} does for a read without regions anyway."
            ),
        }
    ]


def test_v03_a_read_without_regions_that_says_true_opens_without_a_notice(core):
    """V03: nothing to say when the file already outputs the read whole."""
    config = make_config([barcode_list("bc1"), constant_region("con1")], [])
    opened = core.call("config_from_json", config)
    assert opened["model"]["reads"][1]["keep_nonbarcode"] is True
    assert opened["notices"] == []


def test_v03_core_outputs_every_read_without_regions_whole(core):
    """V03: the core helpers treat a read without regions as output whole.

    output_whole_reads_without_blocks changes only reads without blocks and
    returns their indices; remove_block resets a read left without blocks
    and keeps the choice of a read that still has one; describe_outputs
    calls a read without blocks "whole" whatever its flag says.
    """
    result = core.run(
        """(core) => {
            const model = core.new_model();
            model.single_end = false;
            const first = core.add_block(model, 0, 'barcodeList');
            const second = core.add_block(first.model, 0, 'randomBarcode');
            const two_regions = second.model;
            two_regions.reads.forEach((read) => { read.keep_nonbarcode = false; });
            const files = core.describe_outputs(two_regions).output_reads.map((file) => file.kind);
            const changed = core.output_whole_reads_without_blocks(two_regions);
            const keep = two_regions.reads.map((read) => read.keep_nonbarcode);
            const changed_again = core.output_whole_reads_without_blocks(two_regions);
            const one_region = core.remove_block(two_regions, 0, second.block_id);
            const no_region = core.remove_block(one_region, 0, first.block_id);
            return {
                files, changed, keep, changed_again,
                keep_with_one_region: one_region.reads[0].keep_nonbarcode,
                keep_with_no_region: no_region.reads[0].keep_nonbarcode,
                temporary_notes: Object.keys(core.TEMPORARY_OUTPUT_NOTES),
            };
        }"""
    )
    assert result == {
        "files": ["none", "whole"],
        "changed": [1],
        "keep": [False, True],
        "changed_again": [],
        "keep_with_one_region": False,
        "keep_with_no_region": True,
        "temporary_notes": ["ur_without_umi", "empty_alignment_file"],
    }


# ---------------------------------------------------------------------------
# V06: the read lengths a structure needs


def test_v06_the_sdr_seq_rna_caption_gives_minimum_and_recommended_read_lengths(app):
    """V06: SDR-seq RNA's caption gives its length, linkers and read lengths.

    Worked out by hand: cb1 and cb2 9 nt, cso 14–17 nt, commonseq2 23 nt,
    sb 8 nt, umi 8 nt; every region allows 1 error but umi, which allows 0.
    The linker combinations are cso's 4 variants times commonseq2's 1.
    """
    shortest = 9 + 14 + 9 + 23 + 8 + 8
    longest = 9 + 17 + 9 + 23 + 8 + 8
    allowed_errors = 1 + 1 + 1 + 1 + 1 + 0
    open_template(app, "SDR-seq RNA")
    facts = (
        f"{shortest}–{longest} nt · 4 linker combinations · "
        f"reads ≥{NBSP}{shortest - allowed_errors}{NBSP}nt "
        f"(≥{NBSP}{longest + allowed_errors}{NBSP}nt recommended)"
    )
    assert facts.replace(NBSP, " ") == (
        "71–74 nt · 4 linker combinations · reads ≥ 66 nt (≥ 79 nt recommended)"
    )
    assert (
        app.locator("#structure_caption").inner_text() == f"{facts} · {CONTIGUITY_NOTE}"
    )
    assert app.locator("#read_lengths").get_attribute("title") == (
        f"At least the regions' shortest total length ({shortest} nt) minus all "
        f"their allowed errors ({allowed_errors}). Recommended: their longest "
        f"total length ({longest} nt) plus all their allowed errors."
    )


def test_v06_lengths_show_once_every_region_is_filled_in(app):
    """V06: a list, a UMI, a constant region and a second UMI.

    The read lengths wait until the list and the constant region have
    sequences. Then, by hand: bc1 9 nt with 1 allowed error, umi1 8 nt with
    0, con1 21 nt with 1, umi2 8 nt with 0.
    """
    start(app, SINGLE)
    for kind in ("Barcode list", "Random (UMI)", "Constant region", "Random (UMI)"):
        add_region(app, kind)
    assert diagram_regions(app) == ["bc1", "umi1", "con1", "umi2"]
    caption = app.locator("#structure_caption")
    waiting = (
        "1 linker combination · lengths show once every region is filled in · "
        + CONTIGUITY_NOTE
    )
    assert caption.inner_text() == waiting
    assert app.locator("#read_lengths").get_attribute("title") is None

    select_region(app, "bc1")
    assert paste_sequences(app, "\n".join(EIGHT_BARCODES))
    assert app.locator("#sequence_count").inner_text() == "8 barcodes"
    assert caption.inner_text() == waiting

    select_region(app, "con1")
    type_sequence(app, SAMPLE_CONSTANT)
    total = 9 + 8 + 21 + 8
    allowed_errors = 1 + 0 + 1 + 0
    facts = (
        f"{total} nt · 1 linker combination · reads ≥{NBSP}{total - allowed_errors}{NBSP}nt "
        f"(≥{NBSP}{total + allowed_errors}{NBSP}nt recommended)"
    )
    assert facts.replace(NBSP, " ") == (
        "46 nt · 1 linker combination · reads ≥ 44 nt (≥ 48 nt recommended)"
    )
    assert caption.inner_text() == f"{facts} · {CONTIGUITY_NOTE}"


def structure(core, blocks: list[dict]):
    """Return describe_structure for blocks built with gui_support's builders."""
    return core.call("describe_structure", blocks)


def unknown_lengths(linker_combinations: int) -> dict:
    """Return describe_structure's result while some length or error count isn't known."""
    return {
        "lengths_known": False,
        "min_length": None,
        "max_length": None,
        "allowed_errors": None,
        "linker_combinations": linker_combinations,
        "min_read_length": None,
        "recommended_read_length": None,
    }


def test_v06_describe_structure_counts_every_allowed_error(core):
    """V06: shortest total minus all allowed errors, longest total plus them.

    A list's length runs from its shortest to its longest barcode, not its
    first: here the first written is the longest. Empty entries don't count.
    """
    assert structure(core, []) is None
    blocks = [
        barcode_list("bc1", ["CCTTAACAT", "", "ACCTGGA"], maxerrors=1),
        constant_region("con1", ["GGGTCAGTACGTACGAGTCCC", "GGGTCAG"], maxerrors=2),
        random_block("umi1", length=8, maxerrors=0),
    ]
    # bc1 7–9 nt, con1 7–21 nt, umi1 8 nt; 1 + 2 + 0 allowed errors.
    assert structure(core, blocks) == {
        "lengths_known": True,
        "min_length": 7 + 7 + 8,
        "max_length": 9 + 21 + 8,
        "allowed_errors": 3,
        "linker_combinations": 2,
        "min_read_length": 7 + 7 + 8 - 3,
        "recommended_read_length": 9 + 21 + 8 + 3,
    }


def test_v06_the_minimum_read_length_is_at_least_1_nt(core):
    """V06: when the allowed errors reach the structure's length, a read still needs 1 nt."""
    blocks = [constant_region("con1", ["ACGT"], maxerrors=4)]
    result = structure(core, blocks)
    assert result["min_read_length"] == 1
    assert result["recommended_read_length"] == 4 + 4


@pytest.mark.parametrize(
    ("blocks", "linker_combinations"),
    [
        ([barcode_list("bc1", [])], 1),
        ([random_block("umi1", length=0)], 1),
        ([random_block("umi1", length="8")], 1),
        ([barcode_list("bc1", maxerrors=-1)], 1),
        ([barcode_list("bc1", maxerrors=1.5)], 1),
        ([barcode_list("bc1", maxerrors="1")], 1),
        ([constant_region("con1", CSO_VARIANTS), constant_region("con2", [])], 4),
    ],
    ids=[
        "empty list",
        "length 0",
        "length as text",
        "negative errors",
        "fractional errors",
        "errors as text",
        "empty constant region",
    ],
)
def test_v06_lengths_wait_for_every_region(core, blocks, linker_combinations):
    """V06: while any region's length or allowed errors is unknown, no lengths are given.

    The linker combinations are still counted; an empty constant region counts once.
    """
    assert structure(core, blocks) == unknown_lengths(linker_combinations)


# ---------------------------------------------------------------------------
# V07: "only A, C, G and T"


def test_v07_a_column_with_an_n_imports_nothing_and_says_only_acgt(app, tmp_path):
    """V07: a column with an N is no column of DNA, and the message says "only"."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    path = tmp_path / "with_n.txt"
    path.write_text("CCTTAACAT\nAGCGTANAA\nTATAGGTGT\n", encoding="ascii")
    import_file(app, path)
    dialog, heading, text = question(app)
    assert heading == "Nothing imported"
    assert text == "No column of DNA sequences (only A, C, G and T) was found."
    dialog.get_by_role("button", name="OK", exact=True).click()
    wait_ready(app)
    assert app.locator("#sequence_count").inner_text() == "0 barcodes"


# ---------------------------------------------------------------------------
# V09, R01, R06: the Outputs panel's notes


def test_v09_cr_and_ur_say_raw_bases_of_the_cb_regions_and_raw_umi(app):
    """V09: CR is "raw bases of the CB regions" and UR "<UMI regions> raw UMI"."""
    open_template(app, SIMPLE_EXAMPLE)
    assert app.locator("#out_cb").inner_text() == "cb1 corrected"
    assert app.locator("#out_cr").inner_text() == "raw bases of the CB regions"
    assert app.locator("#out_ub").inner_text() == "UR, corrected"
    assert app.locator("#out_ur").inner_text() == "umi raw UMI"
    assert app.locator("#out_ur .note").inner_text() == "raw UMI"


def test_r01_without_a_umi_region_ub_and_ur_say_no_umi_region(app):
    """R01: without a UMI region, UB and UR say "no UMI region", not "block"."""
    open_template(app, "SDR-seq gDNA")
    assert app.locator("#out_ub").inner_text() == "not written (no UMI region)"
    assert app.locator("#out_ur").inner_text() == (
        f"not written (no UMI region); current {tool_name()} writes it empty"
    )


def test_r06_cb_says_which_regions_are_corrected(app):
    """R06: CB's note follows what goes into CB.

    Only random regions: "as read, not corrected"; lists and random regions:
    "barcode lists corrected; <random regions> as read"; only lists: "corrected".
    """
    start(app, SINGLE)
    function = app.get_by_role("combobox", name="Function")
    cb = app.locator("#out_cb")
    add_region(app, "Random (UMI)")
    function.select_option("barcode")
    assert cb.inner_text() == "umi1 as read, not corrected"

    add_region(app, "Barcode list")
    assert cb.inner_text() == "umi1.bc1 barcode lists corrected; umi1 as read"

    add_region(app, "Random (UMI)")
    function.select_option("barcode")
    assert (
        cb.inner_text() == "umi1.bc1.umi2 barcode lists corrected; umi1, umi2 as read"
    )

    for name in ("umi1", "umi2"):
        select_region(app, name)
        function.select_option("UMI")
    assert cb.inner_text() == "bc1 corrected"
    assert app.locator("#out_ur").inner_text() == "umi1.umi2 raw UMI"


def test_r06_describe_outputs_names_corrected_and_uncorrected_cb_regions(core):
    """R06: CB.corrected names the lists in CB, CB.uncorrected its random regions.

    Both follow the file order across the reads. A list set to UMI still goes
    into CB as a list; a region set to discard goes nowhere.
    """
    config = make_config(
        [
            barcode_list("bc1"),
            random_block("rnd1", function="barcode"),
            constant_region("con1"),
            barcode_list("bc2", function="discard"),
        ],
        [random_block("rnd2", function="barcode"), barcode_list("bc3", function="UMI")],
    )
    cb = core.run(
        "(core, config) => core.describe_outputs(core.config_from_json(config).model).CB",
        config,
    )
    assert cb == {
        "pieces": [["bc1", "rnd1"], ["rnd2", "bc3"]],
        "value": "bc1.rnd1.rnd2.bc3",
        "corrected": ["bc1", "bc3"],
        "uncorrected": ["rnd1", "rnd2"],
    }


def test_r01_an_unnamed_region_of_unknown_type_is_named_region1(core):
    """R01: a region of unknown type without a name gets "region1", unique ignoring case."""
    mystery = {
        "maxerrors": 1,
        "blockfunction": "barcode",
        "length": 8,
        "blocktype": "randomBarcod",
        "keepblock": False,
    }
    opened = core.call("config_from_json", make_config([barcode_list("bc1"), mystery]))
    assert opened["model"]["reads"][0]["blocks"][1]["blockname"] == "region1"
    assert {
        "level": "info",
        "text": "A region without a name was named region1.",
    } in opened["notices"]

    taken = make_config([barcode_list("REGION1"), mystery])
    opened = core.call("config_from_json", taken)
    assert opened["model"]["reads"][0]["blocks"][1]["blockname"] == "region2"


# ---------------------------------------------------------------------------
# V10: a random region set to barcode


def test_v10_sdr_seq_rnas_umi_set_to_barcode_goes_into_cb_as_read(app):
    """V10: setting umi's Function to barcode names it "as read" in CB, with W06."""
    open_template(app, "SDR-seq RNA")
    select_region(app, "umi")
    app.get_by_role("combobox", name="Function").select_option("barcode")
    assert (
        app.locator("#out_cb").inner_text()
        == "cb1.cb2.sb.umi barcode lists corrected; umi as read"
    )
    assert (
        "W06",
        "warning",
        (
            "umi is a random region set to barcode, so its bases go into CB "
            "without error correction."
        ),
    ) in checks(app)


# ---------------------------------------------------------------------------
# V11, V15: Cancel is the default and the filled button; its ring always shows


def test_v11_the_replace_question_defaults_to_a_filled_cancel(app):
    """V11, V15: opened by mouse, the Replace question's Cancel is filled and focused.

    Its focus ring shows, and Return keeps the page as it was.
    """
    open_template(app, SIMPLE_EXAMPLE)
    choose_template(app, "SDR-seq RNA")
    _, heading, _ = question(app)
    assert heading == "Replace the current config?"
    assert_cancel_is_the_filled_default(app, ["Replace"])
    app.keyboard.press("Enter")
    wait_closed(app)
    assert diagram_regions(app) == ["cb1", "con1", "umi"]


def test_v11_the_delete_question_defaults_to_a_filled_cancel(app):
    """V11, V13: opened by mouse, the Delete question's Cancel is filled; Return keeps the region."""
    open_template(app, SIMPLE_EXAMPLE)
    app.click("#delete_region_button")
    _, heading, _ = question(app)
    assert heading == "Delete cb1?"
    assert_cancel_is_the_filled_default(app, ["Delete region"])
    app.keyboard.press("Enter")
    wait_closed(app)
    assert diagram_regions(app) == ["cb1", "con1", "umi"]


def test_v11_the_add_or_replace_question_defaults_to_a_filled_cancel(app, tmp_path):
    """V11: the add-or-replace question's Cancel is filled; Return leaves the list as it was."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    for sequence in SAMPLE_BARCODES[:2]:
        type_sequence(app, sequence)
    path = tmp_path / "more.txt"
    path.write_text("\n".join(SAMPLE_BARCODES[2:]) + "\n", encoding="ascii")
    import_file(app, path)
    _, heading, text = question(app)
    assert heading == "Add to the list or replace it?"
    assert text == "The list already has 2 barcodes; the import has 2 barcodes."
    assert_cancel_is_the_filled_default(app, ["Add to the existing 2", "Replace them"])
    app.keyboard.press("Enter")
    wait_closed(app)
    wait_ready(app)
    assert table_sequences(app) == sorted(SAMPLE_BARCODES[:2])


def test_v11_the_header_question_defaults_to_a_filled_cancel(app, tmp_path):
    """V11: the header question's Cancel is filled; Return imports nothing."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    path = tmp_path / "tagged.txt"
    path.write_text("TAG\n" + "\n".join(SAMPLE_BARCODES) + "\n", encoding="ascii")
    import_file(app, path)
    _, heading, _ = question(app)
    assert heading == 'Is "TAG" a header or a barcode?'
    assert_cancel_is_the_filled_default(app, ["A barcode", "A header"])
    app.keyboard.press("Enter")
    wait_closed(app)
    wait_ready(app)
    assert app.locator("#sequence_count").inner_text() == "0 barcodes"


def test_v11_the_unknown_keys_question_defaults_to_a_filled_cancel(app, tmp_path):
    """V11: the unknown-keys question's Cancel is filled; Return opens nothing."""
    config = make_config([barcode_list("bc1"), constant_region("con1")])
    config["comment"] = "made by hand"
    open_file(app, write_json(tmp_path / "commented.json", config))
    _, heading, _ = question(app)
    assert heading == f"Keep the key {tool_name()} doesn't use?"
    assert_cancel_is_the_filled_default(app, ["Leave it out", "Keep it"])
    app.keyboard.press("Enter")
    wait_closed(app)
    wait_ready(app)
    assert on_welcome(app)


# ---------------------------------------------------------------------------
# V12, R02: notices read as full-width list items


def test_v12_the_notice_on_a_config_from_the_page_reads_in_full(app, tmp_path):
    """V12: the one notice on a config this page wrote reads in full, version included."""
    path = tmp_path / "from_the_page.json"
    path.write_text(
        canonical_text(make_config([barcode_list("bc1"), constant_region("con1")])),
        encoding="ascii",
    )
    app.set_input_files("#open_input", str(path))
    wait_ready(app)
    made_for = generated_with()
    items = notice_items(app)
    assert [item["text"] for item in items] == [
        (
            f"This file was made for {made_for['name']} {made_for['version']}, "
            "the version this page writes for."
        )
    ]
    assert items[0]["class"] == "notice-info"
    assert unreadable_notices(app, items) == []
    assert (
        app.locator("#notices header span").inner_text()
        == "Opened from_the_page.json: 1 notice"
    )


def test_v12_several_notices_from_a_legacy_file_read_in_full(app):
    """V12: the three notices on an old master config are list items, not 18-px circles."""
    app.set_input_files("#open_input", str(FIXTURES_DIR / "legacy_master_gDNA.json.gz"))
    app.get_by_role("dialog").get_by_role("button", name="Open", exact=True).click()
    wait_ready(app)
    items = notice_items(app)
    assert [item["class"] for item in items] == [
        "notice-info",
        "notice-info",
        "notice-warning",
    ]
    assert unreadable_notices(app, items) == []
    assert app.locator("#notices li.info, #notices li.warning").count() == 0


# ---------------------------------------------------------------------------
# V13: deleting a region asks first


def test_v13_deleting_a_list_asks_first_and_counts_its_barcodes(app):
    """V13: "Delete cb1?" says the list's 1,536 barcodes go with it; Cancel keeps it."""
    example = read_config(EXAMPLES_DIR / "cDNA_paired_end_config.json")
    barcodes = len(example["barcode_struct_r1"]["blocks"][0]["sequence"])
    open_template(app, "SDR-seq RNA")
    app.click("#delete_region_button")
    dialog, heading, text = question(app)
    assert heading == "Delete cb1?"
    assert text == f"cb1 and its {barcodes:,} barcodes will be removed from Read 1."
    assert text == "cb1 and its 1,536 barcodes will be removed from Read 1."
    buttons = dialog.get_by_role("button")
    assert buttons.all_inner_texts() == ["Cancel", "Delete region"]
    assert classes(buttons.nth(0)) == ["btn", "primary"]
    assert classes(buttons.nth(1)) == ["btn", "danger"]
    buttons.nth(0).click()
    wait_closed(app)
    assert region_tile(app, "cb1").get_attribute("aria-pressed") == "true"


@pytest.mark.parametrize(
    ("name", "text"),
    [
        ("cso", "cso and its 4 variants will be removed from Read 1."),
        ("commonseq2", "commonseq2 and its 1 variant will be removed from Read 1."),
        ("umi", "umi will be removed from Read 1."),
    ],
    ids=["variants", "one variant", "random region"],
)
def test_v13_the_delete_question_says_what_goes_with_the_region(app, name, text):
    """V13: a constant region's variants are counted; a random region has nothing to count."""
    open_template(app, "SDR-seq RNA")
    select_region(app, name)
    app.click("#delete_region_button")
    _, heading, shown = question(app)
    assert (heading, shown) == (f"Delete {name}?", text)


def test_v13_an_empty_list_on_read_2_is_named_without_a_count(app):
    """V13: an empty list has nothing to count, and the question names Read 2."""
    start(app, PAIRED)
    show_read(app, 2)
    add_region(app, "Barcode list")
    app.click("#delete_region_button")
    _, heading, text = question(app)
    assert (heading, text) == ("Delete bc1?", "bc1 will be removed from Read 2.")


@pytest.mark.parametrize("answer", ["Cancel", "Escape"])
def test_v13_cancel_or_escape_keeps_the_region(app, answer):
    """V13: Cancel or Escape keeps the region, still selected."""
    open_template(app, SIMPLE_EXAMPLE)
    select_region(app, "con1")
    delete_selected_region(app, answer)
    assert diagram_regions(app) == ["cb1", "con1", "umi"]
    assert region_tile(app, "con1").get_attribute("aria-pressed") == "true"


def test_v13_delete_region_removes_it_and_focuses_add_barcode_list(app):
    """V13: "Delete region" removes the region, clears the selection and focuses "+ Barcode list"."""
    open_template(app, SIMPLE_EXAMPLE)
    select_region(app, "con1")
    delete_selected_region(app, "Delete region")
    assert diagram_regions(app) == ["cb1", "umi"]
    assert selected_region_card(app).count() == 0
    assert app.locator(".track button.blk[aria-pressed='true']").count() == 0
    assert is_focused(app.locator("#add_list_button"))
    blocks = json_preview(app)["barcode_struct_r1"]["blocks"]
    assert [block["blockname"] for block in blocks] == ["cb1", "umi"]


# ---------------------------------------------------------------------------
# V14: nothing runs past the card at middle widths

BEYOND_MAIN = """() => {
    const main = document.getElementById('main');
    const edge = main.getBoundingClientRect().right;
    const beyond = [];
    for (const element of main.querySelectorAll('*')) {
        const box = element.getBoundingClientRect();
        if (box.width === 0 && box.height === 0) continue;
        if (box.right > edge + 1) {
            const name = [element.tagName.toLowerCase(), element.id ? '#' + element.id : '',
                ...Array.from(element.classList, (item) => '.' + item)].join('');
            beyond.push(name + ' ends at ' + Math.round(box.right) + ' px, #main at ' + Math.round(edge) + ' px');
        }
    }
    return beyond;
}"""


@pytest.mark.parametrize("width", [1040, 1100, 1180, 1260, 1400])
def test_v14_nothing_in_the_card_runs_past_its_right_edge(app, width):
    """V14: with SDR-seq RNA's six regions and cb1 selected, nothing passes #main's right edge."""
    app.set_viewport_size({"width": width, "height": 900})
    open_template(app, "SDR-seq RNA")
    assert region_tile(app, "cb1").get_attribute("aria-pressed") == "true"
    assert app.evaluate(BEYOND_MAIN) == []
    assert app.evaluate(
        "document.documentElement.scrollWidth <= document.documentElement.clientWidth"
    )


def test_v14_at_1100_px_the_rest_tile_wraps_onto_a_second_line(app):
    """V14: at 1100 px the tiles don't fit one line, so the rest of the read wraps below cb1."""
    app.set_viewport_size({"width": 1100, "height": 900})
    open_template(app, "SDR-seq RNA")
    rest = app.locator(".track .blk.rest").bounding_box()
    cb1 = region_tile(app, "cb1").bounding_box()
    assert rest["y"] >= cb1["y"] + cb1["height"]


# ---------------------------------------------------------------------------
# R03: a legacy config that outputs nothing


@pytest.mark.parametrize(
    ("mode", "keep_rest", "expected"),
    [
        ("single", False, [SINGLE_END_KEEPS_NOTHING]),
        ("paired", False, [PAIRED_END_KEEPS_NOTHING]),
        ("single", True, []),
    ],
    ids=["single-end", "paired-end", "single-end, rest output"],
)
def test_r03_a_legacy_config_that_outputs_nothing_says_so_for_its_layout(
    core, mode, keep_rest, expected
):
    """R03: single-end has no other read to use, so its notice says what to set instead."""
    legacy = old_gui_config()
    del legacy["barcode_struct_r2"]
    legacy["barcode_struct_r1"]["keep_nonbarcode"] = keep_rest
    opened = core.call("config_from_json", legacy, {"legacy_mode": mode})
    assert [
        notice
        for notice in opened["notices"]
        if notice["text"].startswith("This config keeps nothing")
    ] == [{"level": "info", "text": text} for text in expected]


def test_r03_opening_a_legacy_file_as_single_end_shows_its_notice(app, tmp_path):
    """R03: answering "Single-end" to the older-format question gives the single-end notice."""
    legacy = old_gui_config()
    del legacy["barcode_struct_r2"]
    app.set_input_files("#open_input", str(write_json(tmp_path / "old.json", legacy)))
    dialog = app.get_by_role("dialog")
    dialog.get_by_label("Single-end", exact=True).check()
    dialog.get_by_role("button", name="Open", exact=True).click()
    wait_ready(app)
    texts = [item["text"] for item in notice_items(app)]
    assert SINGLE_END_KEEPS_NOTHING in texts
    assert PAIRED_END_KEEPS_NOTHING not in texts


# ---------------------------------------------------------------------------
# R04: which regions had no keepblock


@pytest.mark.parametrize(
    ("missing", "expected"),
    [
        (
            ["bc1", "con1", "umi1"],
            (
                '"In output reads" was not set in this file, so every region '
                'is "Leave out".'
            ),
        ),
        (
            ["con1"],
            '"In output reads" was not set for con1, so it is "Leave out".',
        ),
        (
            ["bc1", "umi1"],
            '"In output reads" was not set for bc1, umi1, so they are "Leave out".',
        ),
        ([], None),
    ],
    ids=["every region", "one region", "two regions", "none"],
)
def test_r04_the_keepblock_notice_names_the_regions_that_lacked_it(
    core, missing, expected
):
    """R04: "every region" only when every region lacked keepblock; else their names."""
    blocks = [barcode_list("bc1"), constant_region("con1"), random_block("umi1")]
    for block in blocks:
        if block["blockname"] in missing:
            del block["keepblock"]
    opened = core.call("config_from_json", make_config(blocks))
    assert [
        notice
        for notice in opened["notices"]
        if notice["text"].startswith('"In output reads"')
    ] == ([] if expected is None else [{"level": "info", "text": expected}])
    assert [block["keepblock"] for block in opened["model"]["reads"][0]["blocks"]] == [
        False,
        False,
        False,
    ]


# ---------------------------------------------------------------------------
# R05: a read with regions and no keep_nonbarcode


def test_r05_a_read_with_regions_and_no_keep_nonbarcode_is_left_out(core):
    """R05: the rest of a read with regions is "Leave out"; a read without regions is output."""
    config = make_config([barcode_list("bc1"), constant_region("con1")], [])
    del config["barcode_struct_r1"]["keep_nonbarcode"]
    del config["barcode_struct_r2"]["keep_nonbarcode"]
    opened = core.call("config_from_json", config)
    assert [read["keep_nonbarcode"] for read in opened["model"]["reads"]] == [
        False,
        True,
    ]
    assert opened["notices"] == [
        {
            "level": "warning",
            "text": 'Read 1 had no "keep_nonbarcode", so the rest of Read 1 is "Leave out".',
        },
        {
            "level": "info",
            "text": (
                "Read 2 has no regions, so all of it is kept in the output reads. "
                'The file had no "keep_nonbarcode"; the page writes true, which '
                f"{tool_name()} does for a read without regions anyway."
            ),
        },
    ]


# ---------------------------------------------------------------------------
# R07: an error banner goes with the next change


def test_r07_the_next_change_clears_a_could_not_open_banner(app, tmp_path):
    """R07: the next change, here choosing the read layout, clears "Could not open …"."""
    path = write_json(tmp_path / "not_a_config.json", {"hello": "world"})
    app.set_input_files("#open_input", str(path))
    banner = app.locator("#error_banner")
    banner.wait_for()
    assert banner.inner_text() == (
        f"Could not open not_a_config.json: This is not a {tool_name()} config: "
        'it has no "barcode_struct_r1" section.'
    )
    start(app, PAIRED)
    assert banner.is_hidden()


# ---------------------------------------------------------------------------
# R08: wording


def test_r08_the_decoder_badge_says_what_follows_the_list(app):
    """R08: the decoder badge says "a constant region follows" (or what else follows)."""
    open_template(app, "SDR-seq RNA")
    badge = app.locator("#decoder_badge")
    assert badge.inner_text() == "Levenshtein ≤ 1 · a constant region follows"
    select_region(app, "sb")
    assert badge.inner_text() == "FREE ≤ 1 · a random region follows"


def test_r08_the_replace_question_quotes_the_template_title(app):
    """R08: the Replace question quotes the template's title."""
    open_template(app, SIMPLE_EXAMPLE)
    choose_template(app, "SDR-seq RNA")
    _, _, text = question(app)
    assert text == (
        'Opening the template "SDR-seq RNA" replaces the config on this page.'
    )


def test_r08_the_header_question_for_a_constant_region_says_variant(app, tmp_path):
    """R08: for a constant region, the header question asks "header or a variant?"."""
    start(app, SINGLE)
    add_region(app, "Constant region")
    path = tmp_path / "tagged_variants.txt"
    path.write_text("TAG\n" + "\n".join(OTHER_VARIANTS) + "\n", encoding="ascii")
    import_file(app, path)
    dialog, heading, _ = question(app)
    assert heading == 'Is "TAG" a header or a variant?'
    assert dialog.get_by_role("button").all_inner_texts() == [
        "Cancel",
        "A variant",
        "A header",
    ]
    dialog.get_by_role("button", name="A header", exact=True).click()
    wait_ready(app)
    assert table_sequences(app) == OTHER_VARIANTS


def test_r08_the_add_or_replace_question_for_a_constant_region_says_variants(
    app, tmp_path
):
    """R08: for a constant region the question counts variants, and "1 variant" is singular."""
    start(app, SINGLE)
    add_region(app, "Constant region")
    type_sequence(app, SAMPLE_CONSTANT)
    path = tmp_path / "variants.txt"
    path.write_text("\n".join(OTHER_VARIANTS) + "\n", encoding="ascii")
    import_file(app, path)
    dialog, heading, text = question(app)
    assert heading == "Add to the variants or replace them?"
    assert text == "The region already has 1 variant; the import has 2 variants."
    dialog.get_by_role("button", name="Add to the existing 1", exact=True).click()
    wait_ready(app)
    assert table_sequences(app) == [SAMPLE_CONSTANT, *OTHER_VARIANTS]


def test_r08_the_lowercase_notice_names_regions_after_renaming(core):
    """R08: notices name regions by their names after spaces became "_"."""
    legacy = old_gui_config()
    legacy["barcode_struct_r1"]["blocks"][0]["sequence"] = [
        sequence.lower() for sequence in OLD_GUI_BARCODES
    ]
    texts = [
        notice["text"] for notice in core.call("config_from_json", legacy)["notices"]
    ]
    assert "Lowercase sequences in sequence_0 were converted to uppercase." in texts
    assert (
        'Spaces in region names became "_": sequence 0 → sequence_0, '
        "sequence 1 → sequence_1, sequence 2 → sequence_2, sequence 3 → sequence_3."
    ) in texts


# ---------------------------------------------------------------------------
# R09: "sequence" that is not a list


@pytest.mark.parametrize(
    ("index", "raw", "described", "typed"),
    [
        (0, None, "null", "CCTTAACAT"),
        (0, 7, "7", "CCTTAACAT"),
        (1, SAMPLE_CONSTANT, f'"{SAMPLE_CONSTANT}"', SAMPLE_CONSTANT),
    ],
    ids=["list null", "list number", "constant region string"],
)
def test_r09_sequences_that_are_not_a_list_open_and_can_be_replaced(
    app, tmp_path, index, raw, described, typed
):
    """R09: the page opens the file, shows the region as empty, reports E11, and Enter starts a list."""
    blocks = [barcode_list("bc1"), constant_region("con1")]
    blocks[index]["sequence"] = raw
    name = blocks[index]["blockname"]
    app.set_input_files(
        "#open_input", str(write_json(tmp_path / "odd.json", make_config(blocks)))
    )
    wait_ready(app)
    select_region(app, name)
    assert "sequences not a list" in region_tile(app, name).inner_text().splitlines()
    unit = "barcodes" if index == 0 else "variants"
    assert app.locator("#sequence_count").inner_text() == f"0 {unit}"
    assert app.locator("#list_dropzone").is_visible()
    assert (
        "E11",
        "error",
        f"{name}: the sequences must be a list; the file has {described}.",
    ) in checks(app)

    type_sequence(app, typed)
    assert table_sequences(app) == [typed]
    assert "E11" not in rule_ids(app)
    assert json_preview(app)["barcode_struct_r1"]["blocks"][index]["sequence"] == [
        typed
    ]


# ---------------------------------------------------------------------------
# R10: spaces split a pasted line only in a barcode list


@pytest.mark.parametrize(
    ("kind", "taken", "expected"),
    [
        ("Constant region", False, []),
        ("Barcode list", True, ["CTCGGAGATG", "GTCTCGTGGG"]),
    ],
    ids=["constant region", "barcode list"],
)
def test_r10_a_pasted_line_splits_on_spaces_only_in_a_barcode_list(
    app, kind, taken, expected
):
    """R10: a sequence copied in groups of ten pastes as typed into a constant region."""
    start(app, SINGLE)
    add_region(app, kind)
    assert paste_sequences(app, "gtctcgtggg ctcggagatg") is taken
    assert app.get_by_role("dialog").count() == 0
    assert table_sequences(app) == expected


@pytest.mark.parametrize(
    ("text", "split_on_spaces", "expected"),
    [
        ("GTCTCGTGGG CTCGGAGATG", True, ["GTCTCGTGGG", "CTCGGAGATG"]),
        ("GTCTCGTGGG CTCGGAGATG", False, None),
        ("GTCTCGTGGG\tCTCGGAGATG", False, ["GTCTCGTGGG", "CTCGGAGATG"]),
        ("GTCTCGTGGG, CTCGGAGATG", False, ["GTCTCGTGGG", "CTCGGAGATG"]),
        ("GTCTCGTGGG;CTCGGAGATG", False, ["GTCTCGTGGG", "CTCGGAGATG"]),
        ("GTCTCGTGGG", True, None),
        ("GTCTCGTGGG, hello", True, None),
    ],
    ids=[
        "spaces, list",
        "spaces, constant region",
        "tab",
        "comma",
        "semicolon",
        "one entry",
        "not all DNA",
    ],
)
def test_r10_pasted_list_rows_split_rule(core, text, split_on_spaces, expected):
    """R10: tabs, commas and semicolons always separate; spaces only when asked.

    One entry, or an entry that isn't DNA, is no list (null).
    """
    rows = core.call("pasted_list_rows", text, split_on_spaces)
    if expected is None:
        assert rows is None
    else:
        assert rows == {
            "rows": [
                {"line": number, "cells": [entry]}
                for number, entry in enumerate(expected, start=1)
            ],
            "unit": "item",
        }


def test_r10_pasted_list_rows_splits_on_spaces_by_default(core):
    """R10: without the second argument, spaces separate entries, as for a list."""
    rows = core.call("pasted_list_rows", "GTCTCGTGGG CTCGGAGATG")
    assert [row["cells"] for row in rows["rows"]] == [["GTCTCGTGGG"], ["CTCGGAGATG"]]


# ---------------------------------------------------------------------------
# R11: after the shown read changes, a new region goes at the end


def test_r11_after_a_checks_item_switches_the_read_a_region_goes_at_the_end(app):
    """R11: a Checks item that shows the other read clears the selection there."""
    start(app, PAIRED)
    add_region(app, "Barcode list")
    show_read(app, 2)
    add_region(app, "Constant region")
    app.locator("#checks .issue", has_text="Read 1 has a single region.").click()
    assert (
        app.get_by_role("tab", name="Read 1 · first file").get_attribute(
            "aria-selected"
        )
        == "true"
    )
    assert (
        app.locator("#add_region_where").inner_text()
        == "Add a region at the end of Read 1:"
    )
    add_region(app, "Random (UMI)")
    assert diagram_regions(app) == ["bc1", "umi1"]


def test_r11_after_switching_to_single_end_a_region_goes_at_the_end(app):
    """R11: switching to single-end while a Read 2 region is selected clears the selection."""
    start(app, PAIRED)
    add_region(app, "Barcode list")
    show_read(app, 2)
    add_region(app, "Constant region")
    app.get_by_role("group", name="Read layout").get_by_role(
        "button", name=SINGLE, exact=True
    ).click()
    assert app.get_by_role("tab").count() == 0
    assert (
        app.locator("#add_region_where").inner_text()
        == "Add a region at the end of Read 1:"
    )
    add_region(app, "Random (UMI)")
    assert diagram_regions(app) == ["bc1", "umi1"]


# ---------------------------------------------------------------------------
# R12: menus and help stay inside the window


def assert_inside_window(app, selector: str) -> None:
    """Check that an element keeps 8 px from both window edges (half a pixel of rounding allowed)."""
    left, right, width = horizontal_extent(app, selector)
    assert left >= 8 - 0.5, f"{selector} starts at {left:.1f} px"
    assert right <= width - 8 + 0.5, f"{selector} ends at {right:.1f} px of {width}"


def test_r12_the_list_menu_fits_a_320_px_window(app):
    """R12: the list's "⋯" menu fits a 320-px window."""
    app.set_viewport_size({"width": 320, "height": 800})
    start(app, SINGLE)
    add_region(app, "Barcode list")
    app.click("#list_more_button")
    app.locator("#list_more_menu").wait_for()
    assert_inside_window(app, "#list_more_menu")


@pytest.mark.parametrize("width", [600, 390])
def test_r12_the_orientation_help_fits_the_window(app, width):
    """R12: the orientation help opens fully inside narrow windows."""
    app.set_viewport_size({"width": width, "height": 800})
    start(app, PAIRED)
    app.get_by_role("button", name="About molecule orientation", exact=True).click()
    app.locator("#help_orientation").wait_for()
    assert_inside_window(app, "#help_orientation")


def test_r12_an_open_menu_moves_back_inside_when_the_window_narrows(app):
    """R12: the Templates menu, opened at 1280 px, still fits after a resize to 340 px."""
    app.set_viewport_size({"width": 1280, "height": 800})
    app.click("#templates_button")
    app.locator("#templates_menu").wait_for()
    app.set_viewport_size({"width": 340, "height": 800})
    next_frames(app)
    assert_inside_window(app, "#templates_menu")


# ---------------------------------------------------------------------------
# R13: the keyboard focus is kept or given back


@pytest.mark.parametrize(
    ("layout", "focus_key"), [(PAIRED, "layout_paired"), (SINGLE, "layout_single")]
)
def test_r13_after_the_welcome_choice_the_focus_is_on_the_read_layout_button(
    app, layout, focus_key
):
    """R13: choosing the layout from the keyboard leaves the focus on the matching button."""
    choice = app.locator("#welcome").get_by_role("button", name=layout, exact=True)
    assert choice.get_attribute("data-focus") == focus_key
    choice.focus()
    app.keyboard.press("Enter")
    chosen = app.get_by_role("group", name="Read layout").get_by_role(
        "button", name=layout, exact=True
    )
    assert chosen.get_attribute("aria-pressed") == "true"
    assert is_focused(chosen)


@pytest.mark.parametrize("kind", list(ADD_BUTTONS))
def test_r13_an_add_button_keeps_the_focus(app, kind):
    """R13: an add button keeps the focus, so Return adds another region."""
    start(app, SINGLE)
    button = app.locator(ADD_BUTTONS[kind])
    button.focus()
    app.keyboard.press("Enter")
    assert app.locator(".track button.blk").count() == 1
    assert is_focused(button)
    app.keyboard.press("Enter")
    assert app.locator(".track button.blk").count() == 2


@pytest.mark.parametrize(
    ("button_id", "menu_id"),
    [
        ("templates_button", "templates_menu"),
        ("download_menu_button", "download_menu"),
        ("list_more_button", "list_more_menu"),
    ],
)
def test_r13_escape_from_a_menu_gives_the_focus_back_to_its_button(
    app, button_id, menu_id
):
    """R13: Escape in a menu returns the focus to the button that opened it."""
    open_template(app, SIMPLE_EXAMPLE)
    app.locator(f"#{button_id}").focus()
    app.keyboard.press("Enter")
    menu = app.locator(f"#{menu_id}")
    menu.wait_for()
    menu.get_by_role("menuitem").first.focus()
    app.keyboard.press("Escape")
    menu.wait_for(state="hidden")
    assert is_focused(app.locator(f"#{button_id}"))


def test_r13_escape_from_the_templates_menu_opened_on_the_welcome_panel(app):
    """R13: the welcome link opens the menu at its first item; Escape returns to the link."""
    link = app.locator("#welcome_templates")
    assert link.get_attribute("aria-haspopup") == "menu"
    assert link.get_attribute("aria-expanded") == "false"
    link.focus()
    app.keyboard.press("Enter")
    menu = app.locator("#templates_menu")
    menu.wait_for()
    assert link.get_attribute("aria-expanded") == "true"
    assert is_focused(menu.get_by_role("menuitem").first)
    app.keyboard.press("Escape")
    menu.wait_for(state="hidden")
    assert is_focused(link)
    assert link.get_attribute("aria-expanded") == "false"


@pytest.mark.parametrize("help_id", ["orientation", "name"])
def test_r13_escape_from_a_help_popover_gives_the_focus_back_to_its_button(
    app, help_id
):
    """R13: Escape closes a help popover and focuses its "i" button."""
    start(app, PAIRED)
    add_region(app, "Barcode list")
    app.locator(f"#info_{help_id}").focus()
    app.keyboard.press("Enter")
    help_text = app.locator(f"#help_{help_id}")
    help_text.wait_for()
    app.keyboard.press("Escape")
    help_text.wait_for(state="hidden")
    assert is_focused(app.locator(f"#info_{help_id}"))


@pytest.mark.parametrize("answer", ["Cancel", "Escape"])
def test_r13_closing_a_question_from_the_templates_menu_focuses_the_templates_button(
    app, answer
):
    """R13: after a question that a Templates item opened, the focus is on Templates."""
    open_template(app, SIMPLE_EXAMPLE)
    choose_template(app, "SDR-seq gDNA")
    dialog = app.get_by_role("dialog")
    if answer == "Escape":
        dialog.wait_for()
        app.keyboard.press("Escape")
    else:
        dialog.get_by_role("button", name="Cancel", exact=True).click()
    wait_closed(app)
    assert is_focused(app.locator("#templates_button"))


# ---------------------------------------------------------------------------
# R14: the switch doesn't clip a focus ring

FOCUS_RING = """(button) => {
    const style = getComputedStyle(button);
    const reach = parseFloat(style.outlineOffset) + parseFloat(style.outlineWidth);
    const box = button.getBoundingClientRect();
    const ring = { left: box.left - reach, right: box.right + reach, top: box.top - reach, bottom: box.bottom + reach };
    const owner = button.parentElement;
    const owner_box = owner.getBoundingClientRect();
    const clipped_by = [];
    for (let element = owner; element && element !== document.documentElement; element = element.parentElement) {
        const own = getComputedStyle(element);
        if (own.overflowX === 'visible' && own.overflowY === 'visible') continue;
        const clip = element.getBoundingClientRect();
        if (ring.left < clip.left || ring.right > clip.right || ring.top < clip.top || ring.bottom > clip.bottom) {
            clipped_by.push(element.tagName.toLowerCase() + (element.id ? '#' + element.id : ''));
        }
    }
    return {
        outline: [style.outlineStyle, style.outlineWidth],
        ring,
        switch_box: { left: owner_box.left, right: owner_box.right, top: owner_box.top, bottom: owner_box.bottom },
        switch_overflow: [getComputedStyle(owner).overflowX, getComputedStyle(owner).overflowY],
        clipped_by,
    };
}"""


@pytest.mark.parametrize(
    ("switch_id", "side", "outer_edge"),
    [("region_switch", "Keep", "left"), ("rest_switch", "Leave out", "right")],
)
def test_r14_a_focused_side_of_a_switch_shows_its_whole_ring(
    app, switch_id, side, outer_edge
):
    """R14: the ring of a focused, unpressed side reaches past the switch and isn't clipped.

    The ring lies outside the switch's own box, so it shows only because the
    switch (and every box around it) lets it overflow.
    """
    start(app, SINGLE)
    add_region(app, "Barcode list")
    button = app.locator(f"#{switch_id}").get_by_role("button", name=side, exact=True)
    assert button.get_attribute("aria-pressed") == "false"
    focus_from_keyboard(app, button)
    ring = button.evaluate(FOCUS_RING)
    assert ring["outline"] == ["solid", "2px"]
    assert ring["switch_overflow"] == ["visible", "visible"]
    if outer_edge == "left":
        assert ring["ring"]["left"] < ring["switch_box"]["left"]
    else:
        assert ring["ring"]["right"] > ring["switch_box"]["right"]
    assert ring["ring"]["top"] < ring["switch_box"]["top"]
    assert ring["ring"]["bottom"] > ring["switch_box"]["bottom"]
    assert ring["clipped_by"] == []


# ---------------------------------------------------------------------------
# R15: text contrast of the add buttons


def test_r15_the_add_buttons_text_reaches_4_5_to_1_in_light_mode(app):
    """R15: each add button's text has at least 4.5:1 contrast on its tint (WCAG AA)."""
    assert contrast_ratio([0, 0, 0, 255], [255, 255, 255, 255]) == pytest.approx(21)
    app.emulate_media(color_scheme="light")
    start(app, SINGLE)
    ratios = {}
    for kind, selector in ADD_BUTTONS.items():
        button = app.locator(selector)
        text = computed_rgba(button, "color")
        background = computed_rgba(button, "background-color")
        assert text[3] == 255 and background[3] == 255, f"{kind} is not opaque"
        ratios[kind] = contrast_ratio(text, background)
    assert {
        kind: round(ratio, 2) for kind, ratio in ratios.items() if ratio < 4.5
    } == {}


# ---------------------------------------------------------------------------
# R16: forced colours keep the chosen side visible

FORCES_COLOURS = """() => {
    if (!matchMedia('(forced-colors: active)').matches) return false;
    const probe = document.createElement('span');
    probe.style.backgroundColor = 'rgb(1, 2, 3)';
    document.body.append(probe);
    const forced = getComputedStyle(probe).backgroundColor !== 'rgb(1, 2, 3)';
    probe.remove();
    return forced;
}"""

SYSTEM_HIGHLIGHT = """() => {
    const probe = document.createElement('span');
    probe.style.cssText = 'forced-color-adjust: none; background-color: Highlight';
    document.body.append(probe);
    const value = getComputedStyle(probe).backgroundColor;
    probe.remove();
    const context = document.createElement('canvas').getContext('2d');
    context.fillStyle = value;
    context.fillRect(0, 0, 1, 1);
    return Array.from(context.getImageData(0, 0, 1, 1).data);
}"""

VISIBLE_BACKGROUND = """(element) => {
    const layers = [];
    for (let node = element; node; node = node.parentElement) layers.unshift(getComputedStyle(node).backgroundColor);
    const context = document.createElement('canvas').getContext('2d');
    context.fillStyle = '#ffffff';
    context.fillRect(0, 0, 1, 1);
    for (const colour of layers) {
        context.fillStyle = '#010203';
        context.fillStyle = colour;
        if (context.fillStyle === '#010203') throw new Error('Unreadable colour: ' + colour);
        context.fillRect(0, 0, 1, 1);
    }
    return Array.from(context.getImageData(0, 0, 1, 1).data);
}"""


def test_r16_forced_colours_keep_the_chosen_side_of_a_switch_visible(app, browser_name):
    """R16: in forced colours the pressed side takes the system Highlight.

    Without that, forced colours paint both sides alike. The visible colour
    of each side is its background over every background behind it.
    """
    app.emulate_media(forced_colors="active")
    if not app.evaluate(FORCES_COLOURS):
        pytest.skip(
            f"{browser_name} matches (forced-colors: active) when emulated but "
            "doesn't force colours, so it can't show the difference"
        )
    start(app, SINGLE)
    add_region(app, "Barcode list")
    highlight = app.evaluate(SYSTEM_HIGHLIGHT)
    groups = {
        "rest switch": app.locator("#rest_switch"),
        "Read layout": app.get_by_role("group", name="Read layout"),
    }
    for name, group in groups.items():
        pressed = group.locator("button[aria-pressed='true']")
        unpressed = group.locator("button[aria-pressed='false']").first
        assert computed_rgba(pressed, "background-color") == highlight, name
        assert pressed.evaluate(VISIBLE_BACKGROUND) != unpressed.evaluate(
            VISIBLE_BACKGROUND
        ), name


# ---------------------------------------------------------------------------
# R17: a browser without replaceChildren


def test_r17_without_replace_children_the_page_names_the_browsers_it_needs(app):
    """R17: without Element.replaceChildren the banner names Chrome or Edge 86+ and start-up stops."""
    app.add_init_script("delete Element.prototype.replaceChildren;")
    app.reload()
    banner = app.locator("#feature_banner")
    banner.wait_for()
    assert banner.inner_text() == (
        "This browser lacks replaceChildren. "
        "Use Chrome or Edge 86+, Firefox 113+ or Safari 16.4+."
    )
    assert app.locator("#loading").is_hidden()
    assert app.locator("#layout").is_hidden()


# ---------------------------------------------------------------------------
# R18: colour fallbacks without color-mix()


def test_r18_without_color_mix_the_to_do_bar_still_shows(app, tmp_path):
    """R18: a browser without color-mix() drops those rules and uses the fallbacks.

    A copy of the page whose color-mix() calls are renamed to an unknown
    function stands in for such a browser.
    """
    text = GUI_PATH.read_text(encoding="utf-8")
    assert "color-mix(" in text
    copy = tmp_path / GUI_PATH.name
    copy.write_text(text.replace("color-mix(", "no-such-mix("), encoding="utf-8")
    app.goto(copy.as_uri())
    wait_ready(app)
    start(app, SINGLE)
    add_region(app, "Barcode list")
    tile = region_tile(app, "bc1")
    assert "has-todo" in classes(tile)
    assert computed_style(tile, "box-shadow")["box-shadow"] != "none"


# ---------------------------------------------------------------------------
# R19: text dropped on a field is the browser's


@pytest.mark.parametrize(
    ("field", "taken"), [("Name", False), ("Add a sequence", True)]
)
def test_r19_text_dropped_on_a_field_other_than_the_add_box_is_left_to_the_browser(
    app, field, taken
):
    """R19: a text drop on the Name field isn't taken; on the add box it imports."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    results = drag_from_outside(
        app.get_by_label(field, exact=True), text="AAAACCCC\nCCCCAAAA"
    )
    assert {event: result["prevented"] for event, result in results.items()} == {
        "dragenter": taken,
        "dragover": taken,
        "drop": taken,
    }
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    assert table_sequences(app) == (["AAAACCCC", "CCCCAAAA"] if taken else [])
    assert app.get_by_label("Name", exact=True).input_value() == "bc1"


# ---------------------------------------------------------------------------
# R20: no awkward wrapping, short placeholders


def test_r20_the_welcome_steps_keep_config_on_one_line(app):
    """R20: "--config" in the welcome steps never wraps after its hyphens."""
    code = app.locator("#welcome ol > li").nth(2).locator("code")
    assert code.inner_text() == "--config"
    assert classes(code) == ["nowrap"]
    assert computed_style(code, "white-space")["white-space"] == "nowrap"


@pytest.mark.parametrize(
    ("kind", "placeholder"),
    [
        ("Barcode list", "Type a barcode, or paste a list"),
        ("Constant region", "Type a variant, or paste several"),
    ],
)
def test_r20_the_add_box_has_a_short_placeholder_and_the_full_hint_as_title(
    app, kind, placeholder
):
    """R20: the add box's placeholder is short; its title gives the whole hint."""
    start(app, SINGLE)
    add_region(app, kind)
    box = app.get_by_label("Add a sequence", exact=True)
    assert box.get_attribute("id") == "add_sequence"
    assert box.get_attribute("placeholder") == placeholder
    assert box.get_attribute("title") == (
        "Press Enter to add what you typed. Pasting several lines, or several "
        "sequences on one line, adds them all."
    )


# ---------------------------------------------------------------------------
# R21: a one-line list file


@pytest.mark.parametrize(
    "separator", [", ", "\t", " ", ";"], ids=["comma", "tab", "space", "semicolon"]
)
def test_r21_a_one_line_list_file_imports_every_entry(app, tmp_path, separator):
    """R21: Import file… reads one line of several barcodes as a list, with no question."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    path = tmp_path / "barcodes.txt"
    path.write_text(separator.join(SAMPLE_BARCODES) + "\n", encoding="ascii")
    import_file(app, path)
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    assert table_sequences(app) == sorted(SAMPLE_BARCODES)


# ---------------------------------------------------------------------------
# R22: a constant region's first variant isn't a header for its length alone


def test_r22_four_cso_variants_import_into_a_constant_region_without_a_question(
    app, tmp_path
):
    """R22: variants of four lengths import into a constant region as they are."""
    start(app, SINGLE)
    add_region(app, "Constant region")
    path = tmp_path / "cso.txt"
    path.write_text("\n".join(CSO_VARIANTS) + "\n", encoding="ascii")
    import_file(app, path)
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    assert table_sequences(app) == CSO_VARIANTS


def test_r22_the_same_file_into_a_barcode_list_still_asks(app, tmp_path):
    """R22: in a barcode list, a first entry of a length no other has still asks."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    path = tmp_path / "cso.txt"
    path.write_text("\n".join(CSO_VARIANTS) + "\n", encoding="ascii")
    import_file(app, path)
    dialog, heading, _ = question(app)
    assert heading == f'Is "{CSO_VARIANTS[0]}" a header or a barcode?'
    dialog.get_by_role("button", name="A barcode", exact=True).click()
    wait_ready(app)
    assert app.locator("#sequence_count").inner_text() == "4 barcodes"


def test_r22_detect_header_row_ignores_length_when_told(core):
    """R22: length alone suggests a header only with by_length; "TAG" always does.

    plan_table_import passes by_length = false for a constant region.
    """
    rows = [
        {"line": number, "cells": [variant]}
        for number, variant in enumerate(CSO_VARIANTS, start=1)
    ]
    assert core.call("detect_header_row", rows, 0, True) == "ask"
    assert core.call("detect_header_row", rows, 0, False) is False
    as_list = core.call("plan_table_import", rows, "cso.txt", {}, "line", "barcodeList")
    assert as_list["status"] == "ask_header"
    as_constant = core.call(
        "plan_table_import", rows, "cso.txt", {}, "line", "constantRegion"
    )
    assert as_constant["status"] == "ok"
    assert as_constant["sequences"] == CSO_VARIANTS
    tagged = [{"line": 1, "cells": ["TAG"]}] + [
        {"line": number, "cells": [variant]}
        for number, variant in enumerate(OTHER_VARIANTS, start=2)
    ]
    tagged_plan = core.call(
        "plan_table_import", tagged, "tagged.txt", {}, "line", "constantRegion"
    )
    assert tagged_plan == {"status": "ask_header", "column": 0, "first": "TAG"}


# ---------------------------------------------------------------------------
# R23: unknown keys in a region without a name


def test_r23_unknown_keys_name_a_region_without_a_name_by_its_position(core):
    """R23: '"<key>" in region <k> of Read <n>' for a region without a name."""
    config = make_config(
        [barcode_list("bc1"), constant_region("con1")],
        [random_block("umi1"), barcode_list("bc2")],
    )
    config["barcode_struct_r1"]["blocks"][0]["colour"] = "green"
    del config["barcode_struct_r1"]["blocks"][1]["blockname"]
    config["barcode_struct_r1"]["blocks"][1]["colour"] = "orange"
    del config["barcode_struct_r2"]["blocks"][1]["blockname"]
    config["barcode_struct_r2"]["blocks"][1]["note"] = "plate 3"
    opened = core.call("config_from_json", config)
    keys = [
        '"colour" in region bc1',
        '"colour" in region 2 of Read 1',
        '"note" in region 2 of Read 2',
    ]
    assert opened["unknown_keys"] == keys
    assert {
        "level": "unknown_keys",
        "text": f"{tool_name()} doesn't use " + ", ".join(keys) + ".",
        "keys": keys,
    } in opened["notices"]


def test_r23_the_unknown_keys_notice_is_a_warning(app, tmp_path):
    """R23, R02: the page shows the unknown-keys notice as a warning list item."""
    config = make_config([barcode_list("bc1"), constant_region("con1")])
    config["barcode_struct_r1"]["blocks"][0]["colour"] = "green"
    del config["barcode_struct_r1"]["blocks"][1]["blockname"]
    config["barcode_struct_r1"]["blocks"][1]["colour"] = "orange"
    app.set_input_files("#open_input", str(write_json(tmp_path / "keys.json", config)))
    app.get_by_role("dialog").get_by_role(
        "button", name="Keep them", exact=True
    ).click()
    wait_ready(app)
    warning = app.locator("#notices li.notice-warning")
    assert warning.count() == 1
    assert warning.inner_text().startswith(
        f'{tool_name()} doesn\'t use "colour" in region bc1, '
        '"colour" in region 2 of Read 1.'
    )
    assert unreadable_notices(app, notice_items(app)) == []


# ---------------------------------------------------------------------------
# R24: a dropped line that is no clean list goes into the add box


def test_r24_a_dropped_line_with_an_invisible_character_goes_into_the_add_box(app):
    """R24: a dropped line that is neither a clean list nor one sequence is not imported.

    It goes into the add box, focused, where the user sees all of it.
    """
    zero_width_space = chr(0x200B)
    text = f"CCTTAACAT{zero_width_space} AGCGTAGAA"
    start(app, SINGLE)
    add_region(app, "Barcode list")
    results = drag_from_outside(selected_region_card(app), text=text)
    assert results["drop"]["prevented"]
    wait_ready(app)
    box = app.locator("#add_sequence")
    assert box.input_value() == text
    assert is_focused(box)
    assert app.get_by_role("dialog").count() == 0
    assert app.locator("#sequence_count").inner_text() == "0 barcodes"


# ---------------------------------------------------------------------------
# R25: dragenter is cancelled too


@pytest.mark.parametrize("carries", ["file", "text"])
def test_r25_a_list_card_cancels_dragenter(app, carries):
    """R25: a list's card cancels dragenter as well as dragover, and marks itself."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    card = selected_region_card(app)
    drag = (
        {"files": [("barcodes.txt", BARCODES_TXT)]}
        if carries == "file"
        else {"text": "\n".join(SAMPLE_BARCODES)}
    )
    results = drag_from_outside(card, events=("dragenter",), **drag)
    assert results["dragenter"]["prevented"]
    assert "drop-target" in classes(card)


# ---------------------------------------------------------------------------
# V04: the switch says what changes in the output reads, not what comes next


def test_v04_the_switches_and_outputs_say_in_output_reads(app):
    """V04: "In output reads: Keep | Leave out" on both switches; "Output reads" in Outputs.

    A kept region's tile and the rest tile say "in output reads".
    """
    start(app, SINGLE)
    add_list_and_constant(app)
    assert rest_row_parts(app) == [
        "Rest of Read 1",
        "In output reads:",
        "KeepLeave out",
    ]
    card_row = selected_region_card(app).locator(".align-row")
    assert card_row.locator(".label").inner_text() == "In output reads:"
    assert pressed_state(app.locator("#region_switch")) == {
        "Keep": "false",
        "Leave out": "true",
    }
    set_region_in_output_reads(app, "Keep")
    assert "in output reads" in region_tile(app, "con1").inner_text().splitlines()
    assert app.locator(".track .blk.rest span").inner_text() == "in output reads"
    assert app.locator("#out_files .hint").inner_text() == "Output reads"
    assert output_files(app) == ["Read 1: con1 + rest of read"]


# ---------------------------------------------------------------------------
# R26: the diagram grows when its tiles wrap as the window narrows

TILES_PAST_CAPTION = """() => {
    const tiles = Array.from(document.querySelectorAll('.track > *'));
    const bottom = Math.max(...tiles.map((tile) => tile.getBoundingClientRect().bottom));
    return bottom - document.querySelector('.caption').getBoundingClientRect().top;
}"""


def test_r26_the_caption_stays_below_the_tiles_while_the_window_narrows(app):
    """R26: narrowing the window step by step never draws a tile over the caption.

    WebKit didn't grow a grid row when the diagram's tiles wrapped on a
    resize, so a wrapped tile covered the caption until the next redraw.
    Checked from 1,520 to 640 px in 10-px steps, with SDR-seq RNA open and
    cb1 selected, and no redraw in between.
    """
    app.set_viewport_size({"width": 1520, "height": 900})
    open_template(app, "SDR-seq RNA")
    overlaps = []
    for width in range(1520, 639, -10):
        app.set_viewport_size({"width": width, "height": 900})
        past = app.evaluate(TILES_PAST_CAPTION)
        if past > 0.5:
            overlaps.append(f"{width} px: {past:.0f} px")
    assert overlaps == []


# ---------------------------------------------------------------------------
# R27: E03 names a space or a tab, and keeps it in the quoted sequence


def validate(core, config: dict) -> list[dict]:
    """Return validate_config's issues for a config opened as the page opens a file."""
    return core.run(
        "(core, config) => core.validate_config(core.config_from_json(config).model)",
        config,
    )


@pytest.mark.parametrize(
    ("separator", "named"), [(" ", "a space"), ("\t", "a tab")], ids=["space", "tab"]
)
def test_r27_e03_names_whitespace_and_keeps_it_in_the_quote(core, separator, named):
    """R27: a space or a tab is named, not called invisible, and the quote shows it."""
    sequence = f"GTCTCGTGGG{separator}CTCGGAGATG"
    config = make_config([barcode_list("bc1"), constant_region("con1", [sequence])])
    messages = [
        issue["message"] for issue in validate(core, config) if issue["rule"] == "E03"
    ]
    assert messages == [
        f'con1: "{sequence}" contains {named}; only A, C, G and T are allowed.'
    ]


# ---------------------------------------------------------------------------
# R28: a single-end config's Read 2 without regions needs no notice


def test_r28_a_single_end_config_gets_no_notice_about_read_2_without_regions(core):
    """R28: Read 2 isn't written for single-end, so only the "also describes Read 2" notice comes.

    Its keep_nonbarcode is still set to true, so switching to paired-end keeps it whole.
    """
    config = make_config(
        [barcode_list("bc1"), constant_region("con1")],
        [],
        single_end=True,
        keep_r1=True,
        keep_r2=False,
    )
    opened = core.call("config_from_json", config)
    assert [notice["text"] for notice in opened["notices"]] == [
        (
            "This single-end config also describes Read 2. Read 2 is kept, "
            "but it is not written while single-end is selected."
        )
    ]
    assert opened["model"]["reads"][1]["keep_nonbarcode"] is True


# ---------------------------------------------------------------------------
# R29: an unnamed region is numbered by its place among the regions


def test_r29_an_unnamed_regions_number_skips_items_that_are_not_regions(core):
    """R29: in [cb1, "x", unnamed], the unnamed region is region 2, as the diagram shows it."""
    config = make_config([barcode_list("cb1"), constant_region("con1")])
    unnamed = dict(config["barcode_struct_r1"]["blocks"][1], colour="green")
    del unnamed["blockname"]
    config["barcode_struct_r1"]["blocks"] = [
        config["barcode_struct_r1"]["blocks"][0],
        "x",
        unnamed,
    ]
    opened = core.call("config_from_json", config)
    assert opened["unknown_keys"] == ['"colour" in region 2 of Read 1']


# ---------------------------------------------------------------------------
# R30: keys at the top level say so


def test_r30_unknown_keys_at_the_top_level_say_where_they_are(core):
    """R30: a top-level key no longer reads as if it sat in the region named after it."""
    config = {
        "config_version": 2,
        **make_config([barcode_list("bc1"), constant_region("con1")]),
        "comment": "plate 3",
    }
    config["barcode_struct_r1"]["blocks"][1]["keepBlock"] = True
    opened = core.call("config_from_json", config)
    keys = [
        '"config_version" at the top level',
        '"comment" at the top level',
        '"keepBlock" in region con1',
    ]
    assert opened["unknown_keys"] == keys
    (notice,) = [n for n in opened["notices"] if n["level"] == "unknown_keys"]
    assert notice["text"] == f"{tool_name()} doesn't use " + ", ".join(keys) + "."


# ---------------------------------------------------------------------------
# R31: E11 for a typed value doesn't speak of a file


def test_r31_e11_for_a_typed_value_says_what_the_field_holds(app):
    """R31: clearing or mistyping Length or Allowed errors says "it is empty" or "it is 1.5"."""
    start(app, SINGLE)
    add_region(app, "Random (UMI)")

    def e11() -> list[str]:
        return [message for rule, _, message in checks(app) if rule == "E11"]

    fill(app, "Length (nt)", "")
    assert e11() == ["umi1: the length must be a whole number; it is empty."]
    fill(app, "Length (nt)", "1.5")
    assert e11() == ["umi1: the length must be a whole number; it is 1.5."]
    fill(app, "Length (nt)", "8")
    fill(app, "Allowed errors", "")
    assert e11() == ["umi1: allowed errors must be a whole number; it is empty."]


# ---------------------------------------------------------------------------
# R32: the focus stays on the control used, or moves on to the next one


def press_on(app, locator, key: str = "Enter") -> None:
    """Focus an element and press a key, as a keyboard user does in every engine."""
    locator.focus()
    app.keyboard.press(key)
    wait_ready(app)


def test_r32_a_region_tile_keeps_the_focus(app):
    """R32: Return on a tile selects its region and leaves the focus on the tile."""
    open_template(app, "SDR-seq RNA")
    press_on(app, region_tile(app, "cso"))
    assert region_tile(app, "cso").get_attribute("aria-pressed") == "true"
    assert is_focused(region_tile(app, "cso"))


def test_r32_move_right_keeps_the_focus_and_hands_it_over_at_the_end(app):
    """R32: Move right keeps the focus; at the last place, the focus goes to Move left."""
    open_template(app, "SDR-seq RNA")
    select_region(app, "commonseq2")
    move_right = app.get_by_role("button", name="Move right", exact=True)
    press_on(app, move_right)
    assert diagram_regions(app)[-2:] == ["commonseq2", "umi"]
    assert is_focused(move_right)
    press_on(app, move_right)
    assert diagram_regions(app)[-1] == "commonseq2"
    assert move_right.is_disabled()
    assert is_focused(app.get_by_role("button", name="Move left", exact=True))


@pytest.mark.parametrize(
    "control",
    ["read tab", "orientation", "order", "json"],
)
def test_r32_a_control_that_redraws_the_page_keeps_the_focus(app, control):
    """R32: the Read tabs, the orientation box, the list order and the JSON toggle keep the focus."""
    open_template(app, "SDR-seq RNA")
    if control == "read tab":
        target = app.get_by_role("tab", name="Read 2 · second file")
        press_on(app, target)
        assert target.get_attribute("aria-selected") == "true"
    elif control == "orientation":
        target = app.get_by_role(
            "checkbox", name="Molecules may be oriented either way"
        )
        before = target.is_checked()
        press_on(app, target, "Space")
        assert target.is_checked() != before
    elif control == "order":
        target = app.get_by_role("group", name="Order").get_by_role(
            "button", name="A→Z", exact=True
        )
        press_on(app, target)
        assert target.get_attribute("aria-pressed") == "true"
    else:
        press_on(app, app.get_by_role("button", name="▸ show", exact=True))
        target = app.get_by_role("button", name="▾ hide", exact=True)
    assert is_focused(target)


@pytest.mark.parametrize(
    ("button_id", "menu_id"),
    [
        ("templates_button", "templates_menu"),
        ("download_menu_button", "download_menu"),
        ("list_more_button", "list_more_menu"),
    ],
)
def test_r32_a_menu_opens_with_the_focus_on_its_first_item(app, button_id, menu_id):
    """R32: opening a menu from the keyboard moves the focus to its first item."""
    open_template(app, "SDR-seq RNA")
    press_on(app, app.locator(f"#{button_id}"))
    menu = app.locator(f"#{menu_id}")
    menu.wait_for()
    assert is_focused(menu.get_by_role("menuitem").first)


def test_r32_a_help_button_keeps_the_focus_and_is_described_by_its_text(app):
    """R32: the "i" keeps the focus while its help shows, and the help describes it."""
    open_template(app, "SDR-seq RNA")
    info = app.locator("#info_orientation")
    press_on(app, info)
    app.locator("#help_orientation").wait_for()
    assert is_focused(info)
    assert info.get_attribute("aria-describedby") == "help_orientation"


def test_r32_a_checks_item_moves_the_focus_to_its_regions_tile(app):
    """R32: Return on a Checks item about a region selects it and focuses its tile."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    add_region(app, "Constant region")
    item = app.locator("#checks .issue", has_text="Add barcodes to bc1.")
    press_on(app, item)
    assert region_tile(app, "bc1").get_attribute("aria-pressed") == "true"
    assert is_focused(region_tile(app, "bc1"))


def test_r32_unlink_moves_the_focus_to_add_link(app):
    """R32: after Unlink, the focus is on the "Add link with…" box that replaces the tag."""
    start(app, SINGLE)
    add_list_and_constant(app)
    add_region(app, "Random (UMI)")
    select_region(app, "bc1")
    app.get_by_label("Add link with", exact=True).select_option(label="umi1")
    wait_ready(app)
    press_on(app, app.get_by_role("button", name="Unlink umi1", exact=True))
    assert is_focused(app.get_by_label("Add link with", exact=True))


def test_r32_dismissing_the_notices_moves_the_focus_to_the_card(app):
    """R32: Dismiss removes the notices; the focus goes to the card's first control."""
    app.set_input_files("#open_input", str(FIXTURES_DIR / "legacy_master_gDNA.json.gz"))
    app.get_by_role("dialog").get_by_role("button", name="Open", exact=True).click()
    wait_ready(app)
    press_on(app, app.get_by_role("button", name="Dismiss", exact=True))
    assert app.locator("#notices").is_hidden()
    assert is_focused(
        app.get_by_role("group", name="Read layout").get_by_role(
            "button", name=PAIRED, exact=True
        )
    )


# ---------------------------------------------------------------------------
# R33: a long region name wraps instead of running out of its box

LONG_NAME = ("long_region_name_" * 4)[:64]
PAGE_SCROLLS_SIDEWAYS = (
    "document.documentElement.scrollWidth > document.documentElement.clientWidth"
)


@pytest.mark.parametrize("width", [390, 1280])
def test_r33_a_64_character_name_wraps_in_the_add_row_checks_and_questions(app, width):
    """R33: names of up to 64 characters, with no place to break, still fit the window.

    The add row and Checks name the region; so does the Delete question,
    whose buttons stay inside it.
    """
    assert len(LONG_NAME) == 64
    app.set_viewport_size({"width": width, "height": 900})
    start(app, SINGLE)
    add_region(app, "Barcode list")
    fill(app, "Name", LONG_NAME)
    wait_ready(app)
    assert f"Add barcodes to {LONG_NAME}." in [message for _, _, message in checks(app)]
    assert not app.evaluate(PAGE_SCROLLS_SIDEWAYS)

    app.click("#delete_region_button")
    dialog = app.get_by_role("dialog")
    dialog.wait_for()
    assert dialog.evaluate("(box) => box.scrollWidth <= box.clientWidth")
    dialog_box = dialog.bounding_box()
    button_box = dialog.get_by_role(
        "button", name="Delete region", exact=True
    ).bounding_box()
    assert (
        button_box["x"] + button_box["width"] <= dialog_box["x"] + dialog_box["width"]
    )
    app.keyboard.press("Escape")
    wait_closed(app)


# ---------------------------------------------------------------------------
# R34: Return in a question that starts on a choice


@pytest.mark.parametrize(("arrow", "layout"), [(None, PAIRED), ("ArrowDown", SINGLE)])
def test_r34_return_in_the_layout_question_opens_the_file(app, arrow, layout):
    """R34: in "Single- or paired-end?", Return on the chosen layout presses Open."""
    app.set_input_files("#open_input", str(FIXTURES_DIR / "legacy_master_gDNA.json.gz"))
    dialog = app.get_by_role("dialog")
    dialog.get_by_role("radio", name="Paired-end, barcodes on one read").focus()
    if arrow:
        app.keyboard.press(arrow)
    app.keyboard.press("Enter")
    wait_closed(app)
    wait_ready(app)
    chosen = app.get_by_role("group", name="Read layout").get_by_role(
        "button", name=layout, exact=True
    )
    assert chosen.get_attribute("aria-pressed") == "true"


def test_r34_return_on_a_chosen_column_imports_it(app, tmp_path):
    """R34: in "Which column holds the sequences?", Return on a chosen column presses Import."""
    rows = [
        f"{first},{second}"
        for first, second in zip(SAMPLE_BARCODES, EIGHT_BARCODES[4:])
    ]
    path = tmp_path / "two_columns.csv"
    path.write_text("\n".join(rows) + "\n", encoding="ascii")
    start(app, SINGLE)
    add_region(app, "Barcode list")
    import_file(app, path)
    dialog = app.get_by_role("dialog")
    first_column = dialog.get_by_role("radio").first
    first_column.focus()
    app.keyboard.press("Space")
    app.keyboard.press("Enter")
    wait_closed(app)
    wait_ready(app)
    assert sorted(table_sequences(app)) == sorted(SAMPLE_BARCODES)


# ---------------------------------------------------------------------------
# R35: a warning notice differs by shape as well as colour


def test_r35_a_warning_notice_has_a_marker_of_its_own_shape(app):
    """R35: warnings get a ⚠ marker, info notices the default disc."""
    app.set_input_files("#open_input", str(FIXTURES_DIR / "legacy_master_gDNA.json.gz"))
    app.get_by_role("dialog").get_by_role("button", name="Open", exact=True).click()
    wait_ready(app)
    markers = app.locator("#notices li").evaluate_all(
        "(items) => items.map((item) => [item.className, getComputedStyle(item).listStyleType])"
    )
    assert [name for name, _ in markers] == [
        "notice-info",
        "notice-info",
        "notice-warning",
    ]
    assert [style for name, style in markers if name == "notice-info"] == [
        "disc",
        "disc",
    ]
    (warning_style,) = [style for name, style in markers if name == "notice-warning"]
    assert chr(0x26A0) in warning_style


# ---------------------------------------------------------------------------
# R36: the move buttons and Delete region wrap together


@pytest.mark.parametrize("width", [360, 390, 414])
def test_r36_move_and_delete_buttons_stay_on_one_line(app, width):
    """R36: at phone widths ◀, ▶ and Delete region sit on the same line."""
    app.set_viewport_size({"width": width, "height": 900})
    open_template(app, "SDR-seq RNA")
    tops = [
        app.get_by_role("button", name=name, exact=True).bounding_box()["y"]
        for name in ("Move left", "Move right")
    ]
    tops.append(app.locator("#delete_region_button").bounding_box()["y"])
    assert max(tops) - min(tops) < 1


# ---------------------------------------------------------------------------
# R38: notices name a region without a usable name by its place


def test_r38_notices_name_a_region_without_a_usable_name_by_its_place(core):
    """R38: a region whose name is "" or null is "region <k> of Read <n>" in every notice."""
    config = make_config(
        [barcode_list("cb1"), constant_region("con1"), random_block("umi1")]
    )
    blocks = config["barcode_struct_r1"]["blocks"]
    blocks[1]["blockname"] = ""
    blocks[1]["sequence"] = [SAMPLE_CONSTANT.lower()]
    blocks[1]["name"] = ["linker"]
    blocks[2]["blockname"] = None
    del blocks[1]["keepblock"]
    del blocks[2]["keepblock"]
    texts = [
        notice["text"] for notice in core.call("config_from_json", config)["notices"]
    ]
    assert (
        '"In output reads" was not set for region 2 of Read 1, region 3 of Read 1, '
        'so they are "Leave out".'
    ) in texts
    assert (
        "Per-sequence names were dropped from region 2 of Read 1. "
        f"{tool_name()} uses them only to label simulated reads."
    ) in texts
    assert (
        "Lowercase sequences in region 2 of Read 1 were converted to uppercase."
        in texts
    )


# ---------------------------------------------------------------------------
# R39: a dragged tile is no text


def test_r39_a_tile_dropped_on_the_name_field_leaves_the_name_alone(app):
    """R39: the tile's drag carries no text, so a field it is dropped on doesn't change."""
    open_template(app, "SDR-seq RNA")
    region_tile(app, "cso").drag_to(app.locator("#block_name"))
    wait_ready(app)
    assert app.locator("#block_name").input_value() == "cb1"
    assert diagram_regions(app)[0] == "cb1"


# ---------------------------------------------------------------------------
# R40: a one-line file of a spaced constant sequence


def test_r40_a_one_line_spaced_sequence_file_goes_into_a_constant_regions_add_box(
    app, tmp_path
):
    """R40: spaces don't split a constant region's variant, so the line goes into the add box.

    Import file… used to offer the two halves as columns.
    """
    path = tmp_path / "constant_spaced.txt"
    path.write_text("GTCTCGTGGG CTCGGAGATG\n", encoding="ascii")
    start(app, SINGLE)
    add_region(app, "Barcode list")
    add_region(app, "Constant region")
    import_file(app, path)
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    box = app.locator("#add_sequence")
    assert box.input_value() == "GTCTCGTGGG CTCGGAGATG"
    assert is_focused(box)
    assert app.locator("#sequence_count").inner_text() == "0 variants"


# ---------------------------------------------------------------------------
# R41: a type change clears the row selection


def test_r41_changing_a_lists_type_clears_its_row_selection(app):
    """R41: rows selected in a list don't outlive a type change; Delete then deletes nothing and throws nothing.

    Delete, not Backspace: with nothing selected the page leaves Backspace to
    the browser, and WebKit then goes back a page.
    """
    errors = []
    app.on("pageerror", lambda error: errors.append(str(error)))
    start(app, SINGLE)
    add_region(app, "Barcode list")
    assert paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    app.locator("table.seqs tr").first.click()
    assert app.locator("#selection_count").inner_text() == "1 row selected"
    type_select = selected_region_card(app).locator("select[data-focus=type]")
    # A constant region keeps the sequences, but not the list's selection.
    type_select.select_option("constantRegion")
    wait_ready(app)
    assert app.locator("#sequence_count").inner_text() == "4 variants"
    assert app.locator("#selection_bar").is_hidden()
    app.locator("table.seqs tr").first.click()
    type_select.select_option("randomBarcode")
    wait_ready(app)
    app.locator(".caption").click()
    app.keyboard.press("Delete")
    type_select.select_option("barcodeList")
    wait_ready(app)
    assert app.locator("#selection_bar").is_hidden()
    assert errors == []


# ---------------------------------------------------------------------------
# R42: no import outline over a field


def test_r42_the_cards_drop_outline_goes_away_over_a_field(app):
    """R42: text over the Name field goes into the field, so the card stops offering to import it."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    card = selected_region_card(app)
    drag_from_outside(card, text="CCTTAACAT", events=("dragenter", "dragover"))
    assert "drop-target" in classes(card)
    results = drag_from_outside(
        app.locator("#block_name"), text="CCTTAACAT", events=("dragover",)
    )
    assert not results["dragover"]["prevented"]
    assert "drop-target" not in classes(card)


# ---------------------------------------------------------------------------
# R43: the download caret closes a menu open elsewhere


def test_r43_the_download_caret_closes_the_templates_menu_of_the_welcome_panel(app):
    """R43: the welcome link's menu closes, and the link says so, when the caret opens its menu."""
    link = app.locator("#welcome_templates")
    link.click()
    assert link.get_attribute("aria-expanded") == "true"
    app.click("#download_menu_button")
    assert app.locator("#download_menu").is_visible()
    assert app.locator("#templates_menu").is_hidden()
    assert link.get_attribute("aria-expanded") == "false"


def test_r43_the_download_caret_closes_the_list_menu(app):
    """R43: the list's ⋯ menu closes when the caret opens the download menu."""
    open_template(app, "SDR-seq RNA")
    app.click("#list_more_button")
    assert app.locator("#list_more_menu").is_visible()
    app.click("#download_menu_button")
    assert app.locator("#download_menu").is_visible()
    assert app.locator("#list_more_menu").is_hidden()
    assert app.locator("#list_more_button").get_attribute("aria-expanded") == "false"


# ---------------------------------------------------------------------------
# R44: the add box keeps its text through a redraw


def test_r44_text_in_the_add_box_stays_until_another_region_is_selected(app):
    """R44: opening help or changing the order redraws the card without emptying the add box."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    app.locator("#add_sequence").fill("CCTTAACAG")
    app.locator("#info_name").click()
    app.locator("#help_name").wait_for()
    assert app.locator("#add_sequence").input_value() == "CCTTAACAG"
    app.get_by_role("group", name="Order").get_by_role(
        "button", name="A→Z", exact=True
    ).click()
    assert app.locator("#add_sequence").input_value() == "CCTTAACAG"
    add_region(app, "Constant region")
    assert app.locator("#add_sequence").input_value() == ""


# ---------------------------------------------------------------------------
# R45: FL waits for every constant region


def test_r45_fl_is_unknown_while_a_constant_region_has_no_variant(core):
    """R45: describe_outputs marks FL unknown until every constant region has a valid variant."""
    script = "(core, config) => core.describe_outputs(core.config_from_json(config).model).FL"
    empty = core.run(
        script, make_config([barcode_list("bc1"), constant_region("con1", [])])
    )
    assert empty["known"] is False
    filled = core.run(
        script, make_config([barcode_list("bc1"), constant_region("con1")])
    )
    assert filled["known"] is True
    assert (filled["min"], filled["max"]) == (
        len(SAMPLE_CONSTANT),
        len(SAMPLE_CONSTANT),
    )


def test_r45_outputs_say_fl_shows_once_the_constant_region_is_filled_in(app):
    """R45: Outputs doesn't show FL as 0 for an empty constant region."""
    start(app, SINGLE)
    add_region(app, "Barcode list")
    add_region(app, "Constant region")
    assert (
        app.locator("#out_fl").inner_text()
        == "shows once every constant region is filled in"
    )
    type_sequence(app, SAMPLE_CONSTANT)
    assert app.locator("#out_fl").inner_text() == str(len(SAMPLE_CONSTANT))


# ---------------------------------------------------------------------------
# R46: a region the page named region1 follows its type


def test_r46_a_region_named_region1_takes_its_types_default_name(app, tmp_path):
    """R46: fixing the type of a region the page named region1 renames it as bc1, like any default name."""
    typo = dict(barcode_list("bc1"), blocktype="barcodelist")
    del typo["blockname"]
    config = make_config([typo, constant_region("con1")])
    app.set_input_files("#open_input", str(write_json(tmp_path / "typo.json", config)))
    wait_ready(app)
    assert diagram_regions(app)[0] == "region1"
    select_region(app, "region1")
    selected_region_card(app).locator("select[data-focus=type]").select_option(
        "barcodeList"
    )
    wait_ready(app)
    assert app.locator("#block_name").input_value() == "bc1"
