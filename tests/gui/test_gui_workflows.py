"""Tier B: clicking and typing like a user produces exactly the right file.

Each test drives the real page. Downloads are compared with the canonical
text of bcwithqc's own example configs, computed in Python. The page calls
bcwithqc's blocks "regions"; the config keys keep their names.
"""

import json
from pathlib import Path

import pytest
from gui_support import (
    EXAMPLES_DIR,
    FIXTURES_DIR,
    SAMPLE_BARCODES,
    SAMPLE_CONSTANT,
    SIMPLE_EXAMPLE,
    add_region,
    canonical_text,
    click_button,
    design_doc_example,
    download_text,
    fill,
    generated_with,
    open_config,
    open_template,
    paste_sequences,
    read_config,
    rule_ids,
    select_region,
    set_order,
    set_region_in_output_reads,
    set_rest_in_output_reads,
    show_read,
    table_sequences,
    type_sequence,
    wait_ready,
    without_names,
)

NAME_RULES = (
    'Use letters, digits, "_", "." and "-", up to 64 characters. Two regions in '
    "the same read can't share a name, even in different case: CB1 and cb1 clash."
)
CONTIGUITY_NOTE = (
    "Regions start at base 1 and follow each other without gaps. "
    "Anything after the last region is the rest of the read."
)
# The caption joins each read length to its "≥" and "nt" with a no-break space.
NBSP = chr(0xA0)


# ---------------------------------------------------------------------------
# Helpers for this file


def example(name):
    """Return the canonical text of a bcwithqc example config, without per-sequence names."""
    return canonical_text(without_names(read_config(EXAMPLES_DIR / name)))


def build_read(app, prefix, keep_constant):
    """Add a 4-barcode list, a constant region and a 6-nt UMI to the shown read, named like pe_mini's.

    Args:
        app: The page.
        prefix: Text added to each region's name, such as "_r1" (bc_r1, con_r1, umi_r1).
        keep_constant: Set the constant region to "Keep" in the output reads.
    """
    add_region(app, "Barcode list")
    fill(app, "Name", f"bc{prefix}")
    assert paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    set_order(app, "As entered")
    add_region(app, "Constant region")
    fill(app, "Name", f"con{prefix}")
    type_sequence(app, SAMPLE_CONSTANT)
    if keep_constant:
        set_region_in_output_reads(app, "Keep")
    add_region(app, "Random (UMI)")
    fill(app, "Name", f"umi{prefix}")
    fill(app, "Length (nt)", 6)


# ---------------------------------------------------------------------------
# Building configs from scratch


def test_build_pe_mini_from_scratch(app):
    click_button(app, "Paired-end")
    build_read(app, "_r1", keep_constant=True)
    # pe_mini leaves out the rest of Read 1; both reads output it by default.
    set_rest_in_output_reads(app, "Leave out")
    show_read(app, 2)
    build_read(app, "_r2", keep_constant=False)
    assert rule_ids(app) == []
    assert download_text(app) == example("pe_mini_config.json")


def test_build_design_doc_example_from_scratch(app):
    click_button(app, "Paired-end")
    add_region(app, "Barcode list")
    fill(app, "Name", "cb1")
    assert paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    set_order(app, "As entered")
    add_region(app, "Constant region")
    fill(app, "Name", "con1")
    type_sequence(app, SAMPLE_CONSTANT.lower())
    set_region_in_output_reads(app, "Keep")
    add_region(app, "Random (UMI)")
    fill(app, "Name", "umi")
    fill(app, "Length (nt)", 6)
    # The example leaves out the rest of Read 1 and outputs all of Read 2.
    set_rest_in_output_reads(app, "Leave out")
    assert download_text(app) == canonical_text(design_doc_example())


def test_build_se_mini_by_importing_an_xlsx(app, tmp_path):
    from openpyxl import Workbook

    workbook = Workbook()
    sheet = workbook.active
    sheet.append(["name", "barcode"])
    for index, barcode in enumerate(SAMPLE_BARCODES):
        sheet.append([f"bc{index + 1}", barcode])
    path = tmp_path / "se_mini_barcodes.xlsx"
    workbook.save(path)

    click_button(app, "Single-end")
    add_region(app, "Barcode list")
    fill(app, "Name", "cb1")
    app.set_input_files("#import_input", str(path))
    wait_ready(app)
    set_order(app, "As entered")
    add_region(app, "Constant region")
    type_sequence(app, SAMPLE_CONSTANT)
    set_region_in_output_reads(app, "Keep")
    add_region(app, "Random (UMI)")
    fill(app, "Name", "umi")
    fill(app, "Length (nt)", 6)
    # se_mini outputs the rest of Read 1, which is the default.
    assert download_text(app) == example("se_mini_config.json")


def test_json_gz_download_holds_the_same_text(app):
    open_template(app, "SDR-seq gDNA")
    assert download_text(app, gzipped=True) == example("gDNA_paired_end_config.json")


def test_templates_download_their_sources(app):
    open_template(app, "SDR-seq RNA")
    assert download_text(app) == example("cDNA_paired_end_config.json")


# ---------------------------------------------------------------------------
# Errors, checks and selection


def start_single_end_list(app, sequences):
    """Choose single-end, add list bc1 with these sequences and constant region con1, then select bc1."""
    click_button(app, "Single-end")
    add_region(app, "Barcode list")
    for sequence in sequences:
        type_sequence(app, sequence)
    add_region(app, "Constant region")
    type_sequence(app, SAMPLE_CONSTANT)
    select_region(app, "bc1")


def test_download_is_disabled_while_an_error_exists(app):
    download_button = app.locator("#download_button")
    assert download_button.is_disabled()
    start_single_end_list(app, ["CCTTAACAT", "AGCGTAGAA", "CCTTAACAT"])
    assert "E04" in rule_ids(app)
    assert download_button.is_disabled()
    assert download_button.get_attribute("title") == "Fix the errors in Checks first"
    app.locator("table.seqs tr").first.click()
    click_button(app, "Delete selected")
    assert "E04" not in rule_ids(app)
    assert download_button.is_enabled()


def test_clicking_a_message_selects_its_rows(app):
    start_single_end_list(app, ["CCTTAACAT", "AGCGTAGAA", "CCTTAACAT", "TATAGGTGT"])
    select_region(app, "con1")
    app.locator('#checks .issue[data-rule="E04"]').click()
    selected = app.locator(
        'table.seqs tr[aria-selected="true"] td.seq'
    ).all_inner_texts()
    assert selected == ["CCTTAACAT", "CCTTAACAT"]
    assert app.locator("#main .block-head .name").inner_text() == "bc1"


def test_problem_rows_and_blocks_are_highlighted(app):
    start_single_end_list(app, ["CCTTAACAT", "CCTTAACAT", "ACCTGGA"])
    assert app.locator("table.seqs tr.err").count() == 2
    assert app.locator("table.seqs tr.warn").count() == 1
    region = app.get_by_role("button", name="Region bc1", exact=True)
    assert "has-error" in region.get_attribute("class")


def test_typed_lowercase_is_stored_uppercase(app):
    start_single_end_list(app, ["ccttaacat"])
    assert table_sequences(app) == ["CCTTAACAT"]


def test_contiguity_note_and_structure_caption_are_shown(app):
    """The caption gives Read 1's length, linker combinations and read lengths, then the note.

    SDR-seq RNA's Read 1 holds cb1 (9 nt), cso (4 variants of 14–17 nt), cb2
    (9 nt), commonseq2 (1 variant, 23 nt), sb (8 nt) and umi (8 nt); each
    allows 1 error except umi (0). Each choice of constant-region variants is
    one linker combination. A read must be at least the shortest total minus
    every allowed error, and is recommended to be the longest total plus every
    allowed error.
    """
    open_template(app, "SDR-seq RNA")
    shortest = 9 + 14 + 9 + 23 + 8 + 8  # 71
    longest = 9 + 17 + 9 + 23 + 8 + 8  # 74
    allowed_errors = 1 + 1 + 1 + 1 + 1 + 0  # 5
    linker_combinations = 4 * 1  # cso's variants times commonseq2's
    assert app.locator("#structure_caption").inner_text() == (
        f"{shortest}–{longest} nt · {linker_combinations} linker combinations · "
        f"reads ≥{NBSP}{shortest - allowed_errors}{NBSP}nt "
        f"(≥{NBSP}{longest + allowed_errors}{NBSP}nt recommended) · {CONTIGUITY_NOTE}"
    )


# ---------------------------------------------------------------------------
# Regions: names, order, decoder labels, links


def json_block_names(app):
    """Return Read 1's region names in the order the JSON preview lists them."""
    app.get_by_role("button", name="▸ show").click()
    config = json.loads(app.locator("#json_preview").inner_text())
    app.get_by_role("button", name="▾ hide").click()
    return [block["blockname"] for block in config["barcode_struct_r1"]["blocks"]]


def test_new_blocks_get_default_names(app):
    click_button(app, "Single-end")
    add_region(app, "Barcode list")
    add_region(app, "Constant region")
    add_region(app, "Random (UMI)")
    add_region(app, "Barcode list")
    assert json_block_names(app) == ["bc1", "con1", "umi1", "bc2"]


def test_move_buttons_and_drag_and_drop_reorder_the_json(app):
    click_button(app, "Single-end")
    add_region(app, "Barcode list")
    add_region(app, "Constant region")
    add_region(app, "Random (UMI)")
    click_button(app, "Move left")
    assert json_block_names(app) == ["bc1", "umi1", "con1"]
    app.get_by_role("button", name="Region con1", exact=True).drag_to(
        app.get_by_role("button", name="Region bc1", exact=True)
    )
    assert json_block_names(app) == ["con1", "bc1", "umi1"]


def test_random_block_after_a_list_turns_its_decoder_to_free(app):
    click_button(app, "Single-end")
    add_region(app, "Barcode list")
    add_region(app, "Constant region")
    list_box = app.get_by_role("button", name="Region bc1", exact=True)
    assert "Lev ≤ 1" in list_box.inner_text()
    select_region(app, "bc1")
    add_region(app, "Random (UMI)")
    assert "FREE ≤ 1" in list_box.inner_text()


def link_tags(app):
    """Return the selected region's link tags as shown, each ending in its ✕ button."""
    return app.locator(".links .tag").all_inner_texts()


def test_linking_three_blocks_and_unlinking(app):
    click_button(app, "Single-end")
    for _ in range(3):
        add_region(app, "Barcode list")
    add_region(app, "Constant region")
    select_region(app, "bc1")
    app.get_by_label("Add link with").select_option(label="bc2")
    app.get_by_label("Add link with").select_option(label="bc3")
    assert [tag.replace("✕", "").strip() for tag in link_tags(app)] == ["bc2", "bc3"]
    app.get_by_role("button", name="Unlink bc3").click()
    assert [tag.replace("✕", "").strip() for tag in link_tags(app)] == ["bc2"]
    app.get_by_role("button", name="Unlink bc2").click()
    assert link_tags(app) == []
    assert "E12" not in rule_ids(app)
    for name in ("bc1", "bc2", "bc3"):
        region = app.get_by_role("button", name=f"Region {name}", exact=True)
        assert "link" not in region.inner_text()


# ---------------------------------------------------------------------------
# The sequence table


def test_table_shows_count_and_every_length(app):
    start_single_end_list(app, ["CCTTAACAT", "AGCGTAGAA"])
    assert app.locator("#sequence_count").inner_text() == "2 barcodes"
    assert app.locator("#sequence_lengths").inner_text() == "9 nt"
    type_sequence(app, "ACCTGGA")
    assert app.locator("#sequence_count").inner_text() == "3 barcodes"
    assert app.locator("#sequence_lengths").inner_text() == "9 nt ×2 · 7 nt ×1"
    app.locator("table.seqs tr", has_text="ACCTGGA").click()
    app.keyboard.press("Delete")
    assert app.locator("#sequence_lengths").inner_text() == "9 nt"
    assert app.locator("table.seqs tr").count() == 2


def test_table_has_no_number_column_and_no_row_delete_buttons(app):
    start_single_end_list(app, SAMPLE_BARCODES)
    first_row = app.locator("table.seqs tr").first
    assert first_row.locator("td").count() == 3
    assert app.locator("table.seqs button").count() == 0


def test_shift_click_selects_a_range_and_delete_selected_removes_it(app):
    start_single_end_list(app, SAMPLE_BARCODES)
    set_order(app, "As entered")
    rows = app.locator("table.seqs tr")
    rows.nth(0).click()
    rows.nth(2).click(modifiers=["Shift"])
    assert app.locator(".selbar").inner_text().startswith("3 rows selected")
    click_button(app, "Delete selected")
    assert table_sequences(app) == [SAMPLE_BARCODES[3]]


def test_sort_toggle_controls_display_and_file_order(app):
    entered = ["TTTTCCCC", "AAAACCCC", "GGGGCCCC", "CCCCAAA"]
    start_single_end_list(app, entered)
    select_region(app, "con1")
    type_sequence(app, "TTTGGGAAACCC")
    select_region(app, "bc1")
    # A→Z: the 7-nt row (a warning) comes first, the rest alphabetically; the file is plain A→Z.
    assert table_sequences(app) == ["CCCCAAA", "AAAACCCC", "GGGGCCCC", "TTTTCCCC"]
    config = json.loads(download_text(app))
    assert config["barcode_struct_r1"]["blocks"][0]["sequence"] == sorted(entered)
    set_order(app, "As entered")
    assert table_sequences(app) == entered
    config = json.loads(download_text(app))
    assert config["barcode_struct_r1"]["blocks"][0]["sequence"] == entered
    assert config["barcode_struct_r1"]["blocks"][1]["sequence"] == [
        SAMPLE_CONSTANT,
        "TTTGGGAAACCC",
    ]


def test_double_click_edits_a_sequence(app):
    start_single_end_list(app, ["CCTTAACAT"])
    app.locator("table.seqs td.seq").first.dblclick()
    editor = app.get_by_label("Edit sequence")
    editor.fill("aaccttggt")
    editor.press("Enter")
    assert table_sequences(app) == ["AACCTTGGT"]


# ---------------------------------------------------------------------------
# Reads and modes


def test_barcodes_on_read_2_only_download_with_w12(app):
    click_button(app, "Paired-end")
    show_read(app, 2)
    add_region(app, "Barcode list")
    assert paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    add_region(app, "Constant region")
    type_sequence(app, SAMPLE_CONSTANT)
    assert rule_ids(app) == ["W12"]
    config = json.loads(download_text(app))
    assert config["barcode_struct_r1"]["blocks"] == []


def test_single_block_read_downloads_with_w11(app):
    click_button(app, "Single-end")
    add_region(app, "Barcode list")
    assert paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    assert rule_ids(app) == ["W11"]
    assert app.locator("#checks .chip.temporary").count() == 1
    assert len(json.loads(download_text(app))["barcode_struct_r1"]["blocks"]) == 1


def test_single_end_hides_read_2_and_its_errors(app):
    open_template(app, "SDR-seq RNA")
    show_read(app, 2)
    add_region(app, "Barcode list")
    type_sequence(app, "ACGTN")
    bad_character = app.locator('#checks .issue[data-rule="E03"]')
    assert bad_character.get_attribute("data-level") == "error"
    assert app.locator("#download_button").is_disabled()
    click_button(app, "Single-end")
    assert app.get_by_role("tab").count() == 0
    assert rule_ids(app) == []
    assert app.locator("#download_button").is_enabled()
    config = json.loads(download_text(app))
    assert "barcode_struct_r2" not in config


def test_orientation_choice_survives_a_switch_to_single_end_and_back(app):
    open_template(app, "SDR-seq RNA")
    app.get_by_label("Molecules may be oriented either way").check()
    click_button(app, "Single-end")
    assert json.loads(download_text(app))["unknown_read_orientation"] is False
    click_button(app, "Paired-end")
    assert app.get_by_label("Molecules may be oriented either way").is_checked()
    assert json.loads(download_text(app))["unknown_read_orientation"] is True


# ---------------------------------------------------------------------------
# Opening files and importing


def test_opening_a_legacy_config_asks_the_mode_and_lists_the_changes(app):
    open_config(app, FIXTURES_DIR / "legacy_master_gDNA.json.gz")
    dialog = app.get_by_role("dialog")
    assert dialog.get_by_label("Paired-end, barcodes on one read").is_checked()
    dialog.get_by_role("button", name="Open").click()
    wait_ready(app)
    notices = app.locator("#notices").inner_text()
    assert "older format" in notices
    assert "all of Read 2 is kept in the output reads" in notices
    # No region in the file has keepblock.
    assert (
        '"In output reads" was not set in this file, so every region is '
        '"Leave out".' in notices
    )
    assert "names were dropped" in notices
    assert (
        app.get_by_role("group", name="Order")
        .get_by_role("button", name="As entered")
        .get_attribute("aria-pressed")
        == "true"
    )
    assert download_text(app) == example("gDNA_paired_end_config.json")


def test_opening_a_file_that_is_not_a_config_shows_an_error(app, tmp_path):
    path = tmp_path / "notes.json"
    path.write_text("[1, 2, 3]")
    open_config(app, path)
    wait_ready(app)
    assert "not a bcwithqc config" in app.locator("#error_banner").inner_text()


def test_two_dna_columns_make_the_user_choose(app, tmp_path):
    start_single_end_list(app, [])
    path = tmp_path / "two.csv"
    path.write_text("CCTTAACAT,ATGTTAAGG\nAGCGTAGAA,TTCTACGCT\n")
    app.set_input_files("#import_input", str(path))
    dialog = app.get_by_role("dialog")
    import_button = dialog.get_by_role("button", name="Import")
    assert import_button.is_disabled()
    dialog.get_by_label("Column 2").check()
    import_button.click()
    wait_ready(app)
    assert table_sequences(app) == ["ATGTTAAGG", "TTCTACGCT"]


@pytest.mark.parametrize(
    ("answer", "expected_count"), [("A header", 4), ("A barcode", 5)]
)
def test_tag_in_a_single_column_asks(app, tmp_path, answer, expected_count):
    start_single_end_list(app, [])
    path = tmp_path / "tag.txt"
    path.write_text("TAG\n" + "\n".join(SAMPLE_BARCODES) + "\n")
    app.set_input_files("#import_input", str(path))
    app.get_by_role("dialog").get_by_role("button", name=answer).click()
    wait_ready(app)
    assert app.locator("table.seqs tr").count() == expected_count


@pytest.mark.parametrize(
    ("answer", "expected"), [("Add to the existing 1", 2), ("Replace them", 1)]
)
def test_import_into_a_list_asks_append_or_replace_without_a_default(
    app, tmp_path, answer, expected
):
    start_single_end_list(app, ["CCTTAACAT"])
    path = tmp_path / "one.txt"
    path.write_text("AGCGTAGAA\n")
    app.set_input_files("#import_input", str(path))
    dialog = app.get_by_role("dialog")
    dialog.get_by_role("button", name=answer).wait_for()
    # Neither answer is pre-selected: focus starts on Cancel.
    assert app.evaluate("document.activeElement.textContent") == "Cancel"
    dialog.get_by_role("button", name=answer).click()
    wait_ready(app)
    assert app.locator("table.seqs tr").count() == expected


def test_workbook_with_several_sheets_needs_sheet_and_column(app, tmp_path):
    from openpyxl import Workbook

    workbook = Workbook()
    first = workbook.active
    first.title = "Plate A"
    first.append(["barcode"])
    first.append(["CCTTAACAT"])
    second = workbook.create_sheet("Plate B")
    second.append(["id", "barcode"])
    second.append(["1", "GGGGAAAA"])
    second.append(["2", "TTTTCCCC"])
    path = tmp_path / "plates.xlsx"
    workbook.save(path)

    start_single_end_list(app, [])
    app.set_input_files("#import_input", str(path))
    dialog = app.get_by_role("dialog")
    import_button = dialog.get_by_role("button", name="Import")
    assert import_button.is_disabled()
    dialog.get_by_label("Sheet").select_option("Plate B")
    assert import_button.is_disabled()
    dialog.get_by_label("Column 2").check()
    import_button.click()
    wait_ready(app)
    assert table_sequences(app) == ["GGGGAAAA", "TTTTCCCC"]


# ---------------------------------------------------------------------------
# Leaving and saving


def test_loading_a_template_over_unsaved_changes_asks_first(app):
    """Cancel keeps the page as it is; Replace loads the template."""
    start_single_end_list(app, ["CCTTAACAT"])
    app.click("#templates_button")
    app.get_by_role("menuitem", name="SDR-seq RNA", exact=True).click()
    dialog = app.get_by_role("dialog")
    assert dialog.get_by_role("heading").inner_text() == "Replace the current config?"
    assert dialog.locator("p").inner_text() == (
        'Opening the template "SDR-seq RNA" replaces the config on this page. '
        "It has changes you have not downloaded."
    )
    # Cancel is the default answer: the filled button, with the focus, so
    # Enter alone replaces nothing.
    assert dialog.locator("button.primary").all_inner_texts() == ["Cancel"]
    assert app.evaluate("document.activeElement.textContent") == "Cancel"
    dialog.get_by_role("button", name="Cancel", exact=True).click()
    assert dialog.is_hidden()
    assert app.locator("#main .block-head .name").inner_text() == "bc1"
    assert app.get_by_role("button", name="Region cb1", exact=True).count() == 0
    open_template(app, "SDR-seq RNA")
    assert app.get_by_role("button", name="Region cb1", exact=True).count() == 1
    assert app.get_by_role("button", name="Region bc1", exact=True).count() == 0


def test_control_s_downloads(app):
    open_template(app, SIMPLE_EXAMPLE)
    with app.expect_download() as info:
        app.keyboard.press("Control+s")
    assert info.value.suggested_filename == "simple_example.json"
    assert Path(info.value.path()).read_text() == canonical_text(design_doc_example())


def test_closing_after_an_edit_prompts(app):
    open_template(app, SIMPLE_EXAMPLE)
    fill(app, "Name", "cb_renamed")
    dialogs = []
    app.on("dialog", lambda dialog: (dialogs.append(dialog.type), dialog.dismiss()))
    app.close(run_before_unload=True)
    app.wait_for_timeout(500)
    assert dialogs == ["beforeunload"]


# ---------------------------------------------------------------------------
# Help and notices


def test_orientation_help_opens_without_toggling_the_checkbox(app):
    # The welcome panel has no settings; they appear once a read layout is chosen.
    click_button(app, "Paired-end")
    checkbox = app.get_by_label("Molecules may be oriented either way")
    help_button = app.get_by_role("button", name="About molecule orientation")
    help_button.click()
    assert app.locator("#help_orientation").is_visible()
    assert (
        "Each molecule may be oriented either way"
        in app.locator("#help_orientation").inner_text()
    )
    assert not checkbox.is_checked()
    help_button.click()
    assert app.locator("#help_orientation").is_hidden()
    help_button.click()
    app.keyboard.press("Escape")
    assert app.locator("#help_orientation").is_hidden()


def test_name_help_lists_the_rules(app):
    open_template(app, SIMPLE_EXAMPLE)
    app.get_by_role("button", name="Rules for region names").click()
    assert app.locator("#help_name").inner_text() == NAME_RULES
    app.locator("h1").click()
    assert app.locator("#help_name").is_hidden()


def open_with_extras(app, tmp_path, extras):
    """Open se_mini with extra top-level keys and wait for the unknown-keys question.

    A page that holds a config first asks whether to replace it; this
    answers Replace.

    Returns:
        The config without the extras, and the dialog locator.
    """
    config = without_names(read_config(EXAMPLES_DIR / "se_mini_config.json"))
    path = tmp_path / "with_extras.json"
    path.write_text(json.dumps({**config, **extras}, indent=4))
    open_config(app, path)
    dialog = app.get_by_role("dialog")
    dialog.get_by_role("button", name="Cancel", exact=True).wait_for()
    return config, dialog


def test_opening_a_file_with_an_unknown_key_asks_whether_to_keep_it(app, tmp_path):
    config, dialog = open_with_extras(app, tmp_path, {"config_version": 2})
    assert "Keep the key bcwithqc doesn't use?" in dialog.inner_text()
    assert '"config_version"' in dialog.inner_text()
    dialog.get_by_role("button", name="Keep it", exact=True).click()
    wait_ready(app)
    keep = app.get_by_label("Keep it in the file", exact=True)
    assert keep.is_checked()
    assert download_text(app) == canonical_text({**config, "config_version": 2})
    keep.uncheck()
    assert download_text(app) == canonical_text(config)


def test_unknown_keys_left_out_on_opening_can_be_kept_later(app, tmp_path):
    extras = {"config_version": 2, "comment": "plate 3"}
    config, dialog = open_with_extras(app, tmp_path, extras)
    assert (
        '"config_version" at the top level, "comment" at the top level'
        in dialog.inner_text()
    )
    dialog.get_by_role("button", name="Leave them out", exact=True).click()
    wait_ready(app)
    keep = app.get_by_label("Keep them in the file", exact=True)
    assert not keep.is_checked()
    assert download_text(app) == canonical_text(config)
    keep.check()
    assert download_text(app) == canonical_text({**config, **extras})


@pytest.mark.parametrize("how", ["Escape", "Cancel"])
def test_cancelling_the_unknown_keys_question_opens_nothing(app, tmp_path, how):
    """Answering Replace, then cancelling the unknown-keys question, keeps the template."""
    open_template(app, "SDR-seq RNA")
    _, dialog = open_with_extras(app, tmp_path, {"config_version": 2})
    if how == "Escape":
        app.keyboard.press("Escape")
    else:
        dialog.get_by_role("button", name="Cancel", exact=True).click()
    wait_ready(app)
    assert dialog.is_hidden()
    assert download_text(app) == example("cDNA_paired_end_config.json")


def test_a_legacy_file_with_an_unknown_key_asks_both_questions_in_turn(app, tmp_path):
    legacy = read_config(FIXTURES_DIR / "legacy_master_gDNA.json.gz")
    path = tmp_path / "legacy_with_comment.json"
    path.write_text(json.dumps({**legacy, "comment": "plate 3"}, indent=4))
    open_config(app, path)
    dialog = app.get_by_role("dialog")
    dialog.get_by_role("button", name="Open", exact=True).click()
    dialog.get_by_role("button", name="Leave it out", exact=True).click()
    wait_ready(app)
    assert not app.get_by_label("Keep it in the file", exact=True).is_checked()
    assert download_text(app) == example("gDNA_paired_end_config.json")


def test_a_downloaded_config_reopens_without_questions(app, tmp_path):
    """A file the page wrote opens without the older-format or unknown-keys question.

    The page holds a config, so the replace question comes first and is
    answered with Replace.
    """
    open_template(app, "SDR-seq gDNA")
    path = tmp_path / "gDNA.json"
    path.write_text(download_text(app))
    open_template(app, "SDR-seq RNA")
    open_config(app, path)
    wait_ready(app)
    assert app.get_by_role("dialog").count() == 0
    assert app.locator("#error_banner").is_hidden()
    stamp = generated_with()
    made_for = f"made for {stamp['name']} {stamp['version']}"
    assert made_for in app.locator("#notices").inner_text()
    assert download_text(app) == example("gDNA_paired_end_config.json")


def test_shortcuts_wait_while_a_question_is_open(app, tmp_path):
    """Ctrl/Cmd+S waits while the unknown-keys question is open, after Replace was answered."""
    open_template(app, "SDR-seq gDNA")
    downloads = []
    app.on("download", lambda download: downloads.append(download.suggested_filename))
    _, dialog = open_with_extras(app, tmp_path, {"config_version": 2})
    assert dialog.get_by_role("heading").inner_text() == (
        "Keep the key bcwithqc doesn't use?"
    )
    app.keyboard.press("ControlOrMeta+s")
    assert app.locator("#app").get_attribute("data-busy") == "true"
    dialog.get_by_role("button", name="Keep it", exact=True).click()
    wait_ready(app)
    download_text(app)
    assert downloads == ["with_extras.json"]


def test_outputs_without_a_umi_block_say_ub_and_ur_are_not_written(app):
    open_template(app, "SDR-seq gDNA")
    assert app.locator("#out_ub").inner_text() == "not written (no UMI region)"
    assert app.locator("#out_ur").inner_text().startswith("not written (no UMI region)")
