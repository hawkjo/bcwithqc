"""Tier B: clicking and typing like a user produces exactly the right file.

Each test drives the real page. Downloads are compared with the canonical
text of bcwithqc's own example configs, computed in Python.
"""

import gzip
import json
from pathlib import Path

import pytest
from gui_support import (
    EXAMPLES_DIR,
    FIXTURES_DIR,
    SAMPLE_BARCODES,
    SAMPLE_CONSTANT,
    canonical_text,
    design_doc_example,
    generated_with,
    read_config,
    without_names,
)

READY = "document.getElementById('app').dataset.busy === 'false'"


# ---------------------------------------------------------------------------
# Helpers that act like a user


def wait_ready(app):
    app.wait_for_function(READY)


def click_button(app, name):
    app.get_by_role("button", name=name, exact=True).click()


def add_block(app, kind):
    app.click("#add_block_button")
    app.get_by_role("menuitem", name=kind, exact=True).click()


def select_block(app, name):
    app.get_by_role("button", name=f"Block {name}", exact=True).click()


def fill(app, label, value):
    app.get_by_label(label, exact=True).fill(str(value))


def type_sequence(app, sequence):
    box = app.get_by_label("Add a sequence", exact=True)
    box.fill(sequence)
    box.press("Enter")


def paste_sequences(app, text):
    app.click("#add_sequences_button")
    app.get_by_role("menuitem", name="Paste…").click()
    dialog = app.get_by_role("dialog")
    dialog.get_by_label("Pasted sequences").fill(text)
    dialog.get_by_role("button", name="Import", exact=True).click()
    wait_ready(app)


def set_order(app, order):
    app.get_by_role("group", name="Order").get_by_role(
        "button", name=order, exact=True
    ).click()


def show_read(app, number):
    app.get_by_role(
        "tab", name=f"Read {number} · {'first' if number == 1 else 'second'} file"
    ).click()


def open_template(app, title):
    app.click("#templates_button")
    app.get_by_role("menuitem", name=title, exact=True).click()
    wait_ready(app)


def download_text(app, gzipped=False):
    with app.expect_download() as info:
        if gzipped:
            app.click("#download_menu_button")
            app.click("#download_gz_button")
        else:
            app.click("#download_button")
    data = Path(info.value.path()).read_bytes()
    return gzip.decompress(data).decode("ascii") if gzipped else data.decode("ascii")


def rule_ids(app):
    return app.locator("#checks .issue").evaluate_all(
        "(items) => items.map((item) => item.dataset.rule)"
    )


def table_sequences(app):
    return app.locator("table.seqs td.seq").all_inner_texts()


def example(name):
    return canonical_text(without_names(read_config(EXAMPLES_DIR / name)))


def build_read(app, prefix, keep_constant):
    """Add a 4-barcode list, a constant and a 6-nt UMI, named like pe_mini's."""
    add_block(app, "Barcode list")
    fill(app, "Name", f"bc{prefix}")
    paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    set_order(app, "Entered")
    add_block(app, "Constant region")
    fill(app, "Name", f"con{prefix}")
    type_sequence(app, SAMPLE_CONSTANT)
    if keep_constant:
        app.get_by_label("Keep for alignment", exact=True).check()
    add_block(app, "Random (UMI)")
    fill(app, "Name", f"umi{prefix}")
    fill(app, "Length (nt)", 6)


# ---------------------------------------------------------------------------
# Building configs from scratch


def test_build_pe_mini_from_scratch(app):
    click_button(app, "Paired-end")
    build_read(app, "_r1", keep_constant=True)
    show_read(app, 2)
    build_read(app, "_r2", keep_constant=False)
    assert rule_ids(app) == []
    assert download_text(app) == example("pe_mini_config.json")


def test_build_design_doc_example_from_scratch(app):
    click_button(app, "Paired-end")
    add_block(app, "Barcode list")
    fill(app, "Name", "cb1")
    paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    set_order(app, "Entered")
    add_block(app, "Constant region")
    fill(app, "Name", "con1")
    type_sequence(app, SAMPLE_CONSTANT.lower())
    app.get_by_label("Keep for alignment", exact=True).check()
    add_block(app, "Random (UMI)")
    fill(app, "Name", "umi")
    fill(app, "Length (nt)", 6)
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
    add_block(app, "Barcode list")
    fill(app, "Name", "cb1")
    app.set_input_files("#import_input", str(path))
    wait_ready(app)
    set_order(app, "Entered")
    add_block(app, "Constant region")
    type_sequence(app, SAMPLE_CONSTANT)
    app.get_by_label("Keep for alignment", exact=True).check()
    add_block(app, "Random (UMI)")
    fill(app, "Name", "umi")
    fill(app, "Length (nt)", 6)
    app.get_by_label("Send the rest of Read 1 to alignment").check()
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
    click_button(app, "Single-end")
    add_block(app, "Barcode list")
    for sequence in sequences:
        type_sequence(app, sequence)
    add_block(app, "Constant region")
    type_sequence(app, SAMPLE_CONSTANT)
    select_block(app, "bc1")


def test_download_is_disabled_while_an_error_exists(app):
    assert app.locator("#download_button").is_disabled()
    start_single_end_list(app, ["CCTTAACAT", "AGCGTAGAA", "CCTTAACAT"])
    assert "E04" in rule_ids(app)
    assert app.locator("#download_button").is_disabled()
    app.locator("table.seqs tr").first.click()
    click_button(app, "Delete selected")
    assert "E04" not in rule_ids(app)
    assert app.locator("#download_button").is_enabled()


def test_clicking_a_message_selects_its_rows(app):
    start_single_end_list(app, ["CCTTAACAT", "AGCGTAGAA", "CCTTAACAT", "TATAGGTGT"])
    select_block(app, "con1")
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
    assert "has-error" in app.get_by_role("button", name="Block bc1").get_attribute(
        "class"
    )


def test_typed_lowercase_is_stored_uppercase(app):
    start_single_end_list(app, ["ccttaacat"])
    assert table_sequences(app) == ["CCTTAACAT"]


def test_contiguity_note_and_structure_caption_are_shown(app):
    open_template(app, "SDR-seq RNA")
    caption = app.locator("#structure_caption").inner_text()
    assert caption.startswith(
        "71–74 nt · 4 linker combinations · reads must be ≥ 67 nt"
    )
    assert "Blocks start at base 1 and must be contiguous" in caption


# ---------------------------------------------------------------------------
# Blocks: names, order, decoder labels, links


def json_block_names(app):
    app.get_by_role("button", name="▸ show").click()
    config = json.loads(app.locator("#json_preview").inner_text())
    app.get_by_role("button", name="▾ hide").click()
    return [block["blockname"] for block in config["barcode_struct_r1"]["blocks"]]


def test_new_blocks_get_default_names(app):
    click_button(app, "Single-end")
    add_block(app, "Barcode list")
    add_block(app, "Constant region")
    add_block(app, "Random (UMI)")
    add_block(app, "Barcode list")
    assert json_block_names(app) == ["bc1", "con1", "umi1", "bc2"]


def test_move_buttons_and_drag_and_drop_reorder_the_json(app):
    click_button(app, "Single-end")
    add_block(app, "Barcode list")
    add_block(app, "Constant region")
    add_block(app, "Random (UMI)")
    click_button(app, "Move left")
    assert json_block_names(app) == ["bc1", "umi1", "con1"]
    app.get_by_role("button", name="Block con1").drag_to(
        app.get_by_role("button", name="Block bc1")
    )
    assert json_block_names(app) == ["con1", "bc1", "umi1"]


def test_random_block_after_a_list_turns_its_decoder_to_free(app):
    click_button(app, "Single-end")
    add_block(app, "Barcode list")
    add_block(app, "Constant region")
    list_box = app.get_by_role("button", name="Block bc1")
    assert "Lev ≤ 1" in list_box.inner_text()
    select_block(app, "bc1")
    add_block(app, "Random (UMI)")
    assert "FREE ≤ 1" in app.get_by_role("button", name="Block bc1").inner_text()


def link_tags(app):
    return app.locator(".links .tag").all_inner_texts()


def test_linking_three_blocks_and_unlinking(app):
    click_button(app, "Single-end")
    for _ in range(3):
        add_block(app, "Barcode list")
    add_block(app, "Constant region")
    select_block(app, "bc1")
    app.get_by_label("Add link with").select_option(label="bc2")
    app.get_by_label("Add link with").select_option(label="bc3")
    assert [tag.replace("✕", "").strip() for tag in link_tags(app)] == ["bc2", "bc3"]
    app.get_by_role("button", name="Unlink bc3").click()
    assert [tag.replace("✕", "").strip() for tag in link_tags(app)] == ["bc2"]
    app.get_by_role("button", name="Unlink bc2").click()
    assert link_tags(app) == []
    assert "E12" not in rule_ids(app)
    for name in ("bc1", "bc2", "bc3"):
        assert (
            "link" not in app.get_by_role("button", name=f"Block {name}").inner_text()
        )


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
    set_order(app, "Entered")
    rows = app.locator("table.seqs tr")
    rows.nth(0).click()
    rows.nth(2).click(modifiers=["Shift"])
    assert app.locator(".selbar").inner_text().startswith("3 rows selected")
    click_button(app, "Delete selected")
    assert table_sequences(app) == [SAMPLE_BARCODES[3]]


def test_sort_toggle_controls_display_and_file_order(app):
    entered = ["TTTTCCCC", "AAAACCCC", "GGGGCCCC", "CCCCAAA"]
    start_single_end_list(app, entered)
    select_block(app, "con1")
    type_sequence(app, "TTTGGGAAACCC")
    select_block(app, "bc1")
    # A→Z: the 7-nt row (a warning) comes first, the rest alphabetically; the file is plain A→Z.
    assert table_sequences(app) == ["CCCCAAA", "AAAACCCC", "GGGGCCCC", "TTTTCCCC"]
    config = json.loads(download_text(app))
    assert config["barcode_struct_r1"]["blocks"][0]["sequence"] == sorted(entered)
    set_order(app, "Entered")
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
    add_block(app, "Barcode list")
    paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    add_block(app, "Constant region")
    type_sequence(app, SAMPLE_CONSTANT)
    assert rule_ids(app) == ["W12"]
    config = json.loads(download_text(app))
    assert config["barcode_struct_r1"]["blocks"] == []


def test_single_block_read_downloads_with_w11(app):
    click_button(app, "Single-end")
    add_block(app, "Barcode list")
    paste_sequences(app, "\n".join(SAMPLE_BARCODES))
    assert rule_ids(app) == ["W11"]
    assert app.locator("#checks .chip.temporary").count() == 1
    assert len(json.loads(download_text(app))["barcode_struct_r1"]["blocks"]) == 1


def test_single_end_hides_read_2_and_its_errors(app):
    open_template(app, "SDR-seq RNA")
    show_read(app, 2)
    add_block(app, "Barcode list")
    assert app.locator("#download_button").is_disabled()
    click_button(app, "Single-end")
    assert app.get_by_role("tab").count() == 0
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
    app.set_input_files("#open_input", str(FIXTURES_DIR / "legacy_master_gDNA.json.gz"))
    dialog = app.get_by_role("dialog")
    assert dialog.get_by_label("Paired-end, barcodes on one read").is_checked()
    dialog.get_by_role("button", name="Open").click()
    wait_ready(app)
    notices = app.locator("#notices").inner_text()
    assert "older format" in notices
    assert "names were dropped" in notices
    assert (
        app.get_by_role("group", name="Order")
        .get_by_role("button", name="Entered")
        .get_attribute("aria-pressed")
        == "true"
    )
    assert download_text(app) == example("gDNA_paired_end_config.json")


def test_opening_a_file_that_is_not_a_config_shows_an_error(app, tmp_path):
    path = tmp_path / "notes.json"
    path.write_text("[1, 2, 3]")
    app.set_input_files("#open_input", str(path))
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
    start_single_end_list(app, ["CCTTAACAT"])
    open_template(app, "SDR-seq RNA")
    dialog = app.get_by_role("dialog")
    dialog.get_by_role("button", name="Keep editing").click()
    assert app.locator("#main .block-head .name").inner_text() == "bc1"
    open_template(app, "SDR-seq RNA")
    app.get_by_role("dialog").get_by_role("button", name="Discard changes").click()
    wait_ready(app)
    assert app.get_by_role("button", name="Block cb1").count() == 1


def test_control_s_downloads(app):
    open_template(app, "Design-doc example")
    with app.expect_download() as info:
        app.keyboard.press("Control+s")
    assert Path(info.value.path()).read_text() == canonical_text(design_doc_example())


def test_closing_after_an_edit_prompts(app):
    open_template(app, "Design-doc example")
    fill(app, "Name", "cb_renamed")
    dialogs = []
    app.on("dialog", lambda dialog: (dialogs.append(dialog.type), dialog.dismiss()))
    app.close(run_before_unload=True)
    app.wait_for_timeout(500)
    assert dialogs == ["beforeunload"]


# ---------------------------------------------------------------------------
# Help and notices


def test_orientation_help_opens_without_toggling_the_checkbox(app):
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
    open_template(app, "Design-doc example")
    app.get_by_role("button", name="Rules for block names").click()
    text = app.locator("#help_name").inner_text()
    assert (
        "letters, digits" in text and "64 characters" in text and "CB1 and cb1" in text
    )
    app.locator("h1").click()
    assert app.locator("#help_name").is_hidden()


def open_with_extras(app, tmp_path, extras):
    """Open se_mini with extra top-level keys; return (config, its dialog)."""
    config = without_names(read_config(EXAMPLES_DIR / "se_mini_config.json"))
    path = tmp_path / "with_extras.json"
    path.write_text(json.dumps({**config, **extras}, indent=4))
    app.set_input_files("#open_input", str(path))
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
    assert '"config_version", "comment"' in dialog.inner_text()
    dialog.get_by_role("button", name="Leave them out", exact=True).click()
    wait_ready(app)
    keep = app.get_by_label("Keep them in the file", exact=True)
    assert not keep.is_checked()
    assert download_text(app) == canonical_text(config)
    keep.check()
    assert download_text(app) == canonical_text({**config, **extras})


@pytest.mark.parametrize("how", ["Escape", "Cancel"])
def test_cancelling_the_unknown_keys_question_opens_nothing(app, tmp_path, how):
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
    app.set_input_files("#open_input", str(path))
    dialog = app.get_by_role("dialog")
    dialog.get_by_role("button", name="Open", exact=True).click()
    dialog.get_by_role("button", name="Leave it out", exact=True).click()
    wait_ready(app)
    assert not app.get_by_label("Keep it in the file", exact=True).is_checked()
    assert download_text(app) == example("gDNA_paired_end_config.json")


def test_a_downloaded_config_reopens_without_questions(app, tmp_path):
    open_template(app, "SDR-seq gDNA")
    path = tmp_path / "gDNA.json"
    path.write_text(download_text(app))
    open_template(app, "SDR-seq RNA")
    app.set_input_files("#open_input", str(path))
    wait_ready(app)
    assert app.locator("#error_banner").is_hidden()
    stamp = generated_with()
    made_for = f"made for {stamp['name']} {stamp['version']}"
    assert made_for in app.locator("#notices").inner_text()
    assert download_text(app) == example("gDNA_paired_end_config.json")


def test_shortcuts_wait_while_a_question_is_open(app, tmp_path):
    open_template(app, "SDR-seq gDNA")
    downloads = []
    app.on("download", lambda download: downloads.append(download.suggested_filename))
    _, dialog = open_with_extras(app, tmp_path, {"config_version": 2})
    app.keyboard.press("ControlOrMeta+s")
    assert app.locator("#app").get_attribute("data-busy") == "true"
    dialog.get_by_role("button", name="Keep it", exact=True).click()
    wait_ready(app)
    download_text(app)
    assert downloads == ["with_extras.json"]


def test_outputs_without_a_umi_block_say_ub_and_ur_are_not_written(app):
    open_template(app, "SDR-seq gDNA")
    assert app.locator("#out_ub").inner_text() == "not written (no UMI block)"
    assert app.locator("#out_ur").inner_text().startswith("not written (no UMI block)")
