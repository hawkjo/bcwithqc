"""Tier A: reading and writing configs.

The GUI must reproduce bcwithqc's own configs exactly (without per-sequence
name lists), migrate older formats to them, and write the text Python's
json.dumps(config, indent=4) writes. Expected texts are computed in Python.
"""

import base64
import copy
import gzip
import importlib.util
import json

import pytest
from gui_support import (
    EXAMPLES_DIR,
    FIXTURES_DIR,
    GUI_PATH,
    OLD_GUI_BARCODES,
    OLD_GUI_SAMPLE_BARCODES,
    REPO_ROOT,
    bcwithqc_version,
    canonical_text,
    design_doc_example,
    example_config_names,
    read_bytes_text,
    read_config,
    with_sorted_lists,
    without_names,
    write_old_gui_config,
)

OPEN_AND_EXPORT = """(core, arg) => {
    const opened = core.config_from_json(arg.text, arg.options || {});
    if (!opened.ok) return {error: opened.error};
    return {text: core.config_to_json(opened.model), notices: opened.notices, needs_mode_choice: opened.needs_mode_choice};
}"""

OPEN_SORT_AND_EXPORT = """(core, arg) => {
    const model = core.config_from_json(arg.text).model;
    for (const read of model.reads) for (const block of read.blocks) if (block.blocktype === 'barcodeList') block.sort = 'az';
    return core.config_to_json(model);
}"""


def open_and_export(core, text, **options):
    result = core.run(OPEN_AND_EXPORT, {"text": text, "options": options})
    assert "error" not in result, result.get("error")
    return result


@pytest.mark.parametrize("example", example_config_names())
def test_example_exports_its_canonical_text_in_entered_order(core, example):
    config = read_config(EXAMPLES_DIR / example)
    result = open_and_export(core, (EXAMPLES_DIR / example).read_text(encoding="utf-8"))
    assert result["text"] == canonical_text(without_names(config))
    # The examples carry no bcwithqc_version, so nothing may report one.
    assert not [n for n in result["notices"] if "made for bcwithqc" in n["text"]]


@pytest.mark.parametrize("example", example_config_names())
def test_example_exports_sorted_lists_when_switched_to_a_to_z(core, example):
    config = read_config(EXAMPLES_DIR / example)
    exported = core.run(
        OPEN_SORT_AND_EXPORT,
        {"text": (EXAMPLES_DIR / example).read_text(encoding="utf-8")},
    )
    assert exported == canonical_text(with_sorted_lists(without_names(config)))


def test_examples_cover_unsorted_lists_and_multi_variant_constants():
    """The two tests above would be weak if no example had an unsorted list or several constant variants."""
    config = read_config(EXAMPLES_DIR / "cDNA_paired_end_config.json")
    blocks = config["barcode_struct_r1"]["blocks"]
    assert any(
        block["blocktype"] == "barcodeList"
        and block["sequence"] != sorted(block["sequence"])
        for block in blocks
    )
    assert any(
        block["blocktype"] == "constantRegion"
        and block["sequence"] != sorted(block["sequence"])
        for block in blocks
    )


def test_dropping_names_is_reported(core):
    result = open_and_export(
        core, (EXAMPLES_DIR / "se_mini_config.json").read_text(encoding="utf-8")
    )
    dropped = [
        notice["text"]
        for notice in result["notices"]
        if "names were dropped" in notice["text"]
    ]
    assert dropped and "cb1" in dropped[0]


@pytest.mark.parametrize(
    ("fixture", "example", "r1_keep_nonbarcode"),
    [
        ("legacy_master_cDNA.json", "cDNA_paired_end_config.json", None),
        ("legacy_master_gDNA.json.gz", "gDNA_paired_end_config.json", None),
        ("legacy_sdrranger_gDNA.json.gz", "gDNA_paired_end_config.json", False),
    ],
)
def test_legacy_config_migrates_to_the_paired_end_example(
    core, fixture, example, r1_keep_nonbarcode
):
    expected = without_names(read_config(EXAMPLES_DIR / example))
    if r1_keep_nonbarcode is not None:
        expected["barcode_struct_r1"]["keep_nonbarcode"] = r1_keep_nonbarcode
    result = open_and_export(
        core, read_bytes_text(FIXTURES_DIR / fixture), legacy_mode="paired"
    )
    assert result["text"] == canonical_text(expected)


def test_legacy_config_without_read_2_asks_for_the_mode(core):
    text = read_bytes_text(FIXTURES_DIR / "legacy_master_gDNA.json.gz")
    result = open_and_export(core, text)
    assert result["needs_mode_choice"] is True
    single = open_and_export(core, text, legacy_mode="single")
    expected = without_names(read_config(EXAMPLES_DIR / "gDNA_single_end_config.json"))
    assert single["text"] == canonical_text(expected)
    assert single["needs_mode_choice"] is False


def test_legacy_orientation_true_is_turned_off_with_a_notice(core):
    legacy = read_config(FIXTURES_DIR / "legacy_master_cDNA.json")
    assert legacy["unknown_read_orientation"] is True
    result = open_and_export(core, json.dumps(legacy), legacy_mode="paired")
    assert json.loads(result["text"])["unknown_read_orientation"] is False
    assert any("Orientation was on" in notice["text"] for notice in result["notices"])


def test_legacy_barcode_read_that_keeps_nothing_gets_the_alignment_notice(core):
    result = open_and_export(
        core,
        read_bytes_text(FIXTURES_DIR / "legacy_sdrranger_gDNA.json.gz"),
        legacy_mode="paired",
    )
    assert any(
        "align only the other read's file" in notice["text"]
        for notice in result["notices"]
    )


def test_old_gui_file_migrates_exactly(core, tmp_path):
    path = tmp_path / "old_gui.json.gz"
    write_old_gui_config(path)
    expected = {
        "unknown_read_orientation": False,
        "single_end_reads": False,
        "barcode_struct_r1": {
            "keep_nonbarcode": False,
            "blocks": [
                {
                    "maxerrors": 1,
                    "blockfunction": "barcode",
                    "sequence": OLD_GUI_BARCODES,
                    "blocktype": "barcodeList",
                    "blockname": "sequence_0",
                    "keepblock": False,
                    "linkgroup": 0,
                },
                {
                    "maxerrors": 1,
                    "blockfunction": "barcode",
                    "sequence": ["GGGTCAGTACGTACGAGTCCC"],
                    "blocktype": "constantRegion",
                    "blockname": "sequence_1",
                    "keepblock": False,
                },
                {
                    "maxerrors": 1,
                    "blockfunction": "discard",
                    "sequence": OLD_GUI_SAMPLE_BARCODES,
                    "blocktype": "barcodeList",
                    "blockname": "sequence_2",
                    "keepblock": False,
                    "linkgroup": 0,
                },
                {
                    "maxerrors": 0,
                    "blockfunction": "UMI",
                    "length": 8,
                    "blocktype": "randomBarcode",
                    "blockname": "sequence_3",
                    "keepblock": False,
                },
            ],
        },
        "barcode_struct_r2": {"keep_nonbarcode": True, "blocks": []},
    }
    result = core.run(
        """async (core, arg) => {
            const bytes = Uint8Array.from(atob(arg), (ch) => ch.charCodeAt(0));
            const opened = core.config_from_json(await core.read_config_bytes(bytes));
            return {text: core.config_to_json(opened.model), notices: opened.notices.map((n) => n.text)};
        }""",
        base64.b64encode(path.read_bytes()).decode("ascii"),
    )
    assert result["text"] == canonical_text(expected)
    notices = " | ".join(result["notices"])
    assert "sequence 0 → sequence_0" in notices
    assert "names were dropped from sequence_1, sequence_2" in notices


def test_gdna_single_and_paired_end_convert_into_each_other(core):
    single_text = (EXAMPLES_DIR / "gDNA_single_end_config.json").read_text(
        encoding="utf-8"
    )
    paired_text = (EXAMPLES_DIR / "gDNA_paired_end_config.json").read_text(
        encoding="utf-8"
    )
    switch = """(core, arg) => {
        const model = core.config_from_json(arg.text).model;
        model.single_end = arg.single_end;
        return core.config_to_json(model);
    }"""
    assert core.run(
        switch, {"text": single_text, "single_end": False}
    ) == canonical_text(
        without_names(read_config(EXAMPLES_DIR / "gDNA_paired_end_config.json"))
    )
    assert core.run(
        switch, {"text": paired_text, "single_end": True}
    ) == canonical_text(
        without_names(read_config(EXAMPLES_DIR / "gDNA_single_end_config.json"))
    )


def test_switching_to_single_end_and_back_restores_read_2_and_orientation(core):
    result = core.run(
        """(core, arg) => {
            const model = core.config_from_json(arg).model;
            model.unknown_read_orientation = true;
            model.single_end = true;
            const single = core.config_to_json(model);
            model.single_end = false;
            return {single, paired: core.config_to_json(model)};
        }""",
        (EXAMPLES_DIR / "pe_mini_config.json").read_text(encoding="utf-8"),
    )
    single = json.loads(result["single"])
    assert "barcode_struct_r2" not in single
    assert single["single_end_reads"] is True
    expected_paired = without_names(read_config(EXAMPLES_DIR / "pe_mini_config.json"))
    expected_paired["unknown_read_orientation"] = True
    assert result["paired"] == canonical_text(expected_paired)


def test_single_end_writes_orientation_false_until_se_orientation_is_supported(core):
    """bcwithqc supports orientation only for paired-end reads (notes #se-orientation).

    Delete this test, and write the user's choice, once bcwithqc supports it.
    """
    exported = core.run(
        """(core) => {
            const model = core.new_model();
            model.single_end = true;
            model.unknown_read_orientation = true;
            return core.config_to_json(model);
        }"""
    )
    assert json.loads(exported)["unknown_read_orientation"] is False


def expanded_templates(core):
    return core.run("(core) => window.barcode_templates")


def template_source_config(source):
    if source.endswith(".md"):
        return design_doc_example()
    return read_config(REPO_ROOT / source)


@pytest.mark.parametrize("template_id", ["sdr_seq_rna", "sdr_seq_gdna", "design_doc"])
def test_template_equals_its_source_without_names(core, template_id):
    templates = {template["id"]: template for template in expanded_templates(core)}
    template = templates[template_id]
    exported = core.run(
        "(core, arg) => core.config_to_json(core.config_from_json(arg).model)",
        template["config"],
    )
    assert exported == canonical_text(
        without_names(template_source_config(template["source"]))
    )


def test_design_doc_template_is_read_from_the_doc():
    templates = load_update_script().TEMPLATE_SOURCES
    assert ("design_doc", "Design-doc example", "doc/config_format.md") in templates


def load_update_script():
    spec = importlib.util.spec_from_file_location(
        "update_template_data", REPO_ROOT / "gui" / "update_template_data.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_page_template_data_is_current():
    """Fails when a template source changed; run gui/update_template_data.py to refresh the page."""
    page = GUI_PATH.read_text(encoding="utf-8")
    match = load_update_script().DATA_BLOCK.search(page)
    assert match is not None
    assert json.loads(match.group(2)) == load_update_script().template_data()


def test_json_gz_download_decompresses_to_the_canonical_text(core):
    text = (EXAMPLES_DIR / "pe_mini_config.json").read_text(encoding="utf-8")
    encoded = core.run(
        """async (core, arg) => {
            const bytes = await core.gzip_text(core.config_to_json(core.config_from_json(arg).model));
            let binary = '';
            for (const byte of bytes) binary += String.fromCharCode(byte);
            return btoa(binary);
        }""",
        text,
    )
    assert gzip.decompress(base64.b64decode(encoded)).decode("ascii") == canonical_text(
        without_names(read_config(EXAMPLES_DIR / "pe_mini_config.json"))
    )


def test_python_gzipped_config_opens(core):
    config = without_names(read_config(EXAMPLES_DIR / "se_mini_config.json"))
    payload = base64.b64encode(
        gzip.compress(canonical_text(config).encode("ascii"))
    ).decode("ascii")
    exported = core.run(
        """async (core, arg) => {
            const bytes = Uint8Array.from(atob(arg), (ch) => ch.charCodeAt(0));
            return core.config_to_json(core.config_from_json(await core.read_config_bytes(bytes)).model);
        }""",
        payload,
    )
    assert exported == canonical_text(config)


def linked_read(groups):
    blocks = []
    for index, group in enumerate(groups):
        block = {
            "maxerrors": 1,
            "blockfunction": "barcode",
            "sequence": ["ACGTACGT"],
            "blocktype": "barcodeList",
            "blockname": f"bc{index + 1}",
            "keepblock": False,
        }
        if group is not None:
            block["linkgroup"] = group
        blocks.append(block)
    return {
        "unknown_read_orientation": False,
        "single_end_reads": True,
        "barcode_struct_r1": {"keep_nonbarcode": False, "blocks": blocks},
    }


@pytest.mark.parametrize(
    ("groups", "expected"),
    [
        ([5, None, 9, 5, 9], [0, None, 1, 0, 1]),
        ([7, 7], [0, 0]),
        ([None, 3, 3], [None, 0, 0]),
    ],
)
def test_linkgroups_are_renumbered_from_zero_in_order(core, groups, expected):
    exported = json.loads(
        core.run(
            "(core, arg) => core.config_to_json(core.config_from_json(arg).model)",
            linked_read(groups),
        )
    )
    assert [
        block.get("linkgroup") for block in exported["barcode_struct_r1"]["blocks"]
    ] == expected


@pytest.mark.parametrize(
    "text",
    [
        "plain",
        'quote " and backslash \\',
        "tab\tnewline\ncarriage\rreturn",
        "\b\f\x01\x1f\x7f",
        "é ü ß",
        "€ ✓",
        "😀",
        "\u200b",
    ],
)
def test_strings_are_escaped_as_python_json_dumps_escapes_them(core, text):
    assert core.call("python_json_string", text) == json.dumps(text)


def test_nested_values_match_python_json_dumps(core):
    value = {
        "a": [],
        "b": {},
        "c": [1, -2, [3, []]],
        "d": {"e": True, "f": False, "g": None},
        "h": "x",
    }
    assert core.call("python_json_dumps", value) == json.dumps(value, indent=4)


def test_unknown_keys_are_listed_and_written_last_only_when_kept(core):
    config = without_names(read_config(EXAMPLES_DIR / "se_mini_config.json"))
    first_block = config["barcode_struct_r1"]["blocks"][0]
    # In the file, each unknown key comes first in its section.
    in_file = {"config_version": 1, **copy.deepcopy(config)}
    in_file["barcode_struct_r1"] = {
        "note": "plate 3",
        **copy.deepcopy(config["barcode_struct_r1"]),
    }
    in_file["barcode_struct_r1"]["blocks"][0] = {"keepBlock": True, **first_block}
    # Kept, each is written after the known keys of its section.
    kept = copy.deepcopy(config)
    kept["config_version"] = 1
    kept["barcode_struct_r1"]["note"] = "plate 3"
    kept["barcode_struct_r1"]["blocks"][0]["keepBlock"] = True
    result = core.run(
        """(core, arg) => {
            const opened = core.config_from_json(arg);
            opened.model.keep_unknown_keys = true;
            const kept = core.config_to_json(opened.model);
            opened.model.keep_unknown_keys = false;
            return {kept, left_out: core.config_to_json(opened.model), unknown_keys: opened.unknown_keys, notices: opened.notices};
        }""",
        json.dumps(in_file),
    )
    assert result["kept"] == canonical_text(kept)
    assert result["left_out"] == canonical_text(config)
    expected_keys = ['"config_version"', '"note" in Read 1', '"keepBlock" in block cb1']
    assert result["unknown_keys"] == expected_keys
    (notice,) = [n for n in result["notices"] if n["level"] == "unknown_keys"]
    assert notice["keys"] == expected_keys
    for key in expected_keys:
        assert key in notice["text"]


def test_kept_keys_named_like_numbers_or_proto_are_still_written_last(core):
    """JavaScript lists a key such as "2" first and can't assign "__proto__"; json.dumps writes both where they were added."""
    config = without_names(read_config(EXAMPLES_DIR / "se_mini_config.json"))
    first_block = config["barcode_struct_r1"]["blocks"][0]
    in_file = {
        "2": "plate",
        "__proto__": "x",
        "constructor": 5,
        **copy.deepcopy(config),
    }
    in_file["barcode_struct_r1"]["blocks"][0] = {"0": "x", **first_block}
    kept = copy.deepcopy(config)
    kept["barcode_struct_r1"]["blocks"][0]["0"] = "x"
    kept.update({"2": "plate", "__proto__": "x", "constructor": 5})
    result = core.run(
        """(core, arg) => {
            const opened = core.config_from_json(arg);
            return {text: core.config_to_json(opened.model), unknown_keys: opened.unknown_keys};
        }""",
        json.dumps(in_file),
    )
    assert result["unknown_keys"] == [
        '"2"',
        '"__proto__"',
        '"constructor"',
        '"0" in block cb1',
    ]
    assert result["text"] == canonical_text(kept)


def test_the_page_writes_for_the_bcwithqc_in_this_repo(core):
    assert core.run("(core) => core.BCWITHQC_VERSION") == bcwithqc_version(), (
        "bcwithqc's version changed: set BCWITHQC_VERSION in "
        "gui/barcode_configurator.html to match bcwithqc/__init__.py"
    )


def test_every_download_starts_with_the_bcwithqc_version(core):
    """D15: the first key records the bcwithqc version the config was made for."""
    config = without_names(read_config(EXAMPLES_DIR / "se_mini_config.json"))
    text = open_and_export(core, json.dumps(config))["text"]
    assert text.splitlines()[:2] == [
        "{",
        f'    "bcwithqc_version": "{bcwithqc_version()}",',
    ]
    empty = core.run("(core) => core.config_to_json(core.new_model())")
    assert next(iter(json.loads(empty))) == "bcwithqc_version"


@pytest.mark.parametrize(
    ("version", "shown"), [("0.3.0", "0.3.0"), (0.4, "0.4"), (None, "null")]
)
def test_an_opened_version_is_reported_and_replaced_by_the_pages_own(
    core, version, shown
):
    config = without_names(read_config(EXAMPLES_DIR / "se_mini_config.json"))
    result = open_and_export(core, json.dumps({"bcwithqc_version": version, **config}))
    assert result["text"] == canonical_text(config)
    expected = (
        f"This file was made for bcwithqc {shown}. The file you download will say "
        f"{bcwithqc_version()}, the version this page writes for."
    )
    assert result["notices"] == [{"level": "info", "text": expected}]


def test_a_downloaded_config_reopens_unchanged_and_reports_its_version(core):
    config = without_names(read_config(EXAMPLES_DIR / "pe_mini_config.json"))
    first = open_and_export(core, json.dumps(config))
    assert first["notices"] == [], "a file without bcwithqc_version gets no notice"
    second = open_and_export(core, first["text"])
    assert second["text"] == first["text"]
    expected = (
        f"This file was made for bcwithqc {bcwithqc_version()}, "
        "the version this page writes for."
    )
    assert second["notices"] == [{"level": "info", "text": expected}]


def test_lowercase_sequences_are_written_in_uppercase(core):
    config = without_names(read_config(EXAMPLES_DIR / "se_mini_config.json"))
    lowercase = copy.deepcopy(config)
    lowercase["barcode_struct_r1"]["blocks"][0]["sequence"] = [
        item.lower() for item in config["barcode_struct_r1"]["blocks"][0]["sequence"]
    ]
    result = open_and_export(core, json.dumps(lowercase))
    assert result["text"] == canonical_text(config)
    assert any(
        "Lowercase sequences in cb1" in notice["text"] for notice in result["notices"]
    )


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("{not json", "not valid JSON"),
        ("[1, 2]", "not a bcwithqc config"),
        ('{"single_end_reads": true}', "not a bcwithqc config"),
    ],
)
def test_files_that_are_not_configs_give_a_clear_error(core, text, message):
    result = core.run(OPEN_AND_EXPORT, {"text": text})
    assert message in result["error"]
