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
import re
from pathlib import Path

import pytest
from gui_support import (
    EXAMPLES_DIR,
    FIXTURES_DIR,
    GUI_PATH,
    OLD_GUI_BARCODES,
    OLD_GUI_SAMPLE_BARCODES,
    READ_KEYS,
    REPO_ROOT,
    SAMPLE_BARCODES,
    barcode_list,
    canonical_text,
    constant_region,
    design_doc_example,
    example_config_names,
    generated_with,
    make_config,
    old_gui_config,
    random_block,
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
    # The examples carry no generated_with, so no notice may mention one.
    assert not [
        n
        for n in result["notices"]
        if "made for" in n["text"] or "generated_with" in n["text"]
    ]


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


def test_legacy_barcode_read_that_keeps_nothing_gets_the_keeps_nothing_notice(core):
    result = open_and_export(
        core,
        read_bytes_text(FIXTURES_DIR / "legacy_sdrranger_gDNA.json.gz"),
        legacy_mode="paired",
    )
    assert any(
        "use only the other read's file" in notice["text"]
        for notice in result["notices"]
    )


LEGACY_FIXTURES = (
    "legacy_master_cDNA.json",
    "legacy_master_gDNA.json.gz",
    "legacy_sdrranger_gDNA.json.gz",
)


def without_key(block: dict, key: str) -> dict:
    """Return a copy of a block without one of its keys."""
    return {name: value for name, value in block.items() if name != key}


def as_legacy(config: dict) -> dict:
    """Return a copy of a config without "single_end_reads", as older tools wrote it."""
    return without_key(config, "single_end_reads")


def config_text(source) -> str:
    """Return the text of a config file (a Path) or of a config built here (a dict)."""
    return read_bytes_text(source) if isinstance(source, Path) else json.dumps(source)


# (file or config, a notice it must get on opening), with the page's wording.
OPENING_NOTICES = [
    (
        FIXTURES_DIR / "legacy_master_gDNA.json.gz",
        "This config is in an older format with no Read 2. It opens as paired-end with the barcodes on Read 1, and all of Read 2 is kept in the output reads, as the old pipeline did.",
    ),
    (
        FIXTURES_DIR / "legacy_master_cDNA.json",
        '"In output reads" was not set in this file, so every region is "Leave out".',
    ),
    (
        FIXTURES_DIR / "legacy_sdrranger_gDNA.json.gz",
        "This config keeps nothing from Read 1, the barcode read, in the output reads. As with the old pipeline, use only the other read's file.",
    ),
    (
        make_config(
            [without_key(barcode_list("unused"), "blockname"), constant_region("con1")]
        ),
        "A region without a name was named bc1.",
    ),
    (
        as_legacy(
            make_config([barcode_list("cell barcode"), constant_region("con1")], [])
        ),
        'Spaces in region names became "_": cell barcode → cell_barcode.',
    ),
    (
        make_config([barcode_list("cb1"), "not a region", constant_region("con1")]),
        "Item 2 in Read 1 is not a region and was left out.",
    ),
    (
        dict(
            make_config([barcode_list("cb1"), constant_region("con1")], []),
            barcode_struct_r2={"keep_nonbarcode": True},
        ),
        'Read 2 has no list of regions ("blocks").',
    ),
    (
        make_config(
            [dict(barcode_list("cb1"), colour="green"), constant_region("con1")]
        ),
        'bcwithqc doesn\'t use "colour" in region cb1.',
    ),
]


@pytest.mark.parametrize(
    ("source", "notice"),
    OPENING_NOTICES,
    ids=[
        "legacy_without_read_2",
        "no_keepblock",
        "legacy_barcode_read_outputs_nothing",
        "unnamed_region",
        "spaces_in_names",
        "item_not_a_region",
        "no_blocks_list",
        "unknown_key_in_a_region",
    ],
)
def test_opening_notices_speak_of_regions(core, source, notice):
    result = open_and_export(core, config_text(source), legacy_mode="paired")
    assert notice in [item["text"] for item in result["notices"]]


def test_no_notice_says_block_outside_a_quoted_key(core):
    """Notices call bcwithqc's blocks "regions"; the word stays only in quoted keys such as "blocks"."""
    sources = [EXAMPLES_DIR / name for name in example_config_names()]
    sources += [FIXTURES_DIR / name for name in LEGACY_FIXTURES]
    sources += [old_gui_config()] + [source for source, _ in OPENING_NOTICES]
    texts = []
    for source in sources:
        for legacy_mode in ("paired", "single"):
            result = open_and_export(core, config_text(source), legacy_mode=legacy_mode)
            texts += [notice["text"] for notice in result["notices"]]
    assert texts
    saying_block = [
        text for text in texts if "block" in re.sub(r'"[^"]*"', "", text).lower()
    ]
    assert saying_block == []


LIST_AND_CONSTANT = [barcode_list("cb1"), constant_region("con1")]


@pytest.mark.parametrize(
    ("config", "read_key", "notices"),
    [
        (
            make_config(LIST_AND_CONSTANT, [], keep_r2=False),
            "barcode_struct_r2",
            [
                'Read 2 has no regions, so all of it is kept in the output reads. The file had "keep_nonbarcode": false; the page writes true, which bcwithqc does for a read without regions anyway.'
            ],
        ),
        (
            dict(make_config(LIST_AND_CONSTANT, []), barcode_struct_r2={"blocks": []}),
            "barcode_struct_r2",
            [
                'Read 2 has no regions, so all of it is kept in the output reads. The file had no "keep_nonbarcode"; the page writes true, which bcwithqc does for a read without regions anyway.'
            ],
        ),
        (
            make_config(LIST_AND_CONSTANT, [], keep_r2=None),
            "barcode_struct_r2",
            [
                'Read 2 has no regions, so all of it is kept in the output reads. The file had "keep_nonbarcode": null; the page writes true, which bcwithqc does for a read without regions anyway.'
            ],
        ),
        (make_config(LIST_AND_CONSTANT, [], keep_r2=True), "barcode_struct_r2", []),
        (
            make_config([], LIST_AND_CONSTANT, keep_r1=False),
            "barcode_struct_r1",
            [
                'Read 1 has no regions, so all of it is kept in the output reads. The file had "keep_nonbarcode": false; the page writes true, which bcwithqc does for a read without regions anyway.'
            ],
        ),
    ],
    ids=["false", "missing", "null", "true", "read_1_false"],
)
def test_a_read_without_regions_opens_as_output_whole(core, config, read_key, notices):
    """The notice says what the file had; a file that already says true gets none."""
    result = open_and_export(core, json.dumps(config))
    assert json.loads(result["text"])[read_key]["keep_nonbarcode"] is True
    assert result["notices"] == [{"level": "info", "text": text} for text in notices]


@pytest.mark.parametrize("read_number", [1, 2])
def test_a_read_with_regions_and_no_keep_nonbarcode_leaves_its_rest_out(
    core, read_number
):
    """Both reads say true in the file until the key is removed, so false comes from its absence."""
    read_key = f"barcode_struct_r{read_number}"
    config = make_config(LIST_AND_CONSTANT, LIST_AND_CONSTANT, keep_r1=True)
    config[read_key] = without_key(config[read_key], "keep_nonbarcode")
    result = open_and_export(core, json.dumps(config))
    assert json.loads(result["text"])[read_key]["keep_nonbarcode"] is False
    assert result["notices"] == [
        {
            "level": "warning",
            "text": f'Read {read_number} had no "keep_nonbarcode", so the rest of Read {read_number} is "Leave out".',
        }
    ]


# Every region is set to "Keep", so a region that loses "keepblock" is told apart.
ALL_OUTPUT_ON_TWO_READS = make_config(
    [barcode_list("cb1", keep=True), constant_region("con1", keep=True)],
    [barcode_list("bc2", keep=True), constant_region("con2", keep=True)],
)


def without_keepblock(config: dict, names: set[str]) -> dict:
    """Return a copy of a config in which the regions with these names have no "keepblock"."""
    result = copy.deepcopy(config)
    for key in READ_KEYS:
        result[key]["blocks"] = [
            without_key(block, "keepblock") if block["blockname"] in names else block
            for block in result[key]["blocks"]
        ]
    return result


@pytest.mark.parametrize(
    ("names", "notice"),
    [
        (
            {"cb1", "con1", "bc2", "con2"},
            '"In output reads" was not set in this file, so every region is "Leave out".',
        ),
        (
            {"cb1", "con1"},
            '"In output reads" was not set for cb1, con1, so they are "Leave out".',
        ),
        (
            {"bc2"},
            '"In output reads" was not set for bc2, so it is "Leave out".',
        ),
    ],
    ids=["every_region", "every_region_of_read_1", "one_region"],
)
def test_regions_without_keepblock_are_left_out_and_the_notice_says_which(
    core, names, notice
):
    """The notice says "every region" only when no region on either read has "keepblock"."""
    result = open_and_export(
        core, json.dumps(without_keepblock(ALL_OUTPUT_ON_TWO_READS, names))
    )
    assert result["notices"] == [{"level": "info", "text": notice}]
    exported = json.loads(result["text"])
    keepblock = {
        block["blockname"]: block["keepblock"]
        for key in READ_KEYS
        for block in exported[key]["blocks"]
    }
    assert keepblock == {
        name: name not in names for name in ("cb1", "con1", "bc2", "con2")
    }


@pytest.mark.parametrize(
    ("config", "place"),
    [
        (
            make_config(
                [
                    barcode_list("cb1"),
                    dict(
                        without_key(constant_region("con1"), "blockname"),
                        colour="green",
                    ),
                ]
            ),
            '"colour" in region 2 of Read 1',
        ),
        (
            make_config(
                LIST_AND_CONSTANT,
                [dict(without_key(random_block("umi2"), "blockname"), colour="green")],
            ),
            '"colour" in region 1 of Read 2',
        ),
        (
            make_config(
                [dict(barcode_list(""), colour="green"), constant_region("con1")]
            ),
            '"colour" in region 1 of Read 1',
        ),
    ],
    ids=["read_1", "read_2", "empty_name"],
)
def test_an_unknown_key_in_a_region_without_a_name_gives_the_regions_place(
    core, config, place
):
    result = open_and_export(core, json.dumps(config))
    (notice,) = [n for n in result["notices"] if n["level"] == "unknown_keys"]
    assert notice["keys"] == [place]
    assert notice["text"] == f"bcwithqc doesn't use {place}."


def test_a_region_of_unknown_type_without_a_name_is_named_region1(core):
    unknown_type = without_key(
        dict(constant_region("con1"), blocktype="fancyBarcode"), "blockname"
    )
    result = open_and_export(
        core, json.dumps(make_config([barcode_list("cb1"), unknown_type]))
    )
    assert result["notices"] == [
        {"level": "info", "text": "A region without a name was named region1."}
    ]
    exported = json.loads(result["text"])["barcode_struct_r1"]["blocks"]
    assert [block["blockname"] for block in exported] == ["cb1", "region1"]


def test_the_lowercase_notice_names_a_region_after_spaces_became_underscores(core):
    """A notice names a region as the page writes it: "cell barcode" in an older file is cell_barcode."""
    lowercase = [item.lower() for item in SAMPLE_BARCODES]
    config = as_legacy(
        make_config(
            [barcode_list("cell barcode", lowercase), constant_region("con1")], []
        )
    )
    texts = [
        notice["text"]
        for notice in open_and_export(core, json.dumps(config))["notices"]
    ]
    assert "Lowercase sequences in cell_barcode were converted to uppercase." in texts


@pytest.mark.parametrize(
    ("legacy_mode", "notice"),
    [
        (
            "single",
            'This config keeps nothing from Read 1 in the output reads. With single-end reads they would be empty: set the rest of Read 1, or a region, to "Keep" if you need them.',
        ),
        (
            "paired",
            "This config keeps nothing from Read 1, the barcode read, in the output reads. As with the old pipeline, use only the other read's file.",
        ),
    ],
)
def test_a_legacy_read_1_that_outputs_nothing_gets_the_notice_for_its_layout(
    core, legacy_mode, notice
):
    """Paired-end leaves the other read's file to use; single-end leaves empty output reads."""
    result = open_and_export(
        core,
        read_bytes_text(FIXTURES_DIR / "legacy_sdrranger_gDNA.json.gz"),
        legacy_mode=legacy_mode,
    )
    outputs_nothing = [
        n
        for n in result["notices"]
        if n["text"].startswith("This config keeps nothing")
    ]
    assert outputs_nothing == [{"level": "info", "text": notice}]


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


def test_new_model_keeps_the_rest_of_both_reads(core):
    """Read 1 and Read 2 start with the same default: "Keep" for the rest of the read."""
    reads = core.run("(core) => core.new_model().reads")
    assert [read["keep_nonbarcode"] for read in reads] == [True, True]
    assert [read["blocks"] for read in reads] == [[], []]


def test_paired_end_model_with_one_list_exports_the_rest_of_both_reads(core):
    exported = json.loads(
        core.run(
            """(core) => {
                const model = core.new_model();
                model.single_end = false;
                return core.config_to_json(core.add_block(model, 0, 'barcodeList').model);
            }"""
        )
    )
    assert exported["single_end_reads"] is False
    assert exported["barcode_struct_r1"]["keep_nonbarcode"] is True
    (region,) = exported["barcode_struct_r1"]["blocks"]
    assert region["blocktype"] == "barcodeList"
    assert region["keepblock"] is False
    assert exported["barcode_struct_r2"] == {"keep_nonbarcode": True, "blocks": []}


@pytest.mark.parametrize(
    "blocktype", ["barcodeList", "constantRegion", "randomBarcode"]
)
def test_new_regions_are_left_out_of_the_output_reads(core, blocktype):
    """A new region starts at "Leave out"; only the rest of the read is kept by default."""
    region = core.run(
        "(core, blocktype) => core.add_block(core.new_model(), 0, blocktype).model.reads[0].blocks[0]",
        blocktype,
    )
    assert region["blocktype"] == blocktype
    assert region["keepblock"] is False


def expanded_templates(core):
    return core.run("(core) => window.barcode_templates")


def template_source_config(source):
    if source.endswith(".md"):
        return design_doc_example()
    return read_config(REPO_ROOT / source)


@pytest.mark.parametrize(
    "template_id", ["sdr_seq_rna", "sdr_seq_gdna", "simple_example"]
)
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
    assert (
        "simple_example",
        "Simple example: barcode, linker, UMI",
        "doc/config_format.md",
    ) in templates


def test_page_lists_the_templates_in_order(core):
    """The Templates menu shows these after "New empty config"; "Design-doc example" is gone."""
    listed = [
        (template["id"], template["title"], template["source"])
        for template in expanded_templates(core)
    ]
    assert listed == [
        ("sdr_seq_rna", "SDR-seq RNA", "examples/cDNA_paired_end_config.json"),
        ("sdr_seq_gdna", "SDR-seq gDNA", "examples/gDNA_paired_end_config.json"),
        (
            "simple_example",
            "Simple example: barcode, linker, UMI",
            "doc/config_format.md",
        ),
    ]


def load_update_script():
    spec = importlib.util.spec_from_file_location(
        "update_page", REPO_ROOT / "gui" / "update_page.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_page_template_data_is_current():
    """Fails when a template source changed; run gui/update_page.py to refresh the page."""
    page = GUI_PATH.read_text(encoding="utf-8")
    match = load_update_script().DATA_BLOCK.search(page)
    assert match is not None
    assert json.loads(match.group(2)) == load_update_script().template_data()


def test_page_stamp_is_current():
    """Fails after a version bump or a rename; run gui/update_page.py to restamp the page."""
    script = load_update_script()
    match = script.STAMP_LINE.search(GUI_PATH.read_text(encoding="utf-8"))
    assert match is not None
    assert match.group(0) == script.stamp_line(script.generated_with()), (
        "pyproject.toml's name or version changed: run python gui/update_page.py"
    )


def test_update_script_restamps_a_stale_page(tmp_path, monkeypatch):
    script = load_update_script()
    page = GUI_PATH.read_text(encoding="utf-8")
    old_stamp = script.stamp_line({"name": "old_name", "version": "0.0.1"})
    stale = script.STAMP_LINE.sub(lambda match: old_stamp, page)
    assert stale != page
    stale_copy = tmp_path / GUI_PATH.name
    stale_copy.write_text(stale, encoding="utf-8")
    monkeypatch.setattr(script, "PAGE", stale_copy)
    script.main()
    assert stale_copy.read_text(encoding="utf-8") == page


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
    expected_keys = [
        '"config_version" at the top level',
        '"note" in Read 1',
        '"keepBlock" in region cb1',
    ]
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
        '"2" at the top level',
        '"__proto__" at the top level',
        '"constructor" at the top level',
        '"0" in region cb1',
    ]
    assert result["text"] == canonical_text(kept)


def test_the_page_stamps_the_name_and_version_in_pyproject(core):
    assert core.run("(core) => core.GENERATED_WITH") == generated_with(), (
        "pyproject.toml's name or version changed: run python gui/update_page.py"
    )


def test_every_download_starts_with_generated_with(core):
    """D15: the first key records the tool and version the config was made for."""
    config = without_names(read_config(EXAMPLES_DIR / "se_mini_config.json"))
    text = open_and_export(core, json.dumps(config))["text"]
    stamp = generated_with()
    assert text.splitlines()[:5] == [
        "{",
        '    "generated_with": {',
        f'        "name": "{stamp["name"]}",',
        f'        "version": "{stamp["version"]}"',
        "    },",
    ]
    empty = core.run("(core) => core.config_to_json(core.new_model())")
    assert next(iter(json.loads(empty))) == "generated_with"


@pytest.mark.parametrize(
    ("written", "first_sentence"),
    [
        (
            {"name": "bcwithqc", "version": "0.3.0"},
            "This file was made for bcwithqc 0.3.0.",
        ),
        (
            {"name": "renamed_tool", "version": "1.0"},
            "This file was made for renamed_tool 1.0.",
        ),
        (
            "0.4.1",
            'This file\'s "generated_with" ("0.4.1") doesn\'t give a name and a version.',
        ),
        (
            None,
            "This file's \"generated_with\" (null) doesn't give a name and a version.",
        ),
    ],
)
def test_an_opened_stamp_is_reported_and_replaced_by_the_pages_own(
    core, written, first_sentence
):
    config = without_names(read_config(EXAMPLES_DIR / "se_mini_config.json"))
    result = open_and_export(core, json.dumps({"generated_with": written, **config}))
    assert result["text"] == canonical_text(config)
    stamp = generated_with()
    expected = (
        f"{first_sentence} The file you download will say "
        f"{stamp['name']} {stamp['version']}, the version this page writes for."
    )
    assert result["notices"] == [{"level": "info", "text": expected}]


def test_a_downloaded_config_reopens_unchanged_and_reports_its_stamp(core):
    config = without_names(read_config(EXAMPLES_DIR / "pe_mini_config.json"))
    first = open_and_export(core, json.dumps(config))
    assert first["notices"] == [], "a file without generated_with gets no notice"
    second = open_and_export(core, first["text"])
    assert second["text"] == first["text"]
    stamp = generated_with()
    expected = (
        f"This file was made for {stamp['name']} {stamp['version']}, "
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
