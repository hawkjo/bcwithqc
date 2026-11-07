"""Tier C: importing barcode tables from text.

The importer takes the right column and asks whenever the choice is
ambiguous. Spreadsheets are covered in test_gui_spreadsheets.py.
"""

import pytest
from gui_support import SAMPLE_BARCODES, barcode_list, make_config

ZERO_WIDTH_SPACE = chr(0x200B)
PLAN = """(core, arg) => {
    const parsed = core.parse_delimited_text(arg.text);
    const plan = core.plan_table_import(parsed.rows, arg.name || 'list.csv', arg.choices || {});
    return {delimiter: parsed.delimiter, plan};
}"""


def plan_import(core, text, **choices):
    return core.run(PLAN, {"text": text, "choices": choices})


def table(rows, delimiter):
    return "\n".join(delimiter.join(row) for row in rows) + "\n"


NAMED_ROWS = [["name", "barcode"]] + [
    [f"bc{i + 1}", seq] for i, seq in enumerate(SAMPLE_BARCODES)
]


@pytest.mark.parametrize(
    ("delimiter", "expected"),
    [
        (",", ","),
        (";", ";"),
        ("\t", "\t"),
        (":", ":"),
        (" ", "whitespace"),
        ("   ", "whitespace"),
    ],
)
def test_every_delimiter_gives_the_same_sequences(core, delimiter, expected):
    result = plan_import(core, table(NAMED_ROWS, delimiter))
    assert result["delimiter"] == expected
    assert result["plan"]["status"] == "ok"
    assert result["plan"]["sequences"] == SAMPLE_BARCODES


def test_single_column_without_header_is_taken_directly(core):
    plan = plan_import(core, "\n".join(SAMPLE_BARCODES))["plan"]
    assert plan["status"] == "ok"
    assert plan["header"] is False
    assert plan["sequences"] == SAMPLE_BARCODES


def test_two_dna_columns_make_the_user_choose(core):
    rows = [["CCTTAACAT", "ATGTTAAGG"], ["AGCGTAGAA", "TTCTACGCT"]]
    plan = plan_import(core, table(rows, ","))["plan"]
    assert plan == {"status": "choose_column", "candidates": [0, 1]}
    chosen = plan_import(core, table(rows, ","), column=1)["plan"]
    assert chosen["sequences"] == ["ATGTTAAGG", "TTCTACGCT"]


def test_tag_header_over_named_columns_is_a_header(core):
    rows = [["name", "TAG"]] + [
        [f"bc{i}", seq] for i, seq in enumerate(SAMPLE_BARCODES)
    ]
    plan = plan_import(core, table(rows, ","))["plan"]
    assert plan["header"] is True
    assert plan["sequences"] == SAMPLE_BARCODES


def test_tag_in_a_single_column_asks_whether_it_is_a_header(core):
    text = "TAG\n" + "\n".join(SAMPLE_BARCODES)
    assert plan_import(core, text)["plan"] == {
        "status": "ask_header",
        "column": 0,
        "first": "TAG",
    }
    assert plan_import(core, text, header=True)["plan"]["sequences"] == SAMPLE_BARCODES
    assert (
        plan_import(core, text, header=False)["plan"]["sequences"]
        == ["TAG"] + SAMPLE_BARCODES
    )


def test_first_entry_with_a_length_no_other_entry_has_asks_in_a_single_column(core):
    text = "ACGTA\n" + "\n".join(SAMPLE_BARCODES)
    assert plan_import(core, text)["plan"]["status"] == "ask_header"


def test_text_header_in_a_single_column_is_dropped_without_asking(core):
    plan = plan_import(core, "barcode\n" + "\n".join(SAMPLE_BARCODES))["plan"]
    assert plan["header"] is True
    assert plan["sequences"] == SAMPLE_BARCODES


def test_extra_columns_are_ignored(core):
    rows = [["id", "barcode", "count", "note"]] + [
        [str(i), seq, str(10 * i), "plate A"] for i, seq in enumerate(SAMPLE_BARCODES)
    ]
    plan = plan_import(core, table(rows, "\t"))["plan"]
    assert plan["column"] == 1
    assert plan["sequences"] == SAMPLE_BARCODES


def test_lowercase_crlf_bom_blank_lines_and_quotes_are_handled(core):
    lines = ['"name","barcode"'] + [
        f'"bc{i}","{seq.lower()}"' for i, seq in enumerate(SAMPLE_BARCODES)
    ]
    lines.insert(2, "")
    text = chr(0xFEFF) + "\r\n".join(lines) + "\r\n"
    plan = plan_import(core, text)["plan"]
    assert plan["status"] == "ok"
    assert plan["sequences"] == SAMPLE_BARCODES
    # Line numbers refer to the source, counting the blank line.
    assert [source["line"] for source in plan["sources"]] == [2, 4, 5, 6]


def test_quoted_field_may_contain_the_delimiter(core):
    text = 'name,barcode\n"bc, first",CCTTAACAT\n"bc ""two""",AGCGTAGAA\n'
    parsed = core.call("parse_delimited_text", text)
    assert parsed["rows"][1]["cells"] == ["bc, first", "CCTTAACAT"]
    assert parsed["rows"][2]["cells"] == ['bc "two"', "AGCGTAGAA"]


def test_zero_width_space_is_imported_and_flagged_with_its_source_line(core):
    damaged = SAMPLE_BARCODES[2] + ZERO_WIDTH_SPACE
    lines = [
        "barcode",
        SAMPLE_BARCODES[0],
        SAMPLE_BARCODES[1],
        damaged,
        SAMPLE_BARCODES[3],
    ]
    result = core.run(
        """(core, arg) => {
            const plan = core.plan_table_import(core.parse_delimited_text(arg.text).rows, 'list.csv');
            const opened = core.config_from_json(arg.config).model;
            const block_id = opened.reads[0].blocks[0].id;
            const model = core.apply_import(opened, 0, block_id, plan, 'replace');
            return {plan, issues: core.validate_config(model)};
        }""",
        {"text": "\n".join(lines), "config": make_config([barcode_list("cb1", [])])},
    )
    assert result["plan"]["sequences"][2] == damaged
    e03 = [issue for issue in result["issues"] if issue["rule"] == "E03"]
    assert len(e03) == 1
    assert "line 4 of list.csv" in e03[0]["message"]
    assert "U+200B" in e03[0]["message"]
    assert e03[0]["rows"] == [2]


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("", "the input is empty"),
        ("\n\n  \n", "the input is empty"),
        ("barcode\n", "only a header row"),
        ("name,barcode\n", "only a header row"),
        ("name,count\nbc1,4\nbc2,5\n", "No column of DNA sequences"),
    ],
)
def test_empty_header_only_and_non_dna_input_give_a_clear_message(core, text, message):
    plan = plan_import(core, text)["plan"]
    assert plan["status"] == "error"
    assert message in plan["message"]


APPLY = """(core, arg) => {
    const model = core.config_from_json(arg.config).model;
    const block_id = model.reads[0].blocks[0].id;
    const result = core.apply_import(model, 0, block_id, arg.imported, arg.mode);
    return result.reads[0].blocks[0];
}"""


@pytest.mark.parametrize(
    ("mode", "expected"),
    [
        ("append", ["AAAACCCC", "GGGGTTTT", "ACGTACGT"]),
        ("replace", ["ACGTACGT"]),
    ],
)
def test_import_appends_or_replaces(core, mode, expected):
    config = make_config([barcode_list("cb1", ["AAAACCCC", "GGGGTTTT"])])
    imported = {"sequences": ["ACGTACGT"], "sources": [{"name": "new.csv", "line": 3}]}
    block = core.run(APPLY, {"config": config, "imported": imported, "mode": mode})
    assert block["sequence"] == expected
    assert block["sequence_source"][-1] == {"name": "new.csv", "line": 3}
    assert len(block["sequence_source"]) == len(expected)
