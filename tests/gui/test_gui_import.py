"""Tier C: importing barcode tables from text.

The importer takes the right column and asks whenever the choice is
ambiguous. Spreadsheets are covered in test_gui_spreadsheets.py. A paste of
several sequences into a list's add box is imported the same way: several
lines as a table, one line of sequences as a list of items. Spaces separate
the items of one line only for a barcode list: a constant region's sequence
is often copied in groups of ten. For a constant region, whose variants often
differ in length, a first entry's length alone never suggests a header.
"""

import pytest
from gui_support import (
    EXAMPLES_DIR,
    SAMPLE_BARCODES,
    barcode_list,
    constant_region,
    make_config,
    read_config,
)

ZERO_WIDTH_SPACE = chr(0x200B)
BYTE_ORDER_MARK = chr(0xFEFF)
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


@pytest.mark.parametrize(
    "text",
    ["CCTTAACAT\nCCTTANCAT\nAGCGTAGAA\n", "name,count\nbc1,4\nbc2,5\n"],
    ids=["n_in_a_barcode", "no_dna"],
)
def test_no_dna_column_message_says_which_letters_count(core, text):
    assert plan_import(core, text)["plan"] == {
        "status": "error",
        "message": "No column of DNA sequences (only A, C, G and T) was found.",
    }


def test_a_length_no_other_entry_has_is_suspicious_only_when_judged_by_length(core):
    """ACGTA is the only 5 nt entry: by length (the default) the user is asked; otherwise it is no header."""
    rows = item_rows(["ACGTA"] + SAMPLE_BARCODES)
    assert core.call("detect_header_row", rows, 0) == "ask"
    assert core.call("detect_header_row", rows, 0, True) == "ask"
    assert core.call("detect_header_row", rows, 0, False) is False


@pytest.mark.parametrize("by_length", [True, False])
def test_tag_first_in_a_single_column_asks_either_way(core, by_length):
    rows = item_rows(["TAG"] + SAMPLE_BARCODES)
    assert core.call("detect_header_row", rows, 0, by_length) == "ask"


def sdr_seq_linker_variants() -> list[str]:
    """Return the variants of SDR-seq RNA's constant region cso, in the example config's order."""
    blocks = read_config(EXAMPLES_DIR / "cDNA_paired_end_config.json")[
        "barcode_struct_r1"
    ]["blocks"]
    (cso,) = [block for block in blocks if block["blockname"] == "cso"]
    return cso["sequence"]


PLAN_FOR_BLOCKTYPE = """(core, arg) => core.plan_table_import(
    core.parse_delimited_text(arg.text).rows, 'variants.txt', {}, 'line', arg.blocktype
)"""


def test_variants_of_different_lengths_import_into_a_constant_region_without_a_question(
    core,
):
    """SDR-seq's cso variants are 17, 16, 15 and 14 nt, so the first one's length is its own.

    A constant region takes all four; a barcode list, the default, asks
    whether the first is a header.
    """
    variants = sdr_seq_linker_variants()
    assert [len(variant) for variant in variants] == [17, 16, 15, 14]
    text = "\n".join(variants)
    plan = core.run(PLAN_FOR_BLOCKTYPE, {"text": text, "blocktype": "constantRegion"})
    assert plan["status"] == "ok"
    assert plan["sequences"] == variants
    assert plan["header"] is False
    assert plan_import(core, text)["plan"] == {
        "status": "ask_header",
        "column": 0,
        "first": variants[0],
    }


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


PASTE_AND_PLAN = """(core, text) => {
    const pasted = core.pasted_list_rows(text);
    return {pasted, plan: core.plan_table_import(pasted.rows, 'pasted text', {}, pasted.unit)};
}"""
PASTE_INTO_FIRST_LIST = """(core, arg) => {
    const pasted = core.pasted_list_rows(arg.text);
    const plan = core.plan_table_import(pasted.rows, 'pasted text', {}, pasted.unit);
    const opened = core.config_from_json(arg.config).model;
    const model = core.apply_import(opened, 0, opened.reads[0].blocks[0].id, plan, 'append');
    return {block: model.reads[0].blocks[0], issues: core.validate_config(model)};
}"""


def item_rows(entries: list[str]) -> list[dict]:
    """Return the rows one pasted line of entries should give: one cell each, numbered from 1."""
    return [
        {"line": number, "cells": [entry]} for number, entry in enumerate(entries, 1)
    ]


def item_sources(count: int) -> list[dict]:
    """Return the sources of count pasted items, as messages name them ("item k of pasted text")."""
    return [
        {"name": "pasted text", "line": number, "unit": "item"}
        for number in range(1, count + 1)
    ]


@pytest.mark.parametrize(
    "text",
    [
        "",
        "  \t ",
        "\n\r\n",
        "CCTTAACAT",
        "CCTTAACAT\n",
        "CCTTAACAT\r\n",
        "\n  CCTTAACAT \n\n",
        BYTE_ORDER_MARK + "CCTTAACAT",
        "barcode",
    ],
    ids=[
        "empty",
        "whitespace",
        "blank_lines",
        "one_sequence",
        "trailing_newline",
        "trailing_crlf",
        "blank_lines_around",
        "byte_order_mark",
        "one_word",
    ],
)
def test_nothing_or_a_single_value_is_not_a_pasted_list(core, text):
    """None leaves the paste to the browser, which puts it in the add box as typed."""
    assert core.call("pasted_list_rows", text) is None


@pytest.mark.parametrize(
    ("text", "rows"),
    [
        (
            "CCTTAACAT\n\nAGCGTAGAA\r\nTATAGGTGT\n",
            [
                {"line": 1, "cells": ["CCTTAACAT"]},
                {"line": 3, "cells": ["AGCGTAGAA"]},
                {"line": 4, "cells": ["TATAGGTGT"]},
            ],
        ),
        (
            "name,barcode\nbc1,CCTTAACAT\nbc2,AGCGTAGAA\n",
            [
                {"line": 1, "cells": ["name", "barcode"]},
                {"line": 2, "cells": ["bc1", "CCTTAACAT"]},
                {"line": 3, "cells": ["bc2", "AGCGTAGAA"]},
            ],
        ),
        (
            BYTE_ORDER_MARK + "ccttaacat\nagcgtagaa",
            [
                {"line": 1, "cells": ["ccttaacat"]},
                {"line": 2, "cells": ["agcgtagaa"]},
            ],
        ),
        # Two lines are a table even without DNA; the import then says what is wrong.
        (
            "hello\nworld",
            [{"line": 1, "cells": ["hello"]}, {"line": 2, "cells": ["world"]}],
        ),
    ],
    ids=["blank_line_and_crlf", "csv_table", "byte_order_mark", "not_dna"],
)
def test_several_lines_are_read_as_a_table(core, text, rows):
    """Line numbers count the source's lines, blank ones included, as for a file."""
    assert core.call("pasted_list_rows", text) == {"rows": rows, "unit": "line"}
    assert core.call("parse_delimited_text", text)["rows"] == rows


@pytest.mark.parametrize(
    "text",
    [
        "CCTTAACAT,AGCGTAGAA,TATAGGTGT",
        "CCTTAACAT;AGCGTAGAA;TATAGGTGT",
        "CCTTAACAT\tAGCGTAGAA\tTATAGGTGT",
        "CCTTAACAT AGCGTAGAA   TATAGGTGT",
        "CCTTAACAT, AGCGTAGAA;\tTATAGGTGT",
        " CCTTAACAT,AGCGTAGAA,TATAGGTGT,\r\n",
        BYTE_ORDER_MARK + "CCTTAACAT AGCGTAGAA TATAGGTGT\n",
    ],
    ids=[
        "commas",
        "semicolons",
        "tabs",
        "spaces",
        "mixed",
        "padded_with_trailing_comma",
        "byte_order_mark",
    ],
)
def test_one_line_of_sequences_is_a_list_of_items(core, text):
    assert core.call("pasted_list_rows", text) == {
        "rows": item_rows(["CCTTAACAT", "AGCGTAGAA", "TATAGGTGT"]),
        "unit": "item",
    }


@pytest.mark.parametrize(
    "text",
    ["CCTTAACAT hello", "bc1,CCTTAACAT", "CCTTAACAT;CCTTANCAT"],
    ids=["word", "name_then_barcode", "n_in_a_barcode"],
)
def test_one_line_with_an_entry_that_is_not_dna_is_not_a_pasted_list(core, text):
    assert core.call("pasted_list_rows", text) is None


SPACED_CONSTANT = "GTCTCGTGGG CTCGGAGATG"


def test_spaces_separate_items_unless_split_on_spaces_is_false(core):
    """A constant region's sequence copied in groups of ten stays one value, left to the browser."""
    two_items = {"rows": item_rows(["GTCTCGTGGG", "CTCGGAGATG"]), "unit": "item"}
    assert core.call("pasted_list_rows", SPACED_CONSTANT) == two_items
    assert core.call("pasted_list_rows", SPACED_CONSTANT, True) == two_items
    assert core.call("pasted_list_rows", SPACED_CONSTANT, False) is None


@pytest.mark.parametrize(
    "text",
    ["AAAA,CCCC", "AAAA\tCCCC", "AAAA; CCCC"],
    ids=["comma", "tab", "semicolon_and_space"],
)
def test_tabs_commas_and_semicolons_separate_items_when_spaces_do_not(core, text):
    assert core.call("pasted_list_rows", text, False) == {
        "rows": item_rows(["AAAA", "CCCC"]),
        "unit": "item",
    }


@pytest.mark.parametrize("split_on_spaces", [True, False])
@pytest.mark.parametrize(
    ("text", "rows"),
    [
        (
            "CCTTAACAT\nAGCGTAGAA",
            [
                {"line": 1, "cells": ["CCTTAACAT"]},
                {"line": 2, "cells": ["AGCGTAGAA"]},
            ],
        ),
        (
            "bc1 CCTTAACAT\nbc2 AGCGTAGAA",
            [
                {"line": 1, "cells": ["bc1", "CCTTAACAT"]},
                {"line": 2, "cells": ["bc2", "AGCGTAGAA"]},
            ],
        ),
    ],
    ids=["one_column", "columns_separated_by_spaces"],
)
def test_several_lines_are_a_table_whatever_split_on_spaces_says(
    core, text, rows, split_on_spaces
):
    assert core.call("pasted_list_rows", text, split_on_spaces) == {
        "rows": rows,
        "unit": "line",
    }


def test_pasted_items_keep_their_case_and_all_import_with_item_sources(core):
    """pasted_list_rows keeps the text as typed; plan_table_import uppercases it.

    The repeated entry is imported too, so the duplicate check can flag it.
    """
    result = core.run(PASTE_AND_PLAN, "ccttaacat, AGCGtagaa;TATAGGTGT ccttaacat")
    assert result["pasted"] == {
        "rows": item_rows(["ccttaacat", "AGCGtagaa", "TATAGGTGT", "ccttaacat"]),
        "unit": "item",
    }
    plan = result["plan"]
    assert plan["status"] == "ok"
    assert plan["sequences"] == ["CCTTAACAT", "AGCGTAGAA", "TATAGGTGT", "CCTTAACAT"]
    assert plan["sources"] == item_sources(4)


def test_pasted_items_are_appended_and_a_duplicate_is_flagged(core):
    result = core.run(
        PASTE_INTO_FIRST_LIST,
        {
            "text": "TATAGGTGT AGCGTAGAA GGCGACTCA",
            "config": make_config(
                [
                    barcode_list("cb1", ["CCTTAACAT", "AGCGTAGAA"]),
                    constant_region("con1"),
                ]
            ),
        },
    )
    block = result["block"]
    assert block["sequence"] == [
        "CCTTAACAT",
        "AGCGTAGAA",
        "TATAGGTGT",
        "AGCGTAGAA",
        "GGCGACTCA",
    ]
    assert block["sequence_source"] == [None, None] + item_sources(3)
    (duplicate,) = [issue for issue in result["issues"] if issue["rule"] == "E04"]
    assert "AGCGTAGAA appears twice" in duplicate["message"]
    assert duplicate["rows"] == [1, 3]


def test_a_problem_in_a_pasted_item_names_the_item(core):
    """An item's place reads "(item k of pasted text)", as a file line's reads "(line k of list.csv)".

    E03 is the message that names an entry's place.
    """
    imported = {
        "sequences": ["CCTTAACAT", "AGCGTAGAA" + ZERO_WIDTH_SPACE],
        "sources": item_sources(2),
    }
    issues = core.run(
        """(core, arg) => {
            const opened = core.config_from_json(arg.config).model;
            const model = core.apply_import(opened, 0, opened.reads[0].blocks[0].id, arg.imported, 'replace');
            return core.validate_config(model);
        }""",
        {
            "config": make_config([barcode_list("cb1", []), constant_region("con1")]),
            "imported": imported,
        },
    )
    (e03,) = [issue for issue in issues if issue["rule"] == "E03"]
    assert "(item 2 of pasted text)" in e03["message"]
    assert e03["level"] == "error"
    assert e03["rows"] == [1]
