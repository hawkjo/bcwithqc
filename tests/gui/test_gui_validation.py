"""Tier D: validation rules.

Each rule must fire exactly at its boundary, and bcwithqc's real configs must
pass with no errors and no warnings. Rule ids follow plan §08. Steps not
done yet (no read layout, no regions, an empty list or constant region) are
to-dos rather than errors, and messages call bcwithqc's blocks "regions".
"""

import ast
from pathlib import Path

import pytest
from gui_support import (
    EXAMPLES_DIR,
    FIXTURES_DIR,
    barcode_list,
    constant_region,
    example_config_names,
    make_config,
    random_block,
    read_bytes_text,
)

VALIDATE = "(core, arg) => core.validate_config(core.config_from_json(arg.config, arg.options || {}).model)"
VALIDATE_WITH_IDS = """(core, arg) => {
    const model = core.config_from_json(arg).model;
    return {
        issues: core.validate_config(model),
        block_ids: model.reads.map((read) => read.blocks.map((block) => block.id)),
    };
}"""
# Barcodes of 9, 8 and 9 nt: entry 1 is the one whose length differs.
MIXED_LENGTHS = ["CCTTAACAT", "AGCGTAGA", "TATAGGTGT"]


def issues_for(core, config, **options):
    return core.run(VALIDATE, {"config": config, "options": options})


def rules(issues):
    return sorted(issue["rule"] for issue in issues)


REAL_CONFIG_TEXTS = [(name, EXAMPLES_DIR / name) for name in example_config_names()] + [
    (name, FIXTURES_DIR / name)
    for name in (
        "legacy_master_cDNA.json",
        "legacy_master_gDNA.json.gz",
        "legacy_sdrranger_gDNA.json.gz",
    )
]


@pytest.mark.parametrize(
    ("name", "path"), REAL_CONFIG_TEXTS, ids=[name for name, _ in REAL_CONFIG_TEXTS]
)
def test_real_configs_have_no_errors_and_no_warnings(core, name, path):
    assert len(REAL_CONFIG_TEXTS) == 11
    issues = core.run(
        "(core, arg) => core.validate_config(core.config_from_json(arg, {legacy_mode: 'paired'}).model)",
        read_bytes_text(path),
    )
    assert issues == []


def test_new_model_needs_a_mode_and_a_block(core):
    """A fresh page has nothing wrong yet: both of its steps are to-dos, not errors."""
    issues = core.run("(core) => core.validate_config(core.new_model())")
    assert rules(issues) == ["E01", "E14"]
    assert {issue["rule"]: (issue["level"], issue["message"]) for issue in issues} == {
        "E01": ("todo", "Choose paired-end or single-end reads."),
        "E14": ("todo", "Add at least one region to Read 1 or Read 2."),
    }


@pytest.mark.parametrize(
    ("config", "message"),
    [
        (make_config([]), "Add at least one region to Read 1."),
        (make_config([], []), "Add at least one region to Read 1 or Read 2."),
        # A single-end config doesn't write Read 2, so its regions don't count.
        (
            make_config(
                [], [barcode_list("bc1"), constant_region("con1")], single_end=True
            ),
            "Add at least one region to Read 1.",
        ),
    ],
    ids=["single_end", "paired_end", "single_end_with_regions_on_hidden_read_2"],
)
def test_e14_is_a_to_do_naming_the_reads_the_file_holds(core, config, message):
    issues = issues_for(core, config)
    assert [(issue["rule"], issue["level"], issue["message"]) for issue in issues] == [
        ("E14", "todo", message)
    ]


def test_one_block_anywhere_satisfies_e14(core):
    config = make_config([], [barcode_list("bc1"), constant_region("con1")])
    assert "E14" not in rules(issues_for(core, config))


@pytest.mark.parametrize(
    ("sequences", "expected_rule", "text"),
    [
        (["CCTTAACAT", ""], "E03", "an entry is empty"),
        (["CCTTAACAT", "CCTTANCAT"], "E03", 'contains "N"'),
        (["CCTTAACAT", "CCTTAACAT" + chr(0x200B)], "E03", "an invisible U+200B"),
        (["CCTTAACAT", "AGCGTAGAA", "CCTTAACAT"], "E04", "CCTTAACAT appears twice"),
    ],
)
def test_sequence_rules(core, sequences, expected_rule, text):
    issues = issues_for(
        core, make_config([barcode_list("cb1", sequences), constant_region("con1")])
    )
    matching = [issue for issue in issues if issue["rule"] == expected_rule]
    assert len(matching) == 1
    assert text in matching[0]["message"]
    # A bad entry is a mistake to fix, not a step still to do.
    assert matching[0]["level"] == "error"


@pytest.mark.parametrize(
    ("blocks", "empty_index", "message"),
    [
        (
            [barcode_list("cb1", []), constant_region("con1")],
            0,
            "Add barcodes to cb1.",
        ),
        (
            [barcode_list("cb1"), constant_region("con1", [])],
            1,
            "Add the sequence of con1.",
        ),
    ],
    ids=["barcode_list", "constant_region"],
)
def test_empty_list_or_constant_region_is_a_to_do(core, blocks, empty_index, message):
    """The to-do points at the empty region, so the diagram can mark that region."""
    result = core.run(VALIDATE_WITH_IDS, make_config(blocks))
    assert [
        (issue["rule"], issue["level"], issue["message"]) for issue in result["issues"]
    ] == [("E03", "todo", message)]
    (issue,) = result["issues"]
    assert issue["read"] == 0
    assert issue["block_id"] == result["block_ids"][0][empty_index]


def test_duplicate_message_selects_both_rows(core):
    issues = issues_for(
        core,
        make_config(
            [
                barcode_list("cb1", ["CCTTAACAT", "AGCGTAGAA", "CCTTAACAT"]),
                constant_region("con1"),
            ]
        ),
    )
    (duplicate,) = [issue for issue in issues if issue["rule"] == "E04"]
    assert duplicate["rows"] == [0, 2]


def test_free_list_with_mixed_lengths_is_e05_and_levenshtein_list_is_w10(core):
    mixed = ["CCTTAACAT", "AGCGTAGA", "TATAGGTGT"]
    free = issues_for(
        core, make_config([constant_region("con1"), barcode_list("cb1", mixed)])
    )
    assert "E05" in rules(free) and "W10" not in rules(free)
    levenshtein = issues_for(
        core, make_config([barcode_list("cb1", mixed), constant_region("con1")])
    )
    assert "W10" in rules(levenshtein) and "E05" not in rules(levenshtein)


@pytest.mark.parametrize(
    "blocks",
    [
        [constant_region("con1"), barcode_list("cb1", MIXED_LENGTHS)],
        [
            constant_region("con1"),
            barcode_list("cb1", MIXED_LENGTHS),
            random_block("umi1"),
        ],
    ],
    ids=["last_region", "followed_by_a_random_region"],
)
def test_e05_says_how_to_fix_a_free_list_of_mixed_lengths(core, blocks):
    """E05 names both fixes: a constant region right after the list, or one length."""
    (e05,) = [
        issue
        for issue in issues_for(core, make_config(blocks))
        if issue["rule"] == "E05"
    ]
    prefix = "cb1 has barcodes of different lengths ("
    suffix = (
        "). A list like that needs a constant region right after it, to its right."
        " Add a constant region after cb1, or give all its barcodes the same length."
    )
    assert e05["level"] == "error"
    assert e05["message"].startswith(prefix)
    assert e05["message"].endswith(suffix)
    length_summary = e05["message"][len(prefix) : -len(suffix)]
    assert "9 nt" in length_summary and "8 nt" in length_summary
    assert e05["rows"] == [1]


@pytest.mark.parametrize(
    ("block", "allowed"),
    [
        (barcode_list("cb1", ["CCTTAACAT"]), 4),
        (constant_region("con1", ["GTCAGTACGTACGAGTC", "AGTACGTACGAGTC"]), 7),
        (random_block("umi1", length=8), 4),
    ],
)
def test_maxerrors_cap_is_half_the_shortest_length(core, block, allowed):
    other = (
        constant_region("con9")
        if block["blocktype"] != "constantRegion"
        else barcode_list("cb9")
    )
    for maxerrors, expect_error in [
        (0, False),
        (allowed, False),
        (allowed + 1, True),
        (-1, True),
    ]:
        config = make_config([dict(block, maxerrors=maxerrors), other])
        assert ("E06" in rules(issues_for(core, config))) is expect_error, (
            maxerrors,
            rules(issues_for(core, config)),
        )


@pytest.mark.parametrize(("length", "fires"), [(0, True), (-3, True), (1, False)])
def test_random_length_must_be_at_least_one(core, length, fires):
    config = make_config([barcode_list("cb1"), random_block("umi1", length=length)])
    assert ("E07" in rules(issues_for(core, config))) is fires


@pytest.mark.parametrize(
    "names",
    [
        ["CB1", "cb1"],
        ["cell barcode", "umi"],
        ["cb/1", "umi"],
        ["", "umi"],
        ["x" * 65, "umi"],
    ],
)
def test_block_name_rule_e08(core, names):
    config = make_config([barcode_list(names[0]), random_block(names[1])])
    assert "E08" in rules(issues_for(core, config))


def test_valid_block_names_pass_e08(core):
    config = make_config([barcode_list("cb.1-a_B"), random_block("x" * 64)])
    assert "E08" not in rules(issues_for(core, config))


@pytest.mark.parametrize(
    ("config", "messages"),
    [
        (
            make_config([barcode_list(""), random_block("umi1")]),
            ["Read 1, region 1, needs a name."],
        ),
        (
            make_config(
                [barcode_list("cb1"), constant_region("con1")],
                [random_block("umi2"), barcode_list("")],
            ),
            ["Read 2, region 2, needs a name."],
        ),
        (
            make_config([barcode_list("cell barcode"), random_block("umi1")]),
            [
                'cell barcode: region names may use only letters, digits, "_", "." and "-", up to 64 characters.'
            ],
        ),
        (
            make_config([barcode_list("CB1"), random_block("cb1")]),
            [
                "CB1 and cb1: two regions in Read 1 can't share a name, even in different case."
            ],
        ),
        (
            make_config(
                [barcode_list("cb1"), constant_region("con1")],
                [barcode_list("bc2"), random_block("BC2")],
            ),
            [
                "bc2 and BC2: two regions in Read 2 can't share a name, even in different case."
            ],
        ),
        # Names must differ within a read only, so both reads may use cb1.
        (
            make_config(
                [barcode_list("cb1"), constant_region("con1")],
                [barcode_list("cb1"), constant_region("con2")],
            ),
            [],
        ),
    ],
    ids=[
        "no_name_read_1",
        "no_name_read_2",
        "bad_characters",
        "clash_read_1",
        "clash_read_2",
        "same_name_on_both_reads",
    ],
)
def test_e08_messages_speak_of_regions(core, config, messages):
    issues = issues_for(core, config)
    assert [issue["message"] for issue in issues if issue["rule"] == "E08"] == messages


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("maxerrors", "1"),
        ("maxerrors", 1.5),
        ("blocktype", "fancyBarcode"),
        ("blockfunction", "index"),
        ("keepblock", "yes"),
    ],
)
def test_values_the_app_cannot_represent_are_e11(core, field, value):
    block = barcode_list("cb1")
    block[field] = value
    assert "E11" in rules(
        issues_for(core, make_config([block, constant_region("con1")]))
    )


def test_maxerrors_written_as_1_point_0_is_read_as_1(core):
    text = '{"unknown_read_orientation": false, "single_end_reads": true, "barcode_struct_r1": {"keep_nonbarcode": false, "blocks": [{"maxerrors": 1.0, "blockfunction": "barcode", "sequence": ["CCTTAACAT"], "blocktype": "barcodeList", "blockname": "cb1", "keepblock": false}, {"maxerrors": 0, "blockfunction": "UMI", "length": 8, "blocktype": "randomBarcode", "blockname": "umi", "keepblock": false}]}}'
    result = core.run(
        "(core, arg) => { const model = core.config_from_json(arg).model; return {issues: core.validate_config(model), text: core.config_to_json(model)}; }",
        text,
    )
    assert result["issues"] == []
    assert '"maxerrors": 1,' in result["text"]


def test_one_block_link_group_in_an_opened_file_is_e12(core):
    linked = dict(barcode_list("cb1"), linkgroup=3)
    assert "E12" in rules(
        issues_for(core, make_config([linked, constant_region("con1")]))
    )


def test_removing_the_last_link_of_a_pair_leaves_no_group(core):
    result = core.run(
        """(core, arg) => {
            let model = core.config_from_json(arg).model;
            const [first, second] = model.reads[0].blocks.map((block) => block.id);
            model = core.add_link(model, 0, first, second);
            const linked = model.reads[0].blocks.map((block) => block.linkgroup);
            model = core.remove_link(model, 0, second);
            return {linked, unlinked: model.reads[0].blocks.map((block) => block.linkgroup), issues: core.validate_config(model)};
        }""",
        make_config(
            [barcode_list("cb1"), barcode_list("cb2"), constant_region("con1")]
        ),
    )
    assert result["linked"][0] == result["linked"][1] is not None
    assert result["unlinked"] == [None, None, None]
    assert "E12" not in rules(result["issues"])


def test_linking_three_blocks_makes_one_group(core):
    groups = core.run(
        """(core, arg) => {
            let model = core.config_from_json(arg).model;
            const ids = model.reads[0].blocks.map((block) => block.id);
            model = core.add_link(model, 0, ids[0], ids[1]);
            model = core.add_link(model, 0, ids[0], ids[2]);
            return model.reads[0].blocks.map((block) => block.linkgroup);
        }""",
        make_config(
            [
                barcode_list("cb1"),
                barcode_list("cb2"),
                barcode_list("cb3"),
                constant_region("con1"),
            ]
        ),
    )
    assert groups[0] == groups[1] == groups[2] is not None
    assert groups[3] is None


def test_one_block_read_gives_only_w11(core):
    issues = issues_for(core, make_config([barcode_list("cb1")]))
    assert rules(issues) == ["W11"]
    assert issues[0]["temporary"] is True


def test_barcodes_only_on_read_2_give_only_w12(core):
    issues = issues_for(
        core, make_config([], [barcode_list("bc1"), constant_region("con1")])
    )
    assert rules(issues) == ["W12"]
    assert issues[0]["temporary"] is True


def test_w08_fires_only_without_a_constant_region(core):
    without_constant = make_config([barcode_list("cb1"), random_block("umi1")], [])
    with_constant = make_config([barcode_list("cb1"), constant_region("con1")], [])
    assert rules(issues_for(core, without_constant)) == ["W08"]
    assert rules(issues_for(core, with_constant)) == []


def test_w05_fires_for_paired_end_orientation_only(core):
    paired = make_config(
        [barcode_list("cb1"), constant_region("con1")], [], orientation=True
    )
    issues = issues_for(core, paired)
    assert rules(issues) == ["W05"]
    assert "no effect" not in issues[0]["message"]
    single = make_config(
        [barcode_list("cb1"), constant_region("con1")], orientation=True
    )
    assert rules(issues_for(core, single)) == []


def test_w05_says_the_setting_has_no_effect_without_constants(core):
    config = make_config(
        [barcode_list("cb1"), random_block("umi1")],
        [random_block("umi2"), barcode_list("cb2")],
        orientation=True,
    )
    (w05,) = [issue for issue in issues_for(core, config) if issue["rule"] == "W05"]
    assert "no effect" in w05["message"]


def test_discard_slip_w01_fires_only_when_a_kept_list_follows(core):
    slip = make_config(
        [
            barcode_list("cb1", function="discard"),
            constant_region("con1"),
            barcode_list("cb2"),
        ]
    )
    no_slip = make_config(
        [
            barcode_list("cb1"),
            constant_region("con1"),
            barcode_list("cb2", function="discard"),
        ]
    )
    assert "W01" in rules(issues_for(core, slip))
    assert "W01" not in rules(issues_for(core, no_slip))


def test_discard_slip_on_constant_regions_is_w01(core):
    config = make_config(
        [
            constant_region("con1", function="discard"),
            barcode_list("cb1"),
            constant_region("con2"),
        ]
    )
    (w01,) = [issue for issue in issues_for(core, config) if issue["rule"] == "W01"]
    assert "FL" in w01["message"]


@pytest.mark.parametrize(
    "block",
    [barcode_list("cb1", function="UMI"), constant_region("con1", function="UMI")],
)
def test_umi_on_a_list_or_constant_is_w02(core, block):
    other = (
        constant_region("con9")
        if block["blocktype"] == "barcodeList"
        else barcode_list("cb9")
    )
    assert "W02" in rules(issues_for(core, make_config([block, other])))


def test_random_block_set_to_barcode_is_w06(core):
    config = make_config(
        [barcode_list("cb1"), random_block("rnd1", function="barcode")]
    )
    assert "W06" in rules(issues_for(core, config))


def test_config_without_a_barcode_block_is_w07(core):
    config = make_config([constant_region("con1"), random_block("umi1")])
    assert "W07" in rules(issues_for(core, config))


@pytest.mark.parametrize(
    ("config", "rule", "message"),
    [
        (
            make_config([barcode_list("cb1")]),
            "W11",
            "Read 1 has a single region. Current bcwithqc crashes on a read with only one region.",
        ),
        (
            make_config(
                [barcode_list("cb1"), constant_region("con1")], [random_block("umi2")]
            ),
            "W11",
            "Read 2 has a single region. Current bcwithqc crashes on a read with only one region.",
        ),
        (
            make_config([], [barcode_list("bc1"), constant_region("con1")]),
            "W12",
            "Read 1 has no regions but Read 2 does. Current bcwithqc crashes on this.",
        ),
        (
            make_config([constant_region("con1"), random_block("umi1")]),
            "W07",
            "No region is a barcode, so CB will be empty.",
        ),
        (
            make_config(
                [barcode_list("cb1"), random_block("rnd1", function="barcode")]
            ),
            "W06",
            "rnd1 is a random region set to barcode, so its bases go into CB without error correction.",
        ),
        (
            make_config(
                [dict(barcode_list("cb1"), linkgroup=3), constant_region("con1")]
            ),
            "E12",
            "cb1 is the only region in its link group. Link it with another region or remove the link.",
        ),
    ],
    ids=["w11_read_1", "w11_read_2", "w12", "w07", "w06", "e12"],
)
def test_messages_about_regions_say_region(core, config, rule, message):
    issues = issues_for(core, config)
    assert [issue["message"] for issue in issues if issue["rule"] == rule] == [message]


@pytest.mark.parametrize(("variants", "fires"), [(4, False), (8, False), (9, True)])
def test_more_than_64_linker_combinations_is_w09(core, variants, fires):
    constants = [
        constant_region(f"con{i}", ["ACGTACGT" + "A" * k for k in range(variants)])
        for i in (1, 2)
    ]
    config = make_config([barcode_list("cb1"), *constants])
    assert ("W09" in rules(issues_for(core, config))) is fires


def test_free_table_bytes_follow_freebarcodes_allocation(core):
    # (4^L + 1) entries; 1 byte while barcodes + 1 <= 256, 2 bytes up to 65,536.
    assert core.call("free_table_bytes", 4, 8) == 4**8 + 1
    assert core.call("free_table_bytes", 255, 8) == 4**8 + 1
    assert core.call("free_table_bytes", 256, 8) == (4**8 + 1) * 2
    assert core.call("free_table_bytes", 65535, 8) == (4**8 + 1) * 2
    assert core.call("free_table_bytes", 65536, 8) == (4**8 + 1) * 4


def free_list_config(length, count):
    barcodes = [
        format(i, "b").zfill(length).replace("0", "A").replace("1", "C")
        for i in range(count)
    ]
    return make_config([constant_region("con1"), barcode_list("sb", barcodes)])


@pytest.mark.parametrize(
    ("length", "count", "fires"), [(13, 300, False), (14, 4, True), (13, 4, False)]
)
def test_large_free_tables_are_w04(core, length, count, fires):
    # 4^14 + 1 bytes is just over 256 MiB; 4^13 + 1 entries of 2 bytes is 128 MiB.
    issues = issues_for(core, free_list_config(length, count))
    w04 = [issue["message"] for issue in issues if issue["rule"] == "W04"]
    assert bool(w04) is fires
    if fires:
        # Measured on bcwithqc: --threads 1 builds one table, --threads N >= 2 builds N + 1.
        (message,) = w04
        assert "256 MiB per copy" in message
        assert "one copy with --threads 1" in message
        assert "N + 1 copies with --threads N of 2 or more" in message


def test_hidden_read_2_never_blocks_download(core):
    result = core.run(
        """(core, arg) => {
            const model = core.config_from_json(arg).model;
            model.reads[1].blocks[0].sequence.push(model.reads[1].blocks[0].sequence[0]);
            const paired = core.validate_config(model).map((issue) => issue.rule);
            model.single_end = true;
            return {paired, single: core.validate_config(model).map((issue) => issue.rule)};
        }""",
        (EXAMPLES_DIR / "pe_mini_config.json").read_text(encoding="utf-8"),
    )
    assert "E04" in result["paired"]
    assert result["single"] == []


def test_every_temporary_warning_is_flagged_temporary(core):
    registry = core.run("(core) => core.TEMPORARY_WARNINGS")
    assert sorted(registry) == ["W01", "W05", "W08", "W11", "W12", "W14"]
    for entry in registry.values():
        assert set(entry) == {"bcwithqc", "notes", "tier_f"}


def test_temporary_output_notes_are_the_ur_and_empty_file_notes(core):
    """A read without regions is always kept whole, so it needs no note of its own."""
    registry = core.run("(core) => core.TEMPORARY_OUTPUT_NOTES")
    assert sorted(registry) == ["empty_alignment_file", "ur_without_umi"]
    for entry in registry.values():
        assert set(entry) == {"bcwithqc", "notes", "tier_f"}


def test_messages_name_no_notes_ids_and_no_line_numbers_in_code(core):
    """Rule R4: messages are complete on their own."""
    configs = [
        make_config([barcode_list("cb1")]),
        make_config([], [barcode_list("bc1"), constant_region("con1")]),
        make_config([barcode_list("cb1"), random_block("umi1")], [], orientation=True),
    ]
    for config in configs:
        for issue in issues_for(core, config):
            assert "#" not in issue["message"]
            assert ".py" not in issue["message"]


def test_every_temporary_warning_has_a_strict_xfail_tier_f_case(core):
    """Rule R5: when bcwithqc is fixed, the case passes unexpectedly and the run fails."""
    registry = core.run(
        "(core) => ({...core.TEMPORARY_WARNINGS, ...core.TEMPORARY_OUTPUT_NOTES})"
    )
    source = (Path(__file__).parent / "test_gui_pipeline_contract.py").read_text(
        encoding="utf-8"
    )
    functions = {
        node.name: node
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.FunctionDef)
    }
    for rule, entry in registry.items():
        node = functions.get(entry["tier_f"])
        assert node is not None, f"{rule}: no tier F case named {entry['tier_f']}"
        decorators = [ast.unparse(decorator) for decorator in node.decorator_list]
        assert any("xfail" in text and "strict=True" in text for text in decorators), (
            rule,
            decorators,
        )


@pytest.mark.parametrize(
    ("layout", "maxerrors", "fires"),
    [
        # A list followed by a constant: Levenshtein in bcwithqc today.
        ([barcode_list("cb1"), constant_region("con1")], 2, True),
        ([barcode_list("cb1"), constant_region("con1")], 1, False),
        # The read's last list before a random block: FREE in bcwithqc, unaffected.
        ([constant_region("con1"), barcode_list("cb1"), random_block("umi")], 2, False),
        # A first list before a random block: FREE in the GUI's rule, Levenshtein in bcwithqc.
        (
            [
                barcode_list("cb1"),
                random_block("umi"),
                barcode_list("cb2", maxerrors=1),
            ],
            2,
            True,
        ),
    ],
)
def test_w14_fires_on_lists_bcwithqc_decodes_with_levenshtein(
    core, layout, maxerrors, fires
):
    blocks = [dict(block) for block in layout]
    for block in blocks:
        if block["blockname"] == "cb1":
            block["maxerrors"] = maxerrors
    issues = issues_for(core, make_config(blocks))
    w14 = [issue for issue in issues if issue["rule"] == "W14"]
    assert bool(w14) is fires
    if fires:
        assert w14[0]["temporary"] is True
        assert "farthest" in w14[0]["message"]


# Every rule in plan §08.
ALL_RULES = [
    "E01",
    "E03",
    "E04",
    "E05",
    "E06",
    "E07",
    "E08",
    "E11",
    "E12",
    "E14",
    "W01",
    "W02",
    "W04",
    "W05",
    "W06",
    "W07",
    "W08",
    "W09",
    "W10",
    "W11",
    "W12",
    "W14",
]


def configs_that_trigger_every_rule() -> list[dict]:
    """Return small configs that together trigger every rule, in each of its wordings.

    No region name contains "block", so a message can only say it in its own
    words.
    """
    list_and_constant = [barcode_list("cb1"), constant_region("con1")]

    def with_list_setting(**setting):
        """Return a single-end config of list cb1, with these settings, then con1."""
        return make_config(
            [dict(barcode_list("cb1"), **setting), constant_region("con1")]
        )

    nine_variants_each = [
        constant_region(f"con{i}", ["ACGTACGT" + "A" * k for k in range(9)])
        for i in (1, 2)
    ]
    return [
        # E01 and E14, in its paired-end and single-end wordings.
        dict(make_config([], []), single_end_reads=None),
        make_config([]),
        # E11 on the layout, the orientation, the rest of a read and each region setting.
        dict(make_config(list_and_constant), single_end_reads="yes"),
        dict(make_config(list_and_constant), unknown_read_orientation="no"),
        make_config(list_and_constant, keep_r1="yes"),
        with_list_setting(blocktype="fancyBarcode"),
        with_list_setting(blockfunction="index"),
        with_list_setting(keepblock="yes"),
        with_list_setting(linkgroup=1.5),
        with_list_setting(maxerrors="1"),
        with_list_setting(sequence="CCTTAACAT"),
        with_list_setting(sequence=["CCTTAACAT", 5]),
        make_config([barcode_list("cb1"), dict(random_block("umi1"), length="8")]),
        # E08: no name, characters outside the allowed set, a clash.
        make_config([barcode_list(""), random_block("umi1")]),
        make_config([barcode_list("cell barcode"), random_block("umi1")]),
        make_config([barcode_list("CB1"), random_block("cb1")]),
        # E03: an empty list and constant region, an empty entry, bad characters.
        make_config([barcode_list("cb1", []), constant_region("con1", [])]),
        with_list_setting(sequence=["CCTTAACAT", ""]),
        with_list_setting(
            sequence=["CCTTAACAT", "CCTTANCAT", "CCTTAACAT" + chr(0x200B)]
        ),
        # E04
        with_list_setting(sequence=["CCTTAACAT", "AGCGTAGAA", "CCTTAACAT"]),
        # E05 (decoded with FREE) and W10 (Levenshtein) for mixed lengths.
        make_config([constant_region("con1"), barcode_list("cb1", MIXED_LENGTHS)]),
        make_config([barcode_list("cb1", MIXED_LENGTHS), constant_region("con1")]),
        # E06 on a list and on a random region; E07.
        with_list_setting(maxerrors=5),
        make_config([barcode_list("cb1"), random_block("umi1", maxerrors=5)]),
        make_config([barcode_list("cb1"), random_block("umi1", length=0)]),
        # E12
        with_list_setting(linkgroup=3),
        # W01 on lists and on constant regions.
        make_config(
            [
                barcode_list("cb1", function="discard"),
                constant_region("con1"),
                barcode_list("cb2"),
            ]
        ),
        make_config(
            [
                constant_region("con1", function="discard"),
                barcode_list("cb1"),
                constant_region("con2"),
            ]
        ),
        # W02 on a list and on a constant region.
        with_list_setting(blockfunction="UMI"),
        make_config([barcode_list("cb1"), constant_region("con1", function="UMI")]),
        # W04
        free_list_config(14, 4),
        # W05 with and without a constant region.
        make_config(list_and_constant, [], orientation=True),
        make_config(
            [barcode_list("cb1"), random_block("umi1")],
            [random_block("umi2"), barcode_list("cb2")],
            orientation=True,
        ),
        # W06, W07, W08
        make_config([barcode_list("cb1"), random_block("rnd1", function="barcode")]),
        make_config([constant_region("con1"), random_block("umi1")]),
        make_config([barcode_list("cb1"), random_block("umi1")], []),
        # W09
        make_config([barcode_list("cb1"), *nine_variants_each]),
        # W11, W12
        make_config([barcode_list("cb1")]),
        make_config([], list_and_constant),
        # W14
        with_list_setting(maxerrors=2),
    ]


def issues_from_every_kind_of_config(core) -> list[dict]:
    """Return the issues of the examples, the legacy fixtures and configs_that_trigger_every_rule().

    Each real config is opened both ways a legacy file can open: as
    paired-end and as single-end.
    """
    issues = []
    for _, path in REAL_CONFIG_TEXTS:
        for legacy_mode in ("paired", "single"):
            issues += issues_for(core, read_bytes_text(path), legacy_mode=legacy_mode)
    for config in configs_that_trigger_every_rule():
        issues += issues_for(core, config)
    return issues


def test_no_message_says_block(core):
    """Messages call bcwithqc's blocks "regions"; "block" stays only in the config's keys."""
    issues = issues_from_every_kind_of_config(core)
    assert sorted({issue["rule"] for issue in issues}) == ALL_RULES
    assert [
        issue["message"] for issue in issues if "block" in issue["message"].lower()
    ] == []


def test_only_steps_not_done_yet_are_to_dos(core):
    """E01, E14 and E03 for an empty list or constant region are to-dos.

    Every other E rule is an error and every W rule a warning. The to-do
    E03 messages are the ones that ask to add something.
    """
    issues = issues_from_every_kind_of_config(core)
    for issue in issues:
        if issue["rule"].startswith("W"):
            expected = "warning"
        elif issue["rule"] in ("E01", "E14") or (
            issue["rule"] == "E03" and issue["message"].startswith("Add ")
        ):
            expected = "todo"
        else:
            expected = "error"
        assert issue["level"] == expected, issue
    to_do_rules = {issue["rule"] for issue in issues if issue["level"] == "todo"}
    assert sorted(to_do_rules) == ["E01", "E03", "E14"]
