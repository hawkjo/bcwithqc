"""Tier D: validation rules.

Each rule must fire exactly at its boundary, and bcwithqc's real configs must
pass with no errors and no warnings. Rule ids follow plan §08.
"""

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
    issues = core.run("(core) => core.validate_config(core.new_model())")
    assert rules(issues) == ["E01", "E14"]
    assert all(issue["level"] == "error" for issue in issues)


def test_one_block_anywhere_satisfies_e14(core):
    config = make_config([], [barcode_list("bc1"), constant_region("con1")])
    assert "E14" not in rules(issues_for(core, config))


@pytest.mark.parametrize(
    ("sequences", "expected_rule", "text"),
    [
        ([], "E03", "cb1 has no sequences"),
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
    assert ("W04" in rules(issues_for(core, free_list_config(length, count)))) is fires


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
    assert sorted(registry) == ["W05", "W08", "W11", "W12"]
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
