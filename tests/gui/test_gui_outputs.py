"""Tier D: what the Outputs panel and the structure caption describe (plan §15).

Values are worked out by hand from bcwithqc's tag code: CB, CR and UR join
their pieces with "."; with blocks on both reads, the two reads' values are
joined with "." again; FL sums the constant regions not set to discard. A
read without regions is always kept whole.

The caption's read lengths are worked out by hand from the regions' lengths
and allowed errors: a read must be at least the shortest structure less
every allowed error (and at least 1 nt); the recommended length is the
longest structure plus every allowed error.
"""

import pytest
from gui_support import (
    EXAMPLES_DIR,
    barcode_list,
    constant_region,
    design_doc_example,
    make_config,
    random_block,
)

DESCRIBE = """(core, arg) => {
    const model = core.config_from_json(arg).model;
    return {
        outputs: core.describe_outputs(model),
        structure: model.reads.map((read) => core.describe_structure(read.blocks)),
    };
}"""


def describe(core, config):
    return core.run(DESCRIBE, config)


def example_text(name):
    return (EXAMPLES_DIR / name).read_text(encoding="utf-8")


def test_sdr_seq_rna_outputs(core):
    outputs = describe(core, example_text("cDNA_paired_end_config.json"))["outputs"]
    assert outputs["CB"]["value"] == "cb1.cb2.sb"
    assert outputs["CR"]["value"] == "cb1.cb2.sb"
    assert outputs["UR"]["value"] == "umi"
    assert outputs["UB"]["written"] is True
    assert (outputs["FL"]["min"], outputs["FL"]["max"]) == (37, 40)
    assert outputs["output_reads"] == [
        {"read": 0, "kind": "none", "parts": [], "rest": False},
        {"read": 1, "kind": "whole", "parts": [], "rest": False},
    ]


def test_sdr_seq_rna_structure_line(core):
    """Read 1: cb1 9 nt, cso 14 to 17 nt (4 variants), cb2 9, commonseq2 23, sb 8, umi 8.

    Every region allows 1 error except umi, which allows none. Read 2 has no
    regions, so it has no structure line.
    """
    structure = describe(core, example_text("cDNA_paired_end_config.json"))["structure"]
    shortest = 9 + 14 + 9 + 23 + 8 + 8
    longest = 9 + 17 + 9 + 23 + 8 + 8
    allowed_errors = 1 + 1 + 1 + 1 + 1 + 0
    assert structure[0] == {
        "lengths_known": True,
        "min_length": shortest,  # 71
        "max_length": longest,  # 74
        "allowed_errors": allowed_errors,  # 5
        "linker_combinations": 4,
        "min_read_length": shortest - allowed_errors,  # 66
        "recommended_read_length": longest + allowed_errors,  # 79
    }
    assert structure[1] is None


def test_design_doc_example_structure_line(core):
    """Read 1: cb1 9 nt and con1 21 nt with 1 error each, then umi 6 nt with none."""
    structure = describe(core, design_doc_example())["structure"]
    length = 9 + 21 + 6
    allowed_errors = 1 + 1 + 0
    assert structure[0] == {
        "lengths_known": True,
        "min_length": length,  # 36
        "max_length": length,
        "allowed_errors": allowed_errors,  # 2
        "linker_combinations": 1,
        "min_read_length": length - allowed_errors,  # 34
        "recommended_read_length": length + allowed_errors,  # 38
    }
    assert structure[1] is None


def read_1_structure(core, blocks: list[dict]) -> dict:
    """Return the structure line of Read 1 in a single-end config of these regions."""
    return describe(core, make_config(blocks))["structure"][0]


def test_a_read_without_regions_has_no_structure_line(core):
    assert core.call("describe_structure", []) is None


@pytest.mark.parametrize(
    "barcodes",
    [["AGCGTAGA", "CCTTAACAT"], ["CCTTAACAT", "AGCGTAGA"]],
    ids=["8_nt_first", "9_nt_first"],
)
def test_a_list_of_8_and_9_nt_barcodes_spans_both_lengths_in_either_order(
    core, barcodes
):
    """cb1 is 8 or 9 nt and con1 21 nt, with 1 error each."""
    structure = read_1_structure(
        core, [barcode_list("cb1", barcodes), constant_region("con1")]
    )
    assert structure == {
        "lengths_known": True,
        "min_length": 8 + 21,
        "max_length": 9 + 21,
        "allowed_errors": 2,
        "linker_combinations": 1,
        "min_read_length": 8 + 21 - 2,
        "recommended_read_length": 9 + 21 + 2,
    }


def test_an_empty_entry_does_not_shorten_the_structure(core):
    """An empty entry is an E03 error to fix; the lengths come from cb1's 9 nt barcode."""
    structure = read_1_structure(
        core, [barcode_list("cb1", ["CCTTAACAT", ""]), constant_region("con1")]
    )
    assert (structure["min_length"], structure["max_length"]) == (9 + 21, 9 + 21)


@pytest.mark.parametrize(
    ("maxerrors", "min_read_length"),
    [(1, 1), (2, 1), (3, 1)],
    ids=["1_nt_left", "nothing_left", "less_than_nothing_left"],
)
def test_min_read_length_is_never_below_1_nt(core, maxerrors, min_read_length):
    """A 2 nt random region less its allowed errors leaves 1, 0 or -1 nt: each is shown as 1 nt.

    Allowing that many errors is an E06 error, but the caption still shows lengths.
    """
    structure = read_1_structure(
        core, [random_block("rnd1", length=2, maxerrors=maxerrors)]
    )
    assert structure["lengths_known"] is True
    assert structure["min_read_length"] == min_read_length
    assert structure["recommended_read_length"] == 2 + maxerrors


def test_linker_combinations_multiply_the_variants_of_constant_regions(core):
    """con1's 2 variants times con2's 3; a barcode list adds no combinations."""
    structure = read_1_structure(
        core,
        [
            constant_region("con1", ["GGGTCAGTAC", "GGGTCAGTA"]),
            barcode_list("cb1"),
            constant_region("con2", ["ACGTACGTAC", "ACGTACGTA", "ACGTACGT"]),
        ],
    )
    assert structure["linker_combinations"] == 2 * 3
    assert structure["lengths_known"] is True


TWO_VARIANTS = ["GGGTCAGTACGTACGAGTCCC", "GGTCAGTACGTACGAGTCCC"]


@pytest.mark.parametrize(
    "blocks",
    [
        [
            barcode_list("cb1", []),
            constant_region("con1", TWO_VARIANTS),
            random_block("umi1"),
        ],
        [
            barcode_list("cb1", [""]),
            constant_region("con1", TWO_VARIANTS),
            random_block("umi1"),
        ],
        [
            barcode_list("cb1"),
            constant_region("con1", TWO_VARIANTS),
            constant_region("con2", []),
            random_block("umi1"),
        ],
        [
            barcode_list("cb1"),
            constant_region("con1", TWO_VARIANTS),
            random_block("umi1", length=0),
        ],
        [
            barcode_list("cb1"),
            constant_region("con1", TWO_VARIANTS),
            random_block("umi1", length="8"),
        ],
        [
            barcode_list("cb1", maxerrors="1"),
            constant_region("con1", TWO_VARIANTS),
            random_block("umi1"),
        ],
        [
            barcode_list("cb1", maxerrors=-1),
            constant_region("con1", TWO_VARIANTS),
            random_block("umi1"),
        ],
    ],
    ids=[
        "empty_list",
        "list_of_one_empty_entry",
        "empty_constant_region",
        "random_length_0",
        "random_length_as_text",
        "maxerrors_as_text",
        "negative_maxerrors",
    ],
)
def test_lengths_are_unknown_until_every_region_is_filled_in(core, blocks):
    """The linker combinations are still counted: con1's 2 variants, times 1 for an empty constant region."""
    assert read_1_structure(core, blocks) == {
        "lengths_known": False,
        "min_length": None,
        "max_length": None,
        "allowed_errors": None,
        "linker_combinations": 2,
        "min_read_length": None,
        "recommended_read_length": None,
    }


def test_keeping_the_rest_of_read_1_in_the_output_reads(core):
    text = example_text("cDNA_paired_end_config.json").replace(
        '"keep_nonbarcode": false', '"keep_nonbarcode": true', 1
    )
    outputs = describe(core, text)["outputs"]
    assert outputs["output_reads"][0] == {
        "read": 0,
        "kind": "parts",
        "parts": [],
        "rest": True,
    }


def test_pe_mini_outputs(core):
    outputs = describe(core, example_text("pe_mini_config.json"))["outputs"]
    assert outputs["CB"]["value"] == "bc_r1.bc_r2"
    assert outputs["UR"]["value"] == "umi_r1.umi_r2"
    assert (outputs["FL"]["min"], outputs["FL"]["max"]) == (42, 42)
    assert outputs["output_reads"][0] == {
        "read": 0,
        "kind": "parts",
        "parts": ["con_r1"],
        "rest": False,
    }
    assert outputs["output_reads"][1] == {
        "read": 1,
        "kind": "parts",
        "parts": [],
        "rest": True,
    }


def test_list_on_read_2_only_gives_a_leading_dot(core):
    config = make_config(
        [constant_region("con_r1"), random_block("umi_r1")],
        [barcode_list("bc_r2"), constant_region("con_r2")],
    )
    outputs = describe(core, config)["outputs"]
    assert outputs["CB"]["value"] == ".bc_r2"


def test_umi_on_read_1_only_gives_a_trailing_dot_when_read_2_has_blocks(core):
    config = make_config(
        [barcode_list("bc_r1"), random_block("umi_r1")],
        [barcode_list("bc_r2"), constant_region("con_r2")],
    )
    assert describe(core, config)["outputs"]["UR"]["value"] == "umi_r1."


def test_no_umi_block_means_no_ub_and_an_empty_ur(core):
    outputs = describe(core, example_text("gDNA_paired_end_config.json"))["outputs"]
    assert outputs["UB"]["written"] is False
    # Decision D9: without a UMI block, UR is not written either (bcwithqc today writes it empty).
    assert outputs["UR"]["written"] is False
    assert (outputs["FL"]["min"], outputs["FL"]["max"]) == (29, 32)


def test_discarded_list_is_left_out_of_cb(core):
    config = make_config(
        [
            barcode_list("cb1"),
            constant_region("con1"),
            barcode_list("cb2", function="discard"),
        ]
    )
    assert describe(core, config)["outputs"]["CB"]["value"] == "cb1"


def test_random_block_set_to_barcode_joins_cb(core):
    config = make_config(
        [barcode_list("cb1"), random_block("rnd1", function="barcode")]
    )
    outputs = describe(core, config)["outputs"]
    assert outputs["CB"]["value"] == "cb1.rnd1"
    assert outputs["UR"]["value"] == ""


REGIONS = [barcode_list("cb1"), constant_region("con1")]
DESCRIBE_FILES_WITH_REST_SET = """(core, arg) => {
    const model = core.config_from_json(arg.config).model;
    model.reads[arg.read].keep_nonbarcode = arg.keep_nonbarcode;
    return core.describe_outputs(model).output_reads;
}"""


@pytest.mark.parametrize("empty_read", [0, 1], ids=["read_1", "read_2"])
@pytest.mark.parametrize("keep_nonbarcode", [True, False], ids=["output", "left_out"])
def test_read_without_regions_is_always_output_whole(core, empty_read, keep_nonbarcode):
    """Read 1 and Read 2 behave alike: a read without regions is output whole.

    The model's keep_nonbarcode is set after opening, so "whole" can't come
    from config_from_json having set it to true.
    """
    if empty_read == 0:
        config = make_config([], REGIONS)
    else:
        config = make_config(REGIONS, [])
    files = core.run(
        DESCRIBE_FILES_WITH_REST_SET,
        {"config": config, "read": empty_read, "keep_nonbarcode": keep_nonbarcode},
    )
    assert files[empty_read] == {
        "read": empty_read,
        "kind": "whole",
        "parts": [],
        "rest": False,
    }


OUTPUT_WHOLE_READS = """(core, arg) => {
    const model = core.config_from_json(arg.config).model;
    arg.keep_before.forEach((value, index) => { model.reads[index].keep_nonbarcode = value; });
    const changed = core.output_whole_reads_without_blocks(model);
    return {changed, keep_after: model.reads.map((read) => read.keep_nonbarcode)};
}"""


@pytest.mark.parametrize(
    ("config", "keep_before", "changed", "keep_after"),
    [
        (make_config([], REGIONS), [False, False], [0], [True, False]),
        (make_config(REGIONS, []), [False, False], [1], [False, True]),
        (make_config([], []), [False, False], [0, 1], [True, True]),
        (make_config(REGIONS, []), [False, True], [], [False, True]),
        (make_config(REGIONS, []), [True, "yes"], [1], [True, True]),
        (make_config(REGIONS, REGIONS), [False, False], [], [False, False]),
    ],
    ids=[
        "read_1_without_regions",
        "read_2_without_regions",
        "both_without_regions",
        "already_output",
        "not_a_boolean",
        "both_with_regions",
    ],
)
def test_output_whole_reads_without_blocks_sets_only_reads_without_regions(
    core, config, keep_before, changed, keep_after
):
    """It changes the model it is given and returns the reads whose keep_nonbarcode it changed."""
    result = core.run(
        OUTPUT_WHOLE_READS, {"config": config, "keep_before": keep_before}
    )
    assert result == {"changed": changed, "keep_after": keep_after}


REMOVE_FIRST_REGION = """(core, arg) => {
    const model = core.config_from_json(arg.config).model;
    model.reads[0].keep_nonbarcode = arg.keep_nonbarcode;
    const after = core.remove_block(model, 0, model.reads[0].blocks[0].id);
    return {
        regions_left: after.reads[0].blocks.length,
        keep_nonbarcode: after.reads[0].keep_nonbarcode,
        original_keep_nonbarcode: model.reads[0].keep_nonbarcode,
    };
}"""


def test_removing_the_last_region_of_a_read_outputs_it_whole(core):
    """The rest was "Leave out"; with no regions left, the whole read is output.

    remove_block returns a new model, so the one it was given keeps false.
    """
    result = core.run(
        REMOVE_FIRST_REGION,
        {"config": make_config([barcode_list("cb1")]), "keep_nonbarcode": False},
    )
    assert result == {
        "regions_left": 0,
        "keep_nonbarcode": True,
        "original_keep_nonbarcode": False,
    }


@pytest.mark.parametrize("keep_nonbarcode", [False, True], ids=["left_out", "output"])
def test_removing_one_of_two_regions_leaves_the_rest_of_the_read_alone(
    core, keep_nonbarcode
):
    result = core.run(
        REMOVE_FIRST_REGION,
        {"config": make_config(REGIONS), "keep_nonbarcode": keep_nonbarcode},
    )
    assert result == {
        "regions_left": 1,
        "keep_nonbarcode": keep_nonbarcode,
        "original_keep_nonbarcode": keep_nonbarcode,
    }


@pytest.mark.parametrize(
    ("single_end", "cb_value", "corrected", "uncorrected"),
    [
        (False, "sb1.rnd2.cb3.rnd1.cb2", ["sb1", "cb3", "cb2"], ["rnd2", "rnd1"]),
        (True, "sb1.rnd2.cb3", ["sb1", "cb3"], ["rnd2"]),
    ],
    ids=["paired_end", "single_end_hides_read_2"],
)
def test_corrected_and_uncorrected_name_the_cb_regions_in_cb_order(
    core, single_end, cb_value, corrected, uncorrected
):
    """corrected names the barcode lists in CB, uncorrected its random regions.

    Regions set to discard are in neither. The names put CB order apart
    from alphabetical order.
    """
    config = make_config(
        [
            barcode_list("sb1"),
            random_block("rnd2", function="barcode"),
            barcode_list("dropped_list", function="discard"),
            constant_region("con1"),
            barcode_list("cb3"),
            random_block("umi1"),
        ],
        [
            random_block("rnd1", function="barcode"),
            barcode_list("cb2"),
            random_block("dropped_random", function="discard"),
            constant_region("con2"),
        ],
        single_end=single_end,
    )
    cb = describe(core, config)["outputs"]["CB"]
    assert cb["value"] == cb_value
    assert cb["corrected"] == corrected
    assert cb["uncorrected"] == uncorrected


def test_new_model_outputs_both_reads_whole(core):
    """Both reads start the same way: with no regions, each is output whole."""
    outputs = core.run("(core) => core.describe_outputs(core.new_model())")
    assert outputs["output_reads"] == [
        {"read": 0, "kind": "whole", "parts": [], "rest": False},
        {"read": 1, "kind": "whole", "parts": [], "rest": False},
    ]
    assert outputs["CB"]["value"] == ""
    assert outputs["CB"]["uncorrected"] == []


def test_uncorrected_lists_random_regions_set_to_barcode_in_cb_order(core):
    """CB.uncorrected names the random regions in CB, whose bases go in as read."""
    config = make_config(
        [
            barcode_list("cb1"),
            random_block("rnd1", function="barcode"),
            constant_region("con1"),
            random_block("rnd2", function="barcode"),
            random_block("umi1"),
        ],
        [
            random_block("rnd3", function="barcode"),
            barcode_list("cb2"),
            random_block("dropped", function="discard"),
            constant_region("con2"),
        ],
    )
    cb = describe(core, config)["outputs"]["CB"]
    assert cb["value"] == "cb1.rnd1.rnd2.rnd3.cb2"
    assert cb["uncorrected"] == ["rnd1", "rnd2", "rnd3"]


@pytest.mark.parametrize(
    ("config", "cb_value"),
    [
        (example_text("cDNA_paired_end_config.json"), "cb1.cb2.sb"),
        (example_text("pe_mini_config.json"), "bc_r1.bc_r2"),
        (
            make_config(
                [
                    barcode_list("cb1"),
                    random_block("umi1"),
                    random_block("dropped", function="discard"),
                ]
            ),
            "cb1",
        ),
        # A single-end config doesn't write Read 2, so its random region is not in CB.
        (
            make_config(
                [barcode_list("cb1"), constant_region("con1")],
                [random_block("rnd1", function="barcode"), barcode_list("cb2")],
                single_end=True,
            ),
            "cb1",
        ),
    ],
    ids=["sdr_seq_rna", "pe_mini", "umi_and_discarded_random", "hidden_read_2"],
)
def test_uncorrected_is_empty_without_random_regions_set_to_barcode(
    core, config, cb_value
):
    cb = describe(core, config)["outputs"]["CB"]
    assert cb["value"] == cb_value
    assert cb["uncorrected"] == []
