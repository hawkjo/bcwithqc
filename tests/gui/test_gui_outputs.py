"""Tier D: what the Outputs panel and the structure caption describe (plan §15).

Values are worked out by hand from bcwithqc's tag code: CB, CR and UR join
their pieces with "."; with blocks on both reads, the two reads' values are
joined with "." again; FL sums the constant regions not set to discard.
"""

from gui_support import (
    EXAMPLES_DIR,
    barcode_list,
    constant_region,
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
    assert outputs["alignment_files"] == [
        {"read": 0, "kind": "none", "parts": [], "rest": False},
        {"read": 1, "kind": "whole", "parts": [], "rest": False},
    ]


def test_sdr_seq_rna_structure_line(core):
    structure = describe(core, example_text("cDNA_paired_end_config.json"))["structure"]
    # 9 + (14..17) + 9 + 23 + 8 + 8; reads must hold all but the last block at their longest, plus 1.
    assert structure[0] == {
        "min_length": 71,
        "max_length": 74,
        "linker_combinations": 4,
        "min_read_length": 67,
    }
    assert structure[1] is None


def test_sending_the_rest_of_read_1_to_alignment(core):
    text = example_text("cDNA_paired_end_config.json").replace(
        '"keep_nonbarcode": false', '"keep_nonbarcode": true', 1
    )
    outputs = describe(core, text)["outputs"]
    assert outputs["alignment_files"][0] == {
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
    assert outputs["alignment_files"][0] == {
        "read": 0,
        "kind": "parts",
        "parts": ["con_r1"],
        "rest": False,
    }
    assert outputs["alignment_files"][1] == {
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
