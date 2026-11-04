"""Tier E: the decoder label of each barcode list (decision D1, plan §07).

A list followed by a constant region or another list has its end fixed and
uses Levenshtein. A list that is the read's last block, or is followed by a
random block, uses FREE. The expected labels below are the plan's table.
"""

import pytest
from gui_support import layout_blocks

LEV, FREE = "Levenshtein", "FREE"
LAYOUTS = {
    "LCLCLR": [LEV, LEV, FREE],
    "LCLC": [LEV, LEV],
    "LCLR": [LEV, FREE],
    "LCL": [LEV, FREE],
    "LL": [LEV, FREE],
    "LRC": [FREE],
    "LC": [LEV],
    "LRCL": [FREE, FREE],
    "LRL": [FREE, FREE],
    "LRLC": [FREE, LEV],
    "LRCLC": [FREE, LEV],
    "LLLC": [LEV, LEV, LEV],
    "LLL": [LEV, LEV, FREE],
    "L": [FREE],
    "LR": [FREE],
    "RL": [FREE],
    "CL": [FREE],
    "CRC": [],
}


@pytest.mark.parametrize("layout", sorted(LAYOUTS))
def test_decoder_labels_follow_the_rule(core, layout):
    labels = core.call("decoder_for_each_list", layout_blocks(layout))
    assert [label for label in labels if label is not None] == LAYOUTS[layout]
    assert [label is not None for label in labels] == [
        letter == "L" for letter in layout
    ]


def test_each_read_is_decided_on_its_own(core):
    labels = core.run(
        "(core, arg) => arg.map((blocks) => core.decoder_for_each_list(blocks))",
        [layout_blocks("LC"), layout_blocks("LR")],
    )
    assert labels == [[LEV, None], [FREE, None]]
