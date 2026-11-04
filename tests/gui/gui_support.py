"""Helpers shared by the config builder's browser tests.

Expected values are computed here in Python, independently of the page's
JavaScript: a config's canonical text is Python's own json.dumps output.
"""

import copy
import gzip
import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = REPO_ROOT / "examples"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
GUI_PATH = REPO_ROOT / "gui" / "barcode_configurator.html"
DESIGN_DOC_PATH = REPO_ROOT / "doc" / "config_format.md"
READ_KEYS = ("barcode_struct_r1", "barcode_struct_r2")


def example_config_names() -> list[str]:
    """Return the file names of bcwithqc's example configs, sorted."""
    return sorted(path.name for path in EXAMPLES_DIR.glob("*.json"))


def read_bytes_text(path: Path) -> str:
    """Return a config file's text, decompressing .gz files."""
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode("utf-8")


def read_config(path: Path) -> dict:
    """Return a config file parsed as JSON, decompressing .gz files."""
    return json.loads(read_bytes_text(path))


def canonical_text(config: dict) -> str:
    """Return the exact text the GUI must write: json.dumps(indent=4) plus a newline."""
    return json.dumps(config, indent=4) + "\n"


def without_names(config: dict) -> dict:
    """Return a copy of a config with every block's per-sequence "name" list removed."""
    result = copy.deepcopy(config)
    for key in READ_KEYS:
        for block in result.get(key, {}).get("blocks", []):
            block.pop("name", None)
    return result


def with_sorted_lists(config: dict) -> dict:
    """Return a copy of a config with each barcode list sorted A→Z; constant regions unchanged."""
    result = copy.deepcopy(config)
    for key in READ_KEYS:
        for block in result.get(key, {}).get("blocks", []):
            if block["blocktype"] == "barcodeList":
                block["sequence"] = sorted(block["sequence"])
    return result


def design_doc_example() -> dict:
    """Return the example config under "Example Config.json" in doc/config_format.md."""
    text = DESIGN_DOC_PATH.read_text(encoding="utf-8")
    match = re.search(
        r"^## Example Config\.json\s*```json\n(.*?)\n```",
        text,
        flags=re.DOTALL | re.MULTILINE,
    )
    assert match is not None, "doc/config_format.md has no example config"
    return json.loads(match.group(1))


OLD_GUI_BARCODES = ["CCTTAACAT", "AGCGTAGAA", "TATAGGTGT", "GGCGACTCA"]
OLD_GUI_SAMPLE_BARCODES = ["TCGCCTTA", "CTAGTACG", "TTCTGCCT", "GCTCAGGA"]


def old_gui_config() -> dict:
    """Return a config laid out exactly as barcodinggui 2b473d1 saves one.

    barcodinggui writes maxerrors and blockfunction, then the sequence table
    (a "name" list, left out for barcode lists whose names are all empty, then
    "sequence") or the length, then blocktype, blockname and, for linked
    blocks, linkgroup. It writes no keepblock and no single_end_reads, and
    names new blocks "sequence 0", "sequence 1", ...
    """
    return {
        "unknown_read_orientation": True,
        "barcode_struct_r1": {
            "keep_nonbarcode": False,
            "blocks": [
                {
                    "maxerrors": 1,
                    "blockfunction": "barcode",
                    "sequence": OLD_GUI_BARCODES,
                    "blocktype": "barcodeList",
                    "blockname": "sequence 0",
                    "linkgroup": 0,
                },
                {
                    "maxerrors": 1,
                    "blockfunction": "barcode",
                    "name": ["linker"],
                    "sequence": ["GGGTCAGTACGTACGAGTCCC"],
                    "blocktype": "constantRegion",
                    "blockname": "sequence 1",
                },
                {
                    "maxerrors": 1,
                    "blockfunction": "discard",
                    "name": ["s1", "s2", "s3", "s4"],
                    "sequence": OLD_GUI_SAMPLE_BARCODES,
                    "blocktype": "barcodeList",
                    "blockname": "sequence 2",
                    "linkgroup": 0,
                },
                {
                    "maxerrors": 0,
                    "blockfunction": "UMI",
                    "length": 8,
                    "blocktype": "randomBarcode",
                    "blockname": "sequence 3",
                },
            ],
        },
        "barcode_struct_r2": {"keep_nonbarcode": True, "blocks": []},
    }


def write_old_gui_config(path: Path) -> None:
    """Write old_gui_config() the way barcodinggui saves: gzip level 9, json.dump(indent=4), no final newline."""
    with gzip.open(path, "wt", compresslevel=9, encoding="utf8") as handle:
        json.dump(old_gui_config(), handle, indent=4)
