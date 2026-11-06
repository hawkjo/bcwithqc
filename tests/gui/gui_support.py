"""Helpers shared by the config builder's browser tests.

Expected values are computed here in Python, independently of the page's
JavaScript: a config's canonical text is Python's own json.dumps output, and
the bcwithqc version it records is read from bcwithqc/__init__.py.
"""

import copy
import functools
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


@functools.cache
def bcwithqc_version() -> str:
    """Return bcwithqc's version as bcwithqc/__init__.py sets it.

    Read on first use rather than at import, so a problem here fails the
    tests that need the version instead of stopping pytest's collection.
    """
    text = (REPO_ROOT / "bcwithqc" / "__init__.py").read_text(encoding="utf-8")
    match = re.search(
        r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", text, flags=re.MULTILINE
    )
    assert match is not None, "bcwithqc/__init__.py sets no __version__"
    return match.group(1)


def canonical_text(config: dict) -> str:
    """Return the exact text the GUI must write for a config.

    That is json.dumps(indent=4) plus a newline, with "bcwithqc_version" as
    the first key, set to bcwithqc's version whatever the config says.
    """
    written = {"bcwithqc_version": bcwithqc_version()}
    written.update(
        (key, value) for key, value in config.items() if key != "bcwithqc_version"
    )
    return json.dumps(written, indent=4) + "\n"


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


# Small hand-built configs. Layout letters: L barcode list, C constant region, R random block.
SAMPLE_BARCODES = ["CCTTAACAT", "AGCGTAGAA", "TATAGGTGT", "GGCGACTCA"]
SAMPLE_CONSTANT = "GGGTCAGTACGTACGAGTCCC"


def barcode_list(name, sequences=None, maxerrors=1, function="barcode", keep=False):
    """Return a barcodeList block in bcwithqc's key order."""
    return {
        "maxerrors": maxerrors,
        "blockfunction": function,
        "sequence": list(SAMPLE_BARCODES if sequences is None else sequences),
        "blocktype": "barcodeList",
        "blockname": name,
        "keepblock": keep,
    }


def constant_region(name, sequences=None, maxerrors=1, function="barcode", keep=False):
    """Return a constantRegion block in bcwithqc's key order."""
    return {
        "maxerrors": maxerrors,
        "blockfunction": function,
        "sequence": list([SAMPLE_CONSTANT] if sequences is None else sequences),
        "blocktype": "constantRegion",
        "blockname": name,
        "keepblock": keep,
    }


def random_block(name, length=8, maxerrors=0, function="UMI", keep=False):
    """Return a randomBarcode block in bcwithqc's key order."""
    return {
        "maxerrors": maxerrors,
        "blockfunction": function,
        "length": length,
        "blocktype": "randomBarcode",
        "blockname": name,
        "keepblock": keep,
    }


def make_config(
    read_1,
    read_2=None,
    single_end=None,
    orientation=False,
    keep_r1=False,
    keep_r2=True,
):
    """Return a config: single-end unless read_2 is given (a list, possibly empty)."""
    if single_end is None:
        single_end = read_2 is None
    config = {
        "unknown_read_orientation": orientation,
        "single_end_reads": single_end,
        "barcode_struct_r1": {"keep_nonbarcode": keep_r1, "blocks": list(read_1)},
    }
    if read_2 is not None:
        config["barcode_struct_r2"] = {
            "keep_nonbarcode": keep_r2,
            "blocks": list(read_2),
        }
    return config


def layout_blocks(layout):
    """Return blocks for a layout string such as "LRCL", named L1, R1, C1, L2, ..."""
    counts = {"L": 0, "C": 0, "R": 0}
    blocks = []
    for letter in layout:
        counts[letter] += 1
        name = f"{letter}{counts[letter]}"
        if letter == "L":
            blocks.append(barcode_list(name))
        elif letter == "C":
            blocks.append(constant_region(name))
        else:
            blocks.append(random_block(name))
    return blocks
