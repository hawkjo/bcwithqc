"""Write the config builder's template data from bcwithqc's own configs.

The templates are the shipped example configs without per-sequence name
lists. Barcode lists of 100 or more entries are stored once under "lists" and
referenced as {"$list": key}, so the 1,536-barcode SDR-seq list appears once in
the page.

Run from the repository root after changing a template's source:

    python gui/update_template_data.py

The GUI tests check that the page's template data matches these sources.
"""

import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PAGE = REPO / "gui" / "barcode_configurator.html"

# (template id, title shown in the app, source path relative to the repo root)
TEMPLATE_SOURCES = [
    ("sdr_seq_rna", "SDR-seq RNA", "examples/cDNA_paired_end_config.json"),
    ("sdr_seq_gdna", "SDR-seq gDNA", "examples/gDNA_paired_end_config.json"),
    ("design_doc", "Design-doc example", "doc/config_format.md"),
]
SHARED_LIST_MIN_LENGTH = 100
DATA_BLOCK = re.compile(
    r'(<script type="application/json" id="template_data">\n)(.*?)(\n</script>)',
    re.DOTALL,
)


def design_doc_example(doc_path: Path) -> dict:
    """Return the example config from the design doc.

    Args:
        doc_path: Path to doc/config_format.md.

    Returns:
        The JSON block under the "Example Config.json" heading, parsed.
    """
    text = doc_path.read_text(encoding="utf-8")
    match = re.search(
        r"^## Example Config\.json\s*```json\n(.*?)\n```",
        text,
        flags=re.DOTALL | re.MULTILINE,
    )
    if match is None:
        raise ValueError(f"No example config found in {doc_path}")
    return json.loads(match.group(1))


def template_config(source: str) -> dict:
    """Return a template's config: its source without name lists.

    Args:
        source: Path of the source, relative to the repository root.

    Returns:
        The config as bcwithqc reads it, minus every block's "name" key.
    """
    path = REPO / source
    config = (
        design_doc_example(path)
        if path.suffix == ".md"
        else json.loads(path.read_text(encoding="utf-8"))
    )
    for key in ("barcode_struct_r1", "barcode_struct_r2"):
        for block in config.get(key, {}).get("blocks", []):
            block.pop("name", None)
    return config


def template_data() -> dict:
    """Return the page's template data, with long lists stored once.

    Returns:
        {"lists": {key: [sequences]}, "templates": [{id, title, source, config}]},
        where long sequence lists inside the configs are {"$list": key}.
    """
    shared_lists: dict[str, list[str]] = {}
    templates = []
    for template_id, title, source in TEMPLATE_SOURCES:
        config = template_config(source)
        for key in ("barcode_struct_r1", "barcode_struct_r2"):
            for block in config.get(key, {}).get("blocks", []):
                sequence = block.get("sequence")
                if sequence is None or len(sequence) < SHARED_LIST_MIN_LENGTH:
                    continue
                list_key = next(
                    (
                        name
                        for name, stored in shared_lists.items()
                        if stored == sequence
                    ),
                    None,
                )
                if list_key is None:
                    list_key = f"list_{len(shared_lists) + 1}"
                    shared_lists[list_key] = sequence
                block["sequence"] = {"$list": list_key}
        templates.append(
            {"id": template_id, "title": title, "source": source, "config": config}
        )
    return {"lists": shared_lists, "templates": templates}


def main() -> None:
    """Replace the template data block in the page."""
    page = PAGE.read_text(encoding="utf-8")
    data = json.dumps(template_data(), separators=(",", ":"))
    if DATA_BLOCK.search(page) is None:
        raise ValueError(f"No template data block in {PAGE}")
    PAGE.write_text(
        DATA_BLOCK.sub(lambda match: match.group(1) + data + match.group(3), page),
        encoding="utf-8",
    )
    print(f"Wrote {len(data):,} bytes of template data to {PAGE.relative_to(REPO)}")


if __name__ == "__main__":
    main()
