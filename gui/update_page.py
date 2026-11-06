"""Refresh the parts of the config builder that are copied from the repository.

Two parts of gui/barcode_configurator.html come from elsewhere in the repo:

- The "generated_with" stamp: the tool's name and version, which the page
  writes first in every config. Both are read from pyproject.toml: the name
  from [project] name, and the version from [project] version or, when the
  version is dynamic, from the file [tool.hatch.version] points to.
- The templates: the shipped example configs without per-sequence name lists.
  Barcode lists of 100 or more entries are stored once under "lists" and
  referenced as {"$list": key}, so the 1,536-barcode SDR-seq list appears once
  in the page.

Run from the repository root after a version bump, a rename, or a change to a
template's source:

    python gui/update_page.py

The GUI tests fail while the page is out of date with these sources.
"""

import json
import re
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

REPO = Path(__file__).resolve().parent.parent
PAGE = REPO / "gui" / "barcode_configurator.html"
PYPROJECT = REPO / "pyproject.toml"
STAMP_LINE = re.compile(r"^  const GENERATED_WITH = .*;$", re.MULTILINE)
VERSION_ASSIGNMENT = re.compile(
    r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", flags=re.MULTILINE
)
# Names and versions are written into a JavaScript string, so they are limited
# to the characters package names and PEP 440 versions use.
STAMP_TOKEN = re.compile(r"^[A-Za-z0-9._+!-]+$")

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


def generated_with() -> dict:
    """Return the tool's name and version as pyproject.toml gives them.

    Returns:
        {"name": ..., "version": ...}.

    Raises:
        ValueError: If pyproject.toml gives no version, or a name or version
            has characters other than letters, digits and ._+!-.
    """
    pyproject = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    project = pyproject["project"]
    if "version" in project:
        version = project["version"]
    else:
        version_file = REPO / pyproject["tool"]["hatch"]["version"]["path"]
        match = VERSION_ASSIGNMENT.search(version_file.read_text(encoding="utf-8"))
        if match is None:
            raise ValueError(f"{version_file} sets no __version__")
        version = match.group(1)
    stamp = {"name": project["name"], "version": version}
    for value in stamp.values():
        if not STAMP_TOKEN.match(value):
            raise ValueError(f"Can't write {value!r} into the page's stamp")
    return stamp


def stamp_line(stamp: dict) -> str:
    """Return the line of the page's core script that holds the stamp.

    Args:
        stamp: {"name": ..., "version": ...}, as generated_with() returns it.

    Returns:
        The JavaScript constant GENERATED_WITH, indented as in the page.
    """
    return (
        "  const GENERATED_WITH = Object.freeze("
        f"{{ name: '{stamp['name']}', version: '{stamp['version']}' }});"
    )


def main() -> None:
    """Write the stamp and the template data into the page."""
    page = PAGE.read_text(encoding="utf-8")
    if STAMP_LINE.search(page) is None:
        raise ValueError(f"No GENERATED_WITH line in {PAGE}")
    if DATA_BLOCK.search(page) is None:
        raise ValueError(f"No template data block in {PAGE}")
    stamp = generated_with()
    page = STAMP_LINE.sub(lambda match: stamp_line(stamp), page)
    data = json.dumps(template_data(), separators=(",", ":"))
    page = DATA_BLOCK.sub(lambda match: match.group(1) + data + match.group(3), page)
    PAGE.write_text(page, encoding="utf-8")
    print(
        f"Stamped {stamp['name']} {stamp['version']} and wrote {len(data):,} bytes "
        f"of template data to {PAGE.name}"
    )


if __name__ == "__main__":
    main()
