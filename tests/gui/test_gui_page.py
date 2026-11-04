"""The page as a file: self-contained, UTF-8, and loading cleanly in every engine."""

import re

from gui_support import GUI_PATH

# The vendored SheetJS script is checked by its sha256 instead (tier C9).
VENDOR_SCRIPT = re.compile(r'<script id="vendor_sheetjs">.*?</script>', re.DOTALL)
EXTERNAL_REFERENCES = [
    r"\bsrc\s*=",
    r"<link\b",
    r"href\s*=\s*[\"'](?!#)",
    r"url\(",
    r"\bfetch\(",
    r"\bimport\(",
    r"XMLHttpRequest",
    r"sendBeacon",
    r"WebSocket",
    r"EventSource",
    r"@import",
]


def our_code():
    return VENDOR_SCRIPT.sub("", GUI_PATH.read_text(encoding="utf-8"))


def test_charset_is_declared_within_the_first_1024_bytes():
    head = GUI_PATH.read_bytes()[:1024]
    assert b'<meta charset="utf-8">' in head


def test_page_is_valid_utf8():
    GUI_PATH.read_bytes().decode("utf-8")


def test_our_code_references_nothing_outside_the_file():
    code = our_code()
    found = {pattern: re.findall(pattern, code) for pattern in EXTERNAL_REFERENCES}
    assert not {pattern: hits for pattern, hits in found.items() if hits}


def test_page_loads_without_errors_and_reports_ready(app):
    assert app.evaluate("document.getElementById('app').dataset.busy") == "false"
    assert app.evaluate("typeof window.barcode_core.config_to_json") == "function"
    assert app.locator("#feature_banner").is_hidden()
