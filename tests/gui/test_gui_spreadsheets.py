"""Tier C9: spreadsheet import through the vendored SheetJS.

Workbooks are written by three independent libraries and must give the same
sequences as the equivalent CSV. The vendored SheetJS must be the pinned,
unmodified release.
"""

import base64
import hashlib
import re

import pytest
from gui_support import GUI_PATH, SAMPLE_BARCODES

SHEETJS_SHA256 = "cc015130aa8521e7f088f88898eba949ccdcbfb38df0bd129b44b7273c3a6f41"
SHEETJS_BYTES = 951_904
TABLE = [["name", "barcode"]] + [
    [f"bc{i + 1}", seq] for i, seq in enumerate(SAMPLE_BARCODES)
]


def write_xlsx(path, sheets):
    from openpyxl import Workbook

    workbook = Workbook()
    workbook.remove(workbook.active)
    for title, rows in sheets.items():
        sheet = workbook.create_sheet(title)
        for row in rows:
            sheet.append(row)
    workbook.save(path)


def write_ods(path, sheets):
    from odf.opendocument import OpenDocumentSpreadsheet
    from odf.table import Table, TableCell, TableRow
    from odf.text import P

    document = OpenDocumentSpreadsheet()
    for title, rows in sheets.items():
        table = Table(name=title)
        for row in rows:
            table_row = TableRow()
            for value in row:
                cell = TableCell(valuetype="string")
                cell.addElement(P(text=value))
                table_row.addElement(cell)
            table.addElement(table_row)
        document.spreadsheet.addElement(table)
    document.save(str(path))


def write_xls(path, sheets):
    import xlwt

    workbook = xlwt.Workbook()
    for title, rows in sheets.items():
        sheet = workbook.add_sheet(title)
        for row_index, row in enumerate(rows):
            for column_index, value in enumerate(row):
                sheet.write(row_index, column_index, value)
    workbook.save(str(path))


WRITERS = {"xlsx": write_xlsx, "ods": write_ods, "xls": write_xls}
READ_SHEET = """(core, arg) => {
    const bytes = Uint8Array.from(atob(arg.data), (ch) => ch.charCodeAt(0));
    const sheets = core.list_spreadsheet_sheets(bytes);
    const rows = core.read_spreadsheet_rows(bytes, arg.sheet || sheets[0]);
    return {sheets, rows, plan: core.plan_table_import(rows, arg.name, {}, 'row')};
}"""


def read_workbook(core, path, sheet=None):
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return core.run(READ_SHEET, {"data": data, "name": path.name, "sheet": sheet})


def test_vendored_sheetjs_is_the_pinned_release():
    page = GUI_PATH.read_text(encoding="utf-8")
    match = re.search(
        r'<script id="vendor_sheetjs">\n(.*?)</script>', page, flags=re.DOTALL
    )
    assert match is not None
    vendored = match.group(1).encode("utf-8")
    assert len(vendored) == SHEETJS_BYTES
    assert hashlib.sha256(vendored).hexdigest() == SHEETJS_SHA256
    assert page.rstrip().endswith("</script>\n</body>\n</html>".rstrip()), (
        "SheetJS must be the last script"
    )


def test_sheetjs_reports_its_version(core):
    assert core.run("() => XLSX.version") == "0.20.3"


@pytest.mark.parametrize("extension", sorted(WRITERS))
def test_spreadsheet_gives_the_same_sequences_as_the_csv(core, tmp_path, extension):
    path = tmp_path / f"barcodes.{extension}"
    WRITERS[extension](path, {"Barcodes": TABLE})
    result = read_workbook(core, path)
    assert result["sheets"] == ["Barcodes"]
    csv_plan = core.run(
        "(core, arg) => core.plan_table_import(core.parse_delimited_text(arg).rows, 'barcodes.csv')",
        "\n".join(",".join(row) for row in TABLE),
    )
    assert result["plan"]["status"] == csv_plan["status"] == "ok"
    assert result["plan"]["sequences"] == csv_plan["sequences"] == SAMPLE_BARCODES
    assert [source["line"] for source in result["plan"]["sources"]] == [2, 3, 4, 5]
    assert result["plan"]["sources"][0]["unit"] == "row"


@pytest.mark.parametrize("extension", sorted(WRITERS))
def test_every_sheet_of_a_workbook_is_listed_and_readable(core, tmp_path, extension):
    other = [["barcode"], ["TTTTCCCC"], ["GGGGAAAA"]]
    path = tmp_path / f"two_sheets.{extension}"
    WRITERS[extension](path, {"Plate A": TABLE, "Plate B": other})
    first = read_workbook(core, path)
    assert first["sheets"] == ["Plate A", "Plate B"]
    second = read_workbook(core, path, sheet="Plate B")
    assert second["plan"]["sequences"] == ["TTTTCCCC", "GGGGAAAA"]


def test_spreadsheet_rows_keep_their_row_numbers_past_blank_rows(core, tmp_path):
    path = tmp_path / "gaps.xlsx"
    write_xlsx(path, {"Barcodes": [["barcode"], ["CCTTAACAT"], [], ["AGCGTAGAA"]]})
    rows = read_workbook(core, path)["rows"]
    assert [row["line"] for row in rows] == [1, 2, 4]
