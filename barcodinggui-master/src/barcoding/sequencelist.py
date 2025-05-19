from collections.abc import Sized
from io import StringIO
from typing import Any

import numpy as np
import pandas as pd
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QRegularExpression, Qt, Signal, Slot
from PySide6.QtGui import QBrush, QColor, QIcon, QPalette, QRegularExpressionValidator
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QFileDialog,
    QLineEdit,
    QMessageBox,
    QStyle,
    QStyledItemDelegate,
    QWidget,
)

from .cursormanager import WaitCursor
from .ui.sequencelist import Ui_SequenceList


class countdict(dict):
    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            return 0

    def __setitem__(self, key, value):
        if value == 0:
            try:
                del self[key]
                return
            except KeyError:
                pass
        super().__setitem__(key, value)


class histdict(countdict):
    def __getitem__(self, key):
        if isinstance(key, Sized):
            key = len(key)
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        if isinstance(key, Sized):
            key = len(key)
        super().__setitem__(key, value)

    def mode(self):
        maxkey, maxval = -1, -1
        for k, v in self.items():
            if v > maxval:
                maxkey, maxval = k, v
        return maxkey, maxval


class SequenceListModel(QAbstractTableModel):
    validChanged = Signal(bool, arguments=["valid"])
    variableSequenceLengthChanged = Signal(bool, arguments=["variablelengths"])
    duplicateNamesChanged = Signal(bool, arguments=["duplicates"])
    duplicateSequencesChanged = Signal(bool, arguments=["duplicates"])
    emptySequencesChanged = Signal(bool, arguments=["emptysequences"])
    mixedNamedSquencesChanged = Signal(bool, arguments=["mixednamedsequences"])
    minLengthChanged = Signal(int, arguments=["minlength"])

    lengthBgColor = QColor(166, 166, 89)
    duplicateBgColor = QColor(161, 48, 52)

    def __init__(self, allowVariableLength=False, allowUnnamedSequences=False, parent=None):
        super().__init__(parent)
        self._allowVariableLength = allowVariableLength
        self._allowUnnamedSequences = allowUnnamedSequences
        self.initialize()

    def initialize(self):
        self._tbl = pd.DataFrame(columns=["name", "sequence"])
        self._duplicates = (countdict(), countdict())
        self._sequenceLengthHist = histdict()
        self._mostAbundantLength = -1
        self._minlength = np.inf
        self._valid = False

    @property
    def allowVariableLength(self):
        return self._allowVariableLength

    @property
    def allowUnnamedSequences(self):
        return self._allowUnnamedSequences

    def headerData(self, section, orientation, role):
        if role != Qt.DisplayRole or orientation != Qt.Horizontal:
            return None
        return self._tbl.columns[section]

    def columnCount(self, parent=QModelIndex()):
        return self._tbl.shape[1]

    def rowCount(self, parent=QModelIndex()):
        return self._tbl.shape[0] + 1

    def flags(self, index):
        return Qt.ItemIsSelectable | Qt.ItemIsEditable | Qt.ItemIsEnabled

    def setData(self, index, value, role):
        haveEdit = True
        if role == Qt.EditRole and index.column() < self._tbl.shape[1]:
            if index.column() > 0:
                value = value.upper()
            if index.row() < self._tbl.shape[0]:
                oldval = self._tbl.iloc[index.row(), index.column()]
                self._duplicates[index.column()][oldval] -= 1
                self._tbl.iloc[index.row(), index.column()] = value
                self._duplicates[index.column()][value] += 1
                if index.column() > 0:
                    if not self._allowVariableLength:
                        self._sequenceLengthHist[oldval] -= 1
                        self._sequenceLengthHist[value] += 1
                    if len(value) != len(oldval):
                        if not self._allowVariableLength:
                            if self._sequenceLengthHist[value] > self._sequenceLengthHist[self._mostAbundantLength]:
                                self._mostAbundantLength = len(value)
                            else:
                                self._mostAbundantLength = self._sequenceLengthHist.mode()[0]
                        if len(value) < self._minlength:
                            self._minlength = len(value)
                            self.minLengthChanged.emit(self._minlength)
                        elif self._sequenceLengthHist[self._minlength] == 0:
                            self._minlength = min(self._sequenceLengthHist.keys())
                            self.minLengthChanged.emit(self._minlength)
            elif index.row() == self._tbl.shape[0] and len(value):
                self.beginInsertRows(QModelIndex(), self._tbl.shape[0] + 1, self._tbl.shape[0] + 1)
                newrow = ["", ""]
                newrow[index.column()] = value
                self._tbl.loc[index.row(), :] = newrow
                for col, val in enumerate(newrow):
                    self._duplicates[col][val] += 1
                if not self._allowVariableLength:
                    self._sequenceLengthHist[newrow[1]] += 1
                    if self._sequenceLengthHist[newrow[1]] > self._sequenceLengthHist[self._mostAbundantLength]:
                        self._mostAbundantLength = len(newrow[1])
                if len(value) < self._minlength:
                    self._minlength = len(value)
                    self.minLengthChanged.emit(self._minlength)
                self.endInsertRows()
            else:
                haveEdit = False
        if haveEdit:
            self._updateValid()
            self.dataChanged.emit(index, index)
            return True
        else:
            return False

    def appendDataFrame(self, df):
        if df.shape[1] == 1:
            df.insert(0, None, "")
        df.columns = self._tbl.columns
        self.beginInsertRows(QModelIndex(), self._tbl.shape[0] + 1, self._tbl.shape[0] + df.shape[0])
        self._tbl = pd.concat((self._tbl, df), axis=0, ignore_index=True)

        seqlengths, seqlengthcounts = np.unique(df.iloc[:, 1].str.len().to_numpy(), return_counts=True)
        if not self._allowVariableLength:
            for length, count in zip(seqlengths, seqlengthcounts, strict=True):
                self._sequenceLengthHist[length] += count
                if (
                    length != self._mostAbundantLength
                    and self._sequenceLengthHist[length] > self._sequenceLengthHist[self._mostAbundantLength]
                ):
                    self._mostAbundantLength = length
        minlength = np.min(seqlengths)
        if minlength < self._minlength:
            self._minlength = minlength
            self.minLengthChanged.emit(self._minlength)

        for col in (0, 1):
            for value, count in zip(*np.unique(df.iloc[:, col].to_numpy(), return_counts=True), strict=True):
                self._duplicates[col][value] += count

        self.endInsertRows()
        self._updateValid()

    def removeRowList(self, rows):
        indices = np.empty((len(rows),), int)
        for i, row in enumerate(rows):
            indices[i] = row.row()
        indices = np.sort(indices[indices < self._tbl.shape[0]])
        ranges = np.concatenate(((0,), np.cumsum(np.diff(indices) > 1)))
        uniques, indices2, counts = np.unique(ranges, return_index=True, return_counts=True)
        startindices = indices[indices2]
        startindices[1:] -= np.cumsum(counts)[:-1]
        for startidx, nrows in zip(startindices, counts, strict=True):
            self.beginRemoveRows(QModelIndex(), startidx, startidx + nrows - 1)
        for name, seq in self._tbl.iloc[indices, :].itertuples(index=False):
            self._duplicates[0][name] -= 1
            self._duplicates[1][seq] -= 1
            if not self._allowVariableLength:
                self._sequenceLengthHist[seq] -= 1
        self._tbl.drop(index=indices, inplace=True)
        self._tbl.reset_index(drop=True, inplace=True)
        if not self._allowVariableLength:
            self._mostAbundantLength = self._sequenceLengthHist.mode()[0]
        if not self._sequenceLengthHist[self._minlength]:
            self._minlength = min(self._sequenceLengthHist.keys())
            self.minLengthChanged.emit(self._minlength)
        self.endRemoveRows()
        self._updateValid()

    def data(self, index, role):
        if index.row() < self._tbl.shape[0] and index.column() < self._tbl.shape[1]:
            if role == Qt.DisplayRole or role == Qt.EditRole:
                return self._tbl.iloc[index.row(), index.column()]
            elif role == Qt.BackgroundRole:
                if self._duplicates[index.column()][self._tbl.iloc[index.row(), index.column()]] > 1 and (
                    index.column() > 0
                    or index.column() == 0
                    and (
                        not self._allowUnnamedSequences or len(self._duplicates[0]) > 1 or "" not in self._duplicates[0]
                    )
                ):
                    return QBrush(self.duplicateBgColor)
                elif (
                    not self._allowVariableLength
                    and index.column() > 0
                    and len(self._tbl.iloc[index.row(), index.column()]) != self._mostAbundantLength
                ):
                    return QBrush(self.lengthBgColor)
        return None

    def _updateValid(self):
        nameduplicates = len(self._duplicates[0]) != self._tbl.shape[0] and (
            not self._allowUnnamedSequences or len(self._duplicates[0]) != 1 or "" not in self._duplicates[0]
        )
        sequenceduplicates = len(self._duplicates[1]) != self._tbl.shape[0]
        emptysequences = "" in self._duplicates[1]
        valid = self._tbl.shape[0] > 0 and not nameduplicates and not sequenceduplicates and not emptysequences
        self.duplicateNamesChanged.emit(nameduplicates)
        self.duplicateSequencesChanged.emit(sequenceduplicates)
        self.emptySequencesChanged.emit(emptysequences)
        if not self._allowVariableLength:
            samelengths = len(self._sequenceLengthHist) == 1
            valid = valid and samelengths and list(self._sequenceLengthHist.keys())[0] > 0
            self.variableSequenceLengthChanged.emit(not samelengths)
        if self._allowUnnamedSequences:
            mixednamedsequences = "" in self._duplicates[0] and len(self._duplicates[0]) > 1
            valid = valid and not mixednamedsequences
            self.mixedNamedSquencesChanged.emit(mixednamedsequences)
        else:
            valid = valid and "" not in self._duplicates[0]
        if valid != self._valid:
            self._valid = valid
            self.validChanged.emit(valid)

    def valid(self):
        return self._valid

    def to_dict(self):
        ret = self._tbl.to_dict("list")
        if self._allowUnnamedSequences and len(self._duplicates[0]) == 1 and "" in self._duplicates[0]:
            del ret["name"]
        return ret

    def from_dict(self, d: dict[str, Any]):
        subd = {}
        if "name" in d or not self._allowUnnamedSequences:
            subd["name"] = d["name"]
        subd["sequence"] = d["sequence"]
        df = pd.DataFrame(subd)  # column order follows insertion order

        self.beginResetModel()
        self.initialize()
        self.appendDataFrame(df)
        self.endResetModel()


class SequenceTableDelegate(QStyledItemDelegate):
    def __init__(self, validator=None, parent=None):
        super().__init__(parent)
        self._val = validator

    def createEditor(self, parent, option, index):
        lineedit = QLineEdit(parent)
        lineedit.setValidator(self._val)
        return lineedit

    def setEditorData(self, editor, index):
        editor.setText(index.data(Qt.EditRole))

    def setModelData(self, editor, model, index):
        model.setData(index, editor.text(), Qt.EditRole)


class SequenceDelegate(SequenceTableDelegate):
    def __init__(self, parent=None):
        super().__init__(
            QRegularExpressionValidator(QRegularExpression("[atgc]+", QRegularExpression.CaseInsensitiveOption)), parent
        )


class NameDelegate(SequenceTableDelegate):
    def __init__(self, parent=None):
        super().__init__(QRegularExpressionValidator(QRegularExpression("\\w+")), parent)


class SequenceList(QWidget):
    validChanged = Signal(bool, arguments=["valid"])
    lengthChanged = Signal(int, arguments=["length"])
    edited = Signal()

    def __init__(self, allowVariableLength=False, allowUnnamedSequences=False, parent=None):
        super().__init__(parent)
        self._ui = Ui_SequenceList()
        self._ui.setupUi(self)
        self._ui.openBtn.setIcon(
            QIcon.fromTheme("document-open", QApplication.style().standardIcon(QStyle.SP_DialogOpenButton))
        )
        if QIcon.hasThemeIcon("edit-paste"):
            self._ui.pasteBtn.setIcon(QIcon.fromTheme("edit-paste"))
        else:
            self._ui.pasteBtn.setText("Paste")

        sp = self._ui.sequenceTbl.sizePolicy()
        self._ui.sequenceTbl.setSizePolicy(sp)

        self._sequences = SequenceListModel(allowVariableLength, allowUnnamedSequences, self)

        palette = QPalette()
        palette.setColor(QPalette.Window, self._sequences.lengthBgColor)
        self._ui.seqlenLbl.setPalette(palette)
        palette.setColor(QPalette.Window, self._sequences.duplicateBgColor)
        self._ui.duplicateSeqLbl.setPalette(palette)
        self._ui.duplicateNameLbl.setPalette(palette)
        lbls = [
            self._ui.duplicateSeqLbl,
            self._ui.duplicateNameLbl,
            self._ui.unnamedSeqLbl,
            self._ui.emptySeqLbl,
            self._ui.seqlenLbl,
        ]
        for lbl in lbls:
            lbl.setAutoFillBackground(True)
            lbl.hide()
        if allowVariableLength:
            del lbls[-1]
        if not allowUnnamedSequences:
            del lbls[2]
        for lbl in lbls:
            sp = lbl.sizePolicy()
            sp.setRetainSizeWhenHidden(True)
            lbl.setSizePolicy(sp)

        self._ui.sequenceTbl.setModel(self._sequences)
        self._ui.sequenceTbl.setItemDelegateForColumn(0, NameDelegate(self))
        self._ui.sequenceTbl.setItemDelegateForColumn(1, SequenceDelegate(self))
        self._ui.sequenceTbl.selectionModel().selectionChanged.connect(self.selectionChanged)
        self._ui.deleteBtn.clicked.connect(self.deleteBtnClicked)
        self._ui.openBtn.clicked.connect(self._openBtnClicked)
        self._ui.pasteBtn.clicked.connect(self._pasteBtnClicked)
        self._sequences.validChanged.connect(self.validChanged)
        self._sequences.duplicateNamesChanged.connect(self._ui.duplicateNameLbl.setVisible)
        self._sequences.duplicateSequencesChanged.connect(self._ui.duplicateSeqLbl.setVisible)
        self._sequences.variableSequenceLengthChanged.connect(self._ui.seqlenLbl.setVisible)
        self._sequences.mixedNamedSquencesChanged.connect(self._ui.unnamedSeqLbl.setVisible)
        self._sequences.emptySequencesChanged.connect(self._ui.emptySeqLbl.setVisible)
        self._sequences.dataChanged.connect(self.edited)
        self._sequences.rowsInserted.connect(self.edited)
        self._sequences.modelReset.connect(self.edited)
        self._sequences.minLengthChanged.connect(self.lengthChanged)

    @property
    def valid(self):
        return self._sequences.valid()

    @Slot()
    def selectionChanged(self):
        rows = self._ui.sequenceTbl.selectionModel().selectedRows()
        selection = len(rows) > 0
        if len(rows) == 1 and rows[0].row() == self._sequences.rowCount() - 1:
            selection = False
        self._ui.deleteBtn.setEnabled(selection)

    @Slot()
    def deleteBtnClicked(self):
        self._sequences.removeRowList(self._ui.sequenceTbl.selectionModel().selectedRows())

    def keyPressEvent(self, event):
        if (
            (event.key() == Qt.Key_Enter or event.key() == Qt.Key_Return)
            and self._ui.sequenceTbl.state() != QAbstractItemView.EditingState
            and self._ui.sequenceTbl.currentIndex().isValid()
        ):
            self._ui.sequenceTbl.edit(self._ui.sequenceTbl.currentIndex())
        else:
            super().keyPressEvent(event)

    def _importTable(self, table):
        if table.shape[0] < 1:
            QMessageBox.critical(
                self, "Could not import barcodes", "The selected file is empty and contains no barcodes."
            )
        elif table.shape[1] < 2 and not self._sequences.allowUnnamedSequences:
            QMessageBox.critical(
                self,
                "Could not import barcodes",
                "The selected file contains only one column, at least two columns are required.",
            )
        else:
            barcodecol = None
            namecol = None
            teststring = lambda x: isinstance(x, str) and len(set(x.upper()) - nucleotides) == 0
            candidatecols = []
            nucleotides = {"A", "T", "G", "C"}
            stringcols = np.where(table.dtypes == "O")[0]
            for col, val in zip(range(len(stringcols)), table.iloc[-1, stringcols], strict=True):
                if teststring(val):
                    candidatecols.append(col)
            for scol in candidatecols:
                col = stringcols[scol]
                n_nondna = (~table.iloc[:, col].transform(teststring)).sum()
                if n_nondna == 1 and not teststring(table.iloc[0, col]):
                    barcodecol = scol
                    table.columns = table.iloc[0, :]
                    table = table.iloc[1:, :]
                    break
                elif n_nondna == 0:
                    barcodecol = scol
                    break
            if barcodecol is None:
                QMessageBox.critical(
                    self,
                    "Could not import barcodes",
                    "The selected file does not contain any columns of DNA sequences.",
                )
            else:
                namecol = None
                namecolnotna = np.inf
                if (
                    stringcols[barcodecol] > 0
                    and (notna := pd.isnull(table.iloc[:, stringcols[barcodecol] - 1]).sum()) < table.shape[0]
                ):
                    namecol = stringcols[barcodecol] - 1
                    namecolnotna = notna
                if (
                    stringcols[barcodecol] < table.shape[1] - 1
                    and (notna := pd.isnull(table.iloc[:, stringcols[barcodecol] + 1]).sum()) < table.shape[0]
                    and notna < namecolnotna
                ):
                    namecol = stringcols[barcodecol] + 1
                    namecolnotna = notna
                if namecol is None:
                    if barcodecol > 0:
                        namecol = stringcols[barcodecol - 1]
                    elif len(stringcols) > 1:
                        namecol = stringcols[barcodecol + 1]

                barcodecol = stringcols[barcodecol]
                table.iloc[:, barcodecol] = table.iloc[:, barcodecol].str.upper()
                if namecol is not None:
                    table.iloc[:, namecol].fillna("", inplace=True)
                slice = [namecol, barcodecol] if namecol is not None else [barcodecol]
                self._sequences.appendDataFrame(table.iloc[:, slice])

    @Slot()
    def _openBtnClicked(self):
        fname = QFileDialog.getOpenFileName(
            self, "Import barcodes", filter="Tables (*.csv *.tab *.txt *.xls *.xlsx, *.ods);;All files (*)"
        )[0]
        if not len(fname):
            return
        try:
            with WaitCursor():
                if fname.endswith(".xls") or fname.endswith(".xlsx") or fname.endswith(".ods"):
                    table = pd.read_excel(fname, header=None)
                else:
                    table = pd.read_csv(fname, header=None, index_col=None, sep=",|;|:|\\s+")
            self._importTable(table)
        except Exception as e:  # noqa BLE001
            QMessageBox.critical(self, "Could not import barcodes", str(e))
            return

    @Slot()
    def _pasteBtnClicked(self):
        try:
            with WaitCursor():
                clpbrd = StringIO(QApplication.clipboard().text())
                table = pd.read_csv(clpbrd, header=None, index_col=None, sep=",|;|:|\\s+")
                self._importTable(table)
        except Exception as e:  # noqa BLE001
            QMessageBox.critical(self, "Could not import barcodes", str(e))
            return

    def to_dict(self):
        return self._sequences.to_dict()

    def from_dict(self, d: dict[str, Any]):
        self._sequences.from_dict(d)
