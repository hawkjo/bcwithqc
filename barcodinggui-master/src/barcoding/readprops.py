from enum import IntEnum
from importlib.resources import files
from typing import Any

import numpy as np
import pandas as pd
from PySide6.QtCore import QAbstractListModel, QModelIndex, QObject, Qt, Signal, Slot
from PySide6.QtGui import QAction
from PySide6.QtWidgets import QMenu, QMessageBox, QWidget

from .sequenceprops import sequencePropertiesFactory
from .sequencetypes import SequenceType
from .ui.readprops import Ui_ReadProps


class SequenceBlockRoles(IntEnum):
    stype_name = Qt.UserRole
    stype_display_name = Qt.UserRole + 1
    stype = Qt.UserRole + 2
    blockid = Qt.UserRole + 3
    blockindex = Qt.UserRole + 4
    linkgroup = Qt.UserRole + 5


class SequenceBlockModel(QAbstractListModel):
    linkGroupChanged = Signal(list, arguments=["blockids"])

    def __init__(self, parent=None):
        super().__init__(parent)
        self.clear()

    def clear(self):
        self.beginResetModel()
        self._sequenceBlocks = pd.DataFrame(columns=["name", "type", "linkgroup"])
        self._n_added = 0
        self._n_linkgroups = 0
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return self._sequenceBlocks.shape[0]

    def addBlock(self, blockname: str, blocktype: SequenceType, linkgroup: int = -1):
        blockid = self._n_added + 1
        if blockid in self._sequenceBlocks.index:
            raise ValueError("blockid must be unique")
        self.beginInsertRows(QModelIndex(), self.rowCount(), self.rowCount())
        self._sequenceBlocks.loc[blockid] = blockname, blocktype, linkgroup
        self._n_added = blockid
        self.endInsertRows()
        if linkgroup >= 0:
            linkgrp = np.flatnonzero(self._sequenceBlocks.linkgroup == linkgroup)
            self.linkGroupChanged.emit(linkgrp)
        return blockid

    def removeBlock(self, blockid: int):
        row = self._sequenceBlocks.index.get_loc(blockid)
        self.beginRemoveRows(QModelIndex(), row, row)
        self._sequenceBlocks.drop(index=blockid, inplace=True)
        self.endRemoveRows()
        return True

    def removeRow(self, row: int, parent=QModelIndex()):
        return self.removeBlock(self._sequenceBlocks.index[row])

    def removeRows(self, row: int, count: int, parent=QModelIndex):
        blocks = self._sequenceBlocks.index[row : (row + count)]
        self.beginRemoveRows(QModelIndex(), row, row + count - 1)
        self._sequenceBlocks.drop(index=blocks, inplace=True)
        self.endRemoveRows()
        return True

    def setBlockNameById(self, blockid: int, name: str):
        row = self._sequenceBlocks.index.get_loc(blockid)
        self.setBlockNameByRow(row, name)

    def setBlockNameByRow(self, row: int, name: str):
        self._sequenceBlocks.iloc[row, 0] = name
        self.dataChanged.emit(self.index(row, 0), self.index(row, 0), [Qt.DisplayRole])

    def dataForRow(self, row: int, role: int = Qt.DisplayRole):
        match role:
            case Qt.DisplayRole:
                return self._sequenceBlocks.iloc[row, 0]
            case SequenceBlockRoles.stype_display_name:
                return SequenceType(self._sequenceBlocks.iloc[row, 1]).display_name
            case SequenceBlockRoles.stype_name:
                return SequenceType(self._sequenceBlocks.iloc[row, 1]).name
            case SequenceBlockRoles.stype:
                return SequenceType(self._sequenceBlocks.iloc[row, 1])
            case SequenceBlockRoles.blockid:
                return self._sequenceBlocks.index[row]
            case SequenceBlockRoles.blockindex:
                return row
            case SequenceBlockRoles.linkgroup:
                return self._sequenceBlocks.iloc[row, 2]
            case _:
                return None

    def dataForId(self, blockid: int, role: int = Qt.DisplayRole):
        row = self._sequenceBlocks.index.get_loc(blockid)
        return self.dataForRow(row, role)

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole):
        if not index.isValid():
            return None
        row = index.row()
        return self.dataForRow(row, role)

    def linkBlocksByRow(self, row1, row2):
        linkgroup1 = self._sequenceBlocks.linkgroup.iloc[row1]
        linkgroup2 = self._sequenceBlocks.linkgroup.iloc[row2]

        if linkgroup1 < 0 and linkgroup2 < 0:
            newgroup = self._n_linkgroups
            self._n_linkgroups += 1
            self._sequenceBlocks.iloc[row1, 2] = self._sequenceBlocks.iloc[row2, 2] = newgroup
            self.linkGroupChanged.emit([row1, row2])
        elif linkgroup1 < 0:
            self._sequenceBlocks.iloc[row1, 2] = linkgroup2
            self.linkGroupChanged.emit(np.flatnonzero(self._sequenceBlocks.linkgroup == linkgroup2))
        elif linkgroup2 < 0:
            self._sequenceBlocks.iloc[row2, 2] = linkgroup1
            self.linkGroupChanged.emit(np.flatnonzero(self._sequenceBlocks.linkgroup == linkgroup1))
        else:
            newgroup = max(linkgroup1, linkgroup2)
            mergedidx = np.flatnonzero(
                (self._sequenceBlocks.linkgroup == linkgroup1) | (self._sequenceBlocks.linkgroup == linkgroup2)
            )
            self._sequenceBlocks.iloc[mergedidx, 2] = newgroup
            self.linkGroupChanged.emit(mergedidx)
        self.dataChanged.emit(self.index(row1, 0), self.index(row1, 0), [SequenceBlockRoles.linkgroup])
        self.dataChanged.emit(self.index(row2, 0), self.index(row2, 0), [SequenceBlockRoles.linkgroup])

    def linkBlocksById(self, blockid1, blockid2):
        row1 = self._sequenceBlocks.get_loc(blockid1)
        row2 = self._sequenceBlocks.get_loc(blockid2)
        self.linkBlocksByRow(row1, row2)

    def unlinkBlocksByRow(self, triggerRow, targetRow):
        linkgrp = np.flatnonzero(self._sequenceBlocks.linkgroup == self._sequenceBlocks.linkgroup.iloc[targetRow])
        self._sequenceBlocks.iloc[targetRow, 2] = -1
        self.linkGroupChanged.emit(linkgrp)
        self.dataChanged.emit(self.index(targetRow, 0), self.index(targetRow, 0), [SequenceBlockRoles.linkgroup])

    def unlinkBlocksById(self, triggerBlockId, targetBlockId):
        row1 = self._sequenceBlocks.get_loc(triggerBlockId)
        row2 = self._sequenceBlocks.get_loc(targetBlockId)
        self.unlinkBlocksByRow(row1, row2)

    def rowsOfType(self, stype: SequenceType):
        rows = np.where(self._sequenceBlocks.type == stype)[0]
        return rows, self._sequenceBlocks.index[rows], self._sequenceBlocks.name.iloc[rows]

    def linkGroup(self, groupid: int):
        return self._sequenceBlocks[self._sequenceBlocks.linkgroup == groupid]

    def roleNames(self):
        names = super().roleNames()
        names[SequenceBlockRoles.stype_name] = b"typename"
        return names

    @property
    def addedBlocksCount(self):
        return self._n_added


class SequenceBlockLinkProxyModel(QAbstractListModel):
    linkDescriptionChanged = Signal(str, arguments=["description"])

    def __init__(self, basemodel: SequenceBlockModel, blockid: int, parent=None):
        super().__init__(parent)
        self._basemodel = basemodel
        self._blockid = blockid
        self._row = self._basemodel.dataForId(self._blockid, SequenceBlockRoles.blockindex)

        self._basemodel.rowsInserted.connect(self._baseRowsInserted)
        self._basemodel.rowsRemoved.connect(self._baseRowsRemoved)
        self._basemodel.rowsMoved.connect(self._baseRowsMoved)
        self._basemodel.linkGroupChanged.connect(self._baseModelLinkGroupChanged)
        self._basemodel.modelReset.connect(self._baseModelReset)

    def _rowToBaseRow(self, row):
        return row + 1 if row >= self._row else row

    def _baseRowToRow(self, row):
        return row - 1 if row >= self._row else row

    def _baseRowsInserted(self, parent, first: int, last: int):
        if first < self._row:
            self._row += last - first + 1

    def _baseRowsRemoved(self, parent, first: int, last: int):
        if last < self._row:
            self._row -= last - first + 1
        elif first <= self._row <= last:
            self._row = -1

    def _baseRowsMoved(self, parent, sourceStart: int, sourceEnd: int, destinationParent, destinationRow: int):
        if self._row < 0:
            return
        if sourceEnd < self._row and destinationRow > self._row:
            self._row -= sourceEnd - sourceStart + 1
        elif sourceStart > self._row and destinationRow < self._row:
            self._row += sourceEnd - sourceStart + 1
        elif sourceStart <= self._row <= sourceEnd:
            self._row += destinationRow - sourceStart

    def rowCount(self, parent=QModelIndex()):
        return self._basemodel.rowCount() - 1

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole):
        if not index.isValid() or self._row < 0:
            return None
        row = self._rowToBaseRow(index.row())
        if role == Qt.DisplayRole:
            return self._basemodel.dataForRow(row, role)
        elif role == Qt.CheckStateRole:
            othergroup = self._basemodel.dataForRow(row, SequenceBlockRoles.linkgroup)
            if othergroup >= 0:
                selfgroup = self._basemodel.dataForRow(self._row, SequenceBlockRoles.linkgroup)
                if selfgroup == othergroup:
                    return Qt.Checked
            return Qt.Unchecked

        return None

    def setData(self, index: QModelIndex, value, role: int = Qt.EditRole):
        if not index.isValid() or self._row < 0:
            return False
        if role == Qt.CheckStateRole:
            row = self._rowToBaseRow(index.row())
            if value:
                self._basemodel.linkBlocksByRow(self._row, row)
            else:
                self._basemodel.unlinkBlocksByRow(self._row, row)
            return True
        else:
            return False

    def flags(self, index: QModelIndex):
        return Qt.ItemIsUserCheckable | Qt.ItemIsEnabled

    @Slot()
    def _baseModelLinkGroupChanged(self, rows: list[int]):
        if self._row in rows:
            for row in rows:
                row = self._baseRowToRow(row)
                self.dataChanged.emit(self.index(row, 0), self.index(row, 0), [Qt.CheckStateRole])
            self.linkDescriptionChanged.emit(self.linkDescription)

    @Slot()
    def _baseModelReset(self):
        self._row = -1

    @property
    def linkDescription(self):
        if self._row < 0:
            return ""
        linkgroup = self._basemodel.dataForRow(self._row, SequenceBlockRoles.linkgroup)
        if linkgroup < 0:
            return ""
        else:
            linkgroup = self._basemodel.linkGroup(linkgroup)
            return ", ".join(linkgroup.name[linkgroup.index != self._blockid])


class ReadProps(QWidget):
    _qmlfile = files(__package__).joinpath("SequenceView.qml")

    validChanged = Signal(bool, arguments=["valid"])
    edited = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._ui = Ui_ReadProps()
        self._ui.setupUi(self)

        self._reset()

        self._sequenceModel = SequenceBlockModel(self)

        menu = QMenu(self)
        for stype in SequenceType:
            action = QAction(stype.display_name, menu)
            action.setData(stype.value)
            menu.addAction(action)
        self._ui.addSequenceBtn.setMenu(menu)
        menu.triggered.connect(self._addSequenceClicked)

        self._ui.diagramQuick.rootContext().setContextProperty("myModel", self._sequenceModel)
        self._ui.diagramQuick.setSource(str(self._qmlfile))
        self._listView = self._ui.diagramQuick.rootObject().findChild(QObject, "listView")
        self._listView.currentIndexChanged.connect(self._selectedSequenceChanged)
        self._listView.deleteClicked.connect(self._deleteClicked)
        self._sequenceModel.rowsRemoved.connect(self._sequenceRemoved)
        self._sequenceModel.modelReset.connect(self._reset)

    def _addSequence(self, stype, name, linkgroup: int = -1, state: dict[str, Any] = None):
        blockid = self._sequenceModel.addBlock(name, stype, linkgroup)
        linkmodel = SequenceBlockLinkProxyModel(self._sequenceModel, blockid)
        props = sequencePropertiesFactory(stype, linkmodel, self)
        linkmodel.setParent(props)
        props.setName(name)
        if state is not None:
            props.from_dict(state)
        props.nameEdited.connect(self._sequenceNameEdited)
        props.validChanged.connect(self._validChanged)
        props.edited.connect(self.edited)
        self._ui.sequencePropertiesStck.addWidget(props)
        if props.valid:
            self._nvalid += 1

        return props

    @Slot()
    def _addSequenceClicked(self, action):
        name = f"sequence {self._sequenceModel.addedBlocksCount}"
        self._addSequence(SequenceType(action.data()), name)
        self.edited.emit()
        self._validAction()

    @Slot()
    def _sequenceNameEdited(self, text):
        self._sequenceModel.setBlockNameByRow(self._ui.sequencePropertiesStck.currentIndex(), text)

    @Slot()
    def _selectedSequenceChanged(self):
        self._ui.sequencePropertiesStck.setCurrentIndex(self._listView.property("currentIndex"))

    @Slot()
    def _deleteClicked(self, index):
        sname = self._sequenceModel.dataForRow(index)
        stype = self._sequenceModel.dataForRow(index, role=SequenceBlockRoles.stype_display_name)

        if (
            QMessageBox.question(
                self,
                "Delete sequence?",
                f'Do you really want to delete the {stype} "{sname}"?',
                QMessageBox.Yes | QMessageBox.No,
            )
            == QMessageBox.Yes
        ):
            self._sequenceModel.removeRow(index)

    @Slot()
    def _sequenceRemoved(self, parent, first, last):
        for i in range(first, last + 1):
            widget = self._ui.sequencePropertiesStck.widget(i)
            if widget.valid:
                self._nvalid -= 1
            self._ui.sequencePropertiesStck.removeWidget(widget)
            widget.deleteLater()
            self._validAction()

    @Slot()
    def _reset(self):
        while (widget := self._ui.sequencePropertiesStck.currentWidget()) is not None:
            self._ui.sequencePropertiesStck.removeWidget(widget)
            widget.deleteLater()
        self._nvalid = 0
        self._valid = False

    @Slot()
    def _validChanged(self, valid):
        if valid:
            self._nvalid += 1
        else:
            self._nvalid -= 1
        self._validAction()

    @property
    def valid(self):
        return self._valid

    def _validAction(self):
        valid = self._sequenceModel.rowCount() > 0 and self._nvalid == self._sequenceModel.rowCount()
        if valid != self._valid:
            self._valid = valid
            self.validChanged.emit(valid)

    def to_dict(self):
        blocks = []
        linkgroups = {}
        for i in range(self._sequenceModel.rowCount()):
            ret = self._ui.sequencePropertiesStck.widget(i).to_dict()
            ret["blocktype"] = self._sequenceModel.dataForRow(i, SequenceBlockRoles.stype_name)
            ret["blockname"] = self._sequenceModel.dataForRow(i, Qt.DisplayRole)
            linkgrp = self._sequenceModel.dataForRow(i, SequenceBlockRoles.linkgroup)
            if linkgrp >= 0:
                if linkgrp in linkgroups:
                    linkgroup = linkgroups[linkgrp]
                else:
                    linkgroup = len(linkgroups)
                    linkgroups[linkgrp] = linkgroup
                ret["linkgroup"] = linkgroup
            blocks.append(ret)
        return {"keep_nonbarcode": self._ui.keepNonBrcdChk.isChecked(), "blocks": blocks}

    def from_dict(self, d: dict[str, Any]):
        self._sequenceModel.clear()
        self._ui.keepNonBrcdChk.setChecked(d.get("keep_nonbarcode", False))
        if (blocks := d.get("blocks")) is not None:
            for block in blocks:
                self._addSequence(
                    SequenceType[block["blocktype"]], block["blockname"], block.get("linkgroup", -1), block
                )
        self._validAction()
