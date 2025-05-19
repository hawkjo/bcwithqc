from typing import Any

from PySide6.QtCore import Signal, Slot
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import QSpinBox, QWidget

from .sequencelist import SequenceList
from .sequencetypes import SequenceType
from .ui.sequenceprops import Ui_SequenceProps


class SequenceProps(QWidget):
    nameEdited = Signal(str, arguments=["text"])
    linkedWithChanged = Signal(int, int, arguments=["blockid_self", "blockid_other"])
    validChanged = Signal(bool, arguments=["valid"])
    edited = Signal()

    def __init__(self, linkmodel, parent=None):
        super().__init__(parent)
        self._ui = Ui_SequenceProps()
        self._linkmodel = linkmodel
        self._ui.setupUi(self)
        self._ui.sequenceNameLineEdit.textEdited.connect(self.nameEdited)
        self._ui.sequenceNameLineEdit.textChanged.connect(self._updateValid)
        self._ui.sequenceNameLineEdit.textChanged.connect(self.edited)
        self._ui.linkCmbBx.setModel(self._linkmodel)
        self._ui.linkCmbBx.setPlaceholderText(linkmodel.linkDescription)
        self._ui.linkCmbBx.setCurrentIndex(-1)
        self._linkmodel.linkDescriptionChanged.connect(self._ui.linkCmbBx.setPlaceholderText)
        self._linkmodel.linkDescriptionChanged.connect(lambda x: self._ui.linkCmbBx.setCurrentIndex(-1))
        self._ui.allowedErrorsSpnBx.valueChanged.connect(self.edited)

        pal = self._ui.linkCmbBx.palette()
        pal.setBrush(QPalette.PlaceholderText, pal.text())
        self._ui.linkCmbBx.setPalette(pal)

        self._childWidget = None
        self._valid = False

    @Slot()
    def _updateValid(self):
        valid = len(self._ui.sequenceNameLineEdit.text()) > 0
        if self._childWidget is not None and hasattr(self._childWidget, "valid"):
            valid = valid and self._childWidget.valid
        if valid != self._valid:
            self._valid = valid
            self.validChanged.emit(valid)

    @property
    def valid(self):
        return self._valid

    @property
    def name(self):
        return self._ui.sequenceNameLineEdit.text()

    @Slot()
    def setName(self, text: str):
        self._ui.sequenceNameLineEdit.setText(text)

    @Slot()
    def setAllowedErrorsMax(self, length: int):
        self._ui.allowedErrorsSpnBx.setMaximum(length // 2)

    def setChildWidget(self, widget):
        self.layout().insertWidget(1, widget)
        widget.setParent(self)
        self._childWidget = widget
        if widget.metaObject().indexOfSignal("validChanged(bool)") != -1:
            widget.validChanged.connect(self._updateValid)
        if widget.metaObject().indexOfSignal("edited()") != -1:
            widget.edited.connect(self.edited)
        self._updateValid()

    def to_dict(self):
        state = {
            "maxerrors": self._ui.allowedErrorsSpnBx.value(),
            "blockfunction": self._ui.blockFunctionCmbBx.currentText(),
        }
        if self._childWidget is not None:
            state |= self._childWidget.to_dict()
        return state

    def from_dict(self, d: dict[str, Any]):
        self._ui.allowedErrorsSpnBx.setValue(d.get("maxerrors", 1))
        self._ui.blockFunctionCmbBx.setCurrentText(d.get("blockfunction"))
        if self._childWidget is not None:
            self._childWidget.from_dict(d)


class RandomBarcodeProps(SequenceProps):
    def __init__(self, linkmodel, parent=None):
        super().__init__(linkmodel, parent)
        self._lengthSpnBx = QSpinBox(self)
        self._lengthSpnBx.setMinimum(1)
        self._lengthSpnBx.setSuffix(" bp")
        self._ui.formLayout.addRow("Length:", self._lengthSpnBx)
        self._lengthSpnBx.valueChanged.connect(self.setAllowedErrorsMax)
        self.setAllowedErrorsMax(self._lengthSpnBx.value())

        nitems = self._ui.blockFunctionCmbBx.count()
        self._ui.blockFunctionCmbBx.addItem("UMI")
        self._ui.blockFunctionCmbBx.setCurrentIndex(nitems)

    def to_dict(self):
        d = super().to_dict()
        d["length"] = self._lengthSpnBx.value()
        return d

    def from_dict(self, d: dict[str, Any]):
        self._lengthSpnBx.setValue(d["length"])
        super().from_dict(d)


def sequencePropertiesFactory(stype, linkmodel, parent=None):
    if stype == SequenceType.randomBarcode:
        return RandomBarcodeProps(linkmodel, parent)
    else:
        allowVariableLength = stype == SequenceType.constantRegion
        props = SequenceProps(linkmodel, parent)
        slist = SequenceList(allowVariableLength=allowVariableLength, allowUnnamedSequences=True, parent=props)
        slist.lengthChanged.connect(props.setAllowedErrorsMax)
        props.setChildWidget(slist)
        return props
