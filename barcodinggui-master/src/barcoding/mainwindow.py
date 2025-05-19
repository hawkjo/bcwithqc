import gzip
import json
import sys
import traceback

from PySide6.QtCore import Slot
from PySide6.QtGui import QAction, QIcon, QKeySequence
from PySide6.QtWidgets import QApplication, QDialog, QFileDialog, QMessageBox, QStyle, QWidget

from .cursormanager import WaitCursor
from .readprops import ReadProps
from .sequencetypes import SequencingMode
from .ui.mainwindow import Ui_MainWindow


class MainWindow(QWidget):
    _windowTitle = "barcode configurator"
    _dirtyWindowTitle = _windowTitle + " *"

    def __init__(self, parent=None):
        super().__init__(parent)
        self._ui = Ui_MainWindow()
        self._ui.setupUi(self)
        self._isdirty = False

        self.setWindowTitle(self._windowTitle)

        self._saveAction = QAction(self)
        self._saveAction.setShortcut(QKeySequence.Save)
        self._openAction = QAction(self)
        self._openAction.setShortcut(QKeySequence.Open)
        self._saveAction.triggered.connect(self._save)
        self._openAction.triggered.connect(self._open)
        self.addAction(self._saveAction)
        self.addAction(self._openAction)
        self._ui.saveBtn.setIcon(
            QIcon.fromTheme("document-save", QApplication.style().standardIcon(QStyle.SP_DialogSaveButton))
        )
        self._ui.openBtn.setIcon(
            QIcon.fromTheme("document-open", QApplication.style().standardIcon(QStyle.SP_DialogOpenButton))
        )

        # https://stackoverflow.com/questions/78067789/no-signals-from-qml-if-qquickwidget-not-constructed-with-proper-parent
        self._read1Props = ReadProps(self)
        self._read2Props = ReadProps(self)
        self._read1Props.validChanged.connect(self._validChanged)
        self._read2Props.validChanged.connect(self._validChanged)
        self._read1Props.edited.connect(self._edited)
        self._read2Props.edited.connect(self._edited)
        self._ui.readsTabWdgt.addTab(self._read1Props, "Read 1")
        self._ui.saveBtn.clicked.connect(self._saveAction.trigger)
        self._ui.openBtn.clicked.connect(self._openAction.trigger)
        self._ui.sequencingModeBtnGrp.setId(self._ui.singleEndBtn, SequencingMode.singleEnd)
        self._ui.sequencingModeBtnGrp.setId(self._ui.pairedEndBtn, SequencingMode.pairedEnd)
        self._ui.sequencingModeBtnGrp.idToggled.connect(self._sequencingModeChanged)
        self._ui.singleEndBtn.setChecked(True)
        self._ui.sequencingModeBtnGrp.idToggled.connect(self._edited)
        self._ui.unknownReadOrientationChkBx.toggled.connect(self._edited)
        self._validChanged()

    @Slot()
    def _sequencingModeChanged(self, btnId: int, checked: bool):
        if checked:
            if btnId == SequencingMode.pairedEnd:
                self._ui.readsTabWdgt.addTab(self._read2Props, "Read 2")
                self._read2Props.setVisible(True)
            else:
                self._ui.readsTabWdgt.removeTab(1)
                self._read2Props.setVisible(False)
            self._validChanged()

    @Slot()
    def _validChanged(self, *args):
        enabled = self._read1Props.valid and (self._ui.singleEndBtn.isChecked() or self._read2Props.valid)
        self._ui.saveBtn.setEnabled(enabled)
        self._saveAction.setEnabled(enabled)

    @Slot()
    def _edited(self):
        if not self._isdirty:
            self.setWindowTitle(self._dirtyWindowTitle)
            self._isdirty = True

    def _clean(self):
        self._isdirty = False
        self.setWindowTitle(self._windowTitle)

    def closeEvent(self, e):
        if (
            self._isdirty
            and QMessageBox.warning(
                self,
                "Unsaved configuration",
                "This configuration has not been saved. Do you really want to exit the application?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            != QMessageBox.Yes
        ):
            e.ignore()
        else:
            super().closeEvent(e)

    @Slot()
    def _save(self):
        dlg = QFileDialog(self, caption="Save barcode configuration", filter="Compressed JSON (*.json.gz)")
        dlg.setDefaultSuffix("json.gz")
        dlg.setAcceptMode(QFileDialog.AcceptSave)
        if dlg.exec_() == QDialog.Accepted:
            fname = dlg.selectedFiles()[0]

            with WaitCursor():
                ret = {
                    "unknown_read_orientation": self._ui.unknownReadOrientationChkBx.isChecked(),
                    "barcode_struct_r1": self._read1Props.to_dict(),
                }
                if self._ui.pairedEndBtn.isChecked():
                    ret["barcode_struct_r2"] = self._read2Props.to_dict()

                with gzip.open(fname, "wt", compresslevel=9, encoding="utf8") as f:
                    json.dump(ret, f, indent=4)

                self._clean()

    @Slot()
    def _open(self):
        fname = QFileDialog.getOpenFileName(self, "Open barcode configuration", filter="JSON (*.json *.json.gz)")[0]
        if not len(fname):
            return
        if self._isdirty and (
            QMessageBox.warning(
                self,
                "Unsaved configuration",
                "This configuration has not been saved. Opening a saved configuration will overwrite the current configuration. Do you want to proceed?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            != QMessageBox.Yes
        ):
            return
        try:
            with WaitCursor():
                openf = open
                with open(fname, "rb") as f:
                    if f.read(2) == b"\x1f\x8b":
                        openf = gzip.open
                with openf(fname, "rt", encoding="utf8") as f:
                    d = json.load(f)
                self._ui.unknownReadOrientationChkBx.setChecked(d.get("unknown_read_orientation", False))
                self._read1Props.from_dict(d.get("barcode_struct_r1", {}))
                r2 = d.get("barcode_struct_r2")
                if r2 is not None:
                    self._read2Props.from_dict(r2)
                    self._ui.pairedEndBtn.setChecked(True)
                else:
                    self._read2Props.from_dict({})
                    self._ui.singleEndBtn.setChecked(True)
                self._clean()
                self._validChanged()
        except Exception as e:  # noqa BLE001
            if sys.version_info.major >= 3 and sys.version_info.minor >= 10:
                traceback.print_exception(e)
            QMessageBox.critical(self, "Could not open barcode configuration", str(e))
