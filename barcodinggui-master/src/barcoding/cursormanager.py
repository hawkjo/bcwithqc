from PySide6.QtCore import Qt
from PySide6.QtGui import QCursor
from PySide6.QtWidgets import QApplication


class CursorContext:
    def __init__(self, cursor):
        self._cursor = cursor

    def __enter__(self):
        QApplication.setOverrideCursor(QCursor(self._cursor))

    def __exit__(self, *exc):
        QApplication.restoreOverrideCursor()


class WaitCursor(CursorContext):
    def __init__(self):
        super().__init__(Qt.WaitCursor)
