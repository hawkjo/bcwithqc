################################################################################
## Form generated from reading UI file 'sequencelist.ui'
##
## Created by: Qt User Interface Compiler version 6.8.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (
    QCoreApplication,
    QMetaObject,
    QSize,
)
from PySide6.QtWidgets import (
    QAbstractItemView,
    QAbstractScrollArea,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSpacerItem,
    QTableView,
    QVBoxLayout,
)


class Ui_SequenceList:
    def setupUi(self, SequenceList):
        if not SequenceList.objectName():
            SequenceList.setObjectName("SequenceList")
        SequenceList.resize(615, 308)
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(SequenceList.sizePolicy().hasHeightForWidth())
        SequenceList.setSizePolicy(sizePolicy)
        self.horizontalLayout = QHBoxLayout(SequenceList)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.sequenceTbl = QTableView(SequenceList)
        self.sequenceTbl.setObjectName("sequenceTbl")
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sizePolicy1.setHorizontalStretch(255)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.sequenceTbl.sizePolicy().hasHeightForWidth())
        self.sequenceTbl.setSizePolicy(sizePolicy1)
        self.sequenceTbl.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContentsOnFirstShow)
        self.sequenceTbl.setEditTriggers(QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed)
        self.sequenceTbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.sequenceTbl.horizontalHeader().setVisible(True)
        self.sequenceTbl.horizontalHeader().setStretchLastSection(True)

        self.horizontalLayout.addWidget(self.sequenceTbl)

        self.verticalLayout = QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.seqlenLbl = QLabel(SequenceList)
        self.seqlenLbl.setObjectName("seqlenLbl")
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.seqlenLbl.sizePolicy().hasHeightForWidth())
        self.seqlenLbl.setSizePolicy(sizePolicy2)

        self.verticalLayout.addWidget(self.seqlenLbl)

        self.duplicateSeqLbl = QLabel(SequenceList)
        self.duplicateSeqLbl.setObjectName("duplicateSeqLbl")
        sizePolicy2.setHeightForWidth(self.duplicateSeqLbl.sizePolicy().hasHeightForWidth())
        self.duplicateSeqLbl.setSizePolicy(sizePolicy2)

        self.verticalLayout.addWidget(self.duplicateSeqLbl)

        self.duplicateNameLbl = QLabel(SequenceList)
        self.duplicateNameLbl.setObjectName("duplicateNameLbl")

        self.verticalLayout.addWidget(self.duplicateNameLbl)

        self.unnamedSeqLbl = QLabel(SequenceList)
        self.unnamedSeqLbl.setObjectName("unnamedSeqLbl")

        self.verticalLayout.addWidget(self.unnamedSeqLbl)

        self.emptySeqLbl = QLabel(SequenceList)
        self.emptySeqLbl.setObjectName("emptySeqLbl")

        self.verticalLayout.addWidget(self.emptySeqLbl)

        self.verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout.addItem(self.verticalSpacer)

        self.horizontalLayout_2 = QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.deleteBtn = QPushButton(SequenceList)
        self.deleteBtn.setObjectName("deleteBtn")
        self.deleteBtn.setEnabled(False)
        sizePolicy2.setHeightForWidth(self.deleteBtn.sizePolicy().hasHeightForWidth())
        self.deleteBtn.setSizePolicy(sizePolicy2)
        self.deleteBtn.setMinimumSize(QSize(1, 0))
        self.deleteBtn.setMaximumSize(QSize(16777215, 16777215))
        self.deleteBtn.setStyleSheet("")
        self.deleteBtn.setFlat(True)

        self.verticalLayout_2.addWidget(self.deleteBtn)

        self.line = QFrame(SequenceList)
        self.line.setObjectName("line")
        self.line.setFrameShape(QFrame.Shape.HLine)
        self.line.setFrameShadow(QFrame.Shadow.Sunken)

        self.verticalLayout_2.addWidget(self.line)

        self.openBtn = QPushButton(SequenceList)
        self.openBtn.setObjectName("openBtn")
        sizePolicy2.setHeightForWidth(self.openBtn.sizePolicy().hasHeightForWidth())
        self.openBtn.setSizePolicy(sizePolicy2)
        self.openBtn.setMinimumSize(QSize(1, 0))
        self.openBtn.setStyleSheet("")
        self.openBtn.setFlat(True)

        self.verticalLayout_2.addWidget(self.openBtn)

        self.pasteBtn = QPushButton(SequenceList)
        self.pasteBtn.setObjectName("pasteBtn")
        sizePolicy2.setHeightForWidth(self.pasteBtn.sizePolicy().hasHeightForWidth())
        self.pasteBtn.setSizePolicy(sizePolicy2)
        self.pasteBtn.setMinimumSize(QSize(1, 0))
        self.pasteBtn.setStyleSheet("")
        self.pasteBtn.setFlat(True)

        self.verticalLayout_2.addWidget(self.pasteBtn)

        self.horizontalLayout_2.addLayout(self.verticalLayout_2)

        self.horizontalSpacer = QSpacerItem(31, 0, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_2.addItem(self.horizontalSpacer)

        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.horizontalLayout.addLayout(self.verticalLayout)

        self.retranslateUi(SequenceList)

        QMetaObject.connectSlotsByName(SequenceList)

    # setupUi

    def retranslateUi(self, SequenceList):
        self.seqlenLbl.setText(QCoreApplication.translate("SequenceList", "varying sequence lengths", None))
        self.duplicateSeqLbl.setText(QCoreApplication.translate("SequenceList", "duplicate sequences", None))
        self.duplicateNameLbl.setText(QCoreApplication.translate("SequenceList", "duplicate sequence names", None))
        self.unnamedSeqLbl.setText(QCoreApplication.translate("SequenceList", "unnamed sequences", None))
        self.emptySeqLbl.setText(QCoreApplication.translate("SequenceList", "empty sequences", None))
        # if QT_CONFIG(tooltip)
        self.deleteBtn.setToolTip(QCoreApplication.translate("SequenceList", "Delete selected barcodes", None))
        # endif // QT_CONFIG(tooltip)
        self.deleteBtn.setText(QCoreApplication.translate("SequenceList", "\u274c", None))
        # if QT_CONFIG(tooltip)
        self.openBtn.setToolTip(QCoreApplication.translate("SequenceList", "Import barcodes from file", None))
        # endif // QT_CONFIG(tooltip)
        # if QT_CONFIG(tooltip)
        self.pasteBtn.setToolTip(QCoreApplication.translate("SequenceList", "Import barcodes from clipboard", None))
        # endif // QT_CONFIG(tooltip)
        pass

    # retranslateUi
