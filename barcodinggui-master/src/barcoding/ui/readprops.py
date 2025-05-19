################################################################################
## Form generated from reading UI file 'readprops.ui'
##
## Created by: Qt User Interface Compiler version 6.8.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (
    QCoreApplication,
    QMetaObject,
)
from PySide6.QtQuickWidgets import QQuickWidget
from PySide6.QtWidgets import (
    QCheckBox,
    QGridLayout,
    QPushButton,
    QSizePolicy,
    QSpacerItem,
    QStackedWidget,
    QVBoxLayout,
)


class Ui_ReadProps:
    def setupUi(self, ReadProps):
        if not ReadProps.objectName():
            ReadProps.setObjectName("ReadProps")
        ReadProps.resize(634, 337)
        self.verticalLayout = QVBoxLayout(ReadProps)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.diagramQuick = QQuickWidget(ReadProps)
        self.diagramQuick.setObjectName("diagramQuick")
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.diagramQuick.sizePolicy().hasHeightForWidth())
        self.diagramQuick.setSizePolicy(sizePolicy)
        self.diagramQuick.setResizeMode(QQuickWidget.SizeRootObjectToView)

        self.gridLayout.addWidget(self.diagramQuick, 1, 0, 2, 1)

        self.verticalSpacer = QSpacerItem(20, 15, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)

        self.gridLayout.addItem(self.verticalSpacer, 2, 1, 1, 1)

        self.addSequenceBtn = QPushButton(ReadProps)
        self.addSequenceBtn.setObjectName("addSequenceBtn")
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Fixed)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.addSequenceBtn.sizePolicy().hasHeightForWidth())
        self.addSequenceBtn.setSizePolicy(sizePolicy1)
        self.addSequenceBtn.setStyleSheet("QPushButton::menu-indicator{image:none;width:0px;}")

        self.gridLayout.addWidget(self.addSequenceBtn, 1, 1, 1, 1)

        self.keepNonBrcdChk = QCheckBox(ReadProps)
        self.keepNonBrcdChk.setObjectName("keepNonBrcdChk")
        self.keepNonBrcdChk.setStyleSheet("QCheckBox::indicator {\n" "    subcontrol-position: left top;\n" "}\n" "")

        self.gridLayout.addWidget(self.keepNonBrcdChk, 0, 0, 1, 2)

        self.verticalLayout.addLayout(self.gridLayout)

        self.sequencePropertiesStck = QStackedWidget(ReadProps)
        self.sequencePropertiesStck.setObjectName("sequencePropertiesStck")
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.sequencePropertiesStck.sizePolicy().hasHeightForWidth())
        self.sequencePropertiesStck.setSizePolicy(sizePolicy2)

        self.verticalLayout.addWidget(self.sequencePropertiesStck)

        self.retranslateUi(ReadProps)

        QMetaObject.connectSlotsByName(ReadProps)

    # setupUi

    def retranslateUi(self, ReadProps):
        ReadProps.setWindowTitle(QCoreApplication.translate("ReadProps", "Form", None))
        self.addSequenceBtn.setText(QCoreApplication.translate("ReadProps", "\u2795 Add sequence", None))
        self.keepNonBrcdChk.setText(QCoreApplication.translate("ReadProps", "keep non-barcode sequence", None))

    # retranslateUi
