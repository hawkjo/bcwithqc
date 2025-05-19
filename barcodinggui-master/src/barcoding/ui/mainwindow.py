################################################################################
## Form generated from reading UI file 'mainwindow.ui'
##
## Created by: Qt User Interface Compiler version 6.8.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (
    QCoreApplication,
    QMetaObject,
)
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QHBoxLayout,
    QPushButton,
    QRadioButton,
    QSizePolicy,
    QSpacerItem,
    QTabWidget,
    QVBoxLayout,
)


class Ui_MainWindow:
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName("MainWindow")
        MainWindow.resize(862, 587)
        self.verticalLayout = QVBoxLayout(MainWindow)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.singleEndBtn = QRadioButton(MainWindow)
        self.sequencingModeBtnGrp = QButtonGroup(MainWindow)
        self.sequencingModeBtnGrp.setObjectName("sequencingModeBtnGrp")
        self.sequencingModeBtnGrp.addButton(self.singleEndBtn)
        self.singleEndBtn.setObjectName("singleEndBtn")

        self.horizontalLayout.addWidget(self.singleEndBtn)

        self.pairedEndBtn = QRadioButton(MainWindow)
        self.sequencingModeBtnGrp.addButton(self.pairedEndBtn)
        self.pairedEndBtn.setObjectName("pairedEndBtn")

        self.horizontalLayout.addWidget(self.pairedEndBtn)

        self.horizontalSpacer_2 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout.addItem(self.horizontalSpacer_2)

        self.unknownReadOrientationChkBx = QCheckBox(MainWindow)
        self.unknownReadOrientationChkBx.setObjectName("unknownReadOrientationChkBx")

        self.horizontalLayout.addWidget(self.unknownReadOrientationChkBx)

        self.verticalLayout.addLayout(self.horizontalLayout)

        self.readsTabWdgt = QTabWidget(MainWindow)
        self.readsTabWdgt.setObjectName("readsTabWdgt")
        self.readsTabWdgt.setTabBarAutoHide(True)

        self.verticalLayout.addWidget(self.readsTabWdgt)

        self.horizontalLayout_2 = QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.horizontalSpacer = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_2.addItem(self.horizontalSpacer)

        self.openBtn = QPushButton(MainWindow)
        self.openBtn.setObjectName("openBtn")

        self.horizontalLayout_2.addWidget(self.openBtn)

        self.saveBtn = QPushButton(MainWindow)
        self.saveBtn.setObjectName("saveBtn")

        self.horizontalLayout_2.addWidget(self.saveBtn)

        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.retranslateUi(MainWindow)

        self.readsTabWdgt.setCurrentIndex(-1)

        QMetaObject.connectSlotsByName(MainWindow)

    # setupUi

    def retranslateUi(self, MainWindow):
        self.singleEndBtn.setText(QCoreApplication.translate("MainWindow", "single end", None))
        self.pairedEndBtn.setText(QCoreApplication.translate("MainWindow", "paired end", None))
        self.unknownReadOrientationChkBx.setText(
            QCoreApplication.translate("MainWindow", "unknown read orientation", None)
        )
        self.openBtn.setText(QCoreApplication.translate("MainWindow", "Open", None))
        self.saveBtn.setText(QCoreApplication.translate("MainWindow", "Save", None))
        pass

    # retranslateUi
