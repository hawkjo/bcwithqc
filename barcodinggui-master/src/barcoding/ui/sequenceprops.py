################################################################################
## Form generated from reading UI file 'sequenceprops.ui'
##
## Created by: Qt User Interface Compiler version 6.8.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (
    QCoreApplication,
    QMetaObject,
    QSize,
    Qt,
)
from PySide6.QtWidgets import (
    QComboBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QSizePolicy,
    QSpacerItem,
    QSpinBox,
    QVBoxLayout,
)


class Ui_SequenceProps:
    def setupUi(self, SequenceProps):
        if not SequenceProps.objectName():
            SequenceProps.setObjectName("SequenceProps")
        SequenceProps.resize(400, 338)
        self.verticalLayout = QVBoxLayout(SequenceProps)
        self.verticalLayout.setObjectName("verticalLayout")
        self.verticalLayout.setContentsMargins(0, -1, 0, -1)
        self.formLayout = QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.formLayout.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignTrailing | Qt.AlignmentFlag.AlignVCenter
        )
        self.label = QLabel(SequenceProps)
        self.label.setObjectName("label")

        self.formLayout.setWidget(0, QFormLayout.LabelRole, self.label)

        self.sequenceNameLineEdit = QLineEdit(SequenceProps)
        self.sequenceNameLineEdit.setObjectName("sequenceNameLineEdit")

        self.formLayout.setWidget(0, QFormLayout.FieldRole, self.sequenceNameLineEdit)

        self.label_2 = QLabel(SequenceProps)
        self.label_2.setObjectName("label_2")

        self.formLayout.setWidget(1, QFormLayout.LabelRole, self.label_2)

        self.linkCmbBx = QComboBox(SequenceProps)
        self.linkCmbBx.addItem("")
        self.linkCmbBx.setObjectName("linkCmbBx")
        self.linkCmbBx.setMinimumSize(QSize(128, 0))
        self.linkCmbBx.setStyleSheet("")

        self.formLayout.setWidget(1, QFormLayout.FieldRole, self.linkCmbBx)

        self.label_3 = QLabel(SequenceProps)
        self.label_3.setObjectName("label_3")

        self.formLayout.setWidget(2, QFormLayout.LabelRole, self.label_3)

        self.label_4 = QLabel(SequenceProps)
        self.label_4.setObjectName("label_4")

        self.formLayout.setWidget(3, QFormLayout.LabelRole, self.label_4)

        self.blockFunctionCmbBx = QComboBox(SequenceProps)
        self.blockFunctionCmbBx.addItem("")
        self.blockFunctionCmbBx.addItem("")
        self.blockFunctionCmbBx.setObjectName("blockFunctionCmbBx")
        self.blockFunctionCmbBx.setMinimumSize(QSize(128, 0))

        self.formLayout.setWidget(3, QFormLayout.FieldRole, self.blockFunctionCmbBx)

        self.allowedErrorsSpnBx = QSpinBox(SequenceProps)
        self.allowedErrorsSpnBx.setObjectName("allowedErrorsSpnBx")
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.allowedErrorsSpnBx.sizePolicy().hasHeightForWidth())
        self.allowedErrorsSpnBx.setSizePolicy(sizePolicy)
        self.allowedErrorsSpnBx.setMinimumSize(QSize(128, 0))
        self.allowedErrorsSpnBx.setMaximum(999999999)
        self.allowedErrorsSpnBx.setValue(1)

        self.formLayout.setWidget(2, QFormLayout.FieldRole, self.allowedErrorsSpnBx)

        self.verticalLayout.addLayout(self.formLayout)

        self.verticalSpacer = QSpacerItem(0, 0, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Preferred)

        self.verticalLayout.addItem(self.verticalSpacer)

        # if QT_CONFIG(shortcut)
        self.label.setBuddy(self.sequenceNameLineEdit)
        self.label_2.setBuddy(self.linkCmbBx)
        self.label_3.setBuddy(self.allowedErrorsSpnBx)
        self.label_4.setBuddy(self.blockFunctionCmbBx)
        # endif // QT_CONFIG(shortcut)

        self.retranslateUi(SequenceProps)

        self.linkCmbBx.setCurrentIndex(-1)
        self.blockFunctionCmbBx.setCurrentIndex(1)

        QMetaObject.connectSlotsByName(SequenceProps)

    # setupUi

    def retranslateUi(self, SequenceProps):
        SequenceProps.setWindowTitle(QCoreApplication.translate("SequenceProps", "Form", None))
        self.label.setText(QCoreApplication.translate("SequenceProps", "Name:", None))
        self.label_2.setText(QCoreApplication.translate("SequenceProps", "Linked with:", None))
        self.linkCmbBx.setItemText(0, "")

        self.label_3.setText(QCoreApplication.translate("SequenceProps", "Allowed errors:", None))
        self.label_4.setText(QCoreApplication.translate("SequenceProps", "BAM output:", None))
        self.blockFunctionCmbBx.setItemText(0, QCoreApplication.translate("SequenceProps", "discard", None))
        self.blockFunctionCmbBx.setItemText(1, QCoreApplication.translate("SequenceProps", "barcode", None))

    # retranslateUi
