# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'hypnotoad_preferences.ui'
##
## Created by: Qt User Interface Compiler version 5.14.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from Qt.QtCore import (
    QCoreApplication,
    QDate,
    QDateTime,
    QMetaObject,
    QObject,
    QPoint,
    QRect,
    QSize,
    QTime,
    QUrl,
    Qt,
)
from Qt.QtGui import (
    QBrush,
    QColor,
    QConicalGradient,
    QCursor,
    QFont,
    QFontDatabase,
    QIcon,
    QKeySequence,
    QLinearGradient,
    QPalette,
    QPainter,
    QPixmap,
    QRadialGradient,
)
from Qt.QtWidgets import *


class Ui_Preferences(object):
    def setupUi(self, Preferences):
        if not Preferences.objectName():
            Preferences.setObjectName("Preferences")
        Preferences.resize(400, 300)
        self.verticalLayout = QVBoxLayout(Preferences)
        self.verticalLayout.setObjectName("verticalLayout")
        self.frame = QFrame(Preferences)
        self.frame.setObjectName("frame")
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setFrameShadow(QFrame.Raised)
        self.formLayoutWidget = QWidget(self.frame)
        self.formLayoutWidget.setObjectName("formLayoutWidget")
        self.formLayoutWidget.setGeometry(QRect(10, 10, 361, 211))
        self.formLayout = QFormLayout(self.formLayoutWidget)
        self.formLayout.setObjectName("formLayout")
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.defaultGridFileNameLabel = QLabel(self.formLayoutWidget)
        self.defaultGridFileNameLabel.setObjectName("defaultGridFileNameLabel")

        self.formLayout.setWidget(
            0, QFormLayout.LabelRole, self.defaultGridFileNameLabel
        )

        self.defaultGridFileNameLineEdit = QLineEdit(self.formLayoutWidget)
        self.defaultGridFileNameLineEdit.setObjectName("defaultGridFileNameLineEdit")

        self.formLayout.setWidget(
            0, QFormLayout.FieldRole, self.defaultGridFileNameLineEdit
        )

        self.plotXlowLabel = QLabel(self.formLayoutWidget)
        self.plotXlowLabel.setObjectName("plotXlowLabel")

        self.formLayout.setWidget(1, QFormLayout.LabelRole, self.plotXlowLabel)

        self.plotXlowCheckBox = QCheckBox(self.formLayoutWidget)
        self.plotXlowCheckBox.setObjectName("plotXlowCheckBox")
        self.plotXlowCheckBox.setChecked(True)

        self.formLayout.setWidget(1, QFormLayout.FieldRole, self.plotXlowCheckBox)

        self.plotYlowLabel = QLabel(self.formLayoutWidget)
        self.plotYlowLabel.setObjectName("plotYlowLabel")

        self.formLayout.setWidget(2, QFormLayout.LabelRole, self.plotYlowLabel)

        self.plotYlowCheckBox = QCheckBox(self.formLayoutWidget)
        self.plotYlowCheckBox.setObjectName("plotYlowCheckBox")
        self.plotYlowCheckBox.setChecked(True)

        self.formLayout.setWidget(2, QFormLayout.FieldRole, self.plotYlowCheckBox)

        self.plotCornersLabel = QLabel(self.formLayoutWidget)
        self.plotCornersLabel.setObjectName("plotCornersLabel")

        self.formLayout.setWidget(3, QFormLayout.LabelRole, self.plotCornersLabel)

        self.plotCornersCheckBox = QCheckBox(self.formLayoutWidget)
        self.plotCornersCheckBox.setObjectName("plotCornersCheckBox")
        self.plotCornersCheckBox.setChecked(True)

        self.formLayout.setWidget(3, QFormLayout.FieldRole, self.plotCornersCheckBox)

        self.saveFullYamlLabel = QLabel(self.formLayoutWidget)
        self.saveFullYamlLabel.setObjectName("saveFullYamlLabel")

        self.formLayout.setWidget(4, QFormLayout.LabelRole, self.saveFullYamlLabel)

        self.saveFullYamlCheckBox = QCheckBox(self.formLayoutWidget)
        self.saveFullYamlCheckBox.setObjectName("saveFullYamlCheckBox")
        self.saveFullYamlCheckBox.setChecked(False)

        self.formLayout.setWidget(4, QFormLayout.FieldRole, self.saveFullYamlCheckBox)

        self.verticalLayout.addWidget(self.frame)

        self.buttonBox = QDialogButtonBox(Preferences)
        self.buttonBox.setObjectName("buttonBox")
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel | QDialogButtonBox.Ok)

        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Preferences)
        self.buttonBox.accepted.connect(Preferences.accept)
        self.buttonBox.rejected.connect(Preferences.reject)

        QMetaObject.connectSlotsByName(Preferences)

    # setupUi

    def retranslateUi(self, Preferences):
        Preferences.setWindowTitle(
            QCoreApplication.translate("Preferences", "Hypnotoad Preferences", None)
        )
        self.defaultGridFileNameLabel.setText(
            QCoreApplication.translate("Preferences", "Default grid file name", None)
        )
        self.defaultGridFileNameLineEdit.setPlaceholderText(
            QCoreApplication.translate("Preferences", "bout.grid.nc", None)
        )
        self.plotXlowLabel.setText(
            QCoreApplication.translate("Preferences", "Plot xlow", None)
        )
        self.plotYlowLabel.setText(
            QCoreApplication.translate("Preferences", "Plot ylow", None)
        )
        self.plotCornersLabel.setText(
            QCoreApplication.translate("Preferences", "Plot corners", None)
        )
        self.saveFullYamlLabel.setText(
            QCoreApplication.translate(
                "Preferences",
                "Include options set by\ndefault when saving to YAML",
                None,
            )
        )

    # retranslateUi
