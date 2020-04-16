# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'hypnotoad2_mainWindow.ui'
##
## Created by: Qt User Interface Compiler version 5.14.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from Qt.QtCore import (QCoreApplication, QDate, QDateTime, QMetaObject,
    QObject, QPoint, QRect, QSize, QTime, QUrl, Qt)
from Qt.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont,
    QFontDatabase, QIcon, QKeySequence, QLinearGradient, QPalette, QPainter,
    QPixmap, QRadialGradient)
from Qt.QtWidgets import *


class Ui_Hypnotoad2(object):
    def setupUi(self, Hypnotoad2):
        if not Hypnotoad2.objectName():
            Hypnotoad2.setObjectName(u"Hypnotoad2")
        Hypnotoad2.resize(1215, 863)
        self.actionNew = QAction(Hypnotoad2)
        self.actionNew.setObjectName(u"actionNew")
        icon = QIcon()
        iconThemeName = u"document-new"
        if QIcon.hasThemeIcon(iconThemeName):
            icon = QIcon.fromTheme(iconThemeName)
        else:
            icon.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)
        
        self.actionNew.setIcon(icon)
        self.actionOpen = QAction(Hypnotoad2)
        self.actionOpen.setObjectName(u"actionOpen")
        icon1 = QIcon()
        iconThemeName = u"document-open"
        if QIcon.hasThemeIcon(iconThemeName):
            icon1 = QIcon.fromTheme(iconThemeName)
        else:
            icon1.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)
        
        self.actionOpen.setIcon(icon1)
        self.centralwidget = QWidget(Hypnotoad2)
        self.centralwidget.setObjectName(u"centralwidget")
        self.horizontalLayout_2 = QHBoxLayout(self.centralwidget)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.formLayout = QFormLayout()
        self.formLayout.setObjectName(u"formLayout")
        self.nx_coreLabel = QLabel(self.centralwidget)
        self.nx_coreLabel.setObjectName(u"nx_coreLabel")

        self.formLayout.setWidget(0, QFormLayout.LabelRole, self.nx_coreLabel)

        self.nx_coreSpinBox = QSpinBox(self.centralwidget)
        self.nx_coreSpinBox.setObjectName(u"nx_coreSpinBox")

        self.formLayout.setWidget(0, QFormLayout.FieldRole, self.nx_coreSpinBox)

        self.nx_pfLabel = QLabel(self.centralwidget)
        self.nx_pfLabel.setObjectName(u"nx_pfLabel")

        self.formLayout.setWidget(1, QFormLayout.LabelRole, self.nx_pfLabel)

        self.nx_pfSpinBox = QSpinBox(self.centralwidget)
        self.nx_pfSpinBox.setObjectName(u"nx_pfSpinBox")

        self.formLayout.setWidget(1, QFormLayout.FieldRole, self.nx_pfSpinBox)

        self.nx_pf_lowerLabel = QLabel(self.centralwidget)
        self.nx_pf_lowerLabel.setObjectName(u"nx_pf_lowerLabel")

        self.formLayout.setWidget(2, QFormLayout.LabelRole, self.nx_pf_lowerLabel)

        self.nx_pf_lowerSpinBox = QSpinBox(self.centralwidget)
        self.nx_pf_lowerSpinBox.setObjectName(u"nx_pf_lowerSpinBox")

        self.formLayout.setWidget(2, QFormLayout.FieldRole, self.nx_pf_lowerSpinBox)

        self.nx_pf_upperLabel = QLabel(self.centralwidget)
        self.nx_pf_upperLabel.setObjectName(u"nx_pf_upperLabel")

        self.formLayout.setWidget(3, QFormLayout.LabelRole, self.nx_pf_upperLabel)

        self.nx_pf_upperSpinBox = QSpinBox(self.centralwidget)
        self.nx_pf_upperSpinBox.setObjectName(u"nx_pf_upperSpinBox")

        self.formLayout.setWidget(3, QFormLayout.FieldRole, self.nx_pf_upperSpinBox)

        self.nx_solLabel = QLabel(self.centralwidget)
        self.nx_solLabel.setObjectName(u"nx_solLabel")

        self.formLayout.setWidget(4, QFormLayout.LabelRole, self.nx_solLabel)

        self.nx_solSpinBox = QSpinBox(self.centralwidget)
        self.nx_solSpinBox.setObjectName(u"nx_solSpinBox")

        self.formLayout.setWidget(4, QFormLayout.FieldRole, self.nx_solSpinBox)

        self.nx_sol_lowerLabel = QLabel(self.centralwidget)
        self.nx_sol_lowerLabel.setObjectName(u"nx_sol_lowerLabel")

        self.formLayout.setWidget(5, QFormLayout.LabelRole, self.nx_sol_lowerLabel)

        self.nx_sol_lowerSpinBox = QSpinBox(self.centralwidget)
        self.nx_sol_lowerSpinBox.setObjectName(u"nx_sol_lowerSpinBox")

        self.formLayout.setWidget(5, QFormLayout.FieldRole, self.nx_sol_lowerSpinBox)

        self.nx_sol_upperLabel = QLabel(self.centralwidget)
        self.nx_sol_upperLabel.setObjectName(u"nx_sol_upperLabel")

        self.formLayout.setWidget(6, QFormLayout.LabelRole, self.nx_sol_upperLabel)

        self.nx_sol_upperSpinBox = QSpinBox(self.centralwidget)
        self.nx_sol_upperSpinBox.setObjectName(u"nx_sol_upperSpinBox")

        self.formLayout.setWidget(6, QFormLayout.FieldRole, self.nx_sol_upperSpinBox)

        self.dCTLabel = QLabel(self.centralwidget)
        self.dCTLabel.setObjectName(u"dCTLabel")

        self.formLayout.setWidget(7, QFormLayout.LabelRole, self.dCTLabel)

        self.dCTCheckBox = QCheckBox(self.centralwidget)
        self.dCTCheckBox.setObjectName(u"dCTCheckBox")

        self.formLayout.setWidget(7, QFormLayout.FieldRole, self.dCTCheckBox)


        self.verticalLayout_2.addLayout(self.formLayout)

        self.horizontalLayout_3 = QHBoxLayout()
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.options_file_label = QLabel(self.centralwidget)
        self.options_file_label.setObjectName(u"options_file_label")

        self.horizontalLayout_3.addWidget(self.options_file_label)

        self.options_file_line_edit = QLineEdit(self.centralwidget)
        self.options_file_line_edit.setObjectName(u"options_file_line_edit")

        self.horizontalLayout_3.addWidget(self.options_file_line_edit)

        self.options_file_browse_button = QPushButton(self.centralwidget)
        self.options_file_browse_button.setObjectName(u"options_file_browse_button")

        self.horizontalLayout_3.addWidget(self.options_file_browse_button)

        self.geqdsk_file_label = QLabel(self.centralwidget)
        self.geqdsk_file_label.setObjectName(u"geqdsk_file_label")

        self.horizontalLayout_3.addWidget(self.geqdsk_file_label)

        self.geqdsk_file_line_edit = QLineEdit(self.centralwidget)
        self.geqdsk_file_line_edit.setObjectName(u"geqdsk_file_line_edit")

        self.horizontalLayout_3.addWidget(self.geqdsk_file_line_edit)

        self.geqdsk_file_browse_button = QPushButton(self.centralwidget)
        self.geqdsk_file_browse_button.setObjectName(u"geqdsk_file_browse_button")

        self.horizontalLayout_3.addWidget(self.geqdsk_file_browse_button)

        self.run_button = QPushButton(self.centralwidget)
        self.run_button.setObjectName(u"run_button")

        self.horizontalLayout_3.addWidget(self.run_button)


        self.verticalLayout_2.addLayout(self.horizontalLayout_3)


        self.horizontalLayout.addLayout(self.verticalLayout_2)

        self.plottingArea = QWidget(self.centralwidget)
        self.plottingArea.setObjectName(u"plottingArea")

        self.horizontalLayout.addWidget(self.plottingArea)


        self.horizontalLayout_2.addLayout(self.horizontalLayout)

        Hypnotoad2.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(Hypnotoad2)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1215, 22))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        Hypnotoad2.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(Hypnotoad2)
        self.statusbar.setObjectName(u"statusbar")
        Hypnotoad2.setStatusBar(self.statusbar)
        self.toolBar = QToolBar(Hypnotoad2)
        self.toolBar.setObjectName(u"toolBar")
        Hypnotoad2.addToolBar(Qt.TopToolBarArea, self.toolBar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menuFile.addAction(self.actionNew)
        self.menuFile.addAction(self.actionOpen)
        self.toolBar.addAction(self.actionNew)
        self.toolBar.addAction(self.actionOpen)

        self.retranslateUi(Hypnotoad2)

        QMetaObject.connectSlotsByName(Hypnotoad2)
    # setupUi

    def retranslateUi(self, Hypnotoad2):
        Hypnotoad2.setWindowTitle(QCoreApplication.translate("Hypnotoad2", u"MainWindow", None))
        self.actionNew.setText(QCoreApplication.translate("Hypnotoad2", u"New", None))
        self.actionOpen.setText(QCoreApplication.translate("Hypnotoad2", u"Open", None))
        self.nx_coreLabel.setText(QCoreApplication.translate("Hypnotoad2", u"nx_core", None))
        self.nx_pfLabel.setText(QCoreApplication.translate("Hypnotoad2", u"nx_pf", None))
        self.nx_pf_lowerLabel.setText(QCoreApplication.translate("Hypnotoad2", u"nx_pf_lower", None))
        self.nx_pf_upperLabel.setText(QCoreApplication.translate("Hypnotoad2", u"nx_pf_upper", None))
        self.nx_solLabel.setText(QCoreApplication.translate("Hypnotoad2", u"nx_sol", None))
        self.nx_sol_lowerLabel.setText(QCoreApplication.translate("Hypnotoad2", u"nx_sol_lower", None))
        self.nx_sol_upperLabel.setText(QCoreApplication.translate("Hypnotoad2", u"nx_sol_upper", None))
        self.dCTLabel.setText(QCoreApplication.translate("Hypnotoad2", u"DCT", None))
        self.options_file_label.setText(QCoreApplication.translate("Hypnotoad2", u"Options file", None))
        self.options_file_browse_button.setText(QCoreApplication.translate("Hypnotoad2", u"Browse", None))
        self.geqdsk_file_label.setText(QCoreApplication.translate("Hypnotoad2", u"geqdsk file", None))
        self.geqdsk_file_browse_button.setText(QCoreApplication.translate("Hypnotoad2", u"Browse", None))
        self.run_button.setText(QCoreApplication.translate("Hypnotoad2", u"Run", None))
        self.menuFile.setTitle(QCoreApplication.translate("Hypnotoad2", u"File", None))
        self.toolBar.setWindowTitle(QCoreApplication.translate("Hypnotoad2", u"toolBar", None))
    # retranslateUi

