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
        self.scrollArea = QScrollArea(self.centralwidget)
        self.scrollArea.setObjectName(u"scrollArea")
        self.scrollArea.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName(u"scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 1185, 725))
        self.verticalLayout = QVBoxLayout(self.scrollAreaWidgetContents)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.options_form_layout = QFormLayout()
        self.options_form_layout.setObjectName(u"options_form_layout")

        self.verticalLayout.addLayout(self.options_form_layout)

        self.scrollArea.setWidget(self.scrollAreaWidgetContents)

        self.verticalLayout_2.addWidget(self.scrollArea)

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

        self.write_grid_button = QPushButton(self.centralwidget)
        self.write_grid_button.setObjectName(u"write_grid_button")

        self.horizontalLayout_3.addWidget(self.write_grid_button)


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
        self.options_file_label.setText(QCoreApplication.translate("Hypnotoad2", u"Options file", None))
        self.options_file_browse_button.setText(QCoreApplication.translate("Hypnotoad2", u"Browse", None))
        self.geqdsk_file_label.setText(QCoreApplication.translate("Hypnotoad2", u"geqdsk file", None))
        self.geqdsk_file_browse_button.setText(QCoreApplication.translate("Hypnotoad2", u"Browse", None))
        self.run_button.setText(QCoreApplication.translate("Hypnotoad2", u"Run", None))
        self.write_grid_button.setText(QCoreApplication.translate("Hypnotoad2", u"Write Grid", None))
        self.menuFile.setTitle(QCoreApplication.translate("Hypnotoad2", u"File", None))
        self.toolBar.setWindowTitle(QCoreApplication.translate("Hypnotoad2", u"toolBar", None))
    # retranslateUi

