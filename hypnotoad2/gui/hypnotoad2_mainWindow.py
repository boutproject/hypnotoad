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
        self.action_New = QAction(Hypnotoad2)
        self.action_New.setObjectName(u"action_New")
        icon = QIcon()
        iconThemeName = u"document-new"
        if QIcon.hasThemeIcon(iconThemeName):
            icon = QIcon.fromTheme(iconThemeName)
        else:
            icon.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_New.setIcon(icon)
        self.action_Open = QAction(Hypnotoad2)
        self.action_Open.setObjectName(u"action_Open")
        icon1 = QIcon()
        iconThemeName = u"document-open"
        if QIcon.hasThemeIcon(iconThemeName):
            icon1 = QIcon.fromTheme(iconThemeName)
        else:
            icon1.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Open.setIcon(icon1)
        self.action_Save = QAction(Hypnotoad2)
        self.action_Save.setObjectName(u"action_Save")
        icon2 = QIcon()
        iconThemeName = u"document-save"
        if QIcon.hasThemeIcon(iconThemeName):
            icon2 = QIcon.fromTheme(iconThemeName)
        else:
            icon2.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Save.setIcon(icon2)
        self.action_Save_as = QAction(Hypnotoad2)
        self.action_Save_as.setObjectName(u"action_Save_as")
        icon3 = QIcon()
        iconThemeName = u"document-save-as"
        if QIcon.hasThemeIcon(iconThemeName):
            icon3 = QIcon.fromTheme(iconThemeName)
        else:
            icon3.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Save_as.setIcon(icon3)
        self.action_Quit = QAction(Hypnotoad2)
        self.action_Quit.setObjectName(u"action_Quit")
        icon4 = QIcon()
        iconThemeName = u"application-exit"
        if QIcon.hasThemeIcon(iconThemeName):
            icon4 = QIcon.fromTheme(iconThemeName)
        else:
            icon4.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Quit.setIcon(icon4)
        self.action_About = QAction(Hypnotoad2)
        self.action_About.setObjectName(u"action_About")
        icon5 = QIcon()
        iconThemeName = u"help-about"
        if QIcon.hasThemeIcon(iconThemeName):
            icon5 = QIcon.fromTheme(iconThemeName)
        else:
            icon5.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_About.setIcon(icon5)
        self.action_Run = QAction(Hypnotoad2)
        self.action_Run.setObjectName(u"action_Run")
        icon6 = QIcon()
        iconThemeName = u"system-run"
        if QIcon.hasThemeIcon(iconThemeName):
            icon6 = QIcon.fromTheme(iconThemeName)
        else:
            icon6.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Run.setIcon(icon6)
        self.action_Write_grid = QAction(Hypnotoad2)
        self.action_Write_grid.setObjectName(u"action_Write_grid")
        icon7 = QIcon()
        iconThemeName = u"document-print"
        if QIcon.hasThemeIcon(iconThemeName):
            icon7 = QIcon.fromTheme(iconThemeName)
        else:
            icon7.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Write_grid.setIcon(icon7)
        self.action_Revert = QAction(Hypnotoad2)
        self.action_Revert.setObjectName(u"action_Revert")
        icon8 = QIcon()
        iconThemeName = u"document-revert"
        if QIcon.hasThemeIcon(iconThemeName):
            icon8 = QIcon.fromTheme(iconThemeName)
        else:
            icon8.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Revert.setIcon(icon8)
        self.centralwidget = QWidget(Hypnotoad2)
        self.centralwidget.setObjectName(u"centralwidget")
        self.horizontalLayout_2 = QHBoxLayout(self.centralwidget)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.search_bar = QLineEdit(self.centralwidget)
        self.search_bar.setObjectName(u"search_bar")

        self.verticalLayout_2.addWidget(self.search_bar)

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
        self.menu_File = QMenu(self.menubar)
        self.menu_File.setObjectName(u"menu_File")
        self.menu_Help = QMenu(self.menubar)
        self.menu_Help.setObjectName(u"menu_Help")
        self.menu_Mesh = QMenu(self.menubar)
        self.menu_Mesh.setObjectName(u"menu_Mesh")
        Hypnotoad2.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(Hypnotoad2)
        self.statusbar.setObjectName(u"statusbar")
        Hypnotoad2.setStatusBar(self.statusbar)
        self.toolBar = QToolBar(Hypnotoad2)
        self.toolBar.setObjectName(u"toolBar")
        Hypnotoad2.addToolBar(Qt.TopToolBarArea, self.toolBar)

        self.menubar.addAction(self.menu_File.menuAction())
        self.menubar.addAction(self.menu_Mesh.menuAction())
        self.menubar.addAction(self.menu_Help.menuAction())
        self.menu_File.addAction(self.action_New)
        self.menu_File.addAction(self.action_Open)
        self.menu_File.addAction(self.action_Save)
        self.menu_File.addAction(self.action_Save_as)
        self.menu_File.addAction(self.action_Revert)
        self.menu_File.addSeparator()
        self.menu_File.addAction(self.action_Quit)
        self.menu_Help.addAction(self.action_About)
        self.menu_Mesh.addAction(self.action_Run)
        self.menu_Mesh.addAction(self.action_Write_grid)
        self.toolBar.addAction(self.action_New)
        self.toolBar.addAction(self.action_Open)
        self.toolBar.addAction(self.action_Save)
        self.toolBar.addAction(self.action_Save_as)
        self.toolBar.addAction(self.action_Run)
        self.toolBar.addAction(self.action_Write_grid)

        self.retranslateUi(Hypnotoad2)

        QMetaObject.connectSlotsByName(Hypnotoad2)
    # setupUi

    def retranslateUi(self, Hypnotoad2):
        Hypnotoad2.setWindowTitle(QCoreApplication.translate("Hypnotoad2", u"MainWindow", None))
        self.action_New.setText(QCoreApplication.translate("Hypnotoad2", u"&New", None))
        self.action_Open.setText(QCoreApplication.translate("Hypnotoad2", u"&Open", None))
        self.action_Save.setText(QCoreApplication.translate("Hypnotoad2", u"&Save", None))
#if QT_CONFIG(tooltip)
        self.action_Save.setToolTip(QCoreApplication.translate("Hypnotoad2", u"Save", None))
#endif // QT_CONFIG(tooltip)
        self.action_Save_as.setText(QCoreApplication.translate("Hypnotoad2", u"Save as", None))
        self.action_Quit.setText(QCoreApplication.translate("Hypnotoad2", u"&Quit", None))
        self.action_About.setText(QCoreApplication.translate("Hypnotoad2", u"&About", None))
        self.action_Run.setText(QCoreApplication.translate("Hypnotoad2", u"&Run", None))
        self.action_Write_grid.setText(QCoreApplication.translate("Hypnotoad2", u"&Write grid", None))
        self.action_Revert.setText(QCoreApplication.translate("Hypnotoad2", u"&Revert", None))
        self.options_file_label.setText(QCoreApplication.translate("Hypnotoad2", u"Options file", None))
        self.options_file_browse_button.setText(QCoreApplication.translate("Hypnotoad2", u"Browse", None))
        self.geqdsk_file_label.setText(QCoreApplication.translate("Hypnotoad2", u"geqdsk file", None))
        self.geqdsk_file_browse_button.setText(QCoreApplication.translate("Hypnotoad2", u"Browse", None))
        self.run_button.setText(QCoreApplication.translate("Hypnotoad2", u"Run", None))
        self.write_grid_button.setText(QCoreApplication.translate("Hypnotoad2", u"Write Grid", None))
        self.menu_File.setTitle(QCoreApplication.translate("Hypnotoad2", u"&File", None))
        self.menu_Help.setTitle(QCoreApplication.translate("Hypnotoad2", u"&Help", None))
        self.menu_Mesh.setTitle(QCoreApplication.translate("Hypnotoad2", u"&Mesh", None))
        self.toolBar.setWindowTitle(QCoreApplication.translate("Hypnotoad2", u"toolBar", None))
    # retranslateUi

