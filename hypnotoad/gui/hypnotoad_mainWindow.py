# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'hypnotoad_mainWindow.ui'
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


class Ui_Hypnotoad(object):
    def setupUi(self, Hypnotoad):
        if not Hypnotoad.objectName():
            Hypnotoad.setObjectName(u"Hypnotoad")
        Hypnotoad.resize(1215, 863)
        self.action_New = QAction(Hypnotoad)
        self.action_New.setObjectName(u"action_New")
        icon = QIcon()
        iconThemeName = u"document-new"
        if QIcon.hasThemeIcon(iconThemeName):
            icon = QIcon.fromTheme(iconThemeName)
        else:
            icon.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_New.setIcon(icon)
        self.action_Open = QAction(Hypnotoad)
        self.action_Open.setObjectName(u"action_Open")
        icon1 = QIcon()
        iconThemeName = u"document-open"
        if QIcon.hasThemeIcon(iconThemeName):
            icon1 = QIcon.fromTheme(iconThemeName)
        else:
            icon1.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Open.setIcon(icon1)
        self.action_Save = QAction(Hypnotoad)
        self.action_Save.setObjectName(u"action_Save")
        icon2 = QIcon()
        iconThemeName = u"document-save"
        if QIcon.hasThemeIcon(iconThemeName):
            icon2 = QIcon.fromTheme(iconThemeName)
        else:
            icon2.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Save.setIcon(icon2)
        self.action_Save_as = QAction(Hypnotoad)
        self.action_Save_as.setObjectName(u"action_Save_as")
        icon3 = QIcon()
        iconThemeName = u"document-save-as"
        if QIcon.hasThemeIcon(iconThemeName):
            icon3 = QIcon.fromTheme(iconThemeName)
        else:
            icon3.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Save_as.setIcon(icon3)
        self.action_Quit = QAction(Hypnotoad)
        self.action_Quit.setObjectName(u"action_Quit")
        icon4 = QIcon()
        iconThemeName = u"application-exit"
        if QIcon.hasThemeIcon(iconThemeName):
            icon4 = QIcon.fromTheme(iconThemeName)
        else:
            icon4.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Quit.setIcon(icon4)
        self.action_About = QAction(Hypnotoad)
        self.action_About.setObjectName(u"action_About")
        icon5 = QIcon()
        iconThemeName = u"help-about"
        if QIcon.hasThemeIcon(iconThemeName):
            icon5 = QIcon.fromTheme(iconThemeName)
        else:
            icon5.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_About.setIcon(icon5)
        self.action_Run = QAction(Hypnotoad)
        self.action_Run.setObjectName(u"action_Run")
        icon6 = QIcon()
        iconThemeName = u"system-run"
        if QIcon.hasThemeIcon(iconThemeName):
            icon6 = QIcon.fromTheme(iconThemeName)
        else:
            icon6.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Run.setIcon(icon6)
        self.action_Write_grid = QAction(Hypnotoad)
        self.action_Write_grid.setObjectName(u"action_Write_grid")
        icon7 = QIcon()
        iconThemeName = u"document-print"
        if QIcon.hasThemeIcon(iconThemeName):
            icon7 = QIcon.fromTheme(iconThemeName)
        else:
            icon7.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Write_grid.setIcon(icon7)
        self.action_Revert = QAction(Hypnotoad)
        self.action_Revert.setObjectName(u"action_Revert")
        icon8 = QIcon()
        iconThemeName = u"document-revert"
        if QIcon.hasThemeIcon(iconThemeName):
            icon8 = QIcon.fromTheme(iconThemeName)
        else:
            icon8.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Revert.setIcon(icon8)
        self.action_Preferences = QAction(Hypnotoad)
        self.action_Preferences.setObjectName(u"action_Preferences")
        icon9 = QIcon()
        iconThemeName = u"document-properties"
        if QIcon.hasThemeIcon(iconThemeName):
            icon9 = QIcon.fromTheme(iconThemeName)
        else:
            icon9.addFile(u".", QSize(), QIcon.Normal, QIcon.Off)

        self.action_Preferences.setIcon(icon9)
        self.centralwidget = QWidget(Hypnotoad)
        self.centralwidget.setObjectName(u"centralwidget")
        self.verticalLayout = QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.tabWidget = QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName(u"tabWidget")
        self.equilibrium_tab = QWidget()
        self.equilibrium_tab.setObjectName(u"equilibrium_tab")
        self.horizontalLayout_2 = QHBoxLayout(self.equilibrium_tab)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.verticalLayout_5 = QVBoxLayout()
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.horizontalLayout_4 = QHBoxLayout()
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.eq_geqdsk_label = QLabel(self.equilibrium_tab)
        self.eq_geqdsk_label.setObjectName(u"eq_geqdsk_label")

        self.horizontalLayout_4.addWidget(self.eq_geqdsk_label)

        self.eq_geqdsk_lineedit = QLineEdit(self.equilibrium_tab)
        self.eq_geqdsk_lineedit.setObjectName(u"eq_geqdsk_lineedit")

        self.horizontalLayout_4.addWidget(self.eq_geqdsk_lineedit)

        self.eq_geqdsk_browse = QPushButton(self.equilibrium_tab)
        self.eq_geqdsk_browse.setObjectName(u"eq_geqdsk_browse")

        self.horizontalLayout_4.addWidget(self.eq_geqdsk_browse)


        self.verticalLayout_5.addLayout(self.horizontalLayout_4)

        self.normalised_psi_box = QGroupBox(self.equilibrium_tab)
        self.normalised_psi_box.setObjectName(u"normalised_psi_box")
        self.formLayoutWidget = QWidget(self.normalised_psi_box)
        self.formLayoutWidget.setObjectName(u"formLayoutWidget")
        self.formLayoutWidget.setGeometry(QRect(10, 30, 241, 191))
        self.formLayout_2 = QFormLayout(self.formLayoutWidget)
        self.formLayout_2.setObjectName(u"formLayout_2")
        self.formLayout_2.setSizeConstraint(QLayout.SetDefaultConstraint)
        self.formLayout_2.setContentsMargins(0, 0, 0, 0)
        self.psinorm_coreLabel = QLabel(self.formLayoutWidget)
        self.psinorm_coreLabel.setObjectName(u"psinorm_coreLabel")

        self.formLayout_2.setWidget(1, QFormLayout.LabelRole, self.psinorm_coreLabel)

        self.psinorm_coreDoubleSpinBox = QDoubleSpinBox(self.formLayoutWidget)
        self.psinorm_coreDoubleSpinBox.setObjectName(u"psinorm_coreDoubleSpinBox")

        self.formLayout_2.setWidget(1, QFormLayout.FieldRole, self.psinorm_coreDoubleSpinBox)

        self.psinorm_solLabel = QLabel(self.formLayoutWidget)
        self.psinorm_solLabel.setObjectName(u"psinorm_solLabel")

        self.formLayout_2.setWidget(2, QFormLayout.LabelRole, self.psinorm_solLabel)

        self.psinorm_solDoubleSpinBox = QDoubleSpinBox(self.formLayoutWidget)
        self.psinorm_solDoubleSpinBox.setObjectName(u"psinorm_solDoubleSpinBox")

        self.formLayout_2.setWidget(2, QFormLayout.FieldRole, self.psinorm_solDoubleSpinBox)


        self.verticalLayout_5.addWidget(self.normalised_psi_box)

        self.unnormalised_psi_box = QGroupBox(self.equilibrium_tab)
        self.unnormalised_psi_box.setObjectName(u"unnormalised_psi_box")
        self.formLayoutWidget_2 = QWidget(self.unnormalised_psi_box)
        self.formLayoutWidget_2.setObjectName(u"formLayoutWidget_2")
        self.formLayoutWidget_2.setGeometry(QRect(10, 30, 241, 201))
        self.formLayout_4 = QFormLayout(self.formLayoutWidget_2)
        self.formLayout_4.setObjectName(u"formLayout_4")
        self.formLayout_4.setSizeConstraint(QLayout.SetMaximumSize)
        self.formLayout_4.setContentsMargins(0, 0, 0, 0)
        self.psi_coreLabel = QLabel(self.formLayoutWidget_2)
        self.psi_coreLabel.setObjectName(u"psi_coreLabel")

        self.formLayout_4.setWidget(0, QFormLayout.LabelRole, self.psi_coreLabel)

        self.psi_coreDoubleSpinBox = QDoubleSpinBox(self.formLayoutWidget_2)
        self.psi_coreDoubleSpinBox.setObjectName(u"psi_coreDoubleSpinBox")

        self.formLayout_4.setWidget(0, QFormLayout.FieldRole, self.psi_coreDoubleSpinBox)

        self.psi_solLabel = QLabel(self.formLayoutWidget_2)
        self.psi_solLabel.setObjectName(u"psi_solLabel")

        self.formLayout_4.setWidget(1, QFormLayout.LabelRole, self.psi_solLabel)

        self.psi_solDoubleSpinBox = QDoubleSpinBox(self.formLayoutWidget_2)
        self.psi_solDoubleSpinBox.setObjectName(u"psi_solDoubleSpinBox")

        self.formLayout_4.setWidget(1, QFormLayout.FieldRole, self.psi_solDoubleSpinBox)


        self.verticalLayout_5.addWidget(self.unnormalised_psi_box)


        self.horizontalLayout_2.addLayout(self.verticalLayout_5)

        self.equilibrium_plotting_area = QFrame(self.equilibrium_tab)
        self.equilibrium_plotting_area.setObjectName(u"equilibrium_plotting_area")
        self.equilibrium_plotting_area.setFrameShape(QFrame.StyledPanel)
        self.equilibrium_plotting_area.setFrameShadow(QFrame.Raised)

        self.horizontalLayout_2.addWidget(self.equilibrium_plotting_area)

        self.tabWidget.addTab(self.equilibrium_tab, "")
        self.mesh_tab = QWidget()
        self.mesh_tab.setObjectName(u"mesh_tab")
        self.verticalLayout_3 = QVBoxLayout(self.mesh_tab)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.search_bar = QLineEdit(self.mesh_tab)
        self.search_bar.setObjectName(u"search_bar")

        self.verticalLayout_2.addWidget(self.search_bar)

        self.options_form = QTableWidget(self.mesh_tab)
        if (self.options_form.columnCount() < 2):
            self.options_form.setColumnCount(2)
        __qtablewidgetitem = QTableWidgetItem()
        self.options_form.setHorizontalHeaderItem(0, __qtablewidgetitem)
        __qtablewidgetitem1 = QTableWidgetItem()
        self.options_form.setHorizontalHeaderItem(1, __qtablewidgetitem1)
        self.options_form.setObjectName(u"options_form")

        self.verticalLayout_2.addWidget(self.options_form)

        self.horizontalLayout_3 = QHBoxLayout()
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.options_file_label = QLabel(self.mesh_tab)
        self.options_file_label.setObjectName(u"options_file_label")

        self.horizontalLayout_3.addWidget(self.options_file_label)

        self.options_file_line_edit = QLineEdit(self.mesh_tab)
        self.options_file_line_edit.setObjectName(u"options_file_line_edit")

        self.horizontalLayout_3.addWidget(self.options_file_line_edit)

        self.options_file_browse_button = QPushButton(self.mesh_tab)
        self.options_file_browse_button.setObjectName(u"options_file_browse_button")

        self.horizontalLayout_3.addWidget(self.options_file_browse_button)

        self.geqdsk_file_label = QLabel(self.mesh_tab)
        self.geqdsk_file_label.setObjectName(u"geqdsk_file_label")

        self.horizontalLayout_3.addWidget(self.geqdsk_file_label)

        self.geqdsk_file_line_edit = QLineEdit(self.mesh_tab)
        self.geqdsk_file_line_edit.setObjectName(u"geqdsk_file_line_edit")

        self.horizontalLayout_3.addWidget(self.geqdsk_file_line_edit)

        self.geqdsk_file_browse_button = QPushButton(self.mesh_tab)
        self.geqdsk_file_browse_button.setObjectName(u"geqdsk_file_browse_button")

        self.horizontalLayout_3.addWidget(self.geqdsk_file_browse_button)

        self.run_button = QPushButton(self.mesh_tab)
        self.run_button.setObjectName(u"run_button")

        self.horizontalLayout_3.addWidget(self.run_button)

        self.write_grid_button = QPushButton(self.mesh_tab)
        self.write_grid_button.setObjectName(u"write_grid_button")

        self.horizontalLayout_3.addWidget(self.write_grid_button)


        self.verticalLayout_2.addLayout(self.horizontalLayout_3)


        self.horizontalLayout.addLayout(self.verticalLayout_2)

        self.plottingArea = QFrame(self.mesh_tab)
        self.plottingArea.setObjectName(u"plottingArea")

        self.horizontalLayout.addWidget(self.plottingArea)


        self.verticalLayout_3.addLayout(self.horizontalLayout)

        self.tabWidget.addTab(self.mesh_tab, "")

        self.verticalLayout.addWidget(self.tabWidget)

        Hypnotoad.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(Hypnotoad)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1215, 22))
        self.menu_File = QMenu(self.menubar)
        self.menu_File.setObjectName(u"menu_File")
        self.menu_Help = QMenu(self.menubar)
        self.menu_Help.setObjectName(u"menu_Help")
        self.menu_Mesh = QMenu(self.menubar)
        self.menu_Mesh.setObjectName(u"menu_Mesh")
        Hypnotoad.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(Hypnotoad)
        self.statusbar.setObjectName(u"statusbar")
        Hypnotoad.setStatusBar(self.statusbar)
        self.toolBar = QToolBar(Hypnotoad)
        self.toolBar.setObjectName(u"toolBar")
        Hypnotoad.addToolBar(Qt.TopToolBarArea, self.toolBar)

        self.menubar.addAction(self.menu_File.menuAction())
        self.menubar.addAction(self.menu_Mesh.menuAction())
        self.menubar.addAction(self.menu_Help.menuAction())
        self.menu_File.addAction(self.action_New)
        self.menu_File.addAction(self.action_Open)
        self.menu_File.addAction(self.action_Save)
        self.menu_File.addAction(self.action_Save_as)
        self.menu_File.addAction(self.action_Revert)
        self.menu_File.addSeparator()
        self.menu_File.addAction(self.action_Preferences)
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

        self.retranslateUi(Hypnotoad)

        self.tabWidget.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(Hypnotoad)
    # setupUi

    def retranslateUi(self, Hypnotoad):
        Hypnotoad.setWindowTitle(QCoreApplication.translate("Hypnotoad", u"MainWindow", None))
        self.action_New.setText(QCoreApplication.translate("Hypnotoad", u"&New", None))
#if QT_CONFIG(shortcut)
        self.action_New.setShortcut(QCoreApplication.translate("Hypnotoad", u"Ctrl+N", None))
#endif // QT_CONFIG(shortcut)
        self.action_Open.setText(QCoreApplication.translate("Hypnotoad", u"&Open", None))
#if QT_CONFIG(shortcut)
        self.action_Open.setShortcut(QCoreApplication.translate("Hypnotoad", u"Ctrl+O", None))
#endif // QT_CONFIG(shortcut)
        self.action_Save.setText(QCoreApplication.translate("Hypnotoad", u"&Save", None))
#if QT_CONFIG(tooltip)
        self.action_Save.setToolTip(QCoreApplication.translate("Hypnotoad", u"Save", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(shortcut)
        self.action_Save.setShortcut(QCoreApplication.translate("Hypnotoad", u"Ctrl+S", None))
#endif // QT_CONFIG(shortcut)
        self.action_Save_as.setText(QCoreApplication.translate("Hypnotoad", u"Save as", None))
#if QT_CONFIG(shortcut)
        self.action_Save_as.setShortcut(QCoreApplication.translate("Hypnotoad", u"Ctrl+Shift+S", None))
#endif // QT_CONFIG(shortcut)
        self.action_Quit.setText(QCoreApplication.translate("Hypnotoad", u"&Quit", None))
#if QT_CONFIG(shortcut)
        self.action_Quit.setShortcut(QCoreApplication.translate("Hypnotoad", u"Ctrl+Q", None))
#endif // QT_CONFIG(shortcut)
        self.action_About.setText(QCoreApplication.translate("Hypnotoad", u"&About", None))
        self.action_Run.setText(QCoreApplication.translate("Hypnotoad", u"&Run", None))
#if QT_CONFIG(shortcut)
        self.action_Run.setShortcut(QCoreApplication.translate("Hypnotoad", u"Ctrl+R", None))
#endif // QT_CONFIG(shortcut)
        self.action_Write_grid.setText(QCoreApplication.translate("Hypnotoad", u"&Write grid", None))
#if QT_CONFIG(shortcut)
        self.action_Write_grid.setShortcut(QCoreApplication.translate("Hypnotoad", u"Ctrl+W", None))
#endif // QT_CONFIG(shortcut)
        self.action_Revert.setText(QCoreApplication.translate("Hypnotoad", u"&Revert", None))
        self.action_Preferences.setText(QCoreApplication.translate("Hypnotoad", u"&Preferences...", None))
        self.eq_geqdsk_label.setText(QCoreApplication.translate("Hypnotoad", u"geqdsk file", None))
        self.eq_geqdsk_browse.setText(QCoreApplication.translate("Hypnotoad", u"Browse", None))
        self.normalised_psi_box.setTitle(QCoreApplication.translate("Hypnotoad", u"Normalised Psi", None))
        self.psinorm_coreLabel.setText(QCoreApplication.translate("Hypnotoad", u"core", None))
        self.psinorm_solLabel.setText(QCoreApplication.translate("Hypnotoad", u"sol", None))
        self.unnormalised_psi_box.setTitle(QCoreApplication.translate("Hypnotoad", u"Unormalised Psi", None))
        self.psi_coreLabel.setText(QCoreApplication.translate("Hypnotoad", u"core", None))
        self.psi_solLabel.setText(QCoreApplication.translate("Hypnotoad", u"sol", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.equilibrium_tab), QCoreApplication.translate("Hypnotoad", u"&Equilibrium", None))
        ___qtablewidgetitem = self.options_form.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("Hypnotoad", u"Name", None));
        ___qtablewidgetitem1 = self.options_form.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("Hypnotoad", u"Value", None));
        self.options_file_label.setText(QCoreApplication.translate("Hypnotoad", u"Options file", None))
        self.options_file_browse_button.setText(QCoreApplication.translate("Hypnotoad", u"Browse", None))
        self.geqdsk_file_label.setText(QCoreApplication.translate("Hypnotoad", u"geqdsk file", None))
        self.geqdsk_file_browse_button.setText(QCoreApplication.translate("Hypnotoad", u"Browse", None))
        self.run_button.setText(QCoreApplication.translate("Hypnotoad", u"Run", None))
        self.write_grid_button.setText(QCoreApplication.translate("Hypnotoad", u"Write Grid", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.mesh_tab), QCoreApplication.translate("Hypnotoad", u"&Mesh", None))
        self.menu_File.setTitle(QCoreApplication.translate("Hypnotoad", u"&File", None))
        self.menu_Help.setTitle(QCoreApplication.translate("Hypnotoad", u"&Help", None))
        self.menu_Mesh.setTitle(QCoreApplication.translate("Hypnotoad", u"&Mesh", None))
        self.toolBar.setWindowTitle(QCoreApplication.translate("Hypnotoad", u"toolBar", None))
    # retranslateUi
