from Qt.QtWidgets import (
    QAbstractItemView,
    QAction,
    QFileDialog,
    QMainWindow,
    QMenu,
    QMessageBox,
    QStyle,
    QTableWidgetItem,
    QTreeWidgetItem,
    QWidget,
    QCheckBox,
    QSpinBox,
    QDoubleSpinBox,
    QLineEdit
)

import numbers
import os
import yaml

from .hypnotoad2_mainWindow import Ui_Hypnotoad2
from .matplotlib_widget import MatplotlibWidget
from ..cases import tokamak
from ..core.mesh import BoutMesh


class HypnotoadGui(QMainWindow, Ui_Hypnotoad2):
    def __init__(self):
        super().__init__(None)
        self.setupUi(self)

        self.plot_widget = MatplotlibWidget(self.plottingArea)

        self.geqdsk_file_browse_button.clicked.connect(self.select_geqdsk_file)
        self.geqdsk_file_line_edit.editingFinished.connect(self.read_geqdsk)
        self.options_file_browse_button.clicked.connect(self.select_options_file)
        self.options_file_line_edit.editingFinished.connect(self.read_options)

        self.run_button.clicked.connect(self.run)

        for key, value in tokamak.TokamakEquilibrium.default_options.items():
            if isinstance(value, bool):
                widget = QCheckBox()
                widget.setChecked(value)
            elif isinstance(value, numbers.Integral):
                widget = QSpinBox()
                widget.setValue(value)
            elif isinstance(value, numbers.Real):
                widget = QDoubleSpinBox()
                widget.setValue(value)
            elif isinstance(value, str):
                widget = QLineEdit()
                widget.setText(value)
            else:
                widget = QLineEdit()

            widget.setObjectName(key)
            self.options_form_layout.addRow(key, widget)

    def select_geqdsk_file(self):
        filename, _ = QFileDialog.getOpenFileName(self, "Open geqdsk file", ".")

        if (filename is None) or (filename == ""):
            return  # Cancelled
        if not os.path.exists(filename):
            self.write("Could not find " + filename)
            return

        self.geqdsk_file_line_edit.setText(filename)
        self.read_geqdsk()

    def select_options_file(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open options file", ".", filter="YAML file (*.yml, *.yaml)"
        )

        if (filename is None) or (filename == ""):
            return  # Cancelled
        if not os.path.exists(filename):
            self.write("Could not find " + filename)
            return

        self.options_file_line_edit.setText(filename)
        self.read_options()

    def read_geqdsk(self):
        self.statusbar.showMessage("Reading geqdsk", 2000)

    def read_options(self):
        self.statusbar.showMessage("Reading options", 2000)

    def run(self):
        options_filename = self.options_file_line_edit.text()

        if options_filename:
            with open(options_filename, "r") as f:
                options = yaml.safe_load(f)
        else:
            options = {}

        geqdsk_filename = self.geqdsk_file_line_edit.text()

        if not geqdsk_filename:
            raise ValueError("No geqdsk file given")

        with open(geqdsk_filename, "rt") as fh:
            eq = tokamak.read_geqdsk(fh, options=options)

        mesh = BoutMesh(eq)
        mesh.geometry()

        eq.plotPotential(ncontours=40, axis=self.plot_widget.axes)
        self.plot_widget.axes.plot(*eq.x_points[0], "rx")
        self.plot_widget.canvas.draw()
