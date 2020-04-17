from Qt.QtWidgets import (
    QFileDialog,
    QMainWindow,
    QCheckBox,
    QSpinBox,
    QDoubleSpinBox,
    QLineEdit,
)

import numbers
import os
import yaml

from .hypnotoad2_mainWindow import Ui_Hypnotoad2
from .matplotlib_widget import MatplotlibWidget
from ..cases import tokamak
from ..core.mesh import BoutMesh

def convert_python_type_to_qwidget(value):
    """
    Convert a python type into the appropriate Qt widget
    """
    if isinstance(value, bool):
        return QCheckBox
    if isinstance(value, numbers.Integral):
        return QSpinBox
    if isinstance(value, numbers.Real):
        return QDoubleSpinBox
    if isinstance(value, str):
        return QLineEdit
    return QLineEdit


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

        self.options = dict(tokamak.TokamakEquilibrium.default_options.items())

        for key, value in sorted(self.options.items()):
            self.add_options_widget(key, value)

    def add_options_widget(self, key, value):
        widget_type = convert_python_type_to_qwidget(value)

        widget = widget_type()

        if isinstance(value, bool):
            widget.setChecked(value)
            widget.stateChanged.connect(lambda state: self.options.update(key=value))
        elif isinstance(value, numbers.Integral):
            widget.setMaximum(100000)
            widget.setValue(value)
            widget.valueChanged.connect(lambda value: self.options.update(key=value))
        elif isinstance(value, numbers.Real):
            widget.setDecimals(8)
            widget.setRange(-1e300, 1e300)
            widget.setValue(value)
            widget.valueChanged.connect(lambda value: self.options.update(key=value))
        elif isinstance(value, str):
            widget.setText(value)
            widget.textChanged.connect(lambda text: self.options.update(key=text))
        else:
            widget.textChanged.connect(lambda text: self.options.update(key=text))

        widget.setObjectName(key)
        self.options_form_layout.addRow(key, widget)
        return widget

    def update_options(self):
        for key, value in self.options.items():
            widget_type = convert_python_type_to_qwidget(value)
            widget = self.options_form_layout.findChild(widget_type, key)

            if widget is None:
                widget = self.add_options_widget(key, value)

            if isinstance(widget, QCheckBox):
                widget.setChecked(value)
            elif isinstance(widget, QSpinBox):
                widget.setValue(value)
            elif isinstance(widget, QDoubleSpinBox):
                widget.setValue(value)
            elif isinstance(widget, QLineEdit):
                widget.setText(value)
            else:
                raise RuntimeError(
                    f"Unknown widget when trying to update options ({type(widget)})"
                )

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

    def read_options(self):
        self.statusbar.showMessage("Reading options", 2000)
        options_filename = self.options_file_line_edit.text()

        if options_filename:
            with open(options_filename, "r") as f:
                self.options = yaml.safe_load(f)

        self.update_options()

    def select_geqdsk_file(self):
        filename, _ = QFileDialog.getOpenFileName(self, "Open geqdsk file", ".")

        if (filename is None) or (filename == ""):
            return  # Cancelled
        if not os.path.exists(filename):
            self.write("Could not find " + filename)
            return

        self.geqdsk_file_line_edit.setText(filename)
        self.read_geqdsk()

    def read_geqdsk(self):
        self.statusbar.showMessage("Reading geqdsk", 2000)
        geqdsk_filename = self.geqdsk_file_line_edit.text()

        if not geqdsk_filename:
            raise ValueError("No geqdsk file given")

        with open(geqdsk_filename, "rt") as fh:
            self.eq = tokamak.read_geqdsk(fh, options=self.options)

        self.eq.plotPotential(ncontours=40, axis=self.plot_widget.axes)
        self.plot_widget.axes.plot(*self.eq.x_points[0], "rx")
        self.plot_widget.canvas.draw()

    def run(self):

        if not hasattr(self, "eq"):
            self.statusbar.showMessage("Missing equilibrium!")
            return

        mesh = BoutMesh(self.eq)
        mesh.geometry()
