"""
GUI for Hypnotoad2 using Qt

"""

import copy
import numbers
import os
import yaml

from Qt.QtWidgets import (
    QFileDialog,
    QMainWindow,
    QCheckBox,
    QSpinBox,
    QDoubleSpinBox,
    QLineEdit,
    QMessageBox,
    QCompleter,
)
from Qt.QtCore import Qt

from .hypnotoad2_mainWindow import Ui_Hypnotoad2
from .matplotlib_widget import MatplotlibWidget
from ..cases import tokamak
from ..core.mesh import BoutMesh
from ..__version__ import __version__


colours = {
    "red": "#aa0000",
}


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
    """A graphical interface for Hypnotoad2

    """

    def __init__(self):
        super().__init__(None)
        self.setupUi(self)

        try:
            self.menu_File.setToolTipsVisible(True)
            self.menu_Mesh.setToolTipsVisible(True)
            self.menu_Help.setToolTipsVisible(True)
        except AttributeError:
            pass

        self.plot_widget = MatplotlibWidget(self.plottingArea)

        def set_triggered(widget, function):
            widget.triggered.connect(function)
            if hasattr(function, "__doc__"):
                widget.setToolTip(function.__doc__.strip())

        def set_clicked(widget, function):
            widget.clicked.connect(function)
            if hasattr(function, "__doc__"):
                widget.setToolTip(function.__doc__.strip())

        set_clicked(self.geqdsk_file_browse_button, self.select_geqdsk_file)
        self.geqdsk_file_line_edit.editingFinished.connect(self.read_geqdsk)

        set_clicked(self.options_file_browse_button, self.select_options_file)
        self.options_file_line_edit.editingFinished.connect(self.read_options)

        set_clicked(self.run_button, self.run)
        set_triggered(self.action_Run, self.run)

        set_clicked(self.write_grid_button, self.write_grid)
        set_triggered(self.action_Write_grid, self.write_grid)
        self.write_grid_button.setEnabled(False)

        set_triggered(self.action_Revert, self.revert_options)
        set_triggered(self.action_Save, self.save_options)
        set_triggered(self.action_Save_as, self.save_options_as)
        set_triggered(self.action_New, self.new_options)
        set_triggered(self.action_Open, self.select_options_file)
        set_triggered(self.action_About, self.help_about)

        self.action_Quit.triggered.connect(self.close)

        self.default_options = dict(tokamak.TokamakEquilibrium.default_options.items())
        self.options = copy.deepcopy(self.default_options)
        self.filename = "Untitled.yml"

        for key, value in sorted(self.options.items()):
            self.add_options_widget(key, value)

        self.search_bar.setPlaceholderText("Search options...")
        self.search_bar.textChanged.connect(self.search_options_form)
        self.search_bar.setToolTip(self.search_options_form.__doc__.strip())
        self.search_bar_completer = QCompleter(self.options.keys())
        self.search_bar_completer.setCaseSensitivity(Qt.CaseInsensitive)
        self.search_bar.setCompleter(self.search_bar_completer)

    def help_about(self):
        """About Hypnotoad2

        """

        about_text = __doc__.strip()
        about_text += f"\nVersion : {__version__}"

        about_box = QMessageBox(self)
        about_box.setText(about_text)
        about_box.exec_()

    def revert_options(self):
        """Revert the current options to the loaded file, or defaults if no
        file loaded

        """

        self.statusbar.showMessage("Reverting options", 2000)
        self.options = copy.deepcopy(self.default_options)

        options_filename = self.options_file_line_edit.text()

        if options_filename:
            self.read_options()

        self.update_options_form()

    def new_options(self):
        """New set of options

        """

        self.options = copy.deepcopy(self.default_options)
        self.update_options_form()

    def save_options(self):
        """Save options to file

        """

        self.statusbar.showMessage("Saving...", 2000)

        if not self.filename or self.filename == "Untitled.yml":
            self.save_options_as()

        self.options_file_line_edit.setText(self.filename)

        non_null_options = {k: v for k, v in self.options.items() if v is not None}

        with open(self.filename, "w") as f:
            yaml.dump(non_null_options, f)

    def save_options_as(self):
        """Save options to file with new filename

        """

        self.filename, _ = QFileDialog.getSaveFileName(
            self, "Save grid to file", self.filename, filter="YAML file (*yml *yaml)",
        )

        self.save_options()

    def add_options_widget(self, key, value):
        """Take a key, value pair and add a row with the appropriate widget
        to the options form

        """
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

    def update_options_form(self):
        """Update the widget values in the options form, based on the current
        values in the options dict

        """
        for key, value in sorted(self.options.items()):
            widget_type = convert_python_type_to_qwidget(value)
            widget = self.findChild(widget_type, key)

            # If we didn't already know the type, then it would be a
            # QLineEdit instead of a more specific widget
            if widget is None:
                widget = self.findChild(QLineEdit, key)
                if widget is not None:
                    self.options_form_layout.removeRow(widget)
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

    def search_options_form(self, text):
        """Search for specific options

        """

        for key, value in self.options.items():
            widget_type = convert_python_type_to_qwidget(value)
            widget = self.findChild(widget_type, key)
            if widget is None:
                continue

            label = self.options_form_layout.labelForField(widget)

            if text.lower() in key.lower():
                widget.show()
                label.show()
            else:
                widget.hide()
                label.hide()

    def select_options_file(self):
        """Choose a Hypnotoad2 options file to load

        """

        filename, _ = QFileDialog.getOpenFileName(
            self, "Open options file", ".", filter="YAML file (*.yml *.yaml)"
        )

        if (filename is None) or (filename == ""):
            return  # Cancelled
        if not os.path.exists(filename):
            self.write("Could not find " + filename)
            return

        self.options_file_line_edit.setText(filename)
        self.filename = filename
        self.read_options()

    def read_options(self):
        """Read the options file

        """

        self.statusbar.showMessage("Reading options", 2000)
        options_filename = self.options_file_line_edit.text()

        if options_filename:
            with open(options_filename, "r") as f:
                self.options.update(yaml.safe_load(f))

        self.update_options_form()

    def select_geqdsk_file(self):
        """Choose a "geqdsk" equilibrium file to open

        """

        filename, _ = QFileDialog.getOpenFileName(self, "Open geqdsk file", ".")

        if (filename is None) or (filename == ""):
            return  # Cancelled
        if not os.path.exists(filename):
            self.write("Could not find " + filename)
            self.geqdsk_file_line_edit.setStyleSheet(
                f"QLineEdit {{ background-color: {colours['red']} }}"
            )
            return

        self.geqdsk_file_line_edit.setText(filename)
        self.geqdsk_file_line_edit.setStyleSheet("")

        self.read_geqdsk()

    def read_geqdsk(self):
        """Read the equilibrium file

        """

        self.statusbar.showMessage("Reading geqdsk", 2000)
        geqdsk_filename = self.geqdsk_file_line_edit.text()

        if not os.path.exists(geqdsk_filename):
            self.geqdsk_file_line_edit.setStyleSheet(
                f"QLineEdit {{ background-color : {colours['red']} }}"
            )
            self.statusbar.showMessage(
                f"Could not find equilibrium file '{geqdsk_filename}'"
            )
            return

        with open(geqdsk_filename, "rt") as fh:
            self.eq = tokamak.read_geqdsk(fh, options=self.options)

        self.eq.plotPotential(ncontours=40, axis=self.plot_widget.axes)
        for region in self.eq.regions.values():
            self.plot_widget.axes.plot(
                [p.R for p in region.points], [p.Z for p in region.points], "-o"
            )

        self.plot_widget.axes.plot(*self.eq.x_points[0], "rx")
        self.plot_widget.canvas.draw()

    def run(self):
        """Run Hypnotoad2 and generate the grid

        """

        if not hasattr(self, "eq"):
            self.statusbar.showMessage("Missing equilibrium file!")
            self.geqdsk_file_line_edit.setStyleSheet(
                f"QLineEdit {{ background-color: {colours['red']} }}"
            )
            return

        self.statusbar.showMessage("Running...")
        self.mesh = BoutMesh(self.eq)
        self.mesh.geometry()
        self.statusbar.showMessage("Done!", 2000)

        self.plot_widget._clean_axes()
        self.eq.plotPotential(ncontours=40, axis=self.plot_widget.axes)
        self.mesh.plotPoints(
            xlow=self.options.get("plot_xlow", True),
            ylow=self.options.get("plot_ylow", True),
            corners=self.options.get("plot_corners", True),
            ax=self.plot_widget.axes,
        )
        self.plot_widget.canvas.draw()

        self.write_grid_button.setEnabled(True)

    def write_grid(self):
        """Write generated mesh to file

        """

        if not hasattr(self, "mesh"):
            flags = QMessageBox.StandardButton.Ok
            message_box = QMessageBox.critical(
                self, "Error", "Can't write mesh to file; no mesh found!", flags
            )
            return

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save grid to file",
            self.options.get("grid_file", "bout.grd.nc"),
            filter="NetCDF (*nc)",
        )

        self.mesh.writeGridfile(filename)
