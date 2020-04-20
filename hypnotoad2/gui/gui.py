"""
GUI for Hypnotoad2 using Qt

"""

import ast
import copy
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
    QWidget,
    QTableWidgetItem,
    QHeaderView,
)
from Qt.QtCore import Qt

from .hypnotoad2_mainWindow import Ui_Hypnotoad2
from .matplotlib_widget import MatplotlibWidget
from ..cases import tokamak
from ..core.mesh import BoutMesh
from ..__version__ import __version__


COLOURS = {
    "red": "#aa0000",
}


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

        self.search_bar.setPlaceholderText("Search options...")
        self.search_bar.textChanged.connect(self.search_options_form)
        self.search_bar.setToolTip(self.search_options_form.__doc__.strip())
        self.search_bar_completer = QCompleter(self.options.keys())
        self.search_bar_completer.setCaseSensitivity(Qt.CaseInsensitive)
        self.search_bar.setCompleter(self.search_bar_completer)

        self.options_form.cellChanged.connect(self.options_form_changed)
        self.update_options_form()

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
        else:
            self.options_form.setRowCount(0)
            self.update_options_form()

            if hasattr(self, "eq"):
                self.read_geqdsk()

    def new_options(self):
        """New set of options

        """

        self.options = copy.deepcopy(self.default_options)
        self.options_form.setRowCount(0)
        self.update_options_form()

        if hasattr(self, "eq"):
            self.read_geqdsk()

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

    def update_options_form(self):
        """Update the widget values in the options form, based on the current
        values in the options dict

        """

        self.options_form.setSortingEnabled(False)
        self.options_form.cellChanged.disconnect(self.options_form_changed)
        self.options_form.setRowCount(len(self.options) + 1)

        for row, (key, value) in enumerate(sorted(self.options.items())):
            item = QTableWidgetItem(key)
            item.old_key = key
            self.options_form.setItem(row, 0, item)
            self.options_form.setItem(row, 1, QTableWidgetItem(str(value)))

        self.options_form.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.options_form.setSortingEnabled(True)
        self.options_form.cellChanged.connect(self.options_form_changed)

    def options_form_changed(self, row, column):
        """Change the options form from the widget table

        """

        item = self.options_form.item(row, column)

        if row == self.options_form.rowCount() - 1:
            self.options[item.text()] = None
            self.options_form.setRowCount(len(self.options) + 1)
            return

        if column == 0:
            key = item.text()
            if key == item.old_key:
                return

            self.options[key] = self.options.pop(item.old_key)
            item.old_key = key
        else:
            key = self.options_form.item(row, 0).text()
            self.options[key] = ast.literal_eval(item.text())
            if hasattr(self, "eq"):
                self.read_geqdsk()

    def search_options_form(self, text):
        """Search for specific options

        """

        for i in range(self.options_form.rowCount()):
            row = self.options_form.item(i, 0)

            matches = text.lower() in row.text().lower()
            self.options_form.setRowHidden(i, not matches)

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

        self.options_form.setRowCount(0)
        self.update_options_form()

        if hasattr(self, "eq"):
            self.read_geqdsk()

    def select_geqdsk_file(self):
        """Choose a "geqdsk" equilibrium file to open

        """

        filename, _ = QFileDialog.getOpenFileName(self, "Open geqdsk file", ".")

        if (filename is None) or (filename == ""):
            return  # Cancelled
        if not os.path.exists(filename):
            self.write("Could not find " + filename)
            self.geqdsk_file_line_edit.setStyleSheet(
                f"QLineEdit {{ background-color: {COLOURS['red']} }}"
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
                f"QLineEdit {{ background-color : {COLOURS['red']} }}"
            )
            self.statusbar.showMessage(
                f"Could not find equilibrium file '{geqdsk_filename}'"
            )
            return

        with open(geqdsk_filename, "rt") as f:
            # Need to take a copy so that read_geqdsk doesn't delete used keys
            self.eq = tokamak.read_geqdsk(f, options=copy.deepcopy(self.options))

        self.plot_widget.clear()
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
                f"QLineEdit {{ background-color: {COLOURS['red']} }}"
            )
            return

        self.statusbar.showMessage("Running...")
        self.mesh = BoutMesh(self.eq)
        self.mesh.geometry()
        self.statusbar.showMessage("Done!", 2000)

        self.plot_widget.clear()
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
            QMessageBox.critical(
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