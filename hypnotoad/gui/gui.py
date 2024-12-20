"""
GUI for Hypnotoad using Qt

"""

import ast
import copy
import func_timeout
import gc
import os
import pathlib
import textwrap
import traceback
import yaml

from Qt.QtWidgets import (
    QFileDialog,
    QMainWindow,
    QMessageBox,
    QCompleter,
    QTableWidgetItem,
    QHeaderView,
    QErrorMessage,
    QDialog,
)
from Qt.QtCore import Qt

from .hypnotoad_mainWindow import Ui_Hypnotoad
from .hypnotoad_preferences import Ui_Preferences
from .matplotlib_widget import MatplotlibWidget
from ..cases import tokamak
from ..core.mesh import BoutMesh
from ..core.equilibrium import SolutionError
from ..__init__ import __version__


COLOURS = {
    "red": "#aa0000",
}

DEFAULT_OPTIONS_FILENAME = "Untitled.yml"

# File type filters
YAML_FILTER = "YAML file (*.yml *.yaml)"
NETCDF_FILTER = "NetCDF (*nc)"

DEFAULT_OPTIONS = {
    "orthogonal": tokamak.TokamakEquilibrium.user_options_factory.defaults[
        "orthogonal"
    ].value
}

DEFAULT_GUI_OPTIONS = {
    "grid_file": "bout.grd.nc",
    "plot_flux": True,
    "plot_wall": True,
    "plot_centers": True,
    "plot_xlow": True,
    "plot_ylow": True,
    "plot_corners": True,
    "save_full_yaml": False,
    "plot_legend": True,
    "plot_gridlines": False,
    "plot_celledges": False,
    "plot_penalty": False,
}


def _table_item_edit_display(item):
    """Hide the "(default)" marker on table items in the options form"""
    default_marker = " (default)"
    if item.text().endswith(default_marker):
        item.setText(item.text()[: -len(default_marker)])


class HypnotoadGui(QMainWindow, Ui_Hypnotoad):
    """A graphical interface for Hypnotoad"""

    def __init__(self):
        super().__init__(None)
        self.setupUi(self)

        # Used in file dialogs
        self._current_dir = "."

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

        self.nonorthogonal_box.stateChanged.connect(self.set_nonorthogonal)

        set_clicked(self.regrid_button, self.regrid)
        set_triggered(self.action_Regrid, self.regrid)

        set_triggered(self.action_Revert, self.revert_options)
        set_triggered(self.action_Save, self.save_options)
        set_triggered(self.action_Save_as, self.save_options_as)
        set_triggered(self.action_New, self.new_options)
        set_triggered(self.action_Open, self.select_options_file)
        set_triggered(self.action_About, self.help_about)
        set_triggered(self.action_Preferences, self.open_preferences)

        # View updates trigger a re-plotting
        def trigger_replot():
            """
            Update settings and replot grid
            """
            self.updateGuiOptionsFromMenu()
            # Re-plot, keeping plot limits
            self.plot_grid(keep_limits=True)

        set_triggered(self.action_Flux, trigger_replot)
        set_triggered(self.action_Wall, trigger_replot)
        set_triggered(self.action_Centers, trigger_replot)
        set_triggered(self.action_Corners, trigger_replot)
        set_triggered(self.action_Xlow, trigger_replot)
        set_triggered(self.action_Ylow, trigger_replot)
        set_triggered(self.action_Lines, trigger_replot)
        set_triggered(self.action_Edges, trigger_replot)
        set_triggered(self.action_Legend, trigger_replot)
        set_triggered(self.action_Penalty, trigger_replot)
        set_triggered(self.action_Clear, self.clearPlot)

        self.action_Quit.triggered.connect(self.close)

        self.options = DEFAULT_OPTIONS
        self.gui_options = DEFAULT_GUI_OPTIONS
        self.filename = DEFAULT_OPTIONS_FILENAME

        self.updateMenuFromGuiOptions()

        self.search_bar.setPlaceholderText("Search options...")
        self.search_bar.textChanged.connect(self.search_options_form)
        self.search_bar.setToolTip(self.search_options_form.__doc__.strip())
        option_names = (
            set(BoutMesh.user_options_factory.defaults.keys())
            .union(set(tokamak.TokamakEquilibrium.user_options_factory.defaults.keys()))
            .union(
                set(
                    tokamak.TokamakEquilibrium.nonorthogonal_options_factory.defaults.keys()  # noqa: E501
                )
            )
        )
        self.search_bar_completer = QCompleter(option_names)
        self.search_bar_completer.setCaseSensitivity(Qt.CaseInsensitive)
        self.search_bar.setCompleter(self.search_bar_completer)

        self.options_form.cellChanged.connect(self.options_form_changed)
        self.options_form.itemDoubleClicked.connect(_table_item_edit_display)
        self.update_options_form()

    def updateMenuFromGuiOptions(self):
        """
        Updates menu items from gui_options
        """
        self.action_Flux.setChecked(self.gui_options["plot_flux"])
        self.action_Wall.setChecked(self.gui_options["plot_wall"])
        self.action_Centers.setChecked(self.gui_options["plot_centers"])
        self.action_Xlow.setChecked(self.gui_options["plot_xlow"])
        self.action_Ylow.setChecked(self.gui_options["plot_ylow"])
        self.action_Corners.setChecked(self.gui_options["plot_corners"])
        self.action_Legend.setChecked(self.gui_options["plot_legend"])
        self.action_Lines.setChecked(self.gui_options["plot_gridlines"])
        self.action_Edges.setChecked(self.gui_options["plot_celledges"])
        self.action_Penalty.setChecked(self.gui_options["plot_penalty"])

    def updateGuiOptionsFromMenu(self):
        """
        Update gui_options settings from menu
        """
        self.gui_options["plot_flux"] = self.action_Flux.isChecked()
        self.gui_options["plot_wall"] = self.action_Wall.isChecked()
        self.gui_options["plot_centers"] = self.action_Centers.isChecked()
        self.gui_options["plot_xlow"] = self.action_Xlow.isChecked()
        self.gui_options["plot_ylow"] = self.action_Ylow.isChecked()
        self.gui_options["plot_corners"] = self.action_Corners.isChecked()
        self.gui_options["plot_legend"] = self.action_Legend.isChecked()
        self.gui_options["plot_gridlines"] = self.action_Lines.isChecked()
        self.gui_options["plot_celledges"] = self.action_Edges.isChecked()
        self.gui_options["plot_penalty"] = self.action_Penalty.isChecked()

    def clearPlot(self):
        """
        Clear the grid plot
        """
        self.gui_options["plot_flux"] = False
        self.gui_options["plot_wall"] = False
        self.gui_options["plot_centers"] = False
        self.gui_options["plot_xlow"] = False
        self.gui_options["plot_ylow"] = False
        self.gui_options["plot_corners"] = False
        self.gui_options["plot_legend"] = False
        self.gui_options["plot_gridlines"] = False
        self.gui_options["plot_celledges"] = False
        self.gui_options["plot_penalty"] = False
        self.updateMenuFromGuiOptions()
        self.plot_grid(keep_limits=True)

    def close(self):
        # Delete and garbage-collect hypnotoad objects here so that any ParallelMap
        # instances get deleted if they exists. ParallelMap.__del__() calls
        # terminate on the worker processes. Needs to happen before program exits
        # or program will hang waiting for parallel workers to finish (which they
        # never do because they are waiting for another job in
        # ParallelMap.worker_run()).
        if hasattr(self, "eq"):
            del self.eq
        if hasattr(self, "mesh"):
            del self.mesh
        gc.collect()
        super().close()

    def help_about(self):
        """About Hypnotoad"""

        about_text = __doc__.strip()
        about_text += f"\nVersion : {__version__}"

        about_box = QMessageBox(self)
        about_box.setText(about_text)
        about_box.exec_()

    def open_preferences(self):
        """GUI preferences and settings"""
        preferences_window = Preferences(self)
        preferences_window.exec_()

    def revert_options(self):
        """Revert the current options to the loaded file, or defaults if no
        file loaded

        """

        self.statusbar.showMessage("Reverting options", 2000)
        self.options = DEFAULT_OPTIONS

        options_filename = self.options_file_line_edit.text()

        if options_filename:
            self.read_options()
        else:
            self.options_form.setRowCount(0)
            self.update_options_form()

    def new_options(self):
        """New set of options"""

        self.options = DEFAULT_OPTIONS
        self.options_form.setRowCount(0)
        self.update_options_form()

    def save_options(self):
        """Save options to file"""

        self.statusbar.showMessage("Saving...", 2000)

        if not self.filename or self.filename == DEFAULT_OPTIONS_FILENAME:
            self.save_options_as()

        if not self.filename:
            self.filename = DEFAULT_OPTIONS_FILENAME
            return

        self.options_file_line_edit.setText(self.filename)

        options_to_save = self.options
        if self.gui_options["save_full_yaml"]:
            mesh_options = BoutMesh.user_options_factory.create(self.options)
            eq_options = tokamak.TokamakEquilibrium.user_options_factory.create(
                self.options
            )
            nonorth_options = (
                tokamak.TokamakEquilibrium.nonorthogonal_options_factory.create(
                    self.options
                )
            )

            # This converts any numpy types to native Python using the tolist()
            # method of any numpy objects/types. Note this does return a scalar
            # and not a list for values that aren't arrays. Also remove any
            # private/magic keys
            options_to_save = {}
            for key, value in mesh_options.items():
                options_to_save[key] = getattr(value, "tolist", lambda: value)()
            for key, value in eq_options.items():
                options_to_save[key] = getattr(value, "tolist", lambda: value)()
            for key, value in nonorth_options.items():
                options_to_save[key] = getattr(value, "tolist", lambda: value)()

        with open(self.filename, "w") as f:
            yaml.dump(options_to_save, f)

    def save_options_as(self):
        """Save options to file with new filename"""

        if not self.filename:
            self.filename = DEFAULT_OPTIONS_FILENAME

        self.filename, _ = QFileDialog.getSaveFileName(
            self, "Save grid to file", self.filename, filter=YAML_FILTER
        )

        if not self.filename:
            return

        # If there was no extension, add one, unless the file already exists
        path = pathlib.Path(self.filename)
        if not path.exists() and path.suffix == "":
            self.filename += ".yml"

        self.save_options()

    def update_options_form(self):
        """Update the widget values in the options form, based on the current
        values in self.options

        """

        filtered_options = copy.deepcopy(self.options)

        filtered_defaults = dict(BoutMesh.user_options_factory.defaults)
        filtered_defaults.update(
            tokamak.TokamakEquilibrium.user_options_factory.defaults
        )
        filtered_defaults.update(
            tokamak.TokamakEquilibrium.nonorthogonal_options_factory.defaults
        )

        try:
            # evaluate filtered_defaults using the values in self.options, so that any
            # expressions get evaluated
            filtered_default_values = dict(
                BoutMesh.user_options_factory.create(self.options)
            )

            filtered_default_values.update(
                tokamak.TokamakEquilibrium.user_options_factory.create(self.options)
            )
            if not hasattr(self, "eq"):
                filtered_default_values.update(
                    tokamak.TokamakEquilibrium.nonorthogonal_options_factory.create(
                        self.options
                    )
                )
            else:
                # Use the object if it exists because some defaults are updated when the
                # Equilibrium is created
                filtered_default_values.update(
                    self.eq.nonorthogonal_options_factory.create(self.options)
                )
        except (ValueError, TypeError) as e:
            self._popup_error_message(e)
            return False

        # Skip options handled specially elsewhere
        filtered_options.pop("orthogonal", None)
        del filtered_defaults["orthogonal"]

        self.options_form.setSortingEnabled(False)
        self.options_form.cellChanged.disconnect(self.options_form_changed)
        self.options_form.setRowCount(len(filtered_defaults))

        for row, (key, value) in enumerate(sorted(filtered_defaults.items())):
            item0 = QTableWidgetItem(key)
            item0.setFlags(item0.flags() & ~Qt.ItemIsEditable)
            item0.setToolTip(value.doc)
            self.options_form.setItem(row, 0, item0)
            if key in filtered_options:
                value_to_set = str(filtered_options[key])
            else:
                value_to_set = f"{filtered_default_values[key]} (default)"
            item1 = QTableWidgetItem(value_to_set)
            item1.setToolTip(textwrap.fill(value.doc))
            self.options_form.setItem(row, 1, item1)

        self.options_form.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.options_form.setSortingEnabled(True)
        self.options_form.cellChanged.connect(self.options_form_changed)

        return True

    def options_form_changed(self, row, column):
        """Change the options form from the widget table"""

        item = self.options_form.item(row, column)

        if column == 0:
            # column 0 is not editable, so this should not be possible
            raise ValueError("Not allowed to change option names")
        else:
            key = self.options_form.item(row, 0).text()
            has_old_value = key in self.options
            if has_old_value:
                old_value = self.options[key]

            if item.text() == "":
                # Reset to default
                # Might be better to just keep the old value if nothing is passed, but
                # don't know how to get that
                if key in self.options:
                    del self.options[key]
            else:
                try:
                    self.options[key] = ast.literal_eval(item.text())
                except (ValueError, SyntaxError):
                    # Value might have been string with some illegal character in
                    self.options[key] = item.text()

        if not self.update_options_form():
            # Some error occured, reset the original value
            if has_old_value:
                self.options[key] = old_value
            else:
                del self.options[key]
            self.update_options_form()

    def search_options_form(self, text):
        """Search for specific options"""

        for i in range(self.options_form.rowCount()):
            row = self.options_form.item(i, 0)

            matches = text.lower() in row.text().lower()
            self.options_form.setRowHidden(i, not matches)

    def select_options_file(self):
        """Choose a Hypnotoad options file to load"""

        filename, _ = QFileDialog.getOpenFileName(
            self, "Open options file", self._current_dir, filter=YAML_FILTER
        )

        if (filename is None) or (filename == ""):
            return  # Cancelled
        if not os.path.exists(filename):
            self.write("Could not find " + filename)
            return

        # Record the directory so user doesn't have to navigate again
        self._current_dir = os.path.dirname(filename)

        self.options_file_line_edit.setText(filename)
        self.filename = filename
        self.read_options()
        self.nonorthogonal_box.setChecked(not self.options["orthogonal"])

    def read_options(self):
        """Read the options file"""

        self.statusbar.showMessage("Reading options", 2000)
        options_filename = self.options_file_line_edit.text()

        # Save the existing options in case there is an error loading the options file
        original_options = self.options

        if options_filename:
            with open(options_filename, "r") as f:
                file_options = yaml.safe_load(f)
            # Ensure that all default options keys are in self.options
            # Note: update mutates dict, returns None
            self.options = DEFAULT_OPTIONS.copy()
            self.options.update(file_options)

        possible_options = (
            [opt for opt in tokamak.TokamakEquilibrium.user_options_factory.defaults]
            + [
                opt
                for opt in tokamak.TokamakEquilibrium.nonorthogonal_options_factory.defaults  # noqa: E501
            ]
            + [opt for opt in BoutMesh.user_options_factory.defaults]
        )
        unused_options = [opt for opt in self.options if opt not in possible_options]
        if unused_options != []:
            short_filename = pathlib.Path(options_filename).parts[-1]
            self._popup_error_message(
                f"Error: There were options in the input file that are not used: "
                f"{unused_options}. Cannot load {short_filename}."
            )
            self.options = original_options

        self.options_form.setRowCount(0)
        self.update_options_form()

    def select_geqdsk_file(self):
        """Choose a "geqdsk" equilibrium file to open"""

        filename, _ = QFileDialog.getOpenFileName(
            self, "Open geqdsk file", self._current_dir
        )

        if (filename is None) or (filename == ""):
            return  # Cancelled
        if not os.path.exists(filename):
            self.write("Could not find " + filename)
            self.geqdsk_file_line_edit.setStyleSheet(
                f"QLineEdit {{ background-color: {COLOURS['red']} }}"
            )
            return

        # Record the directory so user doesn't have to navigate again
        self._current_dir = os.path.dirname(filename)

        self.geqdsk_file_line_edit.setText(filename)
        self.geqdsk_file_line_edit.setStyleSheet("")

        self.read_geqdsk()

    def read_geqdsk(self):
        """Read the equilibrium file"""

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

        try:
            with open(geqdsk_filename, "rt") as f:
                self.eq = tokamak.read_geqdsk(
                    f,
                    settings=copy.deepcopy(self.options),
                    nonorthogonal_settings=copy.deepcopy(self.options),
                )
            try:
                # If there was an error in tokamak.read_geqdsk(), it may return both
                # the TokamakEquilibrium and an error, otherwise it would return
                # just a TokamakEquilibrium.
                self.eq, e = self.eq
                self._popup_error_message(e)
            except TypeError:
                # No error, so self.eq is already the TokamakEquilibrium object.
                pass
        except (ValueError, RuntimeError, func_timeout.FunctionTimedOut) as e:
            self._popup_error_message(e)

        self.update_options_form()

        # Delete mesh if it exists, since we have a new self.eq object
        if hasattr(self, "mesh"):
            del self.mesh
        self.regrid_button.setEnabled(False)
        self.action_Regrid.setEnabled(False)

        self.plot_grid()

        self.nonorthogonal_box.setChecked(not self.options["orthogonal"])

    def run(self):
        """Run Hypnotoad and generate the grid"""

        if not hasattr(self, "eq"):
            self.statusbar.showMessage("Missing equilibrium file!")
            self.geqdsk_file_line_edit.setStyleSheet(
                f"QLineEdit {{ background-color: {COLOURS['red']} }}"
            )
            return

        # Call read_geqdsk to recreate self.eq object in case any settings needed in
        # __init__ have been changed
        self.read_geqdsk()

        self.statusbar.showMessage("Running...")
        try:
            self.mesh = BoutMesh(self.eq, self.options)
        except (ValueError, SolutionError, func_timeout.FunctionTimedOut) as e:
            self._popup_error_message(e)
            return

        self.mesh.calculateRZ()
        self.statusbar.showMessage("Done!", 2000)

        self.plot_grid()

        self.write_grid_button.setEnabled(True)
        self.regrid_button.setEnabled(not self.options["orthogonal"])
        self.action_Regrid.setEnabled(not self.options["orthogonal"])

    def set_nonorthogonal(self, state):
        state = bool(state)
        self.options["orthogonal"] = not state
        self.update_options_form()

    def regrid(self):
        """Regrid a nonorthogonal grid after spacing settings are changed"""

        if not hasattr(self, "mesh"):
            self.statusbar.showMessage("Generate grid first!")
            self.geqdsk_file_line_edit.setStyleSheet(
                f"QLineEdit {{ background-color: {COLOURS['red']} }}"
            )
            return

        self.statusbar.showMessage("Running...")

        try:
            self.mesh.redistributePoints(self.options)
            self.mesh.calculateRZ()
        except (ValueError, TypeError, func_timeout.FunctionTimedOut) as e:
            self._popup_error_message(e)
            return

        self.statusbar.showMessage("Done!", 2000)

        self.plot_grid(keep_limits=True)

    def write_grid(self):
        """Write generated mesh to file"""

        # Create all the geometrical quantities
        try:
            self.mesh.geometry()
        except (
            ValueError,
            TypeError,
            func_timeout.FunctionTimedOut,
            SolutionError,
        ) as e:
            self._popup_error_message(e)
            return

        if not hasattr(self, "mesh"):
            flags = QMessageBox.StandardButton.Ok
            QMessageBox.critical(
                self, "Error", "Can't write mesh to file; no mesh found!", flags
            )
            return

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save grid to file",
            os.path.join(self._current_dir, self.gui_options["grid_file"]),
            filter=NETCDF_FILTER,
        )

        if not filename:
            return

        # If there was no extension, add one, unless the file already exists
        path = pathlib.Path(self.filename)
        if not path.exists() and path.suffix == "":
            self.filename += ".nc"

        self.mesh.writeGridfile(filename)

    def plot_grid(self, *, keep_limits: bool = False):
        """
        Re-plot the grid and equilibrium

        # Arguments

        keep_limits: bool
            Keep the axis limits of the plot unchanged?
        """

        if keep_limits:
            xlim = self.plot_widget.axes.get_xlim()
            ylim = self.plot_widget.axes.get_ylim()

        self.plot_widget.clear(keep_limits=keep_limits)

        if hasattr(self, "eq"):
            if self.gui_options["plot_flux"]:
                self.eq.plotPotential(ncontours=40, axis=self.plot_widget.axes)
            if self.gui_options["plot_wall"]:
                self.eq.plotWall(axis=self.plot_widget.axes)

        if hasattr(self, "mesh"):
            # mesh exists, so plot the grid points
            self.mesh.plotPoints(
                centers=self.gui_options["plot_centers"],
                xlow=self.gui_options["plot_xlow"],
                ylow=self.gui_options["plot_ylow"],
                corners=self.gui_options["plot_corners"],
                ax=self.plot_widget.axes,
            )

            if self.gui_options["plot_gridlines"]:
                self.mesh.plotGridLines(ax=self.plot_widget.axes)

            if self.gui_options["plot_celledges"]:
                self.mesh.plotGridCellEdges(ax=self.plot_widget.axes)

            if self.gui_options["plot_penalty"]:
                self.mesh.plotPenaltyMask(ax=self.plot_widget.axes)

        elif hasattr(self, "eq"):
            # no mesh, but do have equilibrium, so plot separatrices
            for region in self.eq.regions.values():
                self.plot_widget.axes.plot(
                    [p.R for p in region.points],
                    [p.Z for p in region.points],
                    marker="o",
                )
            self.plot_widget.axes.plot(*self.eq.x_points[0], "rx")

        if keep_limits:
            self.plot_widget.axes.set_xlim(xlim)
            self.plot_widget.axes.set_ylim(ylim)

        if not self.gui_options["plot_legend"]:
            self.plot_widget.axes.legend().set_visible(False)

        self.plot_widget.canvas.draw()

    def _popup_error_message(self, error):
        error_message = QErrorMessage()
        error_message.showMessage(str(error) + "<br><br>" + traceback.format_exc())
        error_message.exec_()


class Preferences(QDialog, Ui_Preferences):
    """Dialog box for editing Hypnotoad preferences"""

    def __init__(self, parent):
        super().__init__(parent)
        self.setupUi(self)
        self.parent = parent

        self.defaultGridFileNameLineEdit.setText(self.parent.gui_options["grid_file"])
        self.plotXlowCheckBox.setChecked(self.parent.gui_options["plot_xlow"])
        self.plotYlowCheckBox.setChecked(self.parent.gui_options["plot_ylow"])
        self.plotCornersCheckBox.setChecked(self.parent.gui_options["plot_corners"])
        self.saveFullYamlCheckBox.setChecked(self.parent.gui_options["save_full_yaml"])

    def accept(self):
        self.parent.gui_options["grid_file"] = self.defaultGridFileNameLineEdit.text()
        self.parent.gui_options["plot_xlow"] = self.plotXlowCheckBox.isChecked()
        self.parent.gui_options["plot_ylow"] = self.plotYlowCheckBox.isChecked()
        self.parent.gui_options["plot_corners"] = self.plotCornersCheckBox.isChecked()
        self.parent.gui_options["save_full_yaml"] = (
            self.saveFullYamlCheckBox.isChecked()
        )

        self.parent.updateMenuFromGuiOptions()

        self.parent.plot_grid()

        super().accept()
