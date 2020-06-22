# Initial code taken from
# http://stackoverflow.com/questions/6723527/getting-pyside-to-work-with-matplotlib
# Additional bits from https://gist.github.com/jfburkhart/2423179

import matplotlib
from matplotlib.figure import Figure
from Qt.QtWidgets import QVBoxLayout

try:
    matplotlib.use("Qt5Agg")
except ImportError:
    # Continue for now, so that hypnotoad-gui -h works even on systems without a
    # display. Useful for testing existence of the command for conda-forge.
    import warnings

    warnings.warn("Failed to load Qt5Agg backend, plotting widget will fail")
else:
    from matplotlib.backends.backend_qt5agg import (  # noqa: E402
        FigureCanvasQTAgg as FigureCanvas,
        NavigationToolbar2QT as NavigationToolbar,
    )


class MatplotlibWidget:
    def __init__(self, parent):

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(parent)
        self.mpl_toolbar = NavigationToolbar(self.canvas, parent)
        self.axes = self.figure.add_subplot(111)

        self.grid_layout = QVBoxLayout()
        self.grid_layout.addWidget(self.canvas)
        self.grid_layout.addWidget(self.mpl_toolbar)
        parent.setLayout(self.grid_layout)

        self.clear()

    def clear(self, *, keep_limits=False):
        """
        Make sure the figure is in a nice state
        """

        if keep_limits:
            # slightly hacky way to clear axes, but prevents axis limits being reset when
            # we redraw
            for artist in self.axes.lines + self.axes.collections:
                artist.remove()
            self.axes.set_prop_cycle(None)
            return

        self.axes.clear()
        self.figure.clear()
        self.axes = self.figure.add_subplot(111)
        self.axes.grid(True)
        # Reset to some hardcoded default values
        self.figure.subplots_adjust(
            left=0.125, right=0.9, top=0.9, bottom=0.1, wspace=0.2, hspace=0.2
        )

        self.axes.set_xlabel("R", fontsize=16)
        self.axes.set_ylabel("Z", rotation="horizontal", fontsize=16)
        self.canvas.draw()
