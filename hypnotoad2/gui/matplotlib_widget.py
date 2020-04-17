#
# Initial code taken from  http://stackoverflow.com/questions/6723527/getting-pyside-to-work-with-matplotlib
# Additional bits from https://gist.github.com/jfburkhart/2423179
# zplot zoom effect from http://matplotlib.org/examples/pylab_examples/axes_zoom_effect.html
#

import matplotlib

from Qt import __qt_version__

matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
                                                NavigationToolbar2QT as NavigationToolbar)

from matplotlib.figure import Figure
from matplotlib.transforms import Bbox, TransformedBbox, \
    blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
    BboxConnectorPatch

from matplotlib.pyplot import setp

from Qt.QtWidgets import *
from Qt.QtCore import Qt

import numpy as np
import warnings


class MatplotlibWidget():

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

        self.callback_id = None

        warnings.filterwarnings("ignore", "Attempting to set identical left == right.*", UserWarning)

    def _clean_axes(self):
        """
        Make sure the figure is in a nice state
        """
        # Get rid of any extra axes
        if isinstance(self.axes, list):
            for axes in self.axes:
                del axes
        else:
            self.axes.clear()

        self.figure.clear()
        self.axes = self.figure.add_subplot(111)
        self.axes.grid(True)
        # Reset to some hardcoded default values
        self.figure.subplots_adjust(left=0.125, right=0.9, top=0.9, bottom=0.1,
                                    wspace=0.2, hspace=0.2)
        # Remove any event callbacks
        if self.callback_id:
            try:
                self.figure.canvas.mpl_disconnect(self.callback_id)
            except TypeError:
                for callback_id in self.callback_id:
                    self.figure.canvas.mpl_disconnect(callback_id)
                self.callback_id = None

    def plot(self, *args):
        """
        Make multiple plots
        """
        nplots = len(args)
        self._clean_axes()

        for plotnum, p in enumerate(args, 1):
            # For each plot
            if plotnum == 1:
                self.axes = self.figure.add_subplot(nplots, 1, plotnum)
                ax = self.axes
            else:
                self.axes = self.figure.add_subplot(nplots, 1, plotnum, sharex=ax)

            try:
                # Assume each p is a list of data items to be overplotted
                for data in p:
                    label = data.desc
                    if label == "":
                        label = data.name + " (" + data.units + ") " + data.source

                    time = data.time
                    if time is None:
                        if len(data.dim) != 1:
                            print(data.dim)
                            raise ValueError("Cannot plot '"+data.label+"' as it has too many dimensions")
                        time = data.dim[0].data
                    if hasattr(time, "data"):
                        time = time.data

                    self.axes.plot(time, data.data, label=label)

                    loc = "upper right" if len(time) > 10000 else "best"

                # Y label from last plot
                ylabel = data.desc
                if ylabel == "":
                    ylabel = data.label
                    if data.units != "":
                        ylabel += " ("+data.units+")"
                self.axes.set_ylabel(ylabel)

                if plotnum == nplots:
                    x_dim = p[0].dim[p[0].order]
                    xlabel = x_dim.label or x_dim.name
                    if x_dim.units:
                        xlabel = "{} ({})".format(xlabel, x_dim.units)

                    self.axes.set_xlabel(xlabel)

                self.axes.legend(loc=loc)
            except TypeError:
                # p not iterable, so just plot item
                time = p.time.data
                xlabel = p.dim[p.order].label
                if time is None:
                    if len(p.dim) != 1:
                        raise ValueError("Cannot plot '"+p.label+"' as it has too many dimensions")
                    time = p.dim[0].data
                    xlabel = p.dim[0].label

                data = p.data
                # Check the size of the array
                size = len(time)
                if size > 10000:
                    fac = int(len(time) / 10000)
                    time = time[::fac]
                    data = data[::fac]
                    print("Warning: too many samples (%d). Down-sampling to %d points" % (size, len(time)))

                self.axes.plot(time, data)

                ylabel = p.desc
                if ylabel == "":
                    ylabel = p.label
                    if p.units != "":
                        ylabel += " ("+p.units+")"
                self.axes.set_ylabel(ylabel)
                if plotnum == 0:
                    self.axes.set_xlabel(xlabel)

        self.figure.subplots_adjust(hspace=0.001)

        self.canvas.draw()

    def plotxy(self, x, y):
        """
        Plot one variable against another
        """
        self.axes.clear()
        self.figure.clear()
        self.axes = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.08)
        self.axes.plot(x.data, y.data)
        self.axes.set_xlabel(x.name)
        self.axes.set_ylabel(y.name)
        self.canvas.draw()

    def oplot(self, *args):
        """
        Make an overplot from multiple traces
        """
        ntraces = len(args)
        if len(args) == 1:
            self.plot(*args)

        self._clean_axes()

        self.axes = self.figure.add_subplot(1, 1, 1)
        for tracenum, trace in enumerate(args):
            try:
                for data in trace:
                    label = data.desc
                    if label == "":
                        label = data.name + " (" + data.units + ") " + data.source

                    time = data.time
                    if time is None:
                        if len(data.dim) != 1:
                            print(data.dim)
                            raise ValueError("Cannot plot '"+data.label+"' as it has too many dimensions")
                        time = data.dim[0].data
                    if hasattr(time, "data"):
                        time = time.data
                    loc = "upper right" if len(time) > 10000 else "best"

                    self.axes.plot(time, data.data, label=label)
                if tracenum == 0:
                    ylabel = data.desc
                    if ylabel == "":
                        ylabel = data.label
                        if data.units != "":
                            ylabel += " ("+data.units+") "
                    self.axes.set_ylabel(ylabel)
                    x_dim = data.dim[data.order]
                    xlabel = x_dim.label or x_dim.name
                    if x_dim.units:
                        xlabel = "{} ({})".format(xlabel, x_dim.units)

                    self.axes.set_xlabel(xlabel)

                self.axes.legend(loc=loc)
            except TypeError:
                #Trace not iterable
                time = trace.time.data
                xlabel = trace.dim[trace.order].label
                if time is None:
                    if len(trace.dim) != 1:
                        raise ValueError("Cannot plot '"+trace.label+"' as it has too many dimensions")
                    time = trace.dim[0].data
                    xlabel = trace.dim[0].label

                data = trace.data
                # Check array size
                size = len(time)
                if size > 10000:
                    fac = int(len(time) / 10000)
                    time = time[::fac]
                    data = data[::fac]
                    print("Warning: too many samples (%d). Down-sampling to %d points" % (size, len(time)))

                self.axes.plot(time, data)

                ylabel = trace.desc
                if ylabel == "":
                    ylabel = trace.label
                    if trace.units != "":
                        ylabel += " ("+trace.units+") "
                self.axes.set_ylabel(ylabel)

        self.canvas.draw()

    def mplot(self, *args):
        """
        Make a custom number of plots with a custom number of traces on each
        """
        # Create a dialog box to get the plot format
        def getFormat():
            parent = QDialog()
            title = "Plot Format"
            label = 'Enter number of traces in each plot separated by ", ":'
            dialog = QInputDialog.getText(parent, title, label)
            if dialog[1]:
                form = dialog[0].split(',')
                form = [int(i) for i in form]
                return form
            else:
                raise ValueError("No format entered")

        pltform = getFormat()
        nplots = len(pltform)
        ntraces = len(args)

        # Check validity of the plot format
        if ntraces != np.sum(pltform):
            raise ValueError("Number of traces does not equal sum of traces to be plotted in format")

        self._clean_axes()
        # Split data items up into correct plots according to plot format
        traces = []
        maxindx = 0
        for i in range(nplots):
            plotdat = []
            minindx = maxindx
            maxindx += pltform[i]
            for j in range(ntraces):
                if j >= minindx and j < maxindx:
                    plotdat.append(args[j][0])
            traces.append(plotdat)

        # Set axes
        for plotnum in range(nplots):
            if plotnum == 0:
                self.axes = self.figure.add_subplot(nplots, 1, plotnum+1)
                ax = self.axes
            else:
                self.axes = self.figure.add_subplot(nplots, 1, plotnum+1, sharex=ax)
                setp(self.axes.get_xticklabels(), visible=False)

            # Create OPlots for each subfigure
            try:
                for data in traces[plotnum]:
                    label = data.desc
                    if label == "":
                        label = data.name + " (" + data.units + ")"
                    label += " " + data.source

                    time = data.time
                    if time is None:
                        if len(data.dim) != 1:
                            raise ValueError("Cannot plot '"+data.name+"' as it has too many dimensions")
                        time = data.dim[0].data
                    if hasattr(time, "data"):
                        time = time.data

                    self.axes.plot(time, data.data, label=label)
                    loc = "upper right" if len(time) > 10000 else "best"

                # Y label from last plot
                ylabel = data.desc
                if ylabel == "":
                    ylabel = data.label
                    if data.units != "":
                        ylabel += " ("+data.units+")"
                self.axes.set_ylabel(ylabel)

                if plotnum == 0:
                    self.axes.set_xlabel(traces[0][0].dim[traces[0][0].order].label)

                self.axes.legend(loc=loc)

            except TypeError:
                # Traces[plotnum] not iterable
                time = traces[plotnum].time.data
                xlabel = traces[plotnum].dim[traces[plotnum.order]].label
                if time is None:
                    if len(traces[plotnum].dim) != 1:
                        raise ValueError("Cannot plot '"+traces[plotnum].label+"' as it has too many dimensions")
                    time = traces[plotnum].dim[0].data
                    xlabel = traces[plotnum].dim[0].label

                data = traces[plotnum].data
                if size > 10000:
                    fac = int(len(time) / 10000)
                    time = time[::fac]
                    data = data[::fac]
                    print("Warning: too many samples (%d). Down-sampling to %d points" % (size, len(time)))

                self.axes.plot(time, data)

                ylabel = traces[plotnum].desc
                if ylabel == "":
                    ylabel = traces[plotnum].label
                    if traces[plotnum].units != "":
                        ylabel += " ("+traces[plotnum].units+")"
                self.axes.set_ylabel(ylabel)
                if plotnum == 0:
                    self.axes.set_xlabel(xlabel)

        self.figure.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.07, hspace=0.001)
        self.canvas.draw()

    def zplot(self, *args):
        """
        Plot traces and a zoomed in section on subplots
        """

        def connect_bbox(bbox1, bbox2,
                         loc1a, loc2a, loc1b, loc2b,
                         prop_lines, prop_patches=None):
            """
            Connect two bounding boxes
            """
            if prop_patches is None:
                prop_patches = prop_lines.copy()
                prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

            c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
            c1.set_clip_on(False)
            c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
            c2.set_clip_on(False)

            bbox_patch1 = BboxPatch(bbox1, **prop_patches)
            bbox_patch2 = BboxPatch(bbox2, **prop_patches)

            patch = BboxConnectorPatch(bbox1, bbox2,
                                       loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                                       **prop_patches)
            patch.set_clip_on(False)

            return c1, c2, bbox_patch1, bbox_patch2, patch

        def zoom_effect(ax1, ax2, **kwargs):
            """
            ax1 : the main axes
            ax1 : the zoomed axes

            Connect ax1 and ax2.  The xmin & xmax will be taken from the
            ax1.viewLim.
            """

            tt = ax1.transScale + (ax1.transLimits + ax2.transAxes)
            trans = blended_transform_factory(ax2.transData, tt)

            mybbox1 = ax1.bbox
            mybbox2 = TransformedBbox(ax1.viewLim, trans)

            prop_patches = kwargs.copy()
            prop_patches["ec"] = "none"
            prop_patches["alpha"] = 0.2

            c1, c2, bbox_patch1, bbox_patch2, patch = \
                connect_bbox(mybbox1, mybbox2,
                             loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                             prop_lines=kwargs, prop_patches=prop_patches)

            ax1.add_patch(bbox_patch1)
            ax2.add_patch(bbox_patch2)
            ax2.add_patch(c1)
            ax2.add_patch(c2)
            ax2.add_patch(patch)

            return c1, c2, bbox_patch1, bbox_patch2, patch

        # if len(args) == 1:
        #     self.plot(*args)

        self._clean_axes()

        self.ax_top = self.figure.add_subplot(2, 1, 1)
        self.ax_bottom = self.figure.add_subplot(2, 1, 2)
        self.axes = [self.ax_top, self.ax_bottom]

        zoom_effect(self.ax_top, self.ax_bottom)
        self.start_cid = self.ax_bottom.figure.canvas.mpl_connect("button_press_event",
                                                                  self.start_zoom_region)
        self.end_cid = self.ax_bottom.figure.canvas.mpl_connect("button_release_event",
                                                                self.end_zoom_region)
        self.callback_id = [self.start_cid, self.end_cid]

        for tracenum, trace in enumerate(args):
            try:
                for data in trace:
                    label = data.desc
                    if label == "":
                        label = data.name + " (" + data.units + ") " + data.source

                    time = data.time
                    if time is None:
                        if len(data.dim) != 1:
                            print(data.dim)
                            raise ValueError("Cannot plot '"+data.label+"' as it has too many dimensions")
                        time = data.dim[0].data
                    if hasattr(time, "data"):
                        time = time.data

                    loc = "upper right" if len(time) > 10000 else "best"

                    self.ax_top.plot(time, data.data, label=label)
                    self.ax_bottom.plot(time, data.data, label=label)
                if tracenum == 0:
                    x_dim = trace[0].dim[trace[0].order]
                    xlabel = x_dim.label or x_dim.name
                    if x_dim.units:
                        xlabel = "{} ({})".format(xlabel, x_dim.units)

                    self.ax_bottom.set_xlabel(xlabel)

                    ylabel = data.desc
                    if ylabel == "":
                        ylabel = data.label
                        if data.units != "":
                            ylabel += " ("+data.units+") "
                    self.ax_bottom.set_ylabel(ylabel)
                    self.ax_top.set_ylabel(ylabel)
                self.ax_bottom.legend(loc=loc)
            except TypeError:
                # Trace not iterable
                time = trace.time.data
                xlabel = trace.dim[trace.order].label
                if time is None:
                    if len(trace.dim) != 1:
                        raise ValueError("Cannot plot '"+trace.label+"' as it has too many dimensions")
                    time = trace.dim[0].data
                    xlabel = trace.dim[0].label
                self.ax_bottom.set_xlabel(xlabel)

                data = trace.data
                # Check array size
                size = len(time)
                if size > 10000:
                    fac = int(len(time) / 10000)
                    time = time[::fac]
                    data = data[::fac]
                    print("Warning: too many samples (%d). Down-sampling to %d points" % (size, len(time)))

                self.ax_top.plot(time, data)
                self.ax_bottom.plot(time, data)

                ylabel = trace.desc
                if ylabel == "":
                    ylabel = trace.label
                    if trace.units != "":
                        ylabel += " ("+trace.units+") "
                self.ax_top.set_ylabel(ylabel)
                self.ax_bottom.set_ylabel(ylabel)

        self.figure.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.1)
        self.canvas.draw()

    def start_zoom_region(self, event):
        """
        Callback for starting the zoomed region
        """
        if event.inaxes != self.ax_bottom:
            return
        _, right = self.ax_top.get_xlim()
        self.ax_top.set_xlim(event.xdata, right)
        self.ax_bottom.figure.canvas.draw()
        self.drag_cid = self.ax_bottom.figure.canvas.mpl_connect("motion_notify_event",
                                                                 self.middle_zoom_region)

    def middle_zoom_region(self, event):
        """
        Callback for dragging the zoomed region
        """
        if event.inaxes != self.ax_bottom:
            return
        left, _ = self.ax_top.get_xlim()
        self.ax_top.set_xlim(left, event.xdata)
        self.ax_bottom.figure.canvas.draw()

    def end_zoom_region(self, event):
        """
        Callback for ending the zoomed region
        """
        try:
            self.figure.canvas.mpl_disconnect(self.drag_cid)
        except AttributeError:
            pass

    def contour(self, item, levels=10):
        """
        Contour plot of item

        item - Item to plot
        levels - Number of contours to plot [10]
        """

        self._clean_axes()
        self.axes = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.08)

        if not (len(item.dim) == 2 or
                len(item.data) == 2 or
                len(item.data.shape) == 2):
            raise ValueError("Data must be two dimensional")
            return

        xdim = item.dim[item.order]
        ydim = item.dim[item.order+1]

        xaxis = xdim.data
        yaxis = ydim.data
        zaxis = item.data
        self.axes.contour(xaxis, yaxis, zaxis, levels)

        xlabel = xdim.label
        if xlabel == "":
            xlabel = xdim.name
            if xdim.units != "":
                xlabel += " ({})".format(xdim.units)

        self.axes.set_xlabel(xlabel)

        ylabel = ydim.label
        if ylabel == "":
            ylabel = ydim.name
            if ydim.units != "":
                ylabel += " ({})".format(ydim.units)

        self.axes.set_ylabel(ylabel)

        title = item.name
        if title == "":
            title = item.title
        self.axes.set_title(title)

        self.canvas.draw()

    def contourf(self, item, levels=10):
        """
        Filled contour plot of item

        item - Item to plot
        levels - Number of contours to plot [10]
        """

        self._clean_axes()
        self.axes = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.08)

        if not (len(item.dim) == 2 or
                len(item.data) == 2 or
                len(item.data.shape) == 2):
            raise ValueError("Data must be two dimensional")
            return

        xdim = item.dim[item.order]
        ydim = item.dim[item.order+1]

        xaxis = xdim.data
        yaxis = ydim.data
        zaxis = item.data
        self.axes.contourf(xaxis, yaxis, zaxis, levels)

        xlabel = xdim.label
        if xlabel == "":
            xlabel = xdim.name
            if xdim.units != "":
                xlabel += " ({})".format(xdim.units)

        self.axes.set_xlabel(xlabel)

        ylabel = ydim.label
        if ylabel == "":
            ylabel = ydim.name
            if ydim.units != "":
                ylabel += " ({})".format(ydim.units)

        self.axes.set_ylabel(ylabel)

        title = item.name
        if title == "":
            title = item.title
        self.axes.set_title(title)

        self.canvas.draw()

    def clearFig(self):
        """
        Reset the plot widget
        """
        self._clean_axes()
        self.axes = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.07)
        self.canvas.draw()
