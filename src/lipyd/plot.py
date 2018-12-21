#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
# The spectrum plotting based on code by Martin Noehrer
# See original at https://github.com/matrixx567/MassSpectraPlot
#

from future.utils import iteritems
from past.builtins import xrange, range

import imp
import copy
import pandas as pd
import numpy as np
import math as mt

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
from matplotlib import ticker
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class PlotBase(object):
    
    def __init__(
            self,
            fname = None,
            xlab = None,
            ylab = None,
            width = 7,
            height = 7,
            grid_rows = 1,
            grid_cols = 1,
            grid_hratios = None,
            grid_wratios = None,
            font_family = 'Helvetica Neue LT Std',
            font_style = 'normal',
            font_weight = 'normal',
            font_variant = 'normal',
            font_stretch = 'normal',
            font_size = None,
            axis_lab_font = {},
            bar_args = {},
            ticklabel_font = {},
            xticks = None,
            xticklabels = None,
            xticklabel_font = {},
            legend_font = {},
            title_font = {},
            uniform_ylim = False,
            palette = None,
            context = 'poster',
            lab_size = (9, 9),
            axis_lab_size = 10.0,
            lab_angle = 0,
            rc = {},
            title = None,
            maketitle = False,
            fp_title = None,
            title_halign = None,
        ):
        
        for k, v in iteritems(locals()):
            
            if not hasattr(self, k) or getattr(self, k) is None:
                
                setattr(self, k, v)
        
        if type(self.lab_size) is not tuple:
            self.lab_size = (self.lab_size, ) * 2
        if 'axes.labelsize' not in self.rc:
            self.rc['axes.labelsize'] = self.axis_lab_size
        if 'ytick.labelsize' not in self.rc:
            self.rc['ytick.labelsize'] = self.lab_size[0]
        if 'ytick.labelsize' not in self.rc:
            self.rc['ytick.labelsize'] = self.lab_size[1]
        
        self.rc['font.family'] = self.font_family
        self.rc['font.style'] = self.font_style
        self.rc['font.variant'] = self.font_variant
        self.rc['font.weight'] = self.font_weight
        self.rc['font.stretch'] = self.font_stretch
        
        self.fp = mpl.font_manager.FontProperties(
            family=self.font_family,
            style=self.font_style,
            variant=self.font_variant,
            weight=self.font_weight,
            stretch=self.font_stretch
        )
        
        self.grid_hratios = grid_hratios or [1.] * grid_rows
        self.grid_wratios = grid_wratios or [1.] * grid_cols
        
    
    def reload(self):
        """
        Reloads the module and updates the class instance.
        """
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        """
        The total workflow of this class.
        Calls all methods in the correct order.
        """
        self.pre_plot()
        self.do_plot()
        self.post_plot()

    def pre_plot(self):
        """
        Executes all necessary tasks before plotting in the correct order.
        """
        pass

    def do_plot(self):
        """
        Calls the plotting methods in the correct order.
        """
        self.set_figsize()
        self.init_fig()
        self.set_fontproperties()
        self.set_grid()
        self.make_plots()
        self.labels()
        self.set_ylims()
        self.set_title()
        self.set_ticklabels()

    def post_plot(self):
        """
        Saves the plot into file, and closes the figure.
        """
        self.finish()
    
    def set_figsize(self):
        """
        Converts width and height to a tuple so can be used for figsize.
        """
        if hasattr(self, 'width') and hasattr(self, 'height'):
            
            self.figsize = (self.width, self.height)

    def init_fig(self):
        """
        Creates a figure using the object oriented matplotlib interface.
        """
        self.pdf = mpl.backends.backend_pdf.PdfPages(self.fname)
        self.fig = mpl.figure.Figure(figsize = self.figsize)
        self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)

    def set_grid(self):
        """
        Sets up a grid according to the number of subplots,
        with proportions according to the number of elements
        in each subplot.
        """
        self.gs = mpl.gridspec.GridSpec(
            self.grid_rows,
            self.grid_cols,
            height_ratios = self.grid_hratios,
            width_ratios = self.grid_wratios,
        )
        self.axes = [[None] * self.grid_cols] * self.grid_rows

    def get_subplot(self, i, j = 0):
        
        if self.axes[j][i] is None:
            self.axes[j][i] = self.fig.add_subplot(self.gs[j, i])
        self.ax = self.axes[j][i]
    
    def set_ylims(self):
        
        if self.uniform_ylim:
            
            maxy = max(
                ax.get_ylim()[1]
                for ax in self.axes[0][0:]
            )
            
            _ = [None for _ in ax.set_ylim([0, maxy]) for ax in self.axes[0]]

    def set_title(self):
        """
        Sets the main title.
        """
        if self.maketitle:
            self.title_text = self.fig.suptitle(self.title)
            self.title_text.set_fontproperties(self.fp_title)
            self.title_text.set_horizontalalignment(self.title_halign)
    
    def labels(self):
        """
        Sets properties of axis labels and ticklabels.
        """
        if self.xlab is not None:
            self._xlab = self.ax.set_xlabel(self.xlab)
        
        if self.ylab is not None:
            self._ylab = self.ax.set_ylabel(self.ylab)
        
        _ = [
            tick.label.set_fontproperties(self.fp_xticklabel) or (
                self.lab_angle == 0 or self.lab_angle == 90
            ) and (
                tick.label.set_rotation(self.lab_angle) or
                tick.label.set_horizontalalignment('center')
            )
            for tick in self.ax.xaxis.get_major_ticks()
        ]
        
        _ = [
            tick.label.set_fontproperties(self.fp_ticklabel)
            for tick in self.ax.yaxis.get_major_ticks()
        ]
        
        self.ax.set_ylabel(self.ylab, fontproperties = self.fp_axis_lab)
        # self.ax.yaxis.label.set_fontproperties(self)
    
    def set_fontproperties(self):
        
        self.fp_default = (
            mpl.font_manager.FontProperties(
                family = self.font_family,
                style = self.font_style,
                weight = self.font_weight,
                variant = self.font_variant,
                stretch = self.font_stretch,
                size = self.font_size,
            )
        )
        
        self.fp_axis_lab = (
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.axis_lab_font)
            )
            if self.axis_lab_font else
            self.fp_default
        )
        
        self.fp_xticklabel = (
            mpl.font_manager.FontProperties(
                **copy.deepcopy(
                    self.xticklabel_font
                        if self.xticklabel_font else
                    self.ticklabel_font
                )
            )
            if self.xticklabel_font or self.ticklabel_font else
            self.fp_default
        )
        
        self.fp_ticklabel = (
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.ticklabel_font)
            ) if self.ticklabel_font else
            self.fp_default
        )
    
    def set_ticklabels(self):
        
        if self.xticklabels:
            
            self.ax.xaxis.set_ticklabels(
                ('{:,g}'.format(x) for x in self.xticklabels),
            )
        
        _ = [
            tl.set_fontproperties(self.fp_ticklabel)
            for tl in self.ax.get_xticklabels()
        ]
        
        print(self.ax.get_xticklabels()[0].get_fontproperties().get_family())
        
        if self.xticks:
            
            self.ax.xaxis.xticks(self.xticks)
    
    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """
        print(self.ax.get_xticklabels()[0].get_fontproperties().get_family())
        self.fig.tight_layout()
        self.fig.subplots_adjust(top = .92)
        self.cvs.draw()
        print(self.ax.get_xticklabels()[0].get_fontproperties().get_family())
        self.cvs.print_figure(self.pdf)
        print(self.ax.get_xticklabels()[0].get_fontproperties().get_family())
        self.pdf.close()
        print(self.ax.get_xticklabels()[0].get_fontproperties().get_family())
        #self.fig.clf()
        print(self.ax.get_xticklabels()[0].get_fontproperties().get_family())


class TestPlot(PlotBase):
    
    def __init__(self, **kwargs):
        
        PlotBase.__init__(self, **kwargs)
    
    def make_plots(self):
        
        self.get_subplot(0, 0)
        
        x = np.linspace(0, 10, 1000)
        self.ax.plot(x, np.sin(x))


class Profiles(PlotBase):
    
    def __init__(self, **kwargs):
        
        PlotBase.__init__(self, **kwargs)
    
    def make_plots(self, profiles, xlabels):
        
        self.get_subplot(0, 0)


class SpectrumPlotOld(object):
    """The summary line for a class docstring should fit on one line.


    Attributes
    ----------
    mzs : str
        Description of `attr1`.
    intensities :
        Description of `attr2`.
    annotations:

    xmin:

    xmax:

    title:


    """


    filetype = "pdf" # "png" or "pdf"
    """str: Filetype of the outputfile

    Tested with "pdf" and "png"
    
    """

    figdpi = 200
    """int: DPI of the image output file
    
    """


    # correct DPI of figure if using pdf
    if filetype.lower() == "pdf":
        figdpi = 72

    params = {'backend': 'pdf',
            'figure.dpi': figdpi,  # 72 for pdf
            'axes.labelsize': 10,
            'font.size': 10,
            'legend.fontsize': 8,
            'legend.frameon': True,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
            'font.family': 'Helvetica',
            'text.usetex': True,
            'text.latex.unicode': True,
            'axes.linewidth': 0.5,
            'xtick.major.size': 4,  # major tick size in points
            'xtick.minor.size': 2,  # minor tick size in points
            'xtick.direction': 'out',
            'ytick.major.size': 4,  # major tick size in points
            'ytick.minor.size': 2,  # minor tick size in points
            'ytick.direction': 'out',
            }
    
    plt.rcParams.update(params)
    
    def __init__(
            self,
            mzs,
            intensities,
            annotations = None,
            xmin = None,
            xmax = None,
            title = None,
        ):
        
        self.mzs = mzs
        self.intensities = intensities
        self.annotations = annotations
        self.xmin = xmin or int(min(self.mzs) / 10.0) * 10
        self.xmax = xmax or int(max(self.mzs) / 10.0 + 1) * 10
        self.title = title
    
    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def normalize_intensities(self):
        """ """
        
        self.inorm = self.intensities /  np.max(self.intensities) * 100.0
    
    @staticmethod
    def __get_max_peak(x_value, raw_values):
        """Search for the local peak next to the x-value.

        The value searches in the near area of the given x-value for the maximum (local)
        peak, which correspondends to the data.

        Args:
        -----
            x_value (float): m/z value (x-coordinate) of the point, to be annotated
            raw_values (Pandas.DataFrame): Data values, which should be annotated

        Returns:
        --------
            y value of the peak, next to the given x value

        """


        raw_values_index = raw_values[raw_values["m/z"] == x_value].index[0]

        value, index = float(raw_values.loc[raw_values_index - 5, "intensity_normalized"]), raw_values_index - 5

        for z in range(-5, 15):
            if float(raw_values.loc[raw_values_index + z, "intensity_normalized"]) > value:
                value, index = float(raw_values.loc[raw_values_index + z, "intensity_normalized"]), raw_values_index + z
        return value


    def annotate_point(x1, y_pos=0, text='', raw_values=None):
        """Annotate a specific point.
        
        Annotate a point with a label. The label will be placed vertically
        with an additional line.
        
        The function uses the data values and searches for the peak of the
        value to be annotated. Therefore the ``raw_values`` parameter is used

        Parameters
        ----------
        x1 : float
            m/z value (x-coordinate) of the point, to be annotated
        y_pos : `float`, optional
            Position of the label beginning (y-coordinate).
            The value uses the scale of the datapoints.
        text : str
            Label text
        raw_values : Pandas.DataFrame
            Data values, which should be annotated

        """

        delta_point_annotate_line_pixel = 1

        if raw_values is not None:

            if __get_max_peak(x1, raw_values) < (y_pos - delta_point_annotate_line_pixel):
                xy_text = (x1, y_pos)
            else:
                xy_text = (x1, __get_max_peak(x1, raw_values) + delta_point_annotate_line_pixel)

            plt.annotate(
                text, xy=(x1, __get_max_peak(x1, raw_values) + delta_point_annotate_line_pixel), xycoords='data',
                xytext=xy_text, textcoords='data',
                rotation=90, size=8, horizontalalignment='center', verticalalignment='bottom',
                arrowprops=dict(arrowstyle='-', color="#808080", linewidth=0.4, shrinkA=0.05, shrinkB=1))
    
    def annotate_distance(x1=0, x2=0, y_pos=0, text='', raw_values=None, rotate_text=0):
        """Annotate the distance between two peaks
        
        Annotate the distance between two given peaks. The text can be placed with a
        given angle.
        
        Args
        ----
        x1 : float
            m/z value (x-coordinate) of the left point
        x2 : float
            m/z value (x-coordinate) of the right point
        y_pos : float, optional
            Position of the label beginning (y-coordinate).
            The value uses the scale of the datapoints.
        text : str
            Label text
        raw_values : Pandas.DataFrame
            Data values, which should be annotated
        rotate_text : int, optional
            Rotation of the label, should be 0 or 90

        """
        
        delta_point_annotate_line_pixel = 1

        if raw_values is not None:
            if __get_max_peak(x1, raw_values) < (y_pos - delta_point_annotate_line_pixel):
                plt.annotate(
                    '', xy=(x1, __get_max_peak(x1, raw_values) + delta_point_annotate_line_pixel), xycoords='data',
                    xytext=(x1, y_pos), textcoords='data',
                    arrowprops=dict(arrowstyle='-', color="#808080", linewidth=0.4, shrinkA=0.05, shrinkB=0.05))
            
            if __get_max_peak(x2, raw_values) < (y_pos - delta_point_annotate_line_pixel):
                plt.annotate(
                    '', xy=(x2, __get_max_peak(x2, raw_values) + delta_point_annotate_line_pixel), xycoords='data',
                    xytext=(x2, y_pos), textcoords='data',
                    arrowprops=dict(arrowstyle='-', color="#808080", linewidth=0.4, shrinkA=0.05, shrinkB=0.05))
        
        plt.annotate(
            '', xy=(x1, y_pos), xycoords='data',
            xytext=(x2, y_pos), textcoords='data',
            arrowprops=dict(arrowstyle='<|-|>,head_length=0.4,head_width=0.1',
                            color="black", linewidth=0.6, shrinkA=0.05, shrinkB=0.05))
        plt.annotate(
            text, xy=((x1 if x1 <= x2 else x2) + mt.fabs((x1 - x2)) / 2, y_pos), xycoords='data',
            rotation=rotate_text, size=8,
            horizontalalignment='center', verticalalignment='bottom',
            xytext=(0, 2), textcoords='offset points')


    # Helper functions
    @staticmethod
    def __figsize_and_margins(
            plotsize,
            subplots = (1, 1),
            **absolute_margins,
        ):
        """Determine figure size and margins from plot size and absolute margins
        
        Args:
        -----
            plotsize: float
                (width, height) of plot area in inch
            subplots: int
                (nrows, ncols) of subplots
            left, right, top, bottom:
                absolute margins around plot area
            wspace, hspace:
                width and height spacing between subplots
        
        Returns:
        --------
            size
                figure size for figsize argument of figure()
            margins
                relative margins dict suitable for subplots_adjust()
        
        Example:
        --------
            making 2x2 grid of 3" square plots with specific spacings:
        
        sz, rm = figsize_and_margins((3,3), (2,2), left=1, right=.5,
                                                    top=.5, bottom=1,
                                                    wspace=.5, hspace=.5)
        figure(figsize=sz)
        subplots_adjust(**rm)
        subplot(221); subplot(222)
        subplot(223); subplot(224)
        

        """
        
        pw, ph = plotsize
        nr, nc = subplots
        amarg = absolute_margins
        # dictionary for relative margins
        # initialize from rcParams with margins not in amarg
        rmarg = dict((m, rcParams['figure.subplot.' + m])
                    for m in ('left', 'right', 'top', 'bottom', 'wspace', 'hspace')
                    if m not in amarg
                    )
        # subplots_adjust wants wspace and hspace relative to plotsize:
        if 'wspace' in amarg:
            rmarg['wspace'] = float(amarg['wspace']) / pw
        if 'hspace' in amarg:
            rmarg['hspace'] = float(amarg['hspace']) / ph
        # in terms of the relative margins:
        # width  * (right - left)
        #    = ncols * plot_width  + (ncols - 1) * wspace * plot_width
        # height * (top - bottom)
        #    = nrows * plot_height + (nrows - 1) * hspace * plot_height
        # solve for width and height, using absolute margins as necessary:
        width = float((nc + (nc - 1) * rmarg['wspace']) * pw + amarg.get('left', 0) + amarg.get('right', 0)) / (
            rmarg.get('right', 1) - rmarg.get('left', 0))
        height = float((nr + (nr - 1) * rmarg['hspace']) * ph + amarg.get('top', 0) + amarg.get('bottom', 0)) / (
            rmarg.get('top', 1) - rmarg.get('bottom', 0))

        # now we can get any remaining relative margins
        if 'left' in amarg:
            rmarg['left'] = float(amarg['left']) / width
        if 'right' in amarg:
            rmarg['right'] = 1 - float(amarg['right']) / width
        if 'top' in amarg:
            rmarg['top'] = 1 - float(amarg['top']) / height
        if 'bottom' in amarg:
            rmarg['bottom'] = float(amarg['bottom']) / height
        # return figure size and relative margins
        return (width, height), rmarg
    
    @staticmethod
    def __get_ax_size(ax, fig):
        """
        Source:
            http://scipy-central.org/item/65/1/
            absolute-plot-size-and-margins-in-matplotlib
        """
        bbox = ax.get_window_extent().transformed(
            fig.dpi_scale_trans.inverted()
        )
        width, height = bbox.width, bbox.height
        width *= fig.dpi
        height *= fig.dpi
        return width, height
    
    @staticmethod
    def __conv_inch(length_mm):
        """Converts a length from millimeters to inch
        
        Args:
        -----
        :param int,float length_mm :
            Length in millimeters to be converted
        
        Returns:
        --------
        float: The converted length in inch
        """
        return length_mm / 25.4
    
    def generate_massspectra_plot_automatic_labels(input_filename, output_filename):
        """

        Parameters
        ----------
        input_filename :
            
        output_filename :
            

        Returns
        -------

        """
        print("Generate {0} mass spectra plot ({1}.{2}) from {3}.[xy/txt]".format(filetype.upper(), output_filename,
                                                                                filetype, input_filename))

        # min and max value of the x axis
        min_xaxis = 50
        max_xaxis = 425

        # dimensions of the single plot
        figwidth_mm = 141
        left_margin_mm = 14
        right_margin_mm = 2
        top_margin_mm = 26
        bottom_margin_mm = 13

        plotwidth_mm = figwidth_mm - left_margin_mm - right_margin_mm
        plotheight_mm = plotwidth_mm * 0.8

        hspace_mm = 0.1  # horizontal space between two plots

        # read xy input data
        raw_values = pd.read_csv(input_filename + ".xy", sep=" ", names=["m/z", "intensity"],
                                dtype={"m/z": np.float64, "intensity": np.float64})
        # normalize input data
        max_val = np.max(raw_values['intensity'])
        raw_values["intensity_normalized"] = raw_values['intensity'] / max_val * 100.0  # normalize intensity to percent

        # read txt captions
        caption_values = pd.read_csv(input_filename + ".txt", sep=";", names=["m/z", "caption"], dtype={"m/z": np.float64})

        # find maximum value for each annotate peak - necessary for later calculation of annotation lines
        for i in caption_values.index:
            raw_values_index = raw_values[raw_values["m/z"] == caption_values.ix[i]["m/z"]].index[0]
            value, index = float(raw_values.loc[raw_values_index - 5, "intensity_normalized"]), raw_values_index - 5

            for z in range(-5, 15):
                if float(raw_values.loc[raw_values_index + z, "intensity_normalized"]) > value:
                    value, index = float(raw_values.loc[raw_values_index + z, "intensity_normalized"]), raw_values_index + z

            caption_values.loc[i, "m/z"] = float(raw_values.loc[index, "m/z"])
            caption_values.loc[i, "intensity_normalized"] = value  # add intensity to caption table

        # dimension for the annotations in pixel
        delta_x_text_labels_mm = 3.8
        delta_point_annotate_line_mm = 1
        label_annotate_arm_high_mm = 5.3
        delta_y_diagonal_annotate_arm_mm = 3.5
        delta_y_baseline_axistop_mm = 6.3
        delta_label_mznumber_pixel_mm = 7

        fsize, margins = __figsize_and_margins(plotsize=(__conv_inch(plotwidth_mm), __conv_inch(plotheight_mm)),
                                            subplots=(1, 1),
                                            left=__conv_inch(left_margin_mm), right=__conv_inch(right_margin_mm),
                                            top=__conv_inch(top_margin_mm), bottom=__conv_inch(bottom_margin_mm),
                                            wspace=0.01, hspace=__conv_inch(hspace_mm))
        fig = plt.figure(figsize=fsize, )
        plt.subplots_adjust(**margins)
        ax = fig.add_subplot(111)

        # delta_dimension in pixel calculate from mm
        delta_x_text_labels_pixel = __conv_inch(delta_x_text_labels_mm) * fig.dpi
        delta_point_annotate_line_pixel = __conv_inch(delta_point_annotate_line_mm) * fig.dpi
        label_annotate_arm_high_pixel = __conv_inch(label_annotate_arm_high_mm) * fig.dpi
        delta_y_diagonal_annotate_arm = __conv_inch(delta_y_diagonal_annotate_arm_mm) * fig.dpi
        delta_y_baseline_axistop = __conv_inch(delta_y_baseline_axistop_mm) * fig.dpi
        delta_label_mznumber_pixel = __conv_inch(delta_label_mznumber_pixel_mm) * fig.dpi

        # plot spectra line
        ax.plot(raw_values["m/z"], raw_values["intensity_normalized"], color="black", linewidth=0.8)

        # set x axes range
        if min_xaxis is not None and max_xaxis is not None:
            ax.set_xlim([int(min_xaxis), int(max_xaxis)])

        length_x_axis_pixel, length_y_axis_pixel = __get_ax_size(ax, fig)
        x_axis_min, x_axis_max = ax.get_xlim()
        y_axis_min, y_axis_max = ax.get_ylim()

        list_text_pos = list()
        if len(caption_values.index) > 0:

            # annotate from basepeak
            basepeak_index = caption_values[caption_values["intensity_normalized"] >= 100.0].index
            if len(basepeak_index) == 0:
                basepeak_index = [len(caption_values) / 2]

            # annotate peaks left frome basepeak
            for z in caption_values.loc[:basepeak_index[0], "m/z"].sort_index(ascending=False):
                x_value_pixel = (z - x_axis_min) / (x_axis_max - x_axis_min) * length_x_axis_pixel

                if len(list_text_pos) == 0:
                    list_text_pos.append(x_value_pixel)
                else:
                    if x_value_pixel > list_text_pos[-1] + delta_x_text_labels_pixel:
                        list_text_pos.append(x_value_pixel)
                    else:
                        list_text_pos.append(list_text_pos[-1] + delta_x_text_labels_pixel)

            # annotate peaks right from basepeak
            if basepeak_index[0] + 1 in caption_values.index:
                for z in caption_values.loc[basepeak_index[0] + 1:, "m/z"]:
                    x_value_pixel = (z - x_axis_min) / (x_axis_max - x_axis_min) * length_x_axis_pixel

                    if x_value_pixel < list_text_pos[0] - delta_x_text_labels_pixel:
                        list_text_pos.insert(0, x_value_pixel)
                    else:
                        list_text_pos.insert(0, list_text_pos[0] - delta_x_text_labels_pixel)

            # move annotation to left, if the label is out of the axes range
            if len(list_text_pos) > 0 and max(list_text_pos) >= length_x_axis_pixel:
                print("move labels to left")
                calc_delta = max(list_text_pos) - length_x_axis_pixel
                for index in range(len(list_text_pos) - 1, 0, -1):
                    print("change", list_text_pos[index], calc_delta)
                    if index > 1 and list_text_pos[index] - calc_delta > list_text_pos[index - 1] + delta_x_text_labels_pixel:
                        list_text_pos[index] = list_text_pos[index] - calc_delta
                        break
                    else:
                        list_text_pos[index] = list_text_pos[index] - calc_delta

            for i in caption_values.index:

                label_x_pos = list_text_pos.pop()
                x_value = float(caption_values.loc[i, "m/z"])
                y_value = float(caption_values.loc[i, "intensity_normalized"])
                label_text = str(caption_values.loc[i, "caption"])

                x_value_pixel = (x_value - x_axis_min) / (x_axis_max - x_axis_min) * length_x_axis_pixel
                y_value_pixel = (y_value - y_axis_min) / (y_axis_max - y_axis_min) * length_y_axis_pixel
                label_y_pos_pixel = length_y_axis_pixel + delta_y_baseline_axistop
                length_y_annotation_line = label_y_pos_pixel - y_value_pixel - delta_point_annotate_line_pixel
                label_annotate_arm_low_pixel = length_y_annotation_line - label_annotate_arm_high_pixel - delta_y_diagonal_annotate_arm

                if label_annotate_arm_low_pixel < 0:
                    label_annotate_arm_low_pixel = 0

                ax.annotate("{:.0f}".format(x_value),
                            xy=(x_value_pixel, y_value_pixel + delta_point_annotate_line_pixel),
                            xycoords='axes pixels', rotation=90,
                            xytext=(label_x_pos, label_y_pos_pixel),
                            textcoords='axes pixels', size=8,
                            arrowprops=dict(arrowstyle="-",
                                            connectionstyle="arc,angleA=-90,armA=" + str(label_annotate_arm_high_pixel) +
                                                            ",angleB= 90,armB=" + str(label_annotate_arm_low_pixel) +
                                                            ",rad=0", linewidth=0.5, color='#808080'),
                            horizontalalignment='center', verticalalignment='bottom')
                ax.annotate(label_text,
                            xy=(x_value_pixel, y_value_pixel + delta_point_annotate_line_pixel),
                            xycoords='axes pixels', rotation=90,
                            xytext=(label_x_pos, label_y_pos_pixel + delta_label_mznumber_pixel),
                            textcoords='axes pixels', size=8,
                            horizontalalignment='center', verticalalignment='bottom'
                            )

        # remove top and right axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # label axes
        ax.set_xlabel(r"$m/z$")
        ax.set_ylabel(r"$Intensity\,[\%]$")

        # set x labels
        plt.xticks(rotation='vertical')
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end + 1, 25))

        # set y labels
        ax.set_ylim(0, 100)
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
        
        # set grid
        plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)
        
        plt.savefig(output_filename + "." + filetype, dpi=fig.dpi, format=filetype)
        plt.close()
    
    def annotate(self):
        """ """
        
        if self.annotations is None:
            
            return
        
        ann_idx = np.where(np.array([bool(a) for a in self.annotations]))
        
        for i, annot in enumerate(self.annotations):
            
            pass
        
        self.annotate_point(403, 83, r'Peak II', raw_values)
        self.annotate_point(253, 83, r'Peak I', raw_values)
        
        annotate_distance(385, 403, 45, '$-18$', raw_values, rotate_text=90)
        annotate_distance(367, 385, 45, '$-18$', raw_values, rotate_text=90)
        annotate_distance(339, 367, 25, '$xy$', raw_values, rotate_text=90)
        annotate_distance(321, 349, 35, '$xy$', raw_values, rotate_text=90)
        annotate_distance(253, 349, 72, '$-96$', raw_values)
    
    def build(self):
        """ """
        
        # dimensions of the single plot
        figwidth_mm = 141
        left_margin_mm = 14
        right_margin_mm = 2
        top_margin_mm = 5
        bottom_margin_mm = 13
        
        plotwidth_mm = figwidth_mm - left_margin_mm - right_margin_mm
        plotheight_mm = plotwidth_mm * 0.8
        
        hspace_mm = 0.1  # horizontal space between two plots
        
        fsize, margins = self.__figsize_and_margins(
            plotsize = (
                self.__conv_inch(plotwidth_mm),
                self.__conv_inch(plotheight_mm)
            ),
            subplots = (1, 1),
            left   = self.__conv_inch(left_margin_mm),
            right  = self.__conv_inch(right_margin_mm),
            top    = self.__conv_inch(top_margin_mm),
            bottom = self.__conv_inch(bottom_margin_mm),
            wspace = 0.01,
            hspace = self.__conv_inch(hspace_mm)
        )
        self.fig = plt.figure(figsize = fsize, )
        plt.subplots_adjust(**margins)
        self.ax = self.fig.add_subplot(111)
        
        # plot spectra line
        self.ax.plot(
            self.mzs,
            self.inorm,
            color = "black",
            linewidth = 0.8,
            label = self.title,
        )
        
        # set x axes range
        self.ax.set_xlim([self.xmin, self.xmax])
        
        # annotations
        # length in data value (in %)
        self.annotate()
        
        # remove top and right axis
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.get_xaxis().tick_bottom()
        self.ax.get_yaxis().tick_left()
        
        # label axes
        self.ax.set_xlabel(r"$m/z$")
        self.ax.set_ylabel(r"$Intensity\,[\%]$")
        
        # set x labels
        plt.xticks(rotation = 'vertical')
        start, end = self.ax.get_xlim()
        self.ax.xaxis.set_ticks(np.arange(start, end + 1, 25))
        
        # set y labels
        self.ax.set_ylim(0, 100)
        start, end = self.ax.get_ylim()
        self.ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
        
        # set grid
        plt.grid(
            True,
            axis = 'y',
            color = 'black',
            linestyle = ':',
            linewidth = 0.1,
        )
        
        plt.savefig(self.outfile, dpi = self.fig.dpi, format = self.filetype)
        plt.close()
    
    def generate_massspectra_two_plot_manual_annotation(input_filename1, input_filename2, output_filename):
        """

        Parameters
        ----------
        input_filename1 :
            
        input_filename2 :
            
        output_filename :
            

        Returns
        -------

        """
        print("Generate {0} mass spectra plot ({1}.{2}) from {3}.xy and {4}.xy".format(filetype.upper(), output_filename,
                                                                                    filetype, input_filename1,
                                                                                    input_filename2))

        # min and max value of the x axis
        min_xaxis = 50
        max_xaxis = 475

        # labels
        label1 = "Substance 1"
        label2 = "Substance 2"

        # dimensions of the single plot
        figwidth_mm = 141
        left_margin_mm = 14
        right_margin_mm = 2
        top_margin_mm = 5
        bottom_margin_mm = 13

        plotwidth_mm = figwidth_mm - left_margin_mm - right_margin_mm
        plotheight_mm = plotwidth_mm * 0.8

        hspace_mm = 25  # horizontal space between two plots

        # read xy input data
        raw_values1 = pd.read_csv(input_filename1 + ".xy", sep=" ", names=["m/z", "intensity"],
                                dtype={"m/z": np.float64, "intensity": np.float64})
        # normalize input data
        max_val = np.max(raw_values1['intensity'])
        raw_values1["intensity_normalized"] = raw_values1['intensity'] / max_val * 100.0  # normalize intensity to percent

        # read xy input data
        raw_values2 = pd.read_csv(input_filename2 + ".xy", sep=" ", names=["m/z", "intensity"],
                                dtype={"m/z": np.float64, "intensity": np.float64})
        # normalize input data
        max_val = np.max(raw_values2['intensity'])
        raw_values2["intensity_normalized"] = raw_values2['intensity'] / max_val * 100.0  # normalize intensity to percent

        fsize, margins = __figsize_and_margins(plotsize=(__conv_inch(plotwidth_mm), __conv_inch(plotheight_mm)),
                                            subplots=(1, 1),
                                            left=__conv_inch(left_margin_mm), right=__conv_inch(right_margin_mm),
                                            top=__conv_inch(top_margin_mm), bottom=__conv_inch(bottom_margin_mm),
                                            wspace=0.01, hspace=__conv_inch(hspace_mm))
        fig = plt.figure(figsize=fsize, )
        plt.subplots_adjust(**margins)

        ax = fig.add_subplot(211)

        # plot spectra line
        plt1, = ax.plot(raw_values1["m/z"], raw_values1["intensity_normalized"], color="black", linewidth=0.8, label=label1)
        legend = plt.legend(handles=[plt1], loc=2)
        legend.get_frame().set_linewidth(0)

        # set x axes range
        if min_xaxis is not None and max_xaxis is not None:
            ax.set_xlim([int(min_xaxis), int(max_xaxis)])

        # annotate plot1
        annotate_point(253, 60, r'Peak 1', raw_values1)
        annotate_distance(385, 403, 45, r'$-18$', raw_values1, rotate_text=90)
        annotate_distance(367, 385, 45, r'$-18$', raw_values1, rotate_text=90)
        annotate_distance(321, 349, 75, r'$xy$', raw_values1, rotate_text=0)
        annotate_point(403, 60, r'Peak 2', raw_values1)

        # remove top and right axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # label axes
        ax.set_ylabel(r"$Intensity\,[\%]$")

        # set x labels
        plt.xticks(rotation='vertical')
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end + 1, 25))

        # set y labels
        ax.set_ylim(0, 100)
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(np.arange(start, end + 1, 10))

        # set grid
        plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)

        # set x labels
        plt.xticks(rotation='vertical')
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end + 1, 25))
        # set y labels
        ax.set_ylim(0, 100)
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(np.arange(start, end + 1, 10))
        # set grid
        plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)

        # generate plot2
        ax = fig.add_subplot(212)

        plt2, = ax.plot(raw_values2["m/z"], raw_values2["intensity_normalized"], color="#606060", linewidth=0.8,
                        label=label2)
        legend = plt.legend(handles=[plt2], loc=2)
        legend.get_frame().set_linewidth(0)

        # set x axes range
        if min_xaxis is not None and max_xaxis is not None:
            ax.set_xlim([int(min_xaxis), int(max_xaxis)])

        # annotate
        annotate_point(253, 60, r'Peak I', raw_values2)
        annotate_point(189, 60, r'189', raw_values2)
        annotate_point(215, 60, r'215', raw_values2)
        annotate_distance(385, 445, 30, r'Multiple'"\n"r'line'"\n"'annotation', raw_values2)
        annotate_distance(367, 385, 45, r' -18', raw_values2, rotate_text=90)
        annotate_distance(349, 367, 45, r' -18', raw_values2, rotate_text=90)
        annotate_distance(321, 349, 75, r'xy', raw_values2, rotate_text=0)
        annotate_point(445, 60, r'Peak 2', raw_values2)

        # remove top and right axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # label axes
        ax.set_xlabel(r"$m/z$")
        ax.set_ylabel(r"$Intensity\,[\%]$")

        # set x labels
        plt.xticks(rotation='vertical')
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end + 1, 25))

        # set y labels
        ax.set_ylim(0, 100)
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(np.arange(start, end + 1, 10))

        # set grid
        plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.1)

        plt.savefig(output_filename + "." + filetype, dpi=fig.dpi, format=filetype)
        plt.close()


class SpectrumPlot(object):
    """ Class for drawing a plot
        
        Arguments:
        ----------
        title : str
            title of a plot
        result_type : str 
            type of result plot - screen or pdf
        pdf_file_name : str
            filename of your downloading plot (pdf)
        mzs : float
             list of mzs values
        intensities : float
            list of intensity values
        annotations : str
            annotations of peaks
        """

    epsilon = 1e-06
    pdf_type = "pdf"
    screen_type = "screen"

    def __init__(self,
                title = "",
                result_type = "screen",
                pdf_file_name = ". / result.pdf"
                ):
        self.title = title
        self.result_type = result_type
        self.pdf_file_name = pdf_file_name
        self.mzs = []
        self.intensities = []
        self.annotations = []

        
    def build_plot(self,
                    mzs = [],           # the list contains points X
                    intensities = [],           # the list contains points Y
                    annotations = [],        # annotX["X"] - contains annotated points X; annotX["caption"]contains caption of point;
                    title = "",
                    result_type = "screen",
                    pdf_file_name = ". / result.pdf"):
        """ Main method for drawing a plot
        
        Arguments:
        ----------
        li : list
            list of real number
        real_x : float 
            real number
        
        Return:
        -------
        index or None
        
        """

        self.mzs = mzs
        self.intensities = intensities
        max_x = np.max(self.mzs)
        max_y = np.max(self.intensities)
        self.intensities = self.intensities / max_y * 100 #normalize;
        
        self.annotations = annotations          # [(x1,a1),(x2,a2), ...]
        self.title = title
        self.result_type = result_type
        self.pdf_file_name = pdf_file_name


        #For saving plot in pdf format
        #pp = PdfPages(self.pdf_file_name)
        
        fig, ax = plt.subplots(1,1,figsize = (17,9))

        ax.plot(self.mzs, self.intensities, "k")
        ax.set_title(self.title)

        plt.xlabel("m / z")
        plt.ylabel("Relative intensity, %")

        trans = ax.get_xaxis_transform()
        

        vertex_number = 1
        annotate_base_x = 1.1 * max_x  #6% more of max to right after max of X;
        annotate_base_y = 0.07
        annotate_offset_x = 0.
        annotate_offset_y = .06
        annotate_x = annotate_base_x
        annotate_y = annotate_base_y

        for x, c in self.annotations:
            index_x2 = self.BinSearch4Real(self.mzs, x)
            if index_x2 is None:
                print("Value {} not found in annotations.".format(x) )
            else:
                y = self.intensities[index_x2]  # get y - value in point x2;
                ax.annotate("{}".format(vertex_number),
                    xy = (x,y),
                    xycoords = 'data',
                    xytext = (0, 0),
                    textcoords = 'offset points', #"data" 
                    ha = "center",
                    va = "center",
                    bbox = dict(boxstyle = "circle", fc = "w")
                    )

                ax.annotate("{} - {}".format(vertex_number ,c),
                    xy = (annotate_x, annotate_y),
                    xycoords = trans,
                    ha = "center",
                    va = "top",
                    bbox = dict(boxstyle = "round", fc = "w")
                    )

                vertex_number += 1
                annotate_x += annotate_offset_x
                annotate_y += annotate_offset_y
        
        if (self.result_type == Drawer.screen_type):
            plt.show()
        
        #for saving in pdf format
        #pp.savefig(fig)
        #pp.close()

    def BinSearch4Real(self, li, real_x):
        """ Method for annotating peaks using binary search algorithm
        
        Arguments:
        ----------
        li : list
            list of real number
        real_x : float 
            real number
        
        Return:
        -------
            None or value
        
        """

        i = 0
        j = len(li)-1
        m = int(j / 2)
        # look over list while li[m] - real_x greater than Drawer.epsilon
        while abs(li[m] - real_x) > Drawer.epsilon and i <= j:
            if real_x > li[m]:
                i = m + 1
            else:
                j = m-1
            m = int((i + j) / 2)

        if i > j:
            return None
        else:
            return m
        

if __name__ == "__main__":
    

    #Test data
    a=Drawer()
    x=[10., 11.1, 12.2, 13.3, 14.4, 15.5, 16.6, 17.7, 18.8, 19.9, 20.2]
    y=[100., 101., 150., 102., 104., 101., 180., 105., 102., 107., 103.]
    
    # Annotations done by use correlation between x coord and name of peak
    xa=[(12.2, "the first peack"),(16.6, "the second peack")]
    
    a.build_plot(mzs=x,
                intensities=y,
                annotations=xa,
                result_type="screen",
                title="test")
