#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2019 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://denes.omnipathdb.org/
#

#
# The spectrum plotting based on code by Martin Noehrer
# See original at https://github.com/matrixx567/MassSpectraPlot
#

from future.utils import iteritems
from past.builtins import xrange, range

import imp
import copy
import itertools
import pandas as pd
import numpy as np
import math as mt

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
from matplotlib import ticker
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages

import lipyd.settings as settings
import lipyd.session as session
import lipyd.common as common


def is_opentype_cff_font(filename):
    """
    This is necessary to fix a bug in matplotlib:
    https://github.com/matplotlib/matplotlib/pull/6714
    Returns True if the given font is a Postscript Compact Font Format
    Font embedded in an OpenType wrapper.  Used by the PostScript and
    PDF backends that can not subset these fonts.
    """
    
    if os.path.splitext(filename)[1].lower() == '.otf':
        
        result = _is_opentype_cff_font_cache.get(filename)
        
        if result is None:
            
            with open(filename, 'rb') as fd:
                
                tag = fd.read(4)
            
            result = (tag == b'OTTO')
            _is_opentype_cff_font_cache[filename] = result
        
        return result
    
    return False


mpl.font_manager.is_opentype_cff_font = is_opentype_cff_font


class PlotBase(object):
    
    text_elements = (
        'axis_label', 'ticklabel', 'xticklabel',
        'title', 'legend_title', 'legend_label',
        'annotation',
    )
    
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
            font_family = None,
            font_style = 'normal',
            font_weight = 'normal',
            font_variant = 'normal',
            font_stretch = 'normal',
            font_size = None,
            font_sizes = None,
            axis_label_font = None,
            ticklabel_font = None,
            xticklabel_font = None,
            legend_label_font = None,
            legend_title_font = None,
            title_font = None,
            annotation_font = None,
            legend_loc = 2,
            legend_title = None,
            legend = True,
            bar_args = None,
            xticks = None,
            xticklabels = None,
            uniform_ylim = False,
            palette = None,
            lab_size = (9, 9),
            axis_label_size = 10.0,
            lab_angle = 0,
            rc = None,
            title = None,
            maketitle = False,
            title_halign = None,
            usetex = False,
            do_plot = True,
            **kwargs,
        ):
        """
        Base class for various plotting classes. The ``plot`` method is the
        main method in this class, it executes the entire workflow of
        plotting. By default ``__init`` calls this method which means plot is
        made upon instantiation of the class. Plotting consists of three major
        phases done by methods ``pre_plot``, ``make_plot`` and ``post_plot``.
        In this base class ``pre_plot`` does nothing, ``make_plot`` creates
        a objects for figure, axes, grid, etc. It also does the plotting
        itself and applies the additional options like typeface, labels,
        sizes, etc. Then ``post_plot`` applies ``tight_layout`` and closes
        the file.
        
        Parameters
        ----------
        fname : str
            File name for the graphics, e.g. `boxplot01.pdf`.
        xlab : str
            Label for the x axis.
        ylab : str
            Label for the y axis.
        width : float
            Figure width in inches.
        height : float
            Figure height in inches.
        figsizes : tuple
            Figure size as tuple of 2 numbers, in inches.
        grid_rows : int
            Number of rows in the grid. We create a grid even for single plot,
            hence the default is 1.
        grid_cols : int
            Number of columns in the grid. Just like ``grid_rows``.
        grid_hratios : list
            Height ratios of grid columns. List of floats.
        grid_wratios : list
            Width ratios of grid rows. List of floats.
        font_family : str,list
            Font family to use or list of families in order of preference.
            Default is from ``settings.font_family``.
        font_style : str
            Font style, e.g. `bold` or `medium`.
        font_variant : str
            Font variant, e.g. `small-caps`.
        font_stretch : str
            Font stretch, e.g. `condensed` or `expanded`.
        font_size : float
            Base font size.
        font_sizes : dict
            Font size proportions. E.g. ``{'axis_label': 0.8}`` will
            result the axis label to have size of <base size> * 0.8.
        axis_label_font : dict
            Dict with the same font parameters as listed above, specific
            for axis labels.
        ticklabel_font : dict
            Dict with the same font parameters as listed above, specific
            for tick labels.
        xticklabel_font : dict
            Dict with the same font parameters as listed above, specific
            for tick labels of x axis.
        legend_title_font : dict
            Dict with the same font parameters as listed above, specific
            for the legend title.
        legend_label_font : dict
            Dict with the same font parameters as listed above, specific
            for the legend labels.
        title_font : dict
            Dict with the same font parameters as listed above, specific
            for the main title.
        annotation_font : dict
            Dict with the same font parameters as listed above, specific
            for the annotations.
        legend_loc : int
            The location of the legend.
        legend_title : str
            Title for the legend.
        legend : bool
            Create legend for the plot.
        bar_args : dict
            Arguments for barplots ``bar`` method.
        xticks : list
            Locations of ticks on x axis. List of floats.
        xticklabels : list
            Tick labels on x axis. List of strings.
        uniform_ylim : bool
            In case of multiple plots on a grid, the y axis limits should be
            uniform or different for each plot.
        palette : list
            Colours to use.
        lab_size : tuple
            Font size of the tick labels. Tuple of two numbers, for x and y
            axes, respectively.
        axis_label_size : float
            Font size of the axis labels.
        lab_angle : float
            Angle of the x axis labels.
        rc : dict
            Matplotlib rc params.
        title : str
            Main title of the plot.
        maketitle : bool
            Plot with or without main title.
        title_halign : str
            Horizontal alignement of the main title.
        usetex : bool
            Use LaTeX for rendering text elements.
        do_plot : bool
            Execute the plotting workflow upon instatiation.
            This is convenient by default but for resolving issues sometimes
            beneficial to first create the object and call individual methods
            afterwards.
        """
        
        for k, v in itertools.chain(iteritems(locals()), iteritems(kwargs)):
            
            # we check this because derived classes might have set
            # already attributes
            if not hasattr(self, k) or getattr(self, k) is None:
                
                setattr(self, k, v)
        
        if self.do_plot:
            
            self.main()
    
    def reload(self):
        """
        Reloads the module and updates the class instance.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        """
        The total workflow of this class.
        Calls all methods in the correct order.
        """
        
        self.pre_plot()
        self.make_plot()
        self.post_plot()
    
    # a synonym
    main = plot
    
    def pre_plot(self):
        """
        Executes all necessary tasks before plotting in the correct order.
        Derived classes should override this if necessary.
        """
        
        self.set_fonts()
        self.set_bar_args()
        self.set_figsize()
        self.init_fig()
        self.set_grid()

    def make_plot(self):
        """
        Calls the plotting methods in the correct order.
        """
        
        self.make_plots() # this should call post_subplot_hook
                          # after each subplot
    
    def post_subplot_hook(self):
        """
        A method executed each time after creating a subplot.
        """
        
        self.labels()
        self.set_ylims()
        self.set_title()
        self.set_ticklabels()
        self.make_legend()
        self.set_legend_font()

    def post_plot(self):
        """
        Saves the plot into file, and closes the figure.
        """
        
        self.finish()
    
    def set_bar_args(self):
        
        self.bar_args = self.bar_args or {}
    
    def set_fonts(self):
        """
        Sets up everything related to fonts.
        """
        
        # if no font parameters provided by arguments or derived classes
        # we get the defaults from settings
        self.fonts_defaults_from_settings()
        # we create a dict of default font parameters
        self.fonts_default_dict()
        # replace None defaults with dicts
        self.fonts_init_dicts()
        # set font parameters in rc (is this necessary at all?)
        self.fonts_set_rc()
        # create font properties objects from all parameters dicts
        self.fonts_create_fontproperties()
    
    def fonts_init_dicts(self):
        """
        Initializes specific font argument dicts unless they are explicitely
        provided.
        """
        
        font_sizes = copy.deepcopy(settings.get('font_sizes'))
        font_sizes.update(self.font_sizes or {})
        
        for text in self.text_elements:
            
            attr = '%s_font' % text
            this_font = copy.deepcopy(self.font_default)
            specific  = getattr(self, attr) or {}
            
            if 'size' not in specific:
                
                specific['size'] = self.font_size * (
                    font_sizes[text] if text in font_sizes else 1.0
                )
            
            this_font.update(specific)
            
            setattr(self, attr, this_font)
    
    def fonts_set_rc(self):
        """
        Sets up font related settings in matplotlib rc dict.
        """
        
        self.rc = self.rc or {}
        
        if type(self.lab_size) is not tuple:
            self.lab_size = (self.lab_size, ) * 2
        if 'axes.labelsize' not in self.rc:
            self.rc['axes.labelsize'] = self.axis_label_size
        if 'ytick.labelsize' not in self.rc:
            self.rc['ytick.labelsize'] = self.lab_size[0]
        if 'ytick.labelsize' not in self.rc:
            self.rc['ytick.labelsize'] = self.lab_size[1]
        
        self.rc['font.family'] = self.font_family
        self.rc['font.style'] = self.font_style
        self.rc['font.variant'] = self.font_variant
        self.rc['font.weight'] = self.font_weight
        self.rc['font.stretch'] = self.font_stretch
        self.rc['text.usetex'] = self.usetex
    
    def fonts_defaults_from_settings(self):
        """
        Sets default font options from ``settings`` unless they are
        explicitely set already.
        """
        
        self.font_family = self.font_family or settings.get('font_family')
        self.font_style = self.font_style or settings.get('font_style')
        self.font_variant = self.font_variant or settings.get('font_variant')
        self.font_stretch = self.font_stretch or settings.get('font_stretch')
        self.font_size = self.font_size or settings.get('font_size')
    
    def fonts_default_dict(self):
        """
        Creates a dict from default font parameters which can be passed later
        as arguments for ``matplotlib.font_manager.FontProperties``.
        """
        
        # dict for default font parameters
        self.font_default = {
            'family': self.font_family,
            'style': self.font_style,
            'variant': self.font_variant,
            'stretch': self.font_stretch,
            'size': self.font_size,
        }
    
    def fonts_create_fontproperties(self):
        """
        Creates ``matplotlib.font_manager.FontProperties`` objects from
        each font parameter dict.
        """
        
        self.fp = mpl.font_manager.FontProperties(
            **copy.deepcopy(self.font_default)
        )
        
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
        
        for text in self.text_elements:
            
            dictattr = '%s_font' % text
            fpattr   = 'fp_%s' % text
            
            setattr(
                self,
                fpattr,
                mpl.font_manager.FontProperties(
                    **copy.deepcopy(getattr(self, dictattr))
                )
            )
    
    def set_figsize(self):
        """
        Converts width and height to a tuple so can be used for figsize.
        """
        
        if hasattr(self, 'figsize') and self.figsize is not None:
            
            return
            
        elif self.width and self.height:
            
            self.figsize = (self.width, self.height)
            
        else:
            
            self.figsize = settings.get('figsize')

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
        
        self.grid_hratios = self.grid_hratios or [1.] * self.grid_rows
        self.grid_wratios = self.grid_wratios or [1.] * self.grid_cols
        
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
    
    def iter_subplots(self):
        
        for j in xrange(self.grid_rows):
            
            for i in xrange(self.grid_cols):
                
                self.get_subplot(i, j)
                
                yield self.ax
    
    def make_plots(self):
        """
        By default this plots nothing here in the base class.
        Derived classes should override.
        """
        
        self.ax = get_subplot(0, 0)
        
        self.ax.plot(x = [], y = [])
        
        self.post_subplot_hook()
    
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
        
        self.ax.set_ylabel(self.ylab, fontproperties = self.fp_axis_label)
        self.ax.set_xlabel(self.xlab, fontproperties = self.fp_axis_label)
        # self.ax.yaxis.label.set_fontproperties(self)
    
    def set_ticklabels(self):
        
        if self.xticklabels is not None:
            
            self.ax.set_xticklabels(
                ('{0}'.format(x) for x in self.xticklabels),
            )
        
        _ = [
            tl.set_fontproperties(self.fp_ticklabel)
            for tl in self.ax.get_xticklabels()
        ]
        
        if self.xticks is not None:
            
            self.ax.set_xticks(self.xticks)
    
    def make_legend(self):
        
        if self.legend:
            
            self.leg = self.ax.legend(
                loc = self.legend_loc,
                title = self.legend_title,
            )
    
    def set_legend_font(self):
        
        if self.legend:
            
            _ = [
                t.set_fontproperties(self.fp_legend_label)
                for t in self.leg.get_texts()
            ]
            
            self.leg.get_title().set_fontproperties(self.fp_legend_title)
    
    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """
        
        #self.fig.tight_layout()
        self.fig.set_tight_layout(True)
        self.fig.subplots_adjust(top = .92)
        self.cvs.draw()
        self.cvs.print_figure(self.pdf)
        self.pdf.close()
        #self.fig.clf()


class _TestPlot(PlotBase):
    """
    Most minimal class to test the ``PlotBase`` class with plotting only
    a sinus function. Also serves as an example how to build plot classes
    on top of PlotBase.
    """
    
    def __init__(self, **kwargs):
        
        PlotBase.__init__(self, **kwargs)
    
    def make_plots(self):
        
        self.get_subplot(0, 0)
        
        x = np.linspace(0, 10, 1000)
        #self.ax.plot(x, np.sin(x))
        self.post_subplot_hook()


class Profiles(PlotBase):
    
    
    def __init__(
            self,
            profile_data,
            **kwargs,
        ):
        """
        Creates an optionally multi-faceted figure of intensity and
        background variable profiles.
        
        profile_data : list
            A list of ``ProfileData`` objects or one single object.
        """
        
        if (
            'grid_cols' not in kwargs and
            'grid_rows' not in kwargs and
            len(profile_data) > 1
        ):
            
            size = int(np.ceil(np.sqrt(len(profile_data))) ** 2)
            kwargs['grid_cols'] = size
            kwargs['grid_rows'] = size
        
        self.profile_data = profile_data
        
        PlotBase.__init__(self, **kwargs)
    
    
    def make_plots(self):
        
        for param, ax in zip(self.profile_data, self.iter_subplots()):
            
            self.param = param
            self.plot_profiles()
    
    
    def plot_profiles(self):
        
        self.var_barplots()
        self.features()
        self.var_lineplots()
        self.var_ranges()
        self.post_subplot_hook()
    
    
    def features(self):
        
        for i in xrange(self.param.features.shape[0]):
            
            color = (
                '#007B7F'
                    if self.param.colors is None else
                self.param.colors[i]
                    if isinstance(self.param.colors, np.ndarray) else
                self.param.colors
            )
            alpha = (
                1.0
                    if self.param.alphas is None else
                self.param.alphas[i]
                    if isinstance(self.param.alphas, np.ndarray) else
                self.param.alphas
            )
            linewidth = (
                .2
                    if self.param.linewidths is None else
                self.param.linewidths[i]
                    if isinstance(self.param.linewidths, np.ndarray) else
                self.param.linewidths
            )
            
            self.ax.plot(
                self.param.features_x,
                self.param.features[i,:],
                color = color ,
                alpha = alpha,
                linestyle = '-',
                linewidth = linewidth,
            )
        
        self.xticks = self.param.features_x
        self.xticklabels = self.param.labels
        self.title = self.param.title
        self.xlab = self.param.xlab
        self.ylab = self.param.ylab
        self.ylim = [0, 1]
    
    
    def var_barplots(self):
        
        if self.param.variables is None:
            
            return
        
        for var in self.param.variables:
            
            if var.typ != 'bar':
                
                continue
            
            self.ax.bar(
                x = var.x,
                height = var.y,
                color = var.color,
                alpha = var.alpha,
                edgecolor = 'none',
            )
    
    
    def var_lineplots(self):
        
        if self.param.variables is None:
            
            return
        
        for var in self.param.variables:
            
            if var.typ != 'line':
                
                continue
            
            self.ax.plot(
                x = var.x,
                y = var.y,
                color = var.color,
                alpha = var.alpha,
                linestyle = '-',
                marker = var.shape,
            )
    
    
    def var_ranges(self):
        
        if self.param.variables is None:
            
            return
        
        for var in self.param.variables:
            
            if var.typ != 'range':
                
                continue


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
    
    # plt.rcParams.update(params)
    
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


class SpectrumPlot(PlotBase):
    """
    Class for plotting a mass spectrum with annotations and intensities.
    
    Parameters
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
    
    def __init__(
            self,
            mzs,
            intensities = None,
            annotations = None,
            ionmode = None,
            sample_name = None,
            scan_id = None,
            precursor = None,
            intensities_percentage = True,
            annotate = False,
            **kwargs,
        ):
        """
        Creates a plot of a mass spectrum optionally with annotations for
        the peaks.
        
        Parameters
        ----------
        mzs : list,numpy.array
            m/z values of peaks in the spectrum.
        intensities : list,numpy.array
            Intensity values for each peak. If `None` all intensities
            will be plotted as 50%%.
        annotations : list,numpy.array
            String annotations for each peak. `None` values for peaks without
            annotation. If `None` no annotations will be plotted.
        ionmode : str
            Ionmode used in output file name and at fragment database lookup.
        sample_name : str
            Used in output file name.
        scan_id : int
            Used in output file name.
        precursor : float
            Used for fragment database lookup.
        intensities_percentage : bool
            Scale y-axis to 0-100%% or plot the intensity values as they are.
        annotate : bool
            Get annotations by fragment database lookup.
        kwargs
            Passed to ``PlotBase``.
        
        Returns
        -------
        index or None
        
        """
        
        self.ionmode = ionmode
        self.scan_id = scan_id
        self.sample_name = sample_name
        self.kwargs = kwargs
        
        self._set_fname()
        self._set_title()
        self._set_param('figsize')
        self._set_param('xlab')
        self._set_param('legend', False)
        
        self.mzs = np.array(mzs)
        
        self.intensities = intensities
        self.intensities_percentage = intensities_percentage
        self.set_intensities()
        
        self.annotations = annotations
        self.annotate = annotate
        self.set_annotations()
        
        PlotBase.__init__(self, **self.kwargs)
    
    def __len__(self):
        
        return len(self.mzs)
    
    def _set_fname(self):
        
        if 'fname' not in self.kwargs:
            
            self.kwargs['fname'] = (
                'ms2_spectrum%s%s%s.pdf' % (
                    ('__%s' % self.sample_name) if self.sample_name else '',
                    ('__%s' % self.ionmode) if self.ionmode else '',
                    ('__scan_%u' % self.scan_id) if self.scan_id else '',
                )
            )
    
    def _set_title(self):
        
        if 'title' not in self.kwargs:
            
            self.kwargs['title'] = (
                'Sample: %s, ionmode: %s, scan: %s' % (
                    self.sample_name if self.sample_name else 'unknown',
                    self.ionmode if self.ionmode else 'unknown',
                    ('%u' % self.scan_id) if self.scan_id else '?',
                )
            )
    
    def _set_param(self, param, value = None):
        
        if param not in self.kwargs:
            
            self.kwargs[param] = (
                settings.get('spectrum_plot_%s' % param)
                    if value is None else
                value
            )
    
    def get_annotations(self):
        """
        Queries the fragment database for fragment annotations.
        """
        
        if self.ionmode is None:
            
            raise ValueError(
                'SpectrumPlot: Ionmode is necessary for fragment lookup.'
            )
        
        self.annotations = list(fragdb.FragmentAnnotator(
            self.mzs,
            ionmode = self.ionmode,
            precursor = self.precursor,
        ))
    
    @staticmethod
    def _process_annotation(annots):
        
        if not annots:
            
            return None
            
        elif isinstance(annots, common.basestring):
            
            return annots
            
        else:
            
            # only a list of names
            # maybe we should include the fragment m/z
            return '\n'.join(annot.name for annot in annots)
    
    def process_annotations(self):
        """
        Creates strings from annotation arrays if necessary.
        """
        
        self.annotations = np.array([
            self._process_annotation(annots)
            for annots in self.annotations
        ])
    
    def set_annotations(self):
        """
        Sets the annotations for each fragment, optionally attempting
        database lookup.
        """
        
        if self.annotate and self.annotations is None:
            
            try:
                import lipyd.fragdb as fragdb
                self.get_annotations()
            except ImportError:
                session.log.msg(
                    'SpectrumPlot: Could not import `lipyd.fragdb`, '
                    'unable to get peak annotations for MS2 spectrum.'
                )
            except ValueError:
                session.log.msg(
                    'SpectrumPlot: Ionmode is necessary for fragment lookup.'
                )
        
        self.annotations = (
            np.array([None] * len(self))
                if self.annotations is None else
            np.array(self.annotations)
        )
        
        self.process_annotations()
    
    def set_intensities(self):
        """
        Sets the intensity vector: if no intensities provided all values
        will be 50%%, if necessary converts the values to percentage.
        """
        
        self.intensities = (
            np.array([50] * len(self))
                if self.intensities is None else
            np.array(self.intensities)
        )
        self.ylim = None
        self._set_param('ylab', 'Intensity')
        
        if self.intensities_percentage:
            
            self.intensities = (
                self.intensities / np.max(self.intensities) * 100
            )
            self.ylim = (0, 100)
            self._set_param('ylab', r'Relative intensity [%]')
    
    def make_plots(self):
        """
        Main method for drawing the plot.
        """
        
        max_x = np.max(self.mzs)
        max_y = np.max(self.intensities)
        
        #For saving plot in pdf format
        #pp = PdfPages(self.pdf_file_name)
        
        self.get_subplot(0, 0)
        
        self.ax.vlines(
            x = self.mzs,
            ymin = np.array([0.0] * len(self)),
            ymax = self.intensities,
        )
        
        trans = self.ax.get_xaxis_transform()
        
        vertex_number = 1
        annotate_base_x = 1.1 * max(self.mzs) # 10% more of max to right 
                                              # after max of X;
        annotate_base_y = 0.07
        annotate_offset_x = 0.
        annotate_lineheight = .06
        annotate_x = annotate_base_x
        annotate_y = annotate_base_y

        for mz, intens, annot in zip(
            self.mzs, self.intensities, self.annotations
        ):
            
            if not annot:
                
                continue
            
            self.ax.annotate('{}'.format(vertex_number),
                xy = (mz, intens),
                xycoords = 'data',
                xytext = (0, 0),
                textcoords = 'offset points', # "data"
                ha = 'center',
                va = 'center',
                bbox = dict(boxstyle = 'circle', fc = 'w'),
                fontproperties = self.fp_annotation,
            )
            
            self.ax.annotate('{} - {}'.format(vertex_number, annot),
                xy = (annotate_x, annotate_y),
                xycoords = trans,
                #va = 'top',
                bbox = dict(boxstyle = 'round', fc = 'w'),
                fontproperties = self.fp_annotation,
                horizontalalignment = 'left',
            )
            
            vertex_number += 1
            annotate_x += annotate_offset_x
            annotate_y += annotate_lineheight * annot.count('\n')
        
        self.post_subplot_hook()
        self.ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(50))
        self.ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))


if __name__ == "__main__":
    
    #Test data
    a = SpectrumPlot()
    x = [10., 11.1, 12.2, 13.3, 14.4, 15.5, 16.6, 17.7, 18.8, 19.9, 20.2]
    y = [100., 101., 150., 102., 104., 101., 180., 105., 102., 107., 103.]
    
    # Annotations done by use correlation between x coord and name of peak
    xa = ["the first peak", "the second peak"]
    
    a.build_plot(mzs = x,
                intensities = y,
                annotations = xa,
                result_type = "screen",
                title = "test")
