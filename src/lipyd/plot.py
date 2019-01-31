#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2019 - EMBL
#
#  File author(s):
#  Dénes Türei (turei.denes@gmail.com)
#  Igor Bulanov
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
import matplotlib.backends.backend_cairo
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
            format = 'pdf',
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
        
        if self.format == 'pdf':
            
            self.pdf = mpl.backends.backend_pdf.PdfPages(self.fname)
            self.fig = mpl.figure.Figure(figsize = self.figsize)
            self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)
        
        elif self.format == 'png':
            
            self.fig = mpl.figure.Figure(figsize = self.figsize)
            self.cvs = mpl.backends.backend_cairo.FigureCanvasCairo(self.fig)

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
        
        if self.format == 'pdf':
            
            self.cvs.print_figure(self.pdf)
            self.pdf.close()
            
        elif self.format == 'png':
            
            with open(self.fname, 'wb') as fp:
                
                self.cvs.print_png(fp)
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
            np.array( self.annotations )
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
        
        # max_x = np.max(self.mzs)
        # max_y = np.max(self.intensities)
        
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
        max_mz = max(self.mzs)
        annotate_base_x = 1.07 * max_mz # 10% more of max to right 
                                              # after max of X;
        annotate_base_y = 0.99
        detail_annotation = ""
        number_width = 3
        number_separator = 4    #number of space char between: N) and annotation

        annot_list = []
        
        for mz, intens, annot in zip(
                                    self.mzs,
                                    self.intensities,
                                    self.annotations
                                ):
                                    
            if not annot:
                continue
            annot_list.append( (mz, intens, annot) )

        #create layout list:
        import lipyd.plot_layout as plout

        layout_shape = plout.Layout_shape(plot_x_max = max_mz,
                                        annot_list = annot_list
                    )
        
        for mz, intens, annot in annot_list:
            xx, yy = layout_shape.get_free_place(mz, intens)
            
            self.ax.annotate('{}-{:.2f}'.format(vertex_number, mz),
            #self.ax.annotate('{}'.format(vertex_number),
                xy = ( mz, intens ),    #mz, intens
                #xycoords = 'data',
                xytext = (xx, yy),
                #textcoords = 'offset points',
                ha = 'center',
                va = 'center',
                bbox = dict(boxstyle = 'round', fc = 'w'), #circle
                fontproperties = self.fp_annotation,
                fontsize=6,
                arrowprops=dict(arrowstyle="->"),
            )
            # +1 for char ")"
            # replace char new line on new line char 
            # and serial of space char:
            formated_annot = annot.replace("\n", \
                "\n"+(number_width+number_separator+4)*" ")
            #
            detail_annotation += "{:3d}){}{}\n".format(vertex_number, \
                (number_separator+1) * " ", formated_annot)
            
            vertex_number += 1

        detail_annotation=detail_annotation[:-1] #remove last new line char;
        
        self.ax.annotate(detail_annotation,
            xy = (annotate_base_x, annotate_base_y),
            xycoords = trans,
            bbox = dict(boxstyle = 'round', fc = 'w'),
            fontproperties = self.fp_annotation,
            horizontalalignment = 'left',
            verticalalignment = 'top',
            fontsize=6
        )

        self.post_subplot_hook()
        self.ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(50))
        self.ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))

        