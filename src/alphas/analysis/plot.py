import typing
import matplotlib as _mpl
import seaborn as _sns
import gvar as _gvar

_plt = _mpl.pyplot

_plt.rcParams['text.usetex'] = True
_plt.rcParams["mathtext.default"] = "regular"

class Plot(object):
    def __init__(
            self,
            w: float | int = 7.5,
            h: float | int = 4.5,
            grid: bool = True,
            mnrgridlw: float | int = 0.8,
            mjrgridlw: float | int = 1.0,
            mnrgridclr: str = '#EEEEEE',
            mjrgridclr: str = '#DDDDDD',
            xticklblfs: float | int = 17.5,
            yticklblfs: float | int = 17.5,
            xlblfs: float | int = 20.,
            ylblfs: float | int = 20.,
            xoffsetfs: float | int = 17.5,
            yoffsetfs: float | int = 17.5,
            titlefs: float | int = 17.5
    ):
        self._fig = _plt.figure()
        self._fig.set_size_inches(w,h)
        gridspec = self._fig.add_gridspec(nrows=1,ncols=1)
        self._ax = self._fig.add_subplot(gridspec[0,0])

        if grid:
            self._ax.grid(
                which = 'minor',
                color = mnrgridclr,
                linestyle = ':',
                linewidth = mnrgridlw
            )
            self._ax.grid(
                which = 'major',
                color = mjrgridclr,
                linewidth = mjrgridlw
            )
            self._ax.minorticks_on()
            
        self._xticklblfs = xticklblfs
        self._yticklblfs = yticklblfs

        self._xlblfs = xlblfs
        self._ylblfs = ylblfs

        self._xoffsetfs = xoffsetfs
        self._yoffsetfs = yoffsetfs

        self._titlefs = titlefs
        self._nlgs = 0

    def color_palette(self, nm: str, lngth: int, **kwargs):
        return _sns.color_palette(nm, lngth, **kwargs)
        
    def scatter(self, x: list[any], y: list[any],**kwargs):
        return self._ax.scatter(x, y, **kwargs)

    def errorbar(self, x: list[any], y: list[any], **kwargs):
        return self._ax.errorbar(_gvar.mean(x), _gvar.mean(y), yerr=_gvar.sdev(y), **kwargs)

    def line(self, x: list[any], y: list[any], **kwargs):
        return self._ax.plot(x, y, **kwargs)

    def fill_between(
            self,
            x: list[any],
            y: list[any],
            color = 'k',
            alpha = 0.5,
            **kwargs
    ): return self._ax.fill_between(
        x,
        _gvar.mean(y)-_gvar.sdev(y),
        _gvar.mean(y)+_gvar.sdev(y),
        color = color,
        alpha = alpha,
        **kwargs
    )

    def text(
            self,
            x: float | int,
            y: float | int,
            text: str,
            va: str = 'center',
            ha: str = 'center',
            fs: float = 20.,
            transform = 'axes',
            **kwargs
    ):
        if 'transform' in kwargs.keys():
            self._ax.text(
                x, y,
                text,
                va = va,
                ha = ha,
                fontsize = fs,
                **kwargs
            )
        elif transform == 'axes':
            self._ax.text(
                x, y,
                text,
                va = va,
                ha = ha,
                fontsize = fs,
                transform = self._ax.transAxes,
                **kwargs
            )
        else:
            self._ax.text(
                x, y,
                text,
                va = va,
                ha = ha,
                fontsize = fs,
                **kwargs
            )

    def _check_legend_and_add(self, legend):
        #nlgs = len([
        #    None for c in self._ax.get_children()
        #    if isinstance(c, _mpl.legend.Legend)
        #])
        #if nlgs > 1: self._ax.add_artist(legend)
        #if self._nlgs > 0: self._ax.add_artist(legend)
        pass
        
    def legend(
            self,
            x: float | int,
            y: float | int,
            loc = 'center',
            ncol: int = 1,
            framealpha: float | int = 0.,
            fs: float | int = 20.,
            **kwargs
    ):
        self._check_legend_and_add(
            self._ax.legend(
                bbox_to_anchor = (x,y),
                loc = loc,
                ncol = ncol,
                framealpha = framealpha,
                fontsize = fs,
                **kwargs
        ))
        self._nlgs += 1

    def special_legend(
            self,
            x: float | int,
            y: float | int,
            marker_handles: list[any],
            marker_labels: list[any],
            handler_map: dict[any,any] = {
                tuple: _mpl.legend_handler.HandlerTuple(ndivide = None)
            },
            loc = 'center',
            ncol: int = 1,
            framealpha: float | int = 0.,
            fs: float | int = 20.,
            **kwargs
    ):
        _ax = self._ax.twinx()
        self._check_legend_and_add(
            _ax.legend(
                marker_handles,
                marker_labels,
                handler_map = handler_map,
                bbox_to_anchor = (x,y),
                loc = loc,
                ncol = ncol,
                framealpha = framealpha,
                fontsize = fs,
                **kwargs
        ))
        _ax.set_yticks([])
        _ax.get_yaxis().set_visible(False)
        self._nlgs += 1
        
    def decorate(
            self,
            xlim: list[float] | list[int] = None,
            ylim: list[float] | list[int] = None,
            xlabel: list[float] | list[int] = None,
            ylabel: list[float] | list[int] = None,
            
    ):
        if xlim is not None: self._ax.set_xlim(xlim)
        if ylim is not None: self._ax.set_ylim(ylim)    

        if xlabel is not None:
            self._ax.set_xlabel(xlabel, fontsize = self._xlblfs)
        if ylabel is not None:
            self._ax.set_ylabel(ylabel, fontsize = self._ylblfs)

        _plt.setp(self._ax.get_xticklabels(), fontsize = self._xticklblfs)
        _plt.setp(self._ax.get_yticklabels(), fontsize = self._yticklblfs)

        self._ax.xaxis.offsetText.set_fontsize(self._xoffsetfs)
        self._ax.yaxis.offsetText.set_fontsize(self._yoffsetfs)

        self._ax.set_axisbelow(True)

    def set_title(self, title: str, **kwargs):
        if 'fontsize' not in kwargs: kwargs['fontsize'] = self._titlefs
        self._ax.set_title(title,**kwargs)

    def modify_axis(self, axis: str, ticks, labels, **kwargs):
        if axis == 'x': self._ax.set_xticks(ticks,labels,**kwargs)
        elif axis == 'y': self._ax.set_yticks(ticks,labels,**kwargs)

    #def modify_labels(self, axis: str, ticks, labels, **kwargs):
    #    if axis == 'x': self._ax.set_xlabels(labels,**kwargs)
    #    elif axis == 'y': self._ax.set_ylabels(labels,**kwargs)

    @property 
    def axis(self): return self._ax

    @property
    def handle(self): return self._fig
        
    def __call__(
            self,
            save: bool = False,
            show: bool = False,
            fn: str = '',
            **kwargs
    ):
        if save: self._fig.savefig(fn, bbox_inches = 'tight', transparent = True, **kwargs)
        if show: _plt.show()
        _plt.close(self._fig)
