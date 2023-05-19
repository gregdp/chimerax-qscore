
from chimerax.ui.gui import MainToolWindow

from Qt.QtWidgets import (
    QFrame, QLabel,
    QPushButton, QMenu, QRadioButton, QScrollBar,
    QSpinBox, QDoubleSpinBox, QCheckBox,
    QHBoxLayout, QVBoxLayout, QGridLayout
)
from Qt import QtCore
from Qt.QtCore import Qt


class DefaultVLayout(QVBoxLayout):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setContentsMargins(0,2,0,0)
        self.setSpacing(3)

class DefaultHLayout(QHBoxLayout):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setContentsMargins(0,0,0,0)
        self.setSpacing(3)


class QScorePlot(QFrame) :

    MIN_ZOOM = 25
    DEFAULT_ZOOM = 250
    ZOOM_SCALE = 1.1
    MAX_ZOOM_STEPS_PER_FRAME=5

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg as FigureCanvas,
            )
        from matplotlib.widgets import Slider
        from matplotlib.colors import Normalize
        import numpy
        self.setMinimumHeight(150)
        self._slider_blocked = False
        ml = self.main_layout = DefaultVLayout()
        self.setLayout(ml)

        fig = self.plot = Figure()
        axes = self.axes = fig.add_subplot(111)
        axes.autoscale(enable=False)
        axes.set_ylim(0,1)
        axes.set_xlim(0, self.DEFAULT_ZOOM)
        fig.subplots_adjust(bottom=0.25)

        self.volume = None

        fig.tight_layout(rect=(0,0.1,1,1))

        self.residues = []
        resnum = self.residue_numbers = numpy.zeros(2, dtype=numpy.int32)-1
        qscore = self.qscores = numpy.array([0,1], dtype=numpy.float32)

        def _picker(scatter, mouse_event, x_radius=1):
            from matplotlib.backend_bases import MouseButton
            if mouse_event.button != MouseButton.LEFT:
                return False, dict()
            import numpy
            if mouse_event.inaxes != self.axes:
                return False, dict()
            if mouse_event.xdata is None:
                return False, dict()
            x = mouse_event.xdata
            resnums = scatter.get_offsets().T[0]
            d = numpy.sqrt((x - resnums)**2)
            closest = d.min()
            if closest > x_radius:
                return False, dict()
            ind = numpy.argmin(d)
            return True, dict(ind=ind)


        s = axes.scatter(resnum, qscore, s=4, c=qscore, cmap='inferno_r', picker=_picker)
        self._scatter = s

        canvas = self.canvas = FigureCanvas(fig)

        canvas.mpl_connect('scroll_event', self.zoom)
        canvas.mpl_connect('pick_event', self.on_pick)

        ml.addWidget(canvas)

        hpos = self._hpos_slider = QScrollBar(Qt.Orientation.Horizontal)
        hpos.setRange(0,0)
        ml.addWidget(hpos)

        hpos.valueChanged.connect(self._slider_update)

        canvas.draw()

    def initialize_hover_text(self, target):
        self._hover_text_target = target

        def _hover(event, x_radius=1):
            if not len(self.residues):
                target.setText('')
                return
            if event.inaxes != self.axes:
                target.setText('')
                return
            import numpy
            x = event.xdata
            resnums = self._scatter.get_offsets().T[0]
            d = numpy.sqrt((x-resnums)**2)
            closest = d.min()
            if closest > x_radius:
                target.setText('')
                return
            ind = numpy.argmin(d)
            r = self.residues[ind]
            if r.deleted:
                target.setText('{residue deleted}')
                return
            target.setText(f'{r.name} /{r.chain_id}:{r.number}\tQ: {self.qscores[ind]:.3f}')
            return
        self.canvas.mpl_connect('motion_notify_event', _hover)


    def _slider_update(self, val):
        if not len(self.residues) or self._slider_blocked:
            return
        axes = self.axes
        xmin, xmax = axes.get_xlim()
        xrange = xmax-xmin
        axes.set_xlim([val,val+xrange])
        self.canvas.draw_idle()

    def zoom(self, event=None):
        if not len(self.residues):
            return
        axes = self.axes
        hpos = self._hpos_slider
        xmin, xmax = axes.get_xlim()
        xrange = xmax-xmin
        resnum = self.residue_numbers
        if xmax > resnum.max():
            xmax = resnum.max()
            xmin = xmax-xrange
        if xmin < resnum.min():
            xmin = resnum.min()
            xmax = xmin+xrange

        if event is not None:
            if event.inaxes != self.axes:
                return
            xpoint = event.xdata
            xfrac = (xpoint-xmin)/(xmax-xmin)
            if event.button == 'up':
                if xrange <= self.MIN_ZOOM:
                    return
                xrange = int(xrange/(self.ZOOM_SCALE*min(event.step, self.MAX_ZOOM_STEPS_PER_FRAME)))
            else:
                if xrange >= resnum.max()-resnum.min():
                    return
                xrange = int(xrange*self.ZOOM_SCALE*-min(event.step, self.MAX_ZOOM_STEPS_PER_FRAME))
            new_xmin = int(max(min(resnum), min(resnum.max()-xrange, int(xpoint-xrange*xfrac))))
        else:
            new_xmin = max(resnum.min(), xmin)
        new_xmax = min(max(resnum), new_xmin+xrange)
        axes.set_xlim([new_xmin, new_xmax])
        if xrange >= resnum.max()-resnum.min():
            hpos.setRange(resnum.min(), resnum.min())
        else:
            with slot_disconnected(hpos.valueChanged, self._slider_update):
                hpos.setRange(resnum.min(), resnum.max()-xrange)
                hpos.setValue(xmin)
        self.canvas.draw_idle()

    def zoom_extents(self):
        resnum = self.residue_numbers
        self.axes.set_xlim([resnum.min(), resnum.max()])
        self.zoom()

    def on_pick(self, event):
        if not len(self.residues):
            return
        ind = event.ind
        residue = self.residues[ind]
        if residue.deleted:
            return
        session = residue.session
        session.selection.clear()
        residue.atoms.selected = True
        residue.atoms.intra_bonds.selected = True
        atomspec = f'#!{residue.structure.id_string}/{residue.chain_id}:{residue.number}'
        from chimerax.core.commands import run
        from .clipper_compat import model_managed_by_clipper
        m = residue.structure

        if model_managed_by_clipper(m):
            # Just view the model
            run(session, f'view {atomspec}')
        else:
            # TODO: decide what to do here
            from chimerax.atomic import Residues, concise_residue_spec
            neighbors = set([residue])
            # Quick and (very) dirty way to expand the selection. Should probably do something more efficient.
            for _ in range(3):
                new_neighbors = []
                for n in neighbors:
                    for nn in n.neighbors:
                        if nn not in neighbors:
                            new_neighbors.append(nn)
                neighbors.update(new_neighbors)
            residues = Residues(neighbors)
            argspec = concise_residue_spec(session, residues)
            run(session, f'surf zone #{self.volume.id_string} near {argspec} dist 3', log=False)
            run(session, f'~cartoon #{m.id_string}; hide #{m.id_string}; show {argspec}; cartoon {argspec}&~{atomspec}', log=False)
            run(session, f'view {atomspec}', log=False)

    def update_data(self, residues, scores, volume):
        self.volume = volume
        if residues is None or not len(residues):
            self.residues = []
            self._scatter.set_offsets([[0,0]])
        else:
            # Convert to a list so we can gracefully handle residue deletions in callbacks
            import numpy
            update_zoom = False
            if len(residues) != len(self.residues):
                update_zoom = True
            self.residues = list(residues)
            rmin = min(r.number for r in residues)
            rmax = max(r.number for r in residues)
            resnum = self.residue_numbers = numpy.array([r.number for r in residues])
            self.qscores = scores

            clipped_scores = numpy.array(scores, copy=True)
            clipped_scores[clipped_scores<0] = 0
            xy = numpy.array([resnum, clipped_scores]).T

            self._scatter.set_offsets ( xy )
            self._scatter.set_array ( clipped_scores )
            # If the number of plotted residues changes, zoom to the full extent of the data.

            if update_zoom:
                self.zoom_extents()
            else:
                self.zoom(None)
        self.canvas.draw_idle()



from contextlib import contextmanager
@contextmanager
def slot_disconnected(signal, slot):
    '''
    Temporarily disconnect a slot from a signal using

    .. code-block:: python

        with slot_disconnected(signal, slot):
            do_something()

    The signal is guaranteed to be reconnected even if do_something() throws an error.
    '''
    try:
        # disconnect() throws a TypeError if the method is not connected
        signal.disconnect(slot)
        yield
    except TypeError:
        pass
    finally:
        signal.connect(slot)
