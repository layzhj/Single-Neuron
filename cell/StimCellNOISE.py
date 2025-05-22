from neuron import h


class StimCellNOISE(object):
    def __init__(self, vspkNOISE, gid=-1):
        super().__init__()
        self._gid = gid
        self.ncstim = []
        self.vspkNOISE = vspkNOISE
        self._setup_biophysics()
        self.spike_detector = h.NetCon(self.stim, None)
        self.spike_times = h.Vector()
        self.spike_detector.record(self.spike_times)

    def _setup_biophysics(self):
        self.soma = h.Section(name='soma', cell=self)
        self.stim = h.VecStim(self.soma)
        self.stim.play(self.vspkNOISE)
        self.stim.delay = 0

    def connect2target(self, target, delay=1, weight=0.04):
        self.ncstim.append(h.NetCon(self.stim, target))
        self.ncstim[-1].delay = delay
        self.ncstim[-1].weight[0] = weight
        return self.ncstim[-1]