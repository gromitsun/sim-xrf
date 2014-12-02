# detector.py
import numpy as np

try:
    import xraylib as xrl
except ImportError:
    from ..tools import xraylib as xrl


def energy_to_channel(energy, offset=0., gain=10.):
    """
    Convert energy to channel number.
    :param energy:
        The energy in eV.
    :param offset:
        Energy offset in eV. Or the energy of channel 0.
    :param gain:
        Energy gain in eV. The increment of energy in eV from channel to channel.
    :return:
        Channel number.
    """
    return (energy - offset) // gain


class Channel(object):
    """
    Class for detector channels.
    Attributes:
        ev_offset:  Energy offset in eV.
        ev_gain:    Channel gain in eV.
        n_channels: Total number of channels.
    """

    def __init__(self, ev_offset=0, ev_gain=10, n_channels=2048):
        self.ev_offset = ev_offset
        self.ev_gain = ev_gain
        self.n_channels = n_channels
        self.ev_arr = np.arange(ev_offset, ev_offset + ev_gain * n_channels, ev_gain)

    def channel(self, ev):
        """
        Convert energy in eV to channel number.
        :param ev:
            The energy in eV.
        :return:
            Channel number
        """
        return energy_to_channel(ev, self.ev_offset, self.ev_gain)


class Response(object):
    """
    Class for detector response.
    Attributes:
        noise:
        fano:
        fs:
        ft:
        ev_gain:
    """

    def __init__(self,
                 noise=100,
                 fano=0.114,
                 gamma=2.5,
                 fs=0.03,
                 ft=0.02,
                 ev_gain=None):
        self.noise = noise
        self.fano = fano
        self.gamma = gamma
        self.fs = fs
        self.ft = ft
        self.ev_gain = ev_gain

    def FWHM(self, ev, **kwargs):
        """
        FWHM function (of energy in eV).
        :rtype : float
        :param ev: Energy in eV.
        :param kwargs:
        :return:
            FWHM in eV.
        """
        noise = kwargs.get('noise', self.noise)
        fano = kwargs.get('fano', self.fano)
        sigma = np.sqrt((noise / 2.3548) ** 2 + 3.58 * fano * ev)
        return 2.3548 * sigma


class Window(object):
    """
    Class for detector window.
    :Attributes:
        material:
            String. Material of the filtering window on the detector.
        thickness:
            Thickness of the filtering window material. In cm.
        density:
            Mass density of the filtering window material. In g/cm^3.
    """

    def __init__(self,
                 material='Be',
                 thickness=24e-4,
                 density=None):
        self.material = material
        self.thickness = thickness
        if density is None:
            density = xrl.ElementDensity(xrl.SymbolToAtomicNumber(material))
            if not density:
                density = 1
        self.density = density

    def transmission(self, ev):
        """
        The transmission function of the window.
        :param ev:
            The energy in eV of the incident beam.
        :return:
            The ratio of the intensity of the transmitted beam through the window vs that of the incident.
        """
        _mac = xrl.CS_Total_CP(self.material, ev / 1000.)
        # _mac = compound(CP=self.material).mac_total(ev)
        return np.exp(-_mac * self.density * self.thickness)


class Detector(object):
    """
    Class for detector.
    :Attributes:
        channel: A detector channel object.
        response:   A detector response object.
        window: A detector window object.
    """

    def __init__(self, channel=Channel(), response=Response(), window=Window()):
        self.channel = channel
        self.response = response
        response.ev_gain = channel.ev_gain
        self.window = window