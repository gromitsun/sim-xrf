import numpy as np
import snip


def FWHM(x, noise=100, fano=0.114):
    sigma = np.sqrt((noise / 2.3548) ** 2 + 3.58 * fano * x)
    return 2.3548 * sigma


def pnb(ev, y, bg, peak_center, peak_width, bg_range=None, normalize=True):
    """

    :param ev: Energy array in eV.
    :param y: Total spectrum array.
    :param bg: Background spectrum array.
    :param peak_center: peak center in eV.
    :param peak_width: peak width in eV.
    :param bg_range: Range in which background is sampled. List with length of multiples of two.
    :param normalize: Boolean. If true, background is normalized with respect to the width of bg compared to width of y.
    :return: P,B in P/B.
    """

    def ev2ch(x):
        return np.abs(ev - x).argmin()

    p = np.array([ev2ch(peak_center - peak_width / 2.), ev2ch(peak_center + peak_width / 2.)])
    if bg_range is None:
        b = p
    else:
        b = np.array([ev2ch(x) for x in bg_range])

    Nt = np.sum(y[p[0]:p[1]])
    B = 0
    nchB = 0
    for i in range(0, len(b), 2):
        B += np.sum(bg[b[i]:b[i + 1]])
        nchB = b[i + 1] - b[i]
    if normalize:
        B *= 1. * (p[1] - p[0]) / nchB
    return Nt - B, B


def pnb_fixedwidth(ev, y, peak_center, peak_width, bg_range):
    return pnb(ev, y, y, peak_center, peak_width, bg_range)


def pnb_snip(ev, y, peak_center, peak_width=100, FWHM=snip.FWHM, offset=0., gain=10., **kwargs):
    bg = snip.snip(ev, y, FWHM, offset, gain, **kwargs)
    return pnb(ev, y, bg, peak_center, peak_width)


def pnb_ascalc(ev, y, total, peak_center, peak_width):
    return pnb(ev, total, total - y, peak_center, peak_width, normalize=False)