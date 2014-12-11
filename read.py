import numpy as np
import matplotlib.pyplot as plt


def read(fname):
    # Read channel values
    a = np.genfromtxt(fname, invalid_raise=False, skiprows=0)

    # Read detector configurations and line labels
    labels = ['Total']
    f = open(fname)
    for line in f.readlines():
        if 'ev_offset' in line:
            ev_offset = float(line.split('=')[1])
        elif 'ev_gain' in line:
            ev_gain = float(line.split('=')[1])
        elif 'n_channels' in line:
            n_channels = float(line.split('=')[1])
        elif 'Z =' in line:
            labels.append(line.split('(')[-1].strip().strip(')'))
    f.close()
    labels += ['Rayleigh', 'Compton']

    ev = ev_offset + ev_gain * np.arange(n_channels)

    return ev, a, labels

def plot(ev, a, labels, xlim=None, ylim=None, show=True):
    ax = plt.subplot(111)
    plt.plot(ev / 1e3, a.T)
    plt.yscale('log')
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.xlabel('E (KeV)')
    plt.ylabel(r'$I(E)/I_0(E_0)$')
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plt.legend(labels, loc='upper right', ncol=1, bbox_to_anchor=(1.35, 1), borderaxespad=0.)
    if show:
        plt.show()

def dose(fname):
    f = open(fname)
    dose = []
    for line in f.readlines():
        if ('Layer' in line) and ('Gy' in line):
            dose.append(float(line.split(':')[-1].split()[0]))
    f.close()
    return np.array(dose)



if __name__ == '__main__':
    import sys
    from python import plt_kwargs

    try:
        fname = sys.argv[1]
    except IndexError:
        fname = 'out.txt'

    # Read data file
    ev, a, labels = read(fname)

    # Plot spectrum
    plt.figure()
    plt.title('Spectrum read from %s' % fname)
    plot(ev, a, labels, xlim=plt_kwargs['xlim'], ylim=plt_kwargs['ylim'])