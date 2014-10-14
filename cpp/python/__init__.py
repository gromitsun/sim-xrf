__author__ = 'Yue'
__all__ = ["pyapi"]
path = __path__[0]
libpath = path + '/../Lib'

# Load .pyconfig file
configf = open(path + '/../.pyconfig')
config = {}
for line in configf.readlines():
    key, value = line.split('=')
    try:
        config[key.strip()] = [int(x) for x in value.split('#')[0].split(',')]
    except ValueError:
        config[key.strip()] = [float(x) for x in value.split('#')[0].split(',')]

configf.close()

plt_kwargs = {key: value for key, value in config.items()}
del (plt_kwargs['nout'])