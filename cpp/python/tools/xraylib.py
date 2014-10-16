__author__ = 'Yue'
import ctypes
from .. import libpath

try:
    lib = ctypes.cdll.LoadLibrary(libpath + '/../xraylib/Lib/linux/libxrl.so')
except OSError:
    lib = ctypes.cdll.LoadLibrary(libpath + '/../xraylib/Lib/windows/libxrl-7.dll')

lib.CS_Total_CP.restype = ctypes.c_double
lib.CS_Total_CP.argtypes = [ctypes.c_char_p, ctypes.c_double]


def CS_Total_CP(CP, kev):
    assert isinstance(CP, str)
    return lib.CS_Total_CP(CP, kev)


def ElementDensity(*args):
    return lib.ElementDensity(*args)


def SymbolToAtomicNumber(*args):
    return lib.SymbolToAtomicNumber(*args)