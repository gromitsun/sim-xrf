from . import __path__ as datapath

datapath = datapath[0]

def lookup(x):
    """

    :param x:
    :return:
    """
    f = open(datapath + '/periodictable.txt', 'r')  # open periodic table file
    data = f.readlines()
    f.close()

    # dictionaries
    d = {}
    ivd1 = {}
    ivd2 = {}

    # create dictionaries from data in the file
    for line in data:
        Z, symbol, name = (s.strip() for s in line.split('-'))
        Z = int(Z)
        d[Z] = symbol, name
        ivd1[symbol] = Z, name
        ivd2[name] = Z, symbol

    # lookup in the dictionaries
    try:
        x = x.title()
        try:
            return ivd1[x]
        except KeyError:
            return ivd2[x]
    except AttributeError:
        return d[x]


def symbol2number(Z):
    if type(Z) == str:
        try:
            Z = int(Z)
        except ValueError:
            return lookup(Z)[0]
    return Z


def number2symbol(Z):
    if type(Z) == str:
        try:
            Z = int(Z)
        except ValueError:
            return Z
    return lookup(Z)[0]
