from astropy.io import fits


def obtain_coordinates_2(lista, idx):
    x1 = lista[idx, 0]
    x2 = lista[idx, 1]
    y1 = lista[idx, 2]
    y2 = lista[idx, 3]
    output_tuple = (x1, x2, y1, y2)
    return output_tuple


def obtain_coordinates_ccd(image_, mypath_=False):
    """
    Given a string following the CAFOS structure, obtains the coordinates of the ccd (x1, x2, y1 and y2) as a tuple.

    :param image_:
    :param mypath_:
    :return:
    """
    if mypath_:
        str_ccdsec = fits.open(mypath_ + image_)[0].header['CCDSEC']
        longitud = len(str_ccdsec)
    else:
        str_ccdsec = image_
        longitud = len(str_ccdsec)+1

    comma = None
    for x in range(1, longitud):
        if str_ccdsec[x] == ',':
            comma = x
            break
    if comma is None:
        raise ValueError('comma not defined!')

    colon = None
    for x in range(comma + 1, longitud):
        if str_ccdsec[x] == ':':
            colon = x
            break
    if colon is None:
        raise ValueError('colon not defined!')

    coma2 = None
    for x in range(colon + 1, longitud):
        if str_ccdsec[x] == ',':
            coma2 = x
            break
    if coma2 is None:
        raise ValueError('coma2 not defined!')

    x1 = int(str_ccdsec[1:comma])
    y1 = int(str_ccdsec[comma + 1: colon])
    x2 = int(str_ccdsec[colon + 1: coma2])
    y2 = int(str_ccdsec[coma2 + 1: -1])
    tupla_salida = (x1, x2, y1, y2)
    return tupla_salida
