import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from Salida_limpia import mostrarresultados, stdrobusta
from astropy.io import fits
import IMGPlot as ImP


def deshacer_tupla_coord(tupla):
    return tupla[0], tupla[1], tupla[2], tupla[3]


def sacar_coordenadas_ccd(imagen_, mypath_=False):
    if mypath_:
        str_ccdsec = fits.open(mypath_ + imagen_)[0].header['CCDSEC']
        longitud = len(str_ccdsec)
    else:
        str_ccdsec = imagen_
        longitud = len(str_ccdsec)+1
    for x in range(1, longitud):
        if str_ccdsec[x] == ',':
            coma = x
            break
    for x in range(coma + 1, longitud):
        if str_ccdsec[x] == ':':
            puntos = x
            break
    for x in range(puntos + 1, longitud):
        if str_ccdsec[x] == ',':
            coma2 = x
            break

    x1 = int(str_ccdsec[1:coma])
    y1 = int(str_ccdsec[coma + 1: puntos])
    x2 = int(str_ccdsec[puntos + 1: coma2])
    y2 = int(str_ccdsec[coma2 + 1: -1])
    tupla_salida = (x1, x2, y1, y2)
    return tupla_salida


def create_circular_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

#https://stackoverflow.com/questions/44865023/circular-masking-an-image-in-python-using-numpy-arrays
lado = 2048
a = np.arange(lado**2).reshape(lado,lado)
x1 = 1224
x2 = 1624
y1 = 1
y2 = 200
b = np.copy(a[y1:y2, x1:x2])


image_file = 'imagencopiada.fits'
ccdsec = fits.open(image_file)[0].header['CCDSEC']
image_data = fits.getdata(image_file, ext=0)
coordenadas_ = sacar_coordenadas_ccd(ccdsec)
x1, x2, y1, y2 = deshacer_tupla_coord(coordenadas_)


h, w = image_data.shape[:2]
center_real = [1075, 1040]
center = [center_real[0]-x1, center_real[1]-y1]

radius = 809 # con 810 tiene un pixel de borde
mask = create_circular_mask(h, w, center=center, radius=radius)
masked_img = image_data.copy()

masked_img[~mask] = 0
x = ma.masked_array(image_data, ~mask)
mediana = ma.median(x)

mostrarresultados(['ccdsec','h', 'w', 'center_real', 'center', 'radius', 'mediana'],
                  [ccdsec, h, w, center_real, center, radius, mediana])


plt.figure()
plt.imshow(masked_img, origin='low')
plt.title('masked_img - array[~mask]')

plt.figure()
plt.imshow(x, origin='low')
plt.title('x - ma.masked_array(array, ~mask)')

plt.figure()
plt.imshow(mask, origin='low')
plt.title('mask - create_circular_mask(h,w,center,radius)')

plt.show()
