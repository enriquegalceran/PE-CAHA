import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
#https://stackoverflow.com/questions/44865023/circular-masking-an-image-in-python-using-numpy-arrays
lado = 2048
a = np.arange(lado**2).reshape(lado,lado)
x1 = 1224
x2 = 1624
y1 = 1
y2 = 200
b = np.copy(a[y1:y2, x1:x2])
print(a.shape)
print(b.shape)

def create_circular_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

h, w = b.shape[:2]

center_real = [1024,1024]

center = [center_real[0]-x1, center_real[1]-y1]
# Esto habra que cambiarlo, porque las pruebas se han hecho usando zonas de numpy, que empiezan desde arriba. Si despues se usa las coordenadas del CCD, que empiezan desde abajo, habra que cambiar la y1 por y2

print(center)
radius = 1024
mask = create_circular_mask(h, w, center=center, radius=radius)
masked_img = b.copy()

masked_img[~mask] = 0
x = ma.masked_array(b, ~mask)
mediana = ma.median(x)
print(mediana)

#print('mask')
#print(mask)
#print('x')
#print(x)
#plt.figure()
#plt.imshow(a)


plt.figure()
plt.imshow(masked_img)
plt.title('masked_img')

plt.figure()
plt.imshow(x)
plt.title('x')

plt.figure()
plt.imshow(mask)
plt.title('mask')

plt.show()
