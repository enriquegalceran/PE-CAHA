from astropy.io import fits
from Salida_limpia import mostrarresultados
import numpy as np
import IMGPlot as ImP
# import pandas as pd
import matplotlib.pyplot as plt
import os
import csv
from numpy import genfromtxt
import argparse
import time


def deshacer_tupla_coord(tupla):
    return tupla[0], tupla[1], tupla[2], tupla[3]


def listas_archivos(path_):
    lista_cal = []
    lista_sci = []
    lista_misc = []
    listaarchivos = []
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            if '-cal-' in file:
                lista_cal.append(file)
            elif '-sci-' in file:
                lista_sci.append(file)
            else:
                lista_misc.append(file)
            listaarchivos.append(os.path.join(path_, file))
    return lista_cal, lista_sci, lista_misc, listaarchivos


def guardar_listas_csv(csvfile, res):
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in res:
            writer.writerow([val])


def listas_calibracion(path, calib, desc_bias, desc_flats, desc_misc):
    lista_bias = []
    lista_flats = []
    lista_misc = []
    
    for file in calib:
        for texto in desc_bias:
            if texto in fits.open(path + file)[0].header['OBJECT']:
                lista_bias.append(file)
                break
        for texto in desc_flats:
            if texto in fits.open(path + file)[0].header['OBJECT']:
                lista_flats.append(file)
                break
        for text in desc_misc:
            if text in fits.open(path + file)[0].header['OBJECT']:
                lista_misc.append(file)    
    return lista_bias, lista_flats, lista_misc


def crear_listas_cal_y_sci(lista_noches_, dir_listas_, dir_datos_, desc_bias, desc_flats, desc_misc):
    i = 0
    for noche in lista_noches_:
        i += 1
        if noche not in os.listdir(dir_listas_):
            os.mkdir(dir_listas_ + noche + '/')
            # Quiza hay que descomentar para que funcione
            # CAL = []
            # SCI = []
            # MIS = []
            # ARC = []
            path_ = dir_datos_ + noche + '/'
            cal, sci, mis, arc = listas_archivos(path_)
            l_bias, l_flat, l_misc = listas_calibracion(path_, cal, desc_bias, desc_flats, desc_misc)

            mostrarresultados(['Bias', 'Flat', 'Misc'], [len(l_bias), len(l_flat), len(l_misc)], titulo=noche,
                              contador=i, valor_max=len(lista_noches_))

            guardar_listas_csv(dir_listas_ + noche + '/' + 'CAL.csv', cal)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'SCI.csv', sci)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'MIS.csv', mis)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'ARC.csv', arc)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'LBias.csv', l_bias)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'LFlat.csv', l_flat)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'LMisc.csv', l_misc)


#####################################################################
# Separar los Bias por tamaños. Primero definimos una función que devuelva las coordenadas
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


def leer_lista(archivo):
    with open(archivo, 'rt') as f:
        reader = csv.reader(f, delimiter=',')
        your_list = list(reader)
    return your_list


def sacar_coordenadas_2(lista, idx):
    x1 = lista[idx, 0]
    x2 = lista[idx, 1]
    y1 = lista[idx, 2]
    y2 = lista[idx, 3]
    tupla_salida = (x1, x2, y1, y2)
    return tupla_salida


def add_to_file(file, what):
    open(file, 'a+').write(what)


################################################################
# Zonas del CCD mas usadas
def usoccd():
    zona_ccd = np.zeros((2048, 2048), dtype=int)
    my_data = genfromtxt('BiasSecciones.csv', delimiter=';', dtype=str)
    coorden = np.zeros((my_data.shape[0], 5), dtype=int)
    for k in range(my_data.shape[0]):
        x1, x2, y1, y2 = sacar_coordenadas_ccd(my_data[k, 0])
        coorden[k, 0] = x1
        coorden[k, 1] = x2
        coorden[k, 2] = y1
        coorden[k, 3] = y2
        coorden[k, 4] = int(my_data[k, 1])
    for seccion in range(coorden.shape[0]):
        for y in range(coorden[seccion, 2], coorden[seccion, 3]):
            for x in range(coorden[seccion, 0], coorden[seccion, 1]):
                zona_ccd[y, x] += coorden[seccion, 4]

    plt.figure()
    plt.imshow(zona_ccd, cmap='gray')
    plt.colorbar()
    plt.show()


def slicing_data(slicing_push, size_da, size_mb):
    if slicing_push == "NW":
        s1 = size_da[0] - size_mb[0]
        s2 = size_da[0]
        s3 = 0
        s4 = size_mb[1]
    elif slicing_push == "SE":
        s1 = 0
        s2 = size_mb[0]
        s3 = size_da[1] - size_mb[1]
        s4 = size_da[1]
    elif slicing_push == "NE":
        s1 = size_da[0] - size_mb[0]
        s2 = size_da[0]
        s3 = 0
        s4 = size_mb[1]
    else:
        s1 = 0
        s2 = size_mb[0]
        s3 = 0
        s4 = size_mb[1]
    return s1, s2, s3, s4


def juntar_imagenes(noche, secciones_unicas_, coordenadas_secciones_, secciones_count_, indice_seccion_,
                    bin_secciones_, dir_bias_, dir_datos_, lista_bias_, verbose=0, interactive=False, recortar=False):
    for seccion in range(len(secciones_unicas_)):
        print('seccion: ' + str(seccion))
        coordenadas_dibujo = sacar_coordenadas_2(coordenadas_secciones_, seccion)
        x1, x2, y1, y2 = deshacer_tupla_coord(coordenadas_dibujo)

        # Sacar el Binning
        crpix1 = bin_secciones_[seccion, 1]
        crpix2 = bin_secciones_[seccion, 0]

        mostrarresultados(['N', 'Crpix2', 'Crpix1', 'A', 'B'],
                          [len(indice_seccion_[indice_seccion_ == seccion]), crpix2, crpix1,
                           int((coordenadas_dibujo[3] - coordenadas_dibujo[2] + 1) / crpix1),
                           int((coordenadas_dibujo[1] - coordenadas_dibujo[0] + 1) / crpix2)])

        master_biases = np.zeros((secciones_count_[seccion],
                                  int((coordenadas_dibujo[3] - coordenadas_dibujo[2] + 1) / crpix1),
                                  int((coordenadas_dibujo[1] - coordenadas_dibujo[0] + 1) / crpix2)))

        indice0 = 0
        slicing_push = False
        for imagen in range(len(lista_bias_)):
            if indice_seccion_[imagen] == seccion:
                image_file = dir_datos_ + noche + '/' + lista_bias_[imagen]
                image_data = fits.getdata(image_file, ext=0)
                if image_data[:, :].shape == master_biases[indice0, :, :].shape:
                    master_biases[indice0, :, :] = image_data[:, :]                     # Juntar
                else:
                    size_mb = master_biases[indice0, :, :].shape
                    size_da = image_data[:, :].shape
                    if recortar:
                        if not slicing_push:
                            print("Sizes incompatible:")
                            print("Data size: " + str(size_da))
                            print("Master Bias size: " + str(size_mb) + "\n")
                            slicing_push = (input("Slicing fits. "
                                                  "Select side towards to push (SW), SE, NE, NW: ") or "SW")
                        s1, s2, s3, s4 = slicing_data(slicing_push, size_da, size_mb)
                        master_biases[indice0, :, :] = image_data[s1:s2, s3:s4]
                    else:
                        print("Sizes incompatible:")
                        print("Data size: " + str(size_da))
                        print("Master Bias size: " + str(size_mb) + "\n")
                        print("Skipping current Master Bias.")
                        print("Consider using slicing with '--recortar'. ")
                        input("Press Enter to continue")
                        break
                indice0 += 1

        master_bias_colapsado = np.median(master_biases, axis=0)
        nombre_archivo = noche + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}.fits".format(x1, x2, y1, y2)
        masterbias_final = fits.PrimaryHDU(master_bias_colapsado)
        masterbias_header = fits.open(dir_datos_ + noche + '/' + lista_bias_[0])[0].header

        if masterbias_header['BLANK']:
            del masterbias_header['BLANK']
        masterbias_final.header = masterbias_header

        if nombre_archivo in os.listdir(dir_bias_):
            os.remove(dir_bias_ + nombre_archivo)

        masterbias_final.writeto(dir_bias_ + nombre_archivo)

        if verbose >= 1:
            coord_lim = ImP.limites_imagen(*coordenadas_dibujo)
            ImP.imgdibujar(master_bias_colapsado, *coordenadas_dibujo, *coord_lim, verbose_=1)

        if interactive:
            input("Press Enter to continue...")


def realizar_master_biases(lista_noches, dir_listas, dir_datos, dir_bias, verbose, interactive, recortar):
    i_noche = 0
    for noche in lista_noches:
        i_noche += 1
        print('=== NOCHE ' + noche + ' - (' + str(i_noche) + '/' + str(len(lista_noches)) + ') ===')
        secciones = []
        lista_bias = leer_lista(dir_listas + noche + '/' + 'LBias.csv')
        lista_bias = [item for sublist in lista_bias for item in sublist]  # Limpiamos la lista para poder usarla

        for imagen in lista_bias:
            secciones.append(fits.open(dir_datos + noche + '/' + imagen)[0].header['CCDSEC'])
        secciones_unicas, secciones_count = np.unique(secciones, return_counts=True)  # Contamos secciones unicas

        # Variables:
        # secciones_unicas: lista de STR con las diferentes configuraciones de CCD que se usan
        # secciones_count: cuantas veces aparecen estas configuraciones
        # secciones: lista completa de STR de las secciones. se usa de apoyo. Se puede borrar despues
        # indice_seccion: INT con cual de las secciones pertenecen las calibraciones
        # coordenadas_secciones: coordenadas de las dierentes secciones
        # size-secciones: tamanyo de las imagenes en cada seccion

        # print('secciones_unicas')
        # print(secciones_unicas)
        # print('secciones_count')
        # print(secciones_count)
        # print('%%%%%%%%%%%%%%%%')

        bin_secciones = -1 * np.ones((len(secciones_unicas), 2), dtype=int)
        indice_seccion = np.zeros(len(lista_bias), dtype=int)
        for i in range(len(lista_bias)):
            for j in range(len(secciones_unicas)):
                if secciones[i] == secciones_unicas[j]:
                    indice_seccion[i] = j
                    bin_secciones[j, 0] = int(fits.open(dir_datos + noche + '/' + lista_bias[i])[0].header['crpix2'])
                    bin_secciones[j, 1] = int(fits.open(dir_datos + noche + '/' + lista_bias[i])[0].header['crpix1'])

        coordenadas_secciones = np.zeros((len(secciones_unicas), 4), dtype=int)
        for i in range(len(secciones_unicas)):
            coordenadas_unicas = sacar_coordenadas_ccd(secciones_unicas[i])
            coordenadas_secciones[i, :] = [*coordenadas_unicas]

        # for k in range(len(secciones_unicas)):
        #     add_to_file('BiasSecciones.csv', secciones_unicas[k] + ';' + str(secciones_count[k]) + '\n')

        juntar_imagenes(noche, secciones_unicas, coordenadas_secciones, secciones_count, indice_seccion, bin_secciones,
                        dir_bias, dir_datos, lista_bias, verbose=verbose, interactive=interactive, recortar=recortar)


def main():
    # ---------------Valores por defecto-------------------------------------------
    default_dir_datos = 'CAFOS2017/'
    default_dir_bias = 'Biases/'
    default_dir_listas = 'Listas/'
    desc_bias = ['bias', 'Bias', 'BIAS']
    desc_flats = ['flats', 'FLATS', 'Flats', 'Flat', 'flat', 'Skyflat', 'SDkyflat']
    desc_misc = ['arc']
    # -----------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description="Bias and Flat calibration of CAFOS images")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-db", "--dir_bias", default=default_dir_bias, type=str, help='Bias Directory')
    parser.add_argument("-dd", "--dir_datos", default=default_dir_datos, type=str, help='Data Directory')
    parser.add_argument("-dl", "--dir_listas", default=default_dir_listas, type=str, help='Lists Directory')
    parser.add_argument('--cmap', type=str, help="Colormap", default='hot')
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("--recortar", action="store_true")
    args = parser.parse_args()

    lista_noches = os.listdir(args.dir_datos)

    tiempo_inicio_listas = time.time()
    crear_listas_cal_y_sci(lista_noches, args.dir_listas, args.dir_datos, desc_bias, desc_flats, desc_misc)
    tiempo_medio = time.time()
    realizar_master_biases(lista_noches, args.dir_listas, args.dir_datos, args.dir_bias,
                           args.verbose, args.interactive, args.recortar)
    tiempo_final = time.time()

    mostrarresultados(['Tiempo listas', 'Tiempo master bias'],
                      [round(tiempo_medio-tiempo_inicio_listas, 2), round(tiempo_final - tiempo_medio, 2)],
                      titulo='Tiempo Ejecucion')


if __name__ == "__main__":

    main()
    # usoccd()

################################################################
# Dibujar uno de ellos

# print('DIBUJAR')
# image_file = mypath + lista_cal[0]

# hdul = fits.open(image_file)
# hdr = hdul[0].header

# #print(hdr[0:5])
# print(hdr['NAXIS'])
# print(hdr['NAXIS1'])
# print(hdr['NAXIS2'])
# print(hdr['OBJECT'])
# print(hdr['BIASSEC'])
# print(hdr['DATASEC'])
# print(hdr['CCDSEC'])
# print('--...---...---')
# # for imagen in lista_cal:
# #     print(fits.open(mypath + imagen)[0].header['OBJECT'])


# image_data = fits.getdata(image_file, ext=0)
# plt.imshow(image_data)
# plt.show()


# Fijar nuevas variables en el header:
# hdr = hdul[0].header
# >>> hdr['targname'] = ('NGC121-a', 'the observation target')
# >>> hdr['targname']
# 'NGC121-a'
# >>> hdr.comments['targname']
# 'the observation target'

# >>> hdr.set('observer', 'Edwin Hubble')

# History y Comments anyade un valor nuevo
# >>> hdr['history'] = 'I updated this file 2/26/09'
# >>> hdr['comment'] = 'Edwin Hubble really knew his stuff'
# >>> hdr['comment'] = 'I like using HST observations'
# >>> hdr['history']
# I updated this file 2/26/09
# >>> hdr['comment']
# Edwin Hubble really knew his stuff
# I like using HST observations

# Cambiar el valor es muy sencillo, sigue como numpy
# >>> data[30:40, 10:20] = data[1, 4] = 999
# >>> photflam = hdul[1].header['photflam']
# >>> exptime = hdr['exptime']
# >>> data = data * photflam / exptime
# >>> hdul.close()


###################################################################
# Dibujar una imagen
# #image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
# image_file = 'caf-20170224-20:20:04-cal-krek.fits'


# fits.info(image_file)

# hdu = fits.PrimaryHDU()
# print('%&%&%&%&%&%&%&%&%')
# #print(hdu.header)

# hdul = fits.open(image_file)
# hdr = hdul[0].header
# #print(hdr[34])

# #print(hdr)
# print(hdr['OBJECT'])
# print(hdr['BIASSEC'])
# print(hdr['DATASEC'])
# print(hdr['CCDSEC'])


# image_data = fits.getdata(image_file, ext=0)

# print(image_data.shape)

# plt.figure()
# plt.imshow(image_data, cmap='gray')
# plt.colorbar()
# plt.show()


# '#!/bin/bash'
# 'python 3_inicio.py nombre archivo\
# '--parametro1 valor1 \
# '--parametro2 valor2 \
# chmod +x nombre_script.sh

# https://pyformat.info
# "str{0:04d}".format(210)





