from astropy.io import fits
from Salida_limpia import mostrarresultados, stdrobusta
import numpy as np
import IMGPlot as ImP
# import pandas as pd
import matplotlib.pyplot as plt
import os
import csv
from numpy import genfromtxt
import argparse
import time
import warnings


def deshacer_tupla_coord(tupla):
    return tupla[0], tupla[1], tupla[2], tupla[3]


def guardar_listas_csv(csvfile, res):
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in res:
            writer.writerow([val])


def conseguir_listas_archivos(path_):
    listaarchivos = []
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            listaarchivos.append(os.path.join(path_, file))
    return listaarchivos


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


def obtener_bias(dir_bias_, noche, lista_noches, lista_bias, observador, x1, x2, y1, y2, b1, b2):

    bias_asociado_nombre = dir_bias_ + noche + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2,
                                                                                                               y1, y2,
                                                                                                               b1, b2)
    otro_observador = False
    hace_falta_cambio = True
    indice_bias = None
    if bias_asociado_nombre in lista_bias:
        hace_falta_cambio = False
        bias_asociado = fits.getdata(bias_asociado_nombre, ext=0)
    else:
        exitos = []
        observadores = []
        posicion = lista_noches.index(noche)
        busqueda_max = 10  # Cantidad maxima de dias a mirar a izquiera y derecha
        for i in range(1, busqueda_max):
            for mult in [-1, 1]:
                indice = i * mult
                print(indice)
                noche = lista_noches[posicion + indice]
                bias_asociado_nuevo = dir_bias_ + noche + \
                                      "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2,
                                                                                                      y1, y2,
                                                                                                      b1, b2)
                if bias_asociado_nuevo in lista_bias:
                    exitos.append(indice)
                    observadores.append(fits.open(bias_asociado_nuevo)[0].header['OBSERVER'])

        # Ahora tenemos una lista con posibles candidatos
        if all([exitos, observadores]):
            for j in range(len(observadores)):
                if observadores[j] == observador:
                    indice_bias = exitos[j]
                    break

            # Si no ha encontrado nada, cogemos el primero
            if not indice_bias:
                indice_bias = exitos[0]
                otro_observador = True
        else:
            print('No hay biases cerca')
            input('Enter para continuar...')
        if indice_bias:
            print(exitos)
            print(observadores)
            print(observador)
            print(otro_observador)
            print('indice_bias', indice_bias)
            noche = lista_noches[posicion + indice_bias]
            bias_asociado_n = dir_bias_ + noche +\
                              "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2,
                                                                                              y1, y2,
                                                                                              b1, b2)
            bias_asociado = fits.getdata(bias_asociado_n, ext=0)

    if hace_falta_cambio:
        print(exitos)
        print(observadores)
    else:
        print('existe y sin problemas')

    return bias_asociado, noche, otro_observador


def juntar_imagenes(noche, secciones_unicas_, coordenadas_secciones_, secciones_count_, indice_seccion_,
                    bin_secciones_, dir_bias_, dir_datos_, dir_flats_, lista_flats_, lista_noches, lista_bias,
                    verbose=0, interactive=False, recortar=False):

    for seccion in range(len(secciones_unicas_)):
        print('seccion: ' + str(seccion))
        coordenadas_dibujo = sacar_coordenadas_2(coordenadas_secciones_, seccion)
        x1, x2, y1, y2 = deshacer_tupla_coord(coordenadas_dibujo)

        # Creamos una lista con todos aquellos que coinciden el indice
        lista_coincide = []
        for imagen in range(len(lista_flats_)):
            if indice_seccion_[imagen] == seccion:
                lista_coincide.append(lista_flats_[imagen])

        filtros, filtros_unicos, filtros_count, indice_filtro, _ = crear_lista_unicos(dir_datos_, noche,
                                                                                      lista_coincide,
                                                                                      cabecera='INSFLID',
                                                                                      binning=False)

        numero_filtro = np.zeros(len(filtros_unicos), dtype=int)
        p = 0
        for i in filtros_unicos:
            numero_filtro[p] = int(i[-2:].strip())
            p += 1

        for filtro in range(len(filtros_unicos)):
            da_problemas = False

            lista_actual = []
            binning2 = []
            binning1 = []
            for i in range(len(lista_coincide)):
                if indice_filtro[i] == filtro:
                    lista_actual.append(lista_coincide[i])
                    binning2.append(int(fits.open(dir_datos_ + noche + '/' + lista_coincide[i])[0].header['crpix2']))
                    binning1.append(int(fits.open(dir_datos_ + noche + '/' + lista_coincide[i])[0].header['crpix1']))

            bin2_unique, bin2_count = np.unique(binning2, return_counts=True)
            bin1_unique, bin1_count = np.unique(binning1, return_counts=True)
            siempre_coincide_binning = True
            for i in range(len(lista_actual)):
                if not binning2 == binning1:
                    siempre_coincide_binning = False
                    break
            solo_un_caso_de_binning = all([len(bin1_unique) == 1, len(bin2_unique) == 1])
            print(len(lista_actual), len(lista_coincide))
            print(all([len(bin1_unique) == 0, len(bin2_unique) == 0]))
            print('__________')
            print(solo_un_caso_de_binning)
            # if solo_un_caso_de_binning:  ############################################################################

            if siempre_coincide_binning:
                print(bin1_unique, bin1_count, bin2_unique, bin2_count)
                for binning in bin1_unique:
                    lista_actual_bin = []

                    for i in range(len(lista_actual)):
                        if binning1[i] == binning:
                            lista_actual_bin.append(lista_actual[i])

                    # Ahora tenemos una lista con binnings individuales
                    print('binning', binning)
                    print(len(lista_actual_bin), lista_actual_bin)

                    # Como tienen un solo bin no nos preocupamos y generamos las imagenes
                    crpix2 = binning
                    crpix1 = binning
                    naxis1_expected = int((coordenadas_dibujo[3] - coordenadas_dibujo[2] + 1) / crpix1)
                    naxis2_expected = int((coordenadas_dibujo[1] - coordenadas_dibujo[0] + 1) / crpix2)
                    if (coordenadas_dibujo[3] - coordenadas_dibujo[2] + 1) % crpix1 != 0:
                        naxis1_expected += 1
                    if (coordenadas_dibujo[1] - coordenadas_dibujo[0] + 1) % crpix2 != 0:
                        naxis2_expected += 1

                    master_flats = np.zeros((filtros_count[filtro], naxis1_expected, naxis2_expected), dtype=float)

                    mostrarresultados(['N', 'Crpix2', 'Crpix1', 'A', 'B'],
                                      [len(lista_actual_bin), crpix2, crpix1,
                                       naxis1_expected, naxis2_expected],
                                      titulo='Filtro ' + str(filtros_count[filtro]))
                    indice0 = 0
                    slicing_push = False
                    print(len(lista_coincide))
                    print(seccion)
                    print(indice_seccion_)
                    print(filtro)
                    print(indice_filtro)
                    for imagen in range(len(lista_actual_bin)):
                        image_file = dir_datos_ + noche + '/' + lista_actual_bin[imagen]
                        image_data = fits.getdata(image_file, ext=0)
                        bin_prueba2 = int(fits.open(image_file)[0].header['crpix2'])
                        bin_prueba1 = int(fits.open(image_file)[0].header['crpix1'])

                        print(imagen, lista_actual_bin[imagen], lista_actual_bin.index(lista_actual_bin[imagen]), '|',
                              indice_seccion_[lista_flats_.index(lista_actual_bin[imagen])], indice_filtro[imagen],
                              indice_seccion_[lista_flats_.index(lista_actual_bin[imagen])] == seccion,
                              indice_filtro[imagen] == filtro, '|', seccion, filtro, '|', image_data.shape,
                              (crpix1, crpix2), (bin_prueba1, bin_prueba2), '|',
                              bin1_count, bin1_unique, bin2_count, bin2_unique, '|', len(lista_actual_bin),
                              int(fits.open(dir_datos_ + noche + '/' + lista_actual_bin[imagen])[0].header['crpix2']),
                              int(fits.open(dir_datos_ + noche + '/' + lista_actual_bin[imagen])[0].header['crpix1']))

                        if image_data[:, :].shape == master_flats[indice0, :, :].shape:
                            master_flats[indice0, :, :] = image_data[:, :]                  # Juntar
                            if indice0 == 0:
                                cabecera = fits.open(image_file)[0].header
                        else:
                            print('hay un problema con el tamanyo de las imagenes')
                            input('Pausado')
                        # else:
                        #     size_mb = master_flats[indice0, :, :].shape
                        #     size_da = image_data[:, :].shape
                        #     da_problemas = True
                        #     if recortar:
                        #         if not slicing_push:
                        #             warnings.warn("Sizes Incompatible!")
                        #             print("Sizes incompatible:")
                        #             print(image_file)
                        #             print(image_file in lista_coincide, image_file in lista_flats_)
                        #             print("Data size: " + str(size_da))
                        #             print("Master Bias size: " + str(size_mb) + "\n")
                        #             slicing_push = (input("Slicing fits. "
                        #                                   "Select side towards to push (SW), SE, NE, NW: ") or "SW")
                        #         s1, s2, s3, s4 = slicing_data(slicing_push, size_da, size_mb)
                        #         master_flats[indice0, :, :] = image_data[s1:s2, s3:s4]
                        #     else:
                        #         warnings.warn("Sizes Incompatible!")
                        #         print("Sizes incompatible:")
                        #         print("Data size: " + str(size_da))
                        #         print("Master Bias size: " + str(size_mb) + "\n")
                        #         print("Skipping current Master Bias.")
                        #         print("Consider using slicing with '--recortar'. ")
                        #         input("Press Enter to continue...")
                        #         break
                        indice0 += 1

                    # Restamos Bias
                    # Aqui es donde hay que meter la función que busca cual es el bias que va a usar, porque hay que tener
                    # cuidado con que exista ese bias. Si no hay ninguno de la misma noche, buscara en los dias de alrededor
                    # para las imagenes que ha tomado el mismo investigador. Si no estan esos disponibles, buscara de otro

                    observador = cabecera['OBSERVER']
                    bias_asociado, noche_usada, otro_observador = obtener_bias(dir_bias_, noche, lista_noches,
                                                                               lista_bias, observador,
                                                                               x1, x2, y1, y2,
                                                                               crpix1, crpix2)
                    if noche_usada != noche:
                        print('Se ha usado una noche diferente: ', noche_usada)
                    if otro_observador:
                        print('La imagen es de otro observador')

                    for i in range(master_flats.shape[0]):
                        master_flats[i, :, :] = master_flats[i, :, :] - bias_asociado

                    # Normalizamos
                    valor_medio = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):
                        valor_medio[i] = np.mean(master_flats[i, :, :], dtype=float)
                        master_flats[i, :, :] = np.true_divide(master_flats[i, :, :], valor_medio[i])

                    valor_medio2 = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):
                        valor_medio2[i] = np.mean(master_flats[i, :, :], dtype=float)

                    # Colapsamos
                    master_flats_colapsado = np.median(master_flats, axis=0)
                    # plt.imshow(master_flats_colapsado)
                    # plt.show()
                    # ImP.imgdibujar(master_flats_colapsado)

                    nombre_archivo = noche + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:02d}.fits".format(x1, x2, y1, y2, crpix1, crpix2,
                                                                                                     numero_filtro[filtro])

                    masterflats_header = cabecera.copy()
                    if masterflats_header['BLANK']:
                        del masterflats_header['BLANK']

                    masterflats_final = fits.PrimaryHDU(master_flats_colapsado.astype(np.float), masterflats_header)

                    masterflats_final.writeto(dir_flats_ + nombre_archivo, overwrite=True)

            else:
                if siempre_coincide_binning:
                    print(bin1_unique, bin1_count, bin2_unique, bin2_count)
                    for binning in bin1_unique:
                        lista_actual_bin = []
                        for i in range(len(lista_actual)):
                            if binning1[i] == binning:
                                lista_actual_bin.append(lista_actual[i])

                        # Ahora tenemos una lista con binnings individuales
                        print('binning', binning)
                        print(len(lista_actual_bin), lista_actual_bin)
                        print('coinciden entre si')
                else:
                    print('no coincide')
                input('Pausa')

            if verbose >= 1:
                coord_lim = ImP.limites_imagen(*coordenadas_dibujo)
                ImP.imgdibujar(master_flats_colapsado, *coordenadas_dibujo, *coord_lim, verbose_=1)

            if interactive:
                # plt.show()
                input("Press Enter to continue...")


def crear_lista_unicos(dir_datos, noche, lista_cosas, cabecera, binning=False):
    lista = []

    for imagen in lista_cosas:
        lista.append(fits.open(dir_datos + noche + '/' + imagen)[0].header[cabecera])
    lista_unicas, lista_count = np.unique(lista, return_counts=True)  # Contamos secciones unicas

    bin_secciones = -1 * np.ones((len(lista_unicas), 2), dtype=int)
    indice_cosas = np.zeros(len(lista_cosas), dtype=int)
    for i in range(len(lista_cosas)):
        for j in range(len(lista_unicas)):
            if lista[i] == lista_unicas[j]:
                indice_cosas[i] = j
                if binning:
                    bin_secciones[j, 0] = int(fits.open(dir_datos + noche + '/' + lista_cosas[i])[0].header['crpix2'])
                    bin_secciones[j, 1] = int(fits.open(dir_datos + noche + '/' + lista_cosas[i])[0].header['crpix1'])
                break

    return lista, lista_unicas, lista_count, indice_cosas, bin_secciones


def realizar_master_flats(lista_noches, lista_bias, dir_listas, dir_datos, dir_bias, dir_flats,
                          verbose, interactive, recortar):
    i_noche = 0
    for noche in lista_noches[:]:
        i_noche += 1
        print('=== NOCHE ' + noche + ' - (' + str(i_noche) + '/' + str(len(lista_noches)) + ') ===')

        lista_flats = leer_lista(dir_listas + noche + '/' + 'LFlat.csv')
        lista_flats = [item for sublist in lista_flats for item in sublist]  # Limpiamos la lista para poder usarla

        secciones, secciones_unicas, secciones_count, indice_seccion, bin_secciones = crear_lista_unicos(
            dir_datos, noche, lista_flats, cabecera='CCDSEC', binning=True
        )

        # Variables:
        # secciones_unicas: lista de STR con las diferentes configuraciones de CCD que se usan
        # secciones_count: cuantas veces aparecen estas configuraciones
        # secciones: lista completa de STR de las secciones. se usa de apoyo. Se puede borrar despues
        # indice_seccion: INT con cual de las secciones pertenecen las calibraciones
        # coordenadas_secciones: coordenadas de las dierentes secciones
        # size-secciones: tamanyo de las imagenes en cada seccion

        coordenadas_secciones = np.zeros((len(secciones_unicas), 4), dtype=int)
        for i in range(len(secciones_unicas)):
            coordenadas_unicas = sacar_coordenadas_ccd(secciones_unicas[i])
            coordenadas_secciones[i, :] = [*coordenadas_unicas]

        juntar_imagenes(noche, secciones_unicas, coordenadas_secciones, secciones_count, indice_seccion, bin_secciones,
                        dir_bias, dir_datos, dir_flats, lista_flats, lista_noches, lista_bias,
                        verbose=verbose, interactive=interactive, recortar=recortar)


def main():
    # ---------------Valores por defecto-------------------------------------------
    default_dir_datos = 'CAFOS2017/'
    default_dir_bias = 'Biases/'
    default_dir_listas = 'Listas/'
    default_dir_flats = 'Flats/'
    # -----------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description="Bias and Flat calibration of CAFOS images")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-db", "--dir_bias", default=default_dir_bias, type=str, help='Bias Directory')
    parser.add_argument("-df", "--dir_flats", default=default_dir_flats, type=str, help='Flats Directory')
    parser.add_argument("-dd", "--dir_datos", default=default_dir_datos, type=str, help='Data Directory')
    parser.add_argument("-dl", "--dir_listas", default=default_dir_listas, type=str, help='Lists Directory')
    parser.add_argument('--cmap', type=str, help="Colormap", default='hot')
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("--recortar", action="store_true")
    args = parser.parse_args()

    lista_noches = os.listdir(args.dir_datos)
    tiempo_inicio_listas = time.time()

    lista_bias = conseguir_listas_archivos(args.dir_bias)
    realizar_master_flats(lista_noches, lista_bias, args.dir_listas, args.dir_datos, args.dir_bias, args.dir_flats,
                          args.verbose, args.interactive, args.recortar)
    tiempo_final = time.time()

    mostrarresultados(['Tiempo Master Bias', 'Cuantos Flats'],
                      [round(tiempo_final - tiempo_inicio_listas, 2),
                       len(os.listdir(args.dir_flats))],
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
# print(hdr['CCDSEC'])
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


# usar numpy.ma.median para la mascara
# Crear una mascara que calcule el radio con mask=x*x+y*y<=r*r
# Dividir por median
#
# Anadir B{0:02d}-{1:02d}.format(binning1,binning2)
# Hay que comprobar a ver si hay imagenes de ciencia con multiples binning tambien o no
#
# Crear un diccionario de filtros que se van a usar para el nombre. Siempre se usaran los mimso numeros para el mismo
# filtro aunque sean diferentes noches o imagenes. Ese diccionario se ampliara a lo largo del tiempo
