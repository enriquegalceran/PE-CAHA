# -*- coding: utf-8 -*-

"""
    Proyecto de Unidad de Limpieza y Tratamiento de Resultados Observacionales Nativo (Proyecto U.L.T.R.O.N.)

    Recoge los resultados de observaciones del instrumento CAFOS y los reduce correctamente.
    Desarrollado por Enrique Galceran García
"""

from astropy.io import fits
from astropy.time import Time
from Salida_limpia import mostrarresultados, stdrobusta
import numpy as np
import pandas as pd
import numpy.ma as ma
import IMGPlot as ImP
import matplotlib.pyplot as plt
import os
import csv
import argparse
import time
import datetime
import warnings
import json


def deshacer_tupla_coord(tupla):
    return tupla[0], tupla[1], tupla[2], tupla[3]


def str2datetime(string):
    salida = datetime.datetime.strptime(string, '%Y-%m-%dT%H:%M:%S')
    return salida


def guardar_listas_csv(csvfile, res):
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in res:
            writer.writerow([val])


def mostrar_diccionario(nombre_diccionario):
    print("Mostramos el diccionario: ")
    for x, y in nombre_diccionario.items():
        print(x, y)


def guardar_json(variable, filename):
    json_f = json.dumps(variable, indent=2, sort_keys=True)
    f = open(filename, "w")
    f.write(json_f)
    f.close()


def cargar_json(filename='Dic_filtro.json'):
    with open(filename) as json_file:
        data = json.load(json_file)
    return data


def conseguir_listas_archivos(path_):
    listaarchivos = []
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            listaarchivos.append(os.path.join(path_, file))
    return listaarchivos


def leer_diccionario(nombre, diccionario_, filename='Dic_filtro.json'):
    """
    Busca en el diccionario el número del filtro buscado. Si no existe, crea uno nuevo

    :param nombre:
    :param diccionario_:
    :param filename:
    :return:
    """
    if nombre in diccionario_:
        return diccionario_[nombre]
    else:
        len_dic = len(diccionario_)
        diccionario_[nombre] = len_dic
        print('Nueva entrada en el diccionario: Filtro {0} - Indice {1:03d}'.format(nombre, diccionario_[nombre]))
        guardar_json(diccionario_, filename)
        return diccionario_[nombre]


def sacar_coordenadas_ccd(imagen_, mypath_=False):
    """
    Dado un string de texto en formato de CAFOS nos devuelve los valores x1, x2, y1 e y2.

    :param imagen_:
    :param mypath_:
    :return:
    """
    if mypath_:
        str_ccdsec = fits.open(mypath_ + imagen_)[0].header['CCDSEC']
        longitud = len(str_ccdsec)
    else:
        str_ccdsec = imagen_
        longitud = len(str_ccdsec)+1

    coma = None
    for x in range(1, longitud):
        if str_ccdsec[x] == ',':
            coma = x
            break
    if coma is None:
        raise ValueError('coma not defined!')

    puntos = None
    for x in range(coma + 1, longitud):
        if str_ccdsec[x] == ':':
            puntos = x
            break
    if puntos is None:
        raise ValueError('puntos not defined!')

    coma2 = None
    for x in range(puntos + 1, longitud):
        if str_ccdsec[x] == ',':
            coma2 = x
            break
    if coma2 is None:
        raise ValueError('coma2 not defined!')

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
    your_list = [item for sublist in your_list for item in sublist]
    return your_list


def sacar_coordenadas_2(lista, idx):
    x1 = lista[idx, 0]
    x2 = lista[idx, 1]
    y1 = lista[idx, 2]
    y2 = lista[idx, 3]
    tupla_salida = (x1, x2, y1, y2)
    return tupla_salida


def obtener_bias(dir_bias_, noche, lista_noches, lista_bias, x1, x2, y1, y2, b1, b2,
                 busqueda_max=10, relleno=680, verbose=False):
    """
    Busca en la lista de bias disponibles si exite un archivo de bias para la misma noche en la que se hizo la foto.
    Si exite, usará esa directamente. Si no existe, busca una imagen de bias diferente del mismo autor de noches
    cercanas. Si no hay, genera un bias plano.

    :param dir_bias_:
    :param noche:
    :param lista_noches:
    :param lista_bias:
    :param x1:
    :param x2:
    :param y1:
    :param y2:
    :param b1:
    :param b2:
    :param busqueda_max:
    :param relleno:
    :param verbose:
    :return:
    """

    bias_asociado_nombre = dir_bias_ + noche +\
        "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2, y1, y2, b1, b2)
    existe = None
    cogido_otro = False
    if bias_asociado_nombre in lista_bias:
        existe = True
    else:
        posicion = lista_noches.index(noche)
        for i in range(1, busqueda_max):
            for mult in [-1, 1]:
                indice = i * mult
                pos_nueva = posicion + indice
                if pos_nueva >= len(lista_noches):
                    # Se ha llegado al final del año, y no busca para el año siguiente
                    break
                noche = lista_noches[pos_nueva]
                bias_asociado_nuevo = dir_bias_ + noche + \
                    "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits"\
                    .format(x1, x2, y1, y2, b1, b2)

                if bias_asociado_nuevo in lista_bias:
                    bias_asociado_nombre = bias_asociado_nuevo
                    existe = True
                    cogido_otro = True

    if existe:
        if verbose:
            if cogido_otro:
                print('Se ha cogido un bias de un dia diferente. El bias que se ha cogido es:')
                print(bias_asociado_nombre)
            else:
                print('Existe el bias asociado')
        bias_asociado = fits.getdata(bias_asociado_nombre, ext=0)
    else:
        if verbose:
            print('No hay bias cercanos validos. Se ha generado uno.')
        naxis1_expected, naxis2_expected = sacar_naxis((x1, x2, y1, y2), b2, b1)

        bias_asociado = np.full((naxis1_expected, naxis2_expected), relleno, dtype=float)

    return existe, bias_asociado_nombre, bias_asociado


def obtener_flats(dir_flats_, noche_, lista_noches, lista_flats, x1, x2, y1, y2, b1, b2,
                  id_filtro_, busqueda_max=10, relleno=1, verbose=False):
    """
        Busca en la lista de bias disponibles si exite un archivo de bias para la misma noche en la que se hizo la foto.
        Si exite, usará esa directamente. Si no existe, busca una imagen de bias diferente del mismo autor de noches
        cercanas. Si no hay, genera un bias plano.

    :param dir_flats_:
    :param noche_:
    :param lista_noches:
    :param lista_flats:
    :param x1:
    :param x2:
    :param y1:
    :param y2:
    :param b1:
    :param b2:
    :param id_filtro_:
    :param busqueda_max:
    :param relleno:
    :param verbose:
    :return:
    """

    flat_asociado_nombre = dir_flats_ + noche_ + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits" \
        .format(x1, x2, y1, y2, b1, b2, id_filtro_)
    existe = False
    cogido_otro = False

    if flat_asociado_nombre in lista_flats:
        existe = True
    else:
        posicion = lista_noches.index(noche_)
        for i in range(1, busqueda_max):
            if existe:
                break
            for mult in [-1, 1]:
                indice = i * mult
                pos_nueva = posicion + indice
                if pos_nueva >= len(lista_noches):
                    # Se ha llegado al final del año, y no busca para el año siguiente
                    break
                noche = lista_noches[pos_nueva]
                flat_asociado_nuevo = dir_flats_ + noche + \
                    "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits" \
                    .format(x1, x2, y1, y2, b1, b2, id_filtro_)

                if flat_asociado_nuevo in lista_flats:
                    flat_asociado_nombre = flat_asociado_nuevo
                    existe = True
                    cogido_otro = True
                    break

    if existe:
        if verbose:
            if cogido_otro:
                print('Se ha cogido un flat de un dia diferente. El flat que se ha cogido es:')
                print(flat_asociado_nombre)
            else:
                print('Existe el flat asociado')
        flat_asociado = fits.getdata(flat_asociado_nombre, ext=0)
    else:
        if verbose:
            print('No hay flats cercanos validos. Se ha generado uno.')

        # flat_asociado = None
        # if flat_asociado is None:
        #     raise ValueError('No hay definido un valor por defecto para el flat')

        naxis1_expected, naxis2_expected = sacar_naxis((x1, x2, y1, y2), b2, b1)

        flat_asociado = np.full((naxis1_expected, naxis2_expected), relleno, dtype=float)

    return existe, flat_asociado_nombre, flat_asociado


def create_circular_mask(h, w, center=None, radius=None):

    """
    Genera la mascara circular en función de la altura y anchura de la imagen además del centro y radio.

    :param h:
    :param w:
    :param center:
    :param radius:
    :return:
    """

    if center is None:  # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None:  # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    yc, xc = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((xc - center[0])**2 + (yc-center[1])**2)

    mask = dist_from_center <= radius
    return mask


def crear_lista_unicos(dir_datos, noche, lista_cosas, cabecera, binning=False, nombre_filtro=False):

    """
    Dado una lista de parámetros, genera una serie de listas:
        - Lista de archivos
        - Lista de valores únicos (sin repetidos)
        - Para cada valor único, las instancias que se repiten
        - Para cada elemento de la lista, genera un índice de a cual de ĺos valores únicos pertenece
        - Si se pide, la cantidad de secciones que tiene
        - Si se pide, el nombre de los filtros en cuestión

    :param dir_datos:
    :param noche:
    :param lista_cosas:
    :param cabecera:
    :param binning:
    :param nombre_filtro:
    :return:
    """
    lista = []

    for imagen in lista_cosas:
        lista.append(fits.open(dir_datos + noche + '/' + imagen)[0].header[cabecera])
    lista_unicas, lista_count = np.unique(lista, return_counts=True)  # Contamos secciones unicas

    bin_secciones = -1 * np.ones((len(lista_unicas), 2), dtype=int)
    indice_cosas = np.zeros(len(lista_cosas), dtype=int)
    nombres_filtros = []

    for i in range(len(lista_cosas)):
        for j in range(len(lista_unicas)):
            if lista[i] == lista_unicas[j]:
                indice_cosas[i] = j
                if binning:
                    bin_secciones[j, 0] = int(fits.open(dir_datos + noche + '/' + lista_cosas[i])[0].header['ccdbinX'])
                    bin_secciones[j, 1] = int(fits.open(dir_datos + noche + '/' + lista_cosas[i])[0].header['ccdbinY'])
                if nombre_filtro:
                    nombres_filtros.append(fits.open(dir_datos + noche + '/' + lista_cosas[i])[0].header['INSFLNAM'])
                break

    return lista, lista_unicas, lista_count, indice_cosas, bin_secciones, nombres_filtros


def imagen_mas_probable(archivo):
    image_data = fits.getdata(archivo, ext=0)
    v_median = np.median(image_data)
    v_sd = stdrobusta(image_data)
    if v_sd < 50 and v_median < 900:
        probable = 0  # Bias
    elif v_sd < 500:
        probable = 1  # Arc
    else:
        probable = 2  # Flat

    return probable


def comprobar(archivo, descriptor, descriptor2=None, descriptor3=None, verbose=False):
    if any([descriptor2, descriptor3]):
        descriptortemp = descriptor + descriptor2 + descriptor3
        coincide = True
        for texto in descriptortemp:
            if texto in fits.open(archivo)[0].header['OBJECT']:
                # Alguno coincide, no es ciencia
                coincide = False
                break
    else:
        coincide = False
        for texto in descriptor:
            if texto in fits.open(archivo)[0].header['OBJECT']:
                # Sale bien
                coincide = True
                break
    if verbose:
        print(archivo, fits.open(archivo)[0].header['OBJECT'], fits.open(archivo)[0].header['imagetyp'], coincide)

    # Quiza meterle algo de que si object=test que lo ignore

    return coincide


def listas_archivos2(path_, desc_bias, desc_flats, desc_arc, verbose=False, calysci=True):
    lista_cal = []
    lista_sci = []
    lista_misc = []
    lista_bias = []
    lista_flat = []
    lista_ciencia = []
    lista_arc = []
    lista_else = []
    listaarchivos = []
    lista_falla = []

    # Miramos todos los archivos dentro de la carpeta
    for file in os.listdir(path_):
        if file.endswith(".fits"):
            if calysci:
                if '-cal-' in file:
                    lista_cal.append(file)
                elif '-sci-' in file:
                    lista_sci.append(file)
                else:
                    lista_misc.append(file)
            listaarchivos.append(os.path.join(path_, file))

            # Separamos segun IMAGETYP
            tipo = fits.open(path_ + file)[0].header['IMAGETYP'].strip()
            if not fits.open(path_ + file)[0].header['OBJECT'] == 'Test':  # Comprobamos que no es un test
                if tipo == 'BIAS' or tipo == 'bias':
                    coincide = comprobar(path_ + file, desc_bias, verbose=verbose)
                    if coincide:
                        lista_bias.append(file)
                    else:
                        lista_falla.append(file)

                elif tipo == 'flat':
                    coincide = comprobar(path_ + file, desc_flats, verbose=verbose)
                    if coincide:
                        lista_flat.append(file)
                    else:
                        lista_falla.append(file)

                elif tipo == 'arc':
                    coincide = comprobar(path_ + file, desc_arc, verbose=verbose)
                    if coincide:
                        lista_arc.append(file)
                    else:
                        lista_falla.append(file)

                elif tipo == 'science':
                    coincide = comprobar(path_ + file, desc_bias, desc_flats, desc_arc, verbose=verbose)
                    if coincide:
                        lista_ciencia.append(file)
                    else:
                        lista_falla.append(file)

                else:
                    lista_else.append(file)

    for file in lista_falla:

        if calysci:
            if file in lista_sci:
                lista_ciencia.append(file)

        probable = imagen_mas_probable(path_ + file)

        if probable == 0:
            lista_bias.append(file)
        elif probable == 1:
            lista_arc.append(file)
        elif probable == 2:
            lista_flat.append(file)

    return lista_bias, lista_flat, lista_arc, lista_ciencia, listaarchivos, lista_falla


def crear_listas_cal_y_sci(lista_noches_, dir_listas_, dir_datos_, desc_bias, desc_flats, desc_arc, verbose, calysci):
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
            l_bias, l_flat, l_arc, l_ciencia, l_archivos, l_falla = listas_archivos2(path_,
                                                                                     desc_bias,
                                                                                     desc_flats,
                                                                                     desc_arc,
                                                                                     verbose,
                                                                                     calysci)

            mostrarresultados(['Bias', 'Flat', 'Arc', 'Ciencia', 'Falla'],
                              [len(l_bias), len(l_flat), len(l_arc), len(l_ciencia), len(l_falla)],
                              titulo=noche, contador=i, valor_max=len(lista_noches_))

            # guardar_listas_csv(dir_listas_ + noche + '/' + 'CAL.csv', cal)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'SCI.csv', l_ciencia)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'ARC.csv', l_archivos)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'LBias.csv', l_bias)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'LFlat.csv', l_flat)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'LArc.csv', l_arc)
            guardar_listas_csv(dir_listas_ + noche + '/' + 'Errores.csv', l_falla)


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


def obt_naxis(x2, x1, binn):
    naxis = int((x2 - x1 + 1) / binn)
    if (x2 - x1 + 1) % binn != 0:
        naxis += 1
    return naxis


def sacar_naxis(coordenadas_tupla, binn_1, binn_2=None):
    if binn_2 is None:
        binn_2 = binn_1

    x1, x2, y1, y2 = deshacer_tupla_coord(coordenadas_tupla)

    naxis1 = obt_naxis(y2, y1, binn_1)
    naxis2 = obt_naxis(x2, x1, binn_2)

    return naxis1, naxis2
########################################################################################################################


def juntar_imagenes_bias(noche, secciones_unicas_, coordenadas_secciones_, secciones_count_, indice_seccion_,
                         bin_secciones_, dir_bias_, dir_datos_, lista_bias_,
                         interactive=False, recortar=False, verbose_imagen=False):
    elemento_lista = None
    for seccion in range(len(secciones_unicas_)):
        print('seccion: ' + str(seccion))
        coordenadas_dibujo = sacar_coordenadas_2(coordenadas_secciones_, seccion)
        x1, x2, y1, y2 = deshacer_tupla_coord(coordenadas_dibujo)

        # Sacar el Binning
        ccdbinx = bin_secciones_[seccion, 1]
        ccdbiny = bin_secciones_[seccion, 0]

        naxis1_expected, naxis2_expected = sacar_naxis(coordenadas_dibujo, ccdbinx, ccdbiny)

        master_biases = np.zeros((secciones_count_[seccion],
                                  naxis1_expected,
                                  naxis2_expected), dtype=float)

        indice0 = 0
        slicing_push = False
        cabecera = None
        for imagen in range(len(lista_bias_)):
            if indice_seccion_[imagen] == seccion:
                image_file = dir_datos_ + noche + '/' + lista_bias_[imagen]
                image_data = fits.getdata(image_file, ext=0)
                if image_data[:, :].shape == master_biases[indice0, :, :].shape:
                    master_biases[indice0, :, :] = image_data[:, :]                     # Juntar
                    if indice0 == 0:
                        cabecera = fits.open(image_file)[0].header
                else:
                    size_mb = master_biases[indice0, :, :].shape
                    size_da = image_data[:, :].shape
                    if recortar:
                        if not slicing_push:
                            warnings.warn("Sizes Incompatible!")
                            print("Sizes incompatible:")
                            print("Data size: " + str(size_da))
                            print("Master Bias size: " + str(size_mb) + "\n")
                            slicing_push = (input("Slicing fits. "
                                                  "Select side towards to push (SW), SE, NE, NW: ") or "SW")
                        s1, s2, s3, s4 = slicing_data(slicing_push, size_da, size_mb)
                        master_biases[indice0, :, :] = image_data[s1:s2, s3:s4]
                    else:
                        warnings.warn("Sizes Incompatible!")
                        print("Sizes incompatible:")
                        print("Data size: " + str(size_da))
                        print("Master Bias size: " + str(size_mb) + "\n")
                        print("Skipping current Master Bias.")
                        print("Consider using slicing with '--recortar'. ")
                        input("Press Enter to continue...")
                        break
                indice0 += 1

        master_bias_colapsado = np.median(master_biases, axis=0)
        # print(master_bias_colapsado)
        # print(master_bias_colapsado.shape)
        # print(master_bias_colapsado)
        # plt.imshow(master_bias_colapsado)
        if verbose_imagen:
            ImP.imgdibujar(master_bias_colapsado, verbose_=1)
            plt.show()
        nombre_archivo = noche + "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}.fits".format(x1, x2, y1, y2,
                                                                                                 ccdbinx, ccdbiny)

        mostrarresultados(['N', 'ccdbinY', 'ccbinX', 'A', 'B', '-1'],
                          [len(indice_seccion_[indice_seccion_ == seccion]), ccdbiny, ccdbinx,
                           naxis1_expected, naxis2_expected, nombre_archivo],
                          titulo='Bias Realizado')

        if cabecera is None:
            raise ValueError('No se ha creado correctamente la cabecera')

        masterbias_header = cabecera.copy()
        # if masterbias_header['BLANK']:
        #     del masterbias_header['BLANK']

        ahora = datetime.datetime.now()
        ahora_dt = Time(ahora, format='datetime', scale='utc')
        masterbias_header.add_history('Realizado Bias a partir de ' +
                                      str(len(indice_seccion_[indice_seccion_ == seccion])) + ' imagenes. | ' +
                                      str(ahora_dt)[:19])
        numero_bias = 0
        for bias_raw in range(len(lista_bias_)):
            if indice_seccion_[bias_raw] == seccion:
                numero_bias += 1
                masterbias_header.add_history(str(numero_bias) + ': ' + lista_bias_[bias_raw])
                print(str(numero_bias) + ': ' + lista_bias_[bias_raw])

        masterbias_final = fits.PrimaryHDU(master_bias_colapsado.astype(np.float32), masterbias_header)

        masterbias_final.writeto(dir_bias_ + nombre_archivo, overwrite=True)

        if verbose_imagen:
            coord_lim = ImP.limites_imagen(*coordenadas_dibujo)
            ImP.imgdibujar(master_bias_colapsado, *coordenadas_dibujo, *coord_lim, verbose_=1)

        if interactive:
            # plt.show()
            input("Press Enter to continue...")

        # Añadirlo a la tabla
        isot_time = masterbias_header['DATE']
        tiempo = Time(isot_time, format='isot', scale='utc')
        if elemento_lista is None:
            elemento_lista = pd.DataFrame([[naxis1_expected, naxis2_expected,
                                            x1, x2, y1, y2,
                                            ccdbinx, ccdbiny,
                                            nombre_archivo, noche,
                                            tiempo, tiempo.jd]],
                                          columns=['Naxis1', 'Naxis2',
                                                   'x1', 'x2', 'y1', 'y2',
                                                   'Binning1', 'Binning2',
                                                   'nombre_archivo', 'noche',
                                                   'tiempo_astropy', 'julian'])
        else:
            elemento_lista_ = pd.DataFrame([[naxis1_expected, naxis2_expected,
                                            x1, x2, y1, y2,
                                            ccdbinx, ccdbiny,
                                            nombre_archivo, noche,
                                            tiempo, tiempo.jd]],
                                           columns=['Naxis1', 'Naxis2',
                                                    'x1', 'x2', 'y1', 'y2',
                                                    'Binning1', 'Binning2',
                                                    'nombre_archivo', 'noche',
                                                    'tiempo_astropy', 'julian'])
            elemento_lista = pd.concat([elemento_lista, elemento_lista_], ignore_index=True)

    return elemento_lista


def realizar_master_biases(lista_noches, dir_listas, dir_datos, dir_bias, interactive, recortar,
                           verbose_imagen=False):
    i_noche = 0
    df_bias = None
    for noche in lista_noches:
        i_noche += 1
        print('=== NOCHE ' + noche + ' - (' + str(i_noche) + '/' + str(len(lista_noches)) + ') ===')

        lista_bias = leer_lista(dir_listas + noche + '/' + 'LBias.csv')

        secciones, secciones_unicas, secciones_count, indice_seccion, bin_secciones, _ = crear_lista_unicos(
            dir_datos, noche, lista_bias, cabecera='CCDSEC', binning=True
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

        # for k in range(len(secciones_unicas)):
        #     add_to_file('BiasSecciones.csv', secciones_unicas[k] + ';' + str(secciones_count[k]) + '\n')

        df_bias_ = juntar_imagenes_bias(noche, secciones_unicas, coordenadas_secciones, secciones_count,
                                        indice_seccion, bin_secciones, dir_bias, dir_datos, lista_bias,
                                        interactive=interactive, recortar=recortar,
                                        verbose_imagen=verbose_imagen)
        if df_bias is None:
            df_bias = df_bias_
        else:
            df_bias = pd.concat([df_bias, df_bias_], ignore_index=True)

        # print(df_bias)
        # input('listo')

    return df_bias


def juntar_imagenes_flats(noche, secciones_unicas_, coordenadas_secciones_, indice_seccion_,
                          dir_bias_, dir_datos_, dir_flats_, lista_flats_, lista_noches, lista_bias,
                          verbose=0, interactive=False, verbose_imagen=False):

    # Cargamos Diccionario de Filtros
    dic_filtro = cargar_json()
    elemento_lista = None

    # Separamos por secciones para cada noche
    for seccion in range(len(secciones_unicas_)):
        print('seccion: ' + str(seccion))
        coordenadas_dibujo = sacar_coordenadas_2(coordenadas_secciones_, seccion)
        x1, x2, y1, y2 = deshacer_tupla_coord(coordenadas_dibujo)

        # Creamos una lista con todos aquellos que coinciden el indice
        lista_coincide = []
        for imagen in range(len(lista_flats_)):
            if indice_seccion_[imagen] == seccion:
                lista_coincide.append(lista_flats_[imagen])

        filtros, filtros_unicos, filtros_count, indice_filtro, _, nombres_filtros = crear_lista_unicos(
                                                                                      dir_datos_, noche,
                                                                                      lista_coincide,
                                                                                      cabecera='INSFLID',
                                                                                      binning=False, nombre_filtro=True
                                                                                      )

        # Creamos nombres_unicos
        nombres_unicos = []
        for i in range(len(filtros_unicos)):
            for j in range(len(filtros)):
                if filtros_unicos[i] == filtros[j]:
                    nombres_unicos.append(nombres_filtros[j])
                    break

        # Para cada elemento de nombres únicos, buscamos cual es el índice en el diccionario de filtros
        nombre_diccionario = []
        for i in nombres_unicos:
            nombre_diccionario.append(leer_diccionario(i, dic_filtro))

        # Esto sólo será necesario si NO hay grisma, luego miramos cual es el grima que le corresponde
        numero_grisma = int(fits.open(dir_datos_ + noche + '/' + lista_coincide[0])[0]
                            .header['insgrid'].replace(' ', '0')[6:8])

        if numero_grisma == 11:  # Si numero_grisma==11, entonces no hay grisma de por medio
            free_grisma = True
        else:
            free_grisma = False

        # Queremos saber cual es el número del filtro usado (quizá se pueda quitar pues ahora se mira el nombre)
        numero_filtro = np.zeros(len(filtros_unicos), dtype=int)
        p = 0
        for i in filtros_unicos:
            numero_filtro[p] = int(i[-2:].strip())
            p += 1

        # Ahora que hemos separado en función de la sección, Separamos por los filtros
        for filtro in range(len(filtros_unicos)):
            lista_actual = []
            binning2 = []
            binning1 = []

            # Como no siempre mantienen constante los binnings, separamos por binnings
            for i in range(len(lista_coincide)):
                if indice_filtro[i] == filtro:
                    lista_actual.append(lista_coincide[i])
                    binning2.append(int(fits.open(dir_datos_ + noche + '/' + lista_coincide[i])[0].header['ccdbinX']))
                    binning1.append(int(fits.open(dir_datos_ + noche + '/' + lista_coincide[i])[0].header['ccdbinY']))

            # bin2_unique, bin2_count = np.unique(binning2, return_counts=True)
            bin1_unique, bin1_count = np.unique(binning1, return_counts=True)

            siempre_coincide_binning = True
            for i in range(len(lista_actual)):
                if not binning2 == binning1:
                    siempre_coincide_binning = False
                    break
            # solo_un_caso_de_binning = all([len(bin1_unique) == 1, len(bin2_unique) == 1])
            # if solo_un_caso_de_binning:  ############################################################################

            # Si los binnings coinciden
            master_flats_colapsado = None

            if siempre_coincide_binning:
                for binning in bin1_unique:
                    lista_actual_bin = []

                    # Ahora tenemos una lista con binnings individuales
                    for i in range(len(lista_actual)):
                        if binning1[i] == binning:
                            lista_actual_bin.append(lista_actual[i])

                    # Como tienen un solo bin y coinciden en los ejes (no se estira la imagen),
                    # no nos preocupamos y generamos las imagenes
                    ccdbiny = binning
                    ccdbinx = binning

                    naxis1_expected, naxis2_expected = sacar_naxis(coordenadas_dibujo, ccdbinx, ccdbiny)

                    master_flats = np.zeros((filtros_count[filtro], naxis1_expected, naxis2_expected), dtype=float)

                    # Generamos máscara
                    center_real = [1075, 1040]
                    center = [center_real[0] - x1, center_real[1] - y1]
                    radius = 809  # con 810 tiene un pixel de borde
                    mask = create_circular_mask(naxis1_expected, naxis2_expected, center=center, radius=radius)

                    mostrarresultados(['N', 'ccdbiny', 'ccdbinx', 'A', 'B'],
                                      [len(lista_actual_bin), ccdbiny, ccdbinx,
                                       naxis1_expected, naxis2_expected],
                                      titulo='Filtro ' + str(nombre_diccionario[filtro]))

                    # Leemos para cada elemento de imágenes sus valores y los vamos acoplando
                    indice0 = 0
                    cabecera = fits.open(dir_datos_ + noche + '/' + lista_actual_bin[0])[0].header
                    for imagen in range(len(lista_actual_bin)):
                        image_file = dir_datos_ + noche + '/' + lista_actual_bin[imagen]
                        image_data = fits.getdata(image_file, ext=0)

                        if image_data[:, :].shape == master_flats[indice0, :, :].shape:
                            master_flats[indice0, :, :] = image_data[:, :]                  # Juntar
                        else:
                            print('hay un problema con el tamanyo de las imagenes')
                            input('Pausado. Pulsa Enter para continuar...')
                        indice0 += 1

                    existe, bias_asociado_nombre, bias_asociado = obtener_bias(dir_bias_, noche, lista_noches,
                                                                               lista_bias, x1, x2, y1, y2,
                                                                               ccdbinx, ccdbiny)

                    for i in range(master_flats.shape[0]):
                        master_flats[i, :, :] = master_flats[i, :, :] - bias_asociado

                    # Normalizamos
                    valor_medio = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):

                        # Si tiene grisma, usamos una máscara para calcular la mediana
                        if free_grisma:
                            x = ma.masked_array(master_flats[i, :, :], ~mask)
                            valor_medio[i] = ma.median(x)
                        else:
                            valor_medio[i] = np.median(master_flats[i, :, :])

                        if valor_medio[i] <= 0:
                            # raise ValueError('El valor de la mediana no debe ser negativo ni 0.')

                            if valor_medio[i] == 0:
                                valor_medio[i] = 1
                            else:
                                valor_medio[i] = abs(valor_medio[i])

                        # Después ya dividimos
                        master_flats[i, :, :] = np.true_divide(master_flats[i, :, :], valor_medio[i])

                    valor_medio2 = np.zeros(master_flats.shape[0], dtype=float)
                    for i in range(master_flats.shape[0]):
                        valor_medio2[i] = np.mean(master_flats[i, :, :], dtype=float)

                    # Colapsamos
                    master_flats_colapsado = np.median(master_flats, axis=0)

                    # las esquinas las definimos como 1 para evitar problemas de dividir por 0
                    if free_grisma:
                        master_flats_colapsado[~mask] = 1

                    # plt.imshow(master_flats_colapsado)
                    # plt.show()
                    if verbose_imagen:
                        ImP.imgdibujar(master_flats_colapsado)
                        plt.show()

                    # Generamos el nombre del fichero del flat que se va a generar
                    # si tiene numero_filtro[filtro] da la posicion del filtro
                    # si pone nombre_diccionario[filtro] da el numero del diccionario del filtro
                    nombre_archivo = noche +\
                        "-{0:04d}_{1:04d}_{2:04d}_{3:04d}-B{4:02d}_{5:02d}-F{6:03d}.fits"\
                        .format(x1, x2, y1, y2, ccdbinx, ccdbiny, int(nombre_diccionario[filtro]))

                    masterflats_header = cabecera.copy()
                    if masterflats_header['BLANK']:
                        del masterflats_header['BLANK']

                    ahora = datetime.datetime.now()
                    ahora_dt = Time(ahora, format='datetime', scale='utc')
                    masterflats_header.add_history('Realizado Flat a partir de ' +
                                                   str(len(indice_seccion_[indice_seccion_ == seccion])) +
                                                   ' imagenes. | ' + str(ahora_dt)[:19])
                    numero_flats = 0
                    for flat_raw in range(len(lista_flats_)):
                        if indice_seccion_[flat_raw] == seccion:
                            numero_flats += 1
                            masterflats_header.add_history(str(numero_flats) + ': ' + lista_flats_[flat_raw])
                            print(str(numero_flats) + ': ' + lista_flats_[flat_raw])

                    if master_flats_colapsado is None:
                        raise ValueError('No se ha creado correctamente el flat colapsado')
                    masterflats_final = fits.PrimaryHDU(master_flats_colapsado.astype(np.float32), masterflats_header)

                    masterflats_final.writeto(dir_flats_ + nombre_archivo, overwrite=True)

                    # Añadirlo a la tabla
                    isot_time = masterflats_header['DATE']
                    tiempo = Time(isot_time, format='isot', scale='utc')
                    if elemento_lista is None:
                        elemento_lista = pd.DataFrame([[naxis1_expected, naxis2_expected,
                                                        x1, x2, y1, y2,
                                                        ccdbinx, ccdbiny,
                                                        int(nombre_diccionario[filtro]),
                                                        free_grisma, numero_grisma,
                                                        nombre_archivo, noche,
                                                        tiempo, tiempo.jd]],
                                                      columns=['Naxis1', 'Naxis2',
                                                               'x1', 'x2', 'y1', 'y2',
                                                               'Binning1', 'Binning2',
                                                               'filtro',
                                                               'free_grisma', 'num_grisma',
                                                               'nombre_archivo', 'noche',
                                                               'tiempo_astropy', 'julian'])
                    else:
                        elemento_lista_ = pd.DataFrame([[naxis1_expected, naxis2_expected,
                                                        x1, x2, y1, y2,
                                                        ccdbinx, ccdbiny,
                                                        int(nombre_diccionario[filtro]),
                                                        free_grisma, numero_grisma,
                                                        nombre_archivo, noche,
                                                        tiempo, tiempo.jd]],
                                                       columns=['Naxis1', 'Naxis2',
                                                                'x1', 'x2', 'y1', 'y2',
                                                                'Binning1', 'Binning2',
                                                                'filtro',
                                                                'free_grisma', 'num_grisma',
                                                                'nombre_archivo', 'noche',
                                                                'tiempo_astropy', 'julian'])
                        elemento_lista = pd.concat([elemento_lista, elemento_lista_], ignore_index=True)

            else:
                raise ValueError('No Coinciden los binning en ambos ejes')

            if verbose >= 1:
                coord_lim = ImP.limites_imagen(*coordenadas_dibujo)
                ImP.imgdibujar(master_flats_colapsado, *coordenadas_dibujo, *coord_lim, verbose_=1)

            if interactive:
                input("Press Enter to continue...")

    return elemento_lista


def realizar_master_flats(lista_noches, lista_bias, dir_listas, dir_datos, dir_bias, dir_flats,
                          verbose, interactive, verbose_imagen=False):

    """
    Función maestra para generar los master flats.

    :param lista_noches:
    :param lista_bias:
    :param dir_listas:
    :param dir_datos:
    :param dir_bias:
    :param dir_flats:
    :param verbose:
    :param interactive:
    :param verbose_imagen:
    :return:
    """

    i_noche = 0
    df_flat = None
    for noche in lista_noches[:]:
        i_noche += 1
        print('=== NOCHE ' + noche + ' - (' + str(i_noche) + '/' + str(len(lista_noches)) + ') ===')

        lista_flats = leer_lista(dir_listas + noche + '/' + 'LFlat.csv')

        secciones, secciones_unicas, secciones_count, indice_seccion, bin_secciones, _ = crear_lista_unicos(
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

        df_flat_ = juntar_imagenes_flats(noche, secciones_unicas, coordenadas_secciones, indice_seccion,
                                         dir_bias, dir_datos, dir_flats, lista_flats, lista_noches, lista_bias,
                                         verbose=verbose, interactive=interactive, verbose_imagen=verbose_imagen)

        if df_flat is None:
            df_flat = df_flat_
        else:
            df_flat = pd.concat([df_flat, df_flat_], ignore_index=True)

        # print(df_flat)
        # input('listo')
    return df_flat


def realizar_reduccion(lista_noches, dir_listas, dir_datos, dir_bias, dir_flats, dir_reducc,
                       df_bias, df_flat, verbose=2, verbose_imagen=False):
    # Cargar listas
    elemento_lista = None
    dic_filtro = cargar_json()
    no_existen = []
    imagenes_totales_de_ciencia = 0
    imagenes_guardadas = 0

    for noche in lista_noches:
        imagenes_reducidas_noche = 0
        print(noche)
        if noche not in os.listdir(dir_reducc):
            os.mkdir(dir_reducc + noche + '/')

        lista_ciencia = leer_lista(dir_listas + noche + '/' + 'SCI.csv')
        no_existe = (0, 0, len(lista_ciencia))
        secc, secc_unicas, secc_count, indice_secc, bin_secc, nombres_filtros = crear_lista_unicos(
            dir_datos, noche, lista_ciencia, cabecera='CCDSEC', binning=True, nombre_filtro=True
        )
        i_imagen = 0
        for imagen in range(len(lista_ciencia)):
            i_imagen += 1
            # Nombre de la imagen
            nombre_ciencia = dir_datos + noche + '/' + lista_ciencia[imagen]
            cabecera = fits.open(nombre_ciencia)[0].header

            # Coordenadas y binning
            coordenadas = sacar_coordenadas_ccd(secc_unicas[indice_secc[imagen]])
            x1, x2, y1, y2 = deshacer_tupla_coord(coordenadas)
            binning = bin_secc[indice_secc[imagen]]
            naxis1_r = cabecera['Naxis1']
            naxis2_r = cabecera['Naxis2']

            naxis1_ciencia, naxis2_ciencia = sacar_naxis((x1, x2, y1, y2), binning[1], binning[0])

            # Comprobamos si hay overscan
            biassec = cabecera['BIASSEC']
            coordenadas_biassec = sacar_coordenadas_ccd(biassec)
            naxis1_overscan, naxis2_overscan = sacar_naxis(coordenadas_biassec, binning[1], binning[0])

            if coordenadas_biassec[0] == 0 and coordenadas_biassec[1] == 0:
                overscan = False
                x_1, x_2, y_1, y_2 = None, None, None, None
            else:
                print('Hay overscan!')
                overscan = True

                # Hay que recortar. Generamos las secciones de la imagen con la que nos vamos a quedar
                if coordenadas_biassec[0] > x2:
                    print('derecha')
                    print('Zona a recortar:')
                    x_1 = 0
                    x_2 = naxis2_ciencia
                    y_1 = 0
                    y_2 = naxis1_ciencia
                elif coordenadas_biassec[3] > y2:
                    print('encima')
                    x_1 = 0
                    x_2 = naxis2_ciencia
                    y_1 = 0
                    y_2 = naxis1_ciencia
                elif coordenadas_biassec[1] < x1:
                    print('izquierda')
                    x_1 = naxis2_overscan + 1
                    x_2 = naxis2_ciencia + naxis2_overscan
                    y_1 = 0
                    y_2 = naxis1_ciencia
                elif coordenadas_biassec[4] < y1:
                    print('abajo')
                    x_1 = 0
                    x_2 = naxis2_ciencia
                    y_1 = naxis1_overscan + 1
                    y_2 = naxis1_ciencia + naxis2_overscan
                else:
                    raise ValueError('No hay combinacion para donde esta el overscan!')

            # Obtenemos el nombre del filtro y su ID_filtro
            nombrefiltro = nombres_filtros[imagen]
            id_filtro = leer_diccionario(nombrefiltro, dic_filtro)

            # Fecha Juliana
            isot_time = cabecera['DATE']
            fecha = str2datetime(isot_time)
            tiempo = Time(isot_time, format='isot', scale='utc')

            if verbose > 0:
                mostrarresultados(['noche', 'naxis1_r', 'naxis2_r', 'naxis1_c', 'naxis2_c',
                                   'x1', 'x2', 'y1', 'y2',
                                   'binning', 'filtro', 'id_filtro',
                                   'FJM', 'Overscan',
                                   'imagen', 'nombre',
                                   'fecha', 'hora'],
                                  [noche, naxis1_r, naxis2_r, naxis1_ciencia, naxis2_ciencia,
                                   x1, x2, y1, y2,
                                   binning, nombrefiltro, id_filtro,
                                   tiempo.jd, overscan,
                                   imagen, lista_ciencia[imagen],
                                   fecha.date(), fecha.time()],
                                  contador=i_imagen, valor_max=len(lista_ciencia))

            # Buscamos el Bias
            # Seleccionamos los que coinciden con el tamaño teorico
            datat = df_bias[df_bias.Naxis1 == naxis1_ciencia]
            datat = datat[datat.Naxis2 == naxis2_ciencia]
            datat = datat[datat.Binning1 == binning[0]]
            datat = datat[datat.Binning2 == binning[1]]
            fechasjd_b = abs(datat.julian.values - tiempo.jd)
            if len(fechasjd_b) > 0:
                pos_min = np.argmin(fechasjd_b)
                nombre_bias_buscado = datat.iloc[pos_min]['nombre_archivo']
                if datat.iloc[pos_min]['noche'] != noche:
                    if verbose > 1:
                        print('El bias mas cercano no es de esta noche.')

                bias_asociado = fits.getdata(dir_bias + nombre_bias_buscado, ext=0)
            else:
                nombre_bias_buscado = 'Ausente'
                relleno_b = 680
                if verbose > 1:
                    print('No se han encontrado bias de esa forma, se genera uno artificial.')
                bias_asociado = np.full((naxis1_ciencia, naxis2_ciencia), relleno_b, dtype=float)

            # Buscamos el Flat
            # Seleccionamos los que coinciden con el tamaño teorico
            dataf = df_flat[df_flat.Naxis1 == naxis1_ciencia]
            dataf = dataf[dataf.Naxis2 == naxis2_ciencia]
            dataf = dataf[dataf.Binning1 == binning[0]]
            dataf = dataf[dataf.Binning2 == binning[1]]
            dataf = dataf[dataf.filtro == id_filtro]

            # Lista de los bias en función de la distancia en tiempo a la imagen de ciencia
            fechasjd_f = abs(dataf.julian.values - tiempo.jd)
            if len(fechasjd_f) > 0:
                pos_min = np.argmin(fechasjd_f)
                nombre_flat_buscado = dataf.iloc[pos_min]['nombre_archivo']
                if dataf.iloc[pos_min]['noche'] != noche:
                    if verbose > 1:
                        print('Existe un Flat, pero no es de esta noche.')
                flat_asociado = fits.getdata(dir_flats + nombre_flat_buscado, ext=0)

            else:
                nombre_flat_buscado = 'Ausente'
                relleno_f = 1
                print('No se han encontrado flats, se genera uno artificialmente.')
                flat_asociado = np.full((naxis1_ciencia, naxis2_ciencia), relleno_f, dtype=float)

            # Obtenemos la información de la imagen de ciencia
            image_data = fits.getdata(nombre_ciencia, ext=0)
            reducido_header = cabecera.copy()

            # Comprobamos si hay o no overscan que haya que tener en cuenta
            if overscan:
                print('Recortamos los datos')
                image_data = image_data[y_1:y_2, x_1:x_2]

            ##############################################################################
            reducido_datos = (image_data - bias_asociado) / flat_asociado
            ##############################################################################

            if verbose_imagen:
                ImP.imgdibujar(reducido_datos)
                plt.show()

            if reducido_header['BLANK']:
                del reducido_header['BLANK']

            # Sacar si tiene o no un grisma
            numero_grisma = reducido_header['insgrid'].replace(' ', '0')[6:8]
            if numero_grisma == 11:  # Si numero_grisma==11, entonces no hay grisma de por medio
                free_grisma = True
            else:
                free_grisma = False

            if free_grisma:
                mask = create_circular_mask(naxis1_ciencia, naxis2_ciencia)  # Se puede cambiar el centro y radio
                reducido_datos[~mask] = 1

            ahora = datetime.datetime.now()
            ahora_dt = Time(ahora, format='datetime', scale='utc')
            reducido_header.add_history('Realizada la reduccion a partir de las siguientes imagenes: | ' +
                                        str(ahora_dt)[:19])
            reducido_header.add_history('Bias: ' + nombre_bias_buscado)
            reducido_header.add_history('Flat: ' + nombre_flat_buscado)

            # Guardamos la imagen
            reducido_final = fits.PrimaryHDU(reducido_datos.astype(np.float32), reducido_header)

            reducido_final.writeto(dir_reducc + noche + '/' + 'r_' + lista_ciencia[imagen], overwrite=True)

            imagenes_reducidas_noche += 1

            ahora = datetime.datetime.now()
            tiempo = Time(ahora, format='datetime', scale='utc')
            if elemento_lista is None:
                elemento_lista = pd.DataFrame([[naxis1_ciencia, naxis2_ciencia,
                                                x1, x2, y1, y2,
                                                binning[0], binning[1],
                                                id_filtro,
                                                free_grisma, numero_grisma,
                                                'r_' + lista_ciencia[imagen],
                                                nombre_bias_buscado, nombre_flat_buscado,
                                                noche, tiempo, tiempo.jd]],
                                              columns=['Naxis1', 'Naxis2',
                                                       'x1', 'x2', 'y1', 'y2',
                                                       'Binning1', 'Binning2',
                                                       'filtro',
                                                       'free_grisma', 'num_grisma',
                                                       'nombre_archivo',
                                                       'nombre bias', 'nombre flat',
                                                       'noche', 'Fecha realizacion', 'julian'])
            else:
                elemento_lista_ = pd.DataFrame([[naxis1_ciencia, naxis2_ciencia,
                                                 x1, x2, y1, y2,
                                                 binning[0], binning[1],
                                                 id_filtro,
                                                 free_grisma, numero_grisma,
                                                 'r_' + lista_ciencia[imagen],
                                                 nombre_bias_buscado, nombre_flat_buscado,
                                                 noche, tiempo, tiempo.jd]],
                                               columns=['Naxis1', 'Naxis2',
                                                        'x1', 'x2', 'y1', 'y2',
                                                        'Binning1', 'Binning2',
                                                        'filtro',
                                                        'free_grisma', 'num_grisma',
                                                        'nombre_archivo',
                                                        'nombre bias', 'nombre flat',
                                                        'noche', 'Fecha realizacion', 'julian'])
                elemento_lista = pd.concat([elemento_lista, elemento_lista_], ignore_index=True)

        # Al final de cada noche se hace el recuento
        no_existen.append(no_existe)
        imagenes_guardadas += imagenes_reducidas_noche
        print('verbosidad reduccion', verbose)
        if verbose > 0:
            mostrarresultados(['Imagenes reducidas', 'Acumulado'], [imagenes_reducidas_noche, imagenes_guardadas],
                              titulo="Reduccion de Imagenes",
                              contador=lista_noches.index(noche), valor_max=len(lista_noches))

    mostrarresultados(lista_noches, no_existen)
    no_existen_2 = [0, 0]
    _ = elemento_lista.to_csv('df_reducido.csv', index=None, header=True)
    for i in range(len(no_existen)):
        imagenes_totales_de_ciencia += no_existen[i][2]
        no_existen_2[0] += no_existen[i][0]
        no_existen_2[1] += no_existen[i][1]

    print('Biases que no ha funcionado: ', no_existen_2[0], '| Flats que no ha funcionado: ', no_existen_2[1],
          '| Imagenes en total: ', imagenes_totales_de_ciencia, '| Imagenes reducidas: ', imagenes_guardadas)

    return imagenes_guardadas


def pruebas_pandas(data):
    print(data.head())
    # print(data[data.noche == '170225_t2_CAFOS'])
    datat = data[data.Naxis1 == 2048]
    datat = datat[datat.Naxis2 == 1000]
    print(datat)
    print(data.nombre_archivo.values[0])


def decidir_repetir_calculos(norealizar, sirealizar, sujeto, dir_df, dir_sujetos):
    existe_bias = os.path.exists(dir_df + 'df_' + sujeto + '.csv')
    lista_existentes = os.listdir(dir_sujetos)
    if not any([norealizar, sirealizar]):
        if existe_bias:
            print('Se va a coger una version ya existente de los ' + sujeto)
            importado = pd.read_csv(dir_df + 'df_' + sujeto + '.csv')
            lista_teorica = importado.nombre_archivo.values.tolist()

            contador_t = 0
            contador_e = 0
            for i in range(len(lista_existentes)):
                if lista_teorica[i] in lista_existentes:
                    contador_t += 1
                if lista_existentes[i] in lista_teorica:
                    contador_e += 1

            if contador_e == len(lista_existentes) and contador_t == len(lista_teorica):
                print('Estan todos contabilizados, no hace falta volver a calcular los ' + sujeto)
                realizar = False
            else:
                print('No estan todos los ' + sujeto + ', se vuelven a calcular')
                realizar = True
        else:
            print('No existe el dataframe, se calcula el dataframe de los ' + sujeto)
            realizar = True
    elif norealizar:
        if existe_bias:
            print('Existe el dataframe de ' + sujeto + '. No se repite')
            realizar = False
        else:
            print('No existe el dataframe de ' + sujeto + ', se fuerza que se vuelva a hacer')
            realizar = True
    elif sirealizar:
        if existe_bias:
            print('Existe el dataframe de ' + sujeto + '. Se fuerza calcularlo de nuevo')
            realizar = True
        else:
            print('No existe el dataframe de ' + sujeto + ', Se fuerza que se haga')
            realizar = True
    else:
        raise ValueError('No cuadran las condiciones para los ' + sujeto + '!')

    return realizar


def main():

    # ---------------Valores por defecto-------------------------------------------
    default_dir_datos = '/media/enrique/TOSHIBA EXT/DataCAHA/CAFOS2017/'
    default_dir_bias = '/media/enrique/TOSHIBA EXT/CAHA/Biases/'
    default_dir_listas = '/media/enrique/TOSHIBA EXT/CAHA/Listas/'
    default_dir_flats = '/media/enrique/TOSHIBA EXT/CAHA/Flats/'
    default_dir_reduccion = '/media/enrique/TOSHIBA EXT/CAHA/Reduccion/'
    default_dir_dataframe = ''
    desc_bias = ['bias', 'Bias', 'BIAS']
    desc_flats = ['flats', 'FLATS', 'FLAT', 'Flats', 'Flat', 'flat', 'Skyflat', 'SDkyflat']
    desc_arc = ['arc', 'ARC']
    # -----------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description="Bias and Flat calibration of CAFOS images")
    group = parser.add_mutually_exclusive_group()
    grupo2 = parser.add_mutually_exclusive_group()
    grupo3 = parser.add_mutually_exclusive_group()

    group.add_argument("-v", "--verbose", action="store_true")
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-db", "--dir_bias", default=default_dir_bias, type=str, help='Bias Directory')
    parser.add_argument("-df", "--dir_flats", default=default_dir_flats, type=str, help='Flats Directory')
    parser.add_argument("-dd", "--dir_datos", default=default_dir_datos, type=str, help='Data Directory')
    parser.add_argument("-dl", "--dir_listas", default=default_dir_listas, type=str, help='Lists Directory')
    parser.add_argument("-de", "--dir_reducc", default=default_dir_reduccion, type=str, help='Reducction Directory')
    parser.add_argument("-ddf", "--dir_dataf", default=default_dir_dataframe, type=str, help='DataFrame Directory')
    parser.add_argument("-vi", "--verboseimage", action="store_true", help="Mostrar Imagenes")
    parser.add_argument('--cmap', type=str, help="Colormap", default='hot')
    parser.add_argument("-i", "--interactive", action="store_true")
    parser.add_argument("--recortar", action="store_true", help="Activar el recorte de imagenes")
    parser.add_argument("--calysci", action="store_false",
                        help="Usar cuando los archivos no tienen '-cal-' y '-sci-'"
                             + "en el nombre para diferenciar entre calibración y ciencia.")
    parser.add_argument("-nr", "--noreducc", action="store_false", help="No realizar la reduccion.")

    grupo2.add_argument("-nb", "--nobias", action="store_true",
                        help="No realizar los Master Bias.")
    grupo2.add_argument("-sb", "--sibias", action="store_true",
                        help="Fuerza realizar los Master Bias.")

    grupo3.add_argument("-nf", "--noflat", action="store_true",
                        help="No realizar los Master Flat.")
    grupo3.add_argument("-sf", "--siflat", action="store_true",
                        help="Fuerza realizar los Master Flat.")

    args = parser.parse_args()

    if args.verbose:
        verbosidad = 2
    elif args.quiet:
        verbosidad = 0
    else:
        verbosidad = 1
    print('verbosidad', verbosidad)

    # Comprobamos si queremos/hace falta calcular los bias/flats
    print('bias:', args.nobias, args.sibias)
    realizarbias = decidir_repetir_calculos(args.nobias, args.sibias, 'bias', args.dir_dataf, args.dir_bias)
    print('flat:', args.noflat, args.siflat)
    realizarflat = decidir_repetir_calculos(args.noflat, args.siflat, 'flat', args.dir_dataf, args.dir_flats)

    # Creamos una lista de las noches disponibles
    lista_noches = os.listdir(args.dir_datos)
    lista_noches.sort()
    tiempo_inicio = time.time()

    # Separamos entre calibración y ciencia
    crear_listas_cal_y_sci(lista_noches, args.dir_listas, args.dir_datos, desc_bias, desc_flats, desc_arc,
                           verbosidad, args.calysci)
    tiempo_listas = time.time()

    print(args.nobias, args.noflat, args.noreducc)

    # importado_b = pd.read_csv('df_bias.csv')
    # importado_f = pd.read_csv('df_flat.csv')
    # pruebas_pandas(importado_f)

    # Creamos los Master Biases
    if realizarbias:
        df_bias = realizar_master_biases(lista_noches, args.dir_listas, args.dir_datos, args.dir_bias,
                                         args.interactive, args.recortar, verbose_imagen=args.verboseimage)
        numero_bias = len(os.listdir(args.dir_bias))
        print(df_bias)
        _ = df_bias.to_csv('df_bias.csv', index=None, header=True)
    else:
        df_bias = pd.read_csv('df_bias.csv')
        print('Se han importado los bias')
        numero_bias = '-'

    lista_bias = conseguir_listas_archivos(args.dir_bias)
    tiempo_biases = time.time()

    # Creamos los Master Flats
    if realizarflat:
        df_flat = realizar_master_flats(lista_noches, lista_bias,
                                        args.dir_listas, args.dir_datos, args.dir_bias, args.dir_flats,
                                        verbosidad, args.interactive, verbose_imagen=args.verboseimage)
        numero_flats = len(os.listdir(args.dir_flats))
        print(df_flat)
        _ = df_flat.to_csv('df_flat.csv', index=None, header=True)
    else:
        df_flat = pd.read_csv('df_flat.csv')
        print('Se han importado los flats')
        numero_flats = '-'

    tiempo_flats = time.time()

    # Juntamos todos los procesos y relizamos la reducción
    if args.noreducc:
        numeros_reducidos = realizar_reduccion(lista_noches, args.dir_listas, args.dir_datos, args.dir_bias,
                                               args.dir_flats, args.dir_reducc, df_bias, df_flat,
                                               verbosidad, verbose_imagen=args.verboseimage)
    else:
        numeros_reducidos = '-'
    tiempo_reducc = time.time()

    # Mostramos resultados de ambos procesos
    mostrarresultados(['Tiempo Listas', 'Tiempo Master Bias', 'Tiempo Master Flats', 'Tiempo Reduccion', 'Tiempo Total',
                       'Cuantos Biases', 'Cuantos Flats', 'Cuantos Reducidos'],
                      [round(tiempo_listas - tiempo_inicio, 2), round(tiempo_biases - tiempo_listas, 2),
                       round(tiempo_flats - tiempo_biases, 2), round(tiempo_reducc - tiempo_flats, 2),
                       round(tiempo_reducc - tiempo_listas, 2),
                       numero_bias, numero_flats, numeros_reducidos],
                      titulo='Tiempo Ejecucion')


if __name__ == "__main__":

    main()
